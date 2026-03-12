%% =============================================================================
% Description:
%{
    TITLE: PHASE FRACTION CALCULATION (FUNCTION)
    AUTHOR: DERIAN R. MORPHEW / PATRICK J. SHAMBERGER
    VERSION: 1 (2026)  

    SUMMARY:
        [Calculates phase fraction values for the points in each path.]

    INPUTS:
        - Paths structure containing all paths.
        - MatProps structure containing phase boundary points.
        - Plotting structure containing grid points.
        - Boolean logic to display plots (1 - Yes, 0 - No).

    OUTPUTS:
        - Updated FullTable in the Paths structure for each path, populated
          with phase fraction values.

    NOTES:
        - Austenite phase = 1
        - Martensite phase = 0
%}
%  =============================================================================
function Paths = phasepop(Paths, MatProps, Plotting, Plot)

    % === Extract Required Inputs ===
    TG = Plotting.TempGrid;
    StG = Plotting.StressGrid; 
    T_ax = TG(1,:);
    St_ax = StG(:,1);

    PB = MatProps.Phaseboundary;
    val = [PB.Stress] > 0;
    PBSt = [PB(val).Stress]';
    PBAs = [PB(val).As]';
    PBAf = [PB(val).Af]';
    PBMs = [PB(val).Ms]';
    PBMf = [PB(val).Mf]';

    % === Loop over each path ===
    pns = fieldnames(Paths);
    pns = pns(contains(pns, 'Path_'));

    for i = 1:length(pns)
        pn = pns{i};
        FT = Paths.(pn).FullTable;

        T = FT.Temperature;
        St = FT.Stress;
        PF = NaN(height(FT), 1);

        %% First Point
        T0 = T(1);
        St0 = St(1);
        As_0 = interp1(PBSt, PBAs, St0, 'linear');
        Af_0 = interp1(PBSt, PBAf, St0, 'linear');

        if T0 > Af_0
            PF(1) = 1; % Austenite

        elseif T0 > As_0 && T0 <= Af_0
            StrainGridMA = path.StrainShifted_MA;

            [~, St_i0] = min(abs(St_az - St0));
            val_i = find(~isnan(StrainGridMA(St_io, :)));
            [~, T_i0] = min(abs(T_ax(val_i) - T0));
            
            % === Strains ===
            eps_0 = StrainGridMA(St_i0, T_i0);
            epsAs_0 = StrainGridMA(St_i0, val_i(1));
            epsAf_0 = StrainGridMA(St_i0, val_i(end));

            % === Compute Phase Fraction ===
            PF(1) = (eps_0 - epsAs_0) / (epsAf_0 - epsAs_0);
        end

        % === Interpolate Phase Boundaries for all Path Points ===
        Ms_k = interp1(PBSt, PBMs, St, 'linear');
        Mf_k = interp1(PBSt, PBMf, St, 'linear');
        As_k = interp1(PBSt, PBAs, St, 'linear');
        Af_k = interp1(PBSt, PBAf, St, 'linear');

        %% CASE 1: Austenite region
        %  =====================================================================
        %{
             Conditions:
               - T > Ms
               - (σ increasing) or (σ constant & T decreasing)
             Populate:
               - Repeat initial phase fraction value
        %}
        %  =====================================================================

        % === Conditions ===
        St_cons = abs([0; diff(St)]) < 1e-8;
        T_dec = [false; diff(T) < 0];
        St_inc = [false; diff(St) > 0];
        motion1 = St_inc | (St_cons & T_dec);

        % === Populate ===
        c1_i = motion1 & (T >= Ms_k);
        PF(c1_i) = PF(1);

        %% CASE 2: Austenite to Martensite Region
        %  =====================================================================
        %{
             Conditions:
               - Ms >= T >= Mf
               - (σ increasing) or (σ constant before unloading begins in
                 an eCE type cycle)
               - If not an eCE type cycle (ISHC prediction), then (σ
                 constant & T decreasing)
             Populate:
               - Calculate phase fraction value using StrainGridAM
        %}
        %  =====================================================================
        
        % === Conditions ===
        Tcon_c2  = (T <= Ms_k) & (T >= Mf_k);
        dS = [0; diff(St)];
        dT = [0; diff(T)];
        St_cons = abs(dS) < 1e-8;
        St_inc  = dS > 0;
        T_dec = dT < 0;
        i_drop = find(dS < 0, 1, 'first'); % First stress decrease
        k_i = (1:length(St))';
        if isempty(i_drop)
            motion2 = St_inc | (St_cons & T_dec);
        else
            motion2 = St_inc | (St_cons & (k_i <= i_drop));
        end

        % === Populate ===
        c2_i = find(Tcon_c2 & motion2);
        StrainGridAM = Paths.(pn).StrainShifted_AM;

        T_c2 = T(c2_i);
        St_c2 = St(c2_i);

        eps_k_c2 = NaN(length(c2_i), 1);
        eps_Ms_c2 = NaN(length(c2_i), 1);
        eps_Mf_c2 = NaN(length(c2_i), 1);

        for ii = 1:length(c2_i)
            [~, St_i2] = min(abs(St_ax - St_c2(ii)));
            val_i = find(~isnan(StrainGridAM(St_i2, :)));
            [~, T_i2] = min(abs(T_ax(val_i) - T_c2(ii)));
            eps_k_c2(ii) = StrainGridAM(St_i2, val_i(T_i2));
            eps_Ms_c2(ii) = StrainGridAM(St_i2, val_i(end));
            eps_Mf_c2(ii) = StrainGridAM(St_i2, val_i(1));
        end

        PF_calc = 1 - (eps_k_c2 - eps_Ms_c2) ./ (eps_Mf_c2 - eps_Ms_c2);

        for ii = 1:length(c2_i)
            k = c2_i(ii);
                if k > 1 && ~isnan(PF(k-1)) && PF_calc(ii) <= PF(k-1)
                    PF(k) = PF_calc(ii);
                elseif k == 1
                    PF(k) = PF_calc(ii);
                else
                    PF(k) = PF(k-1);
                end
        end

        %% CASE 3: M Detwinning Region (Forward)
        %  =====================================================================
        %{
             Conditions:
               - T <= Mf
               - (σ increasing) or (σ constant & T decreasing)
             Populate:
               - Phase fraction is 0
        %}
        %  =====================================================================
        
        % === Conditions ===
        Tcon_c3 = (T <= Mf_k);
        St_cons = abs([0; diff(St)]) < 1e-8;
        St_inc = [false; diff(St) > 0];
        temp_decr = [false; diff(T) < 0];
        motion3 = St_inc | (St_cons & temp_decr);

        % === Populate ===
        c3_i = Tcon_c3 & motion3;
        PF(c3_i) = 0;

        %% CASE 4: M Detwinned Region (Reverse)
        %  =====================================================================
        %{
             Conditions:
               - T <= As
               - (σ decreasing) or (σ constant & T increasing)
             Populate:
               - Repeat previous phase fraction value
        %}
        %  =====================================================================
        
        % === Conditions ===
        Tcon_c4 = (T <= As_k);
        St_cons = abs([0; diff(St)]) < 1e-8;
        St_dec = [false; diff(St) < 0];
        T_inc = [false; diff(T) > 0];
        motion4 = St_dec | (St_cons & T_inc);

        % === Populate ===
        c4_i = find(Tcon_c4 & motion4);
        for ii = 1:length(c4_i)
            k = c4_i(ii);
            if k > 1 && ~isnan(PF(k-1))
                PF(k) = PF(k-1);
            end
        end

        %% === CASE 5: M → A Transition Region ===
        %  =====================================================================
        %{
             Conditions:
               - Af >= T >= As
               - (σ decreasing) or (σ constant after unloading begins in
                 an eCE type cycle)
               - If not an eCE type cycle (ISHC prediction), then (σ
                 constant & T increasing)
             Populate:
               - Calculate phase fraction value using StrainGridMA
        %}
        %  =====================================================================

        % === Conditions ===
        Tcon_c5 = (T >= As_k) & (T <= Af_k);
        dS = [0; diff(St)];
        dT = [0; diff(T)];
        St_cons = abs(dS) < 1e-8;
        St_dec  = dS < 0;
        T_inc = dT > 0;
        i_drop = find(dS < 0, 1, 'first'); % First stress decrease
        k_i = (1:length(St))';
        if isempty(i_drop)
            motion5 = St_dec | (St_cons & T_inc);
        else
            motion5 = St_dec | (St_cons & (k_i > i_drop));
        end

        % === Populate ===
        c5_i = find(Tcon_c5 & motion5);
        StrainGridMA = Paths.(pn).StrainShifted_MA;

        T_c5  = T(c5_i);
        St_c5  = St(c5_i);

        % Linear interpolation
        eps_k_c5  = NaN(length(c5_i), 1);
        eps_As_c5 = NaN(length(c5_i), 1);
        eps_Af_c5 = NaN(length(c5_i), 1);

        for ii = 1:length(c5_i)
            [~, St_i5] = min(abs(St_ax - St_c5(ii)));
            val_i = find(~isnan(StrainGridMA(St_i5, :)));
            [~, T_i5] = min(abs(T_ax(val_i) - T_c5(ii)));
            eps_k_c5(ii) = StrainGridMA(St_i5, val_i(T_i5));
            eps_As_c5(ii) = StrainGridMA(St_i5, val_i(1));
            eps_Af_c5(ii) = StrainGridMA(St_i5, val_i(end));
        end

        PF_calc = (eps_k_c5 - eps_As_c5) ./ (eps_Af_c5 - eps_As_c5);

        for ii = 1:length(c5_i)
            k = c5_i(ii);
                if k > 1 && ~isnan(PF(k-1)) && PF_calc(ii) >= PF(k-1)
                    PF(k) = PF_calc(ii);
                elseif k == 1
                    PF(k) = PF_calc(ii);
                else
                    PF(k) = PF(k-1);
                end
        end

        %% CASE 6: Final Austenite Region 
        %  =====================================================================
        %{
             Conditions:
               - T > Af
               - Phase fraction is 1
        %}
        %  =====================================================================
        c6_i = T > Af_k;
        PF(c6_i) = 1;

        %% Store
        FT.PhaseFraction = PF;
        Paths.(pn).FullTable = FT;
    end

    %% === Optional Plots ===
    if Plot == 1

        colormap = [
            161 218 180
            65 182 196
            44 127 184
            37  52 148
            ] / 255;

        nP = length(pns);
        colors = interp1(linspace(0,1,4), colormap, linspace(0,1,nP));

        % === Phase Fraction vs Temperature ===
        figure; 
        hold on;
        xlabel('Temperature (K)');
        ylabel('Phase Fraction');
        for i = 1:length(pns)
            pn = pns{i};
            FT = Paths.(pn).FullTable;
            T = FT.Temperature;
            PF = FT.PhaseFraction;
            plot(T, PF, '-', 'LineWidth', 1, 'Color', colors(i,:));
        end
        hold off;

        % === Phase Fraction vs Stress ===
        figure; 
        hold on;
        xlabel('Stress (MPa)');
        ylabel('Phase Fraction');
        for i = 1:length(pns)
            pn = pns{i};
            FT = Paths.(pn).FullTable;
            S  = FT.Stress;
            PF = FT.PhaseFraction;
            plot(S, PF, '-', 'LineWidth', 1, 'Color', colors(i,:));
        end
        hold off;
    end
end
