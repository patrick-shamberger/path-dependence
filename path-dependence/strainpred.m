%% =============================================================================
% Description:
%{
    TITLE: STRAIN GRID PREDICTION (FUNCTION)
    AUTHOR: DERIAN R. MORPHEW / PATRICK J. SHAMBERGER
    VERSION: 1 (2026)

    SUMMARY:
        [Reconstructs and applies path-dependence corrections to the strain 
         grids].

    INPUTS:
        - Paths structure containing all paths.
        - Plotting structure that contains the strain grid of the
          martensite phase.
        - Boolean logic to apply path-dependence or not (1 - Yes, 0 - No).
        - Boolean logic to display plots (1 - Yes, 0 - No).

    OUTPUTS:
        - All strain grids for each path. 

    NOTES:
%}
%  =============================================================================
function Paths = strainpred(Paths, Plotting, PD, Plot)

    % === Extract Required Inputs ===
    StrainGrid_A = Plotting.StrainGrid_A;
    StrainGrid_M = Plotting.StrainGrid_M;
    MsMfSig = Plotting.MsMfSigmoids;
    AsAfSig = Plotting.AsAfSigmoids;

    Tax = Plotting.TempAxis;
    Stax = Plotting.StressAxis;
    TG = Plotting.TempGrid;

    % === Path-Specific Corrections ===
    pns = fieldnames(Paths);
    pns = pns(contains(pns, 'Path_'));

    for i = 1:length(pns)
        pn = pns{i};

        % === Extract Path-Specific Inputs ===
        FT = Paths.(pn).FullTable;
        St = FT.Stress;
        T = FT.Temperature;
        ef = Paths.(pn).FinalStrain;

        % === Path-Dependence Flag ===
        noPD = (PD == 0) || isnan(ef);

        % === Find First Decrease in Stress ===
        dS = diff(St);
        s_i = find(dS < 0, 1, 'first');
        Ts_i = T(s_i);
        Sts_i = St(s_i);

        % === Closest Grid Point ===
        [~, StG_i] = min(abs(Stax - Sts_i));
        [~, TG_i] = min(abs(Tax - Ts_i));

        %% Path-Dependent Martensite Grid Shift
        eM = StrainGrid_M(StG_i, TG_i);

        if isnan(eM) % Occurs if path doesn't unload in pure Martensite
            MTs = StrainGrid_M(StG_i, :);
            TG_i = find(~isnan(MTs), 1, 'last');
            eM = MTs(TG_i);
        end

        de_ref = ef - eM;

        % === Ignore shifts if no Path Dependence ===
        if noPD
            de_ref = 0;
        end

        St_pmax = max(St);
        St_pmin = min(St);

        Mgrid_shift = zeros(size(StrainGrid_M));

        for k = 1:length(Stax)
            scale = (Stax(k) - St_pmin) / (St_pmax - St_pmin);
            scale = max(0, min(1, scale));
            de_M = scale * de_ref;
            Mgrid_shift(k, :) = StrainGrid_M(k, :) + de_M;
        end
        Paths.(pn).StrainShifted_M = Mgrid_shift;


        %% Path-Dependent Detwinned Martensite Grid Shift
        StrainGrid_DM = Paths.(pn).StrainGrid_DM;
        eDM = StrainGrid_DM(StG_i, TG_i);

        de_DM = ef - eDM;

        % === Ignore shifts if no Path Dependence ===
        if noPD
            de_DM = 0;
        end
        DMgrid_shift = StrainGrid_DM + de_DM;
        Paths.(pn).StrainShifted_DM = DMgrid_shift;

        %% Path-Dependent Austenite to Martensite Grid Shift
        for k = 1:length(MsMfSig)
            Mf_T = MsMfSig(k).Mf_T;
            [~, T_k] = min(abs(Tax - Mf_T));
            % === Update Shifted Mf Value ===
            eps_Mf = Mgrid_shift(k, T_k);
            MsMfSig(k).eps_Mf = eps_Mf;
        end

        %% Path-Dependent Martensite to Austenite Grid Shift
        for k = 1:length(AsAfSig)
            DMTs = DMgrid_shift(k, :);
            T_k = find(~isnan(DMTs), 1, 'last');
            AsAfSig(k).eps_As = DMTs(T_k);
        end

        %% Construct Strain Surfaces
        Strain_MsMf = NaN(size(TG));
        Strain_AsAf = NaN(size(TG));

        for k = 1:length(MsMfSig)
            Mf_T = MsMfSig(k).Mf_T;
            Ms_T = MsMfSig(k).Ms_T;
            eps_Mf = MsMfSig(k).eps_Mf;
            eps_Ms = MsMfSig(k).eps_Ms;
            ks = MsMfSig(k).k;
            T0 = MsMfSig(k).T0;
            range = TG(k,:) >= Mf_T & TG(k,:) <= Ms_T;
            Tval = TG(k, range);
            Strain_MsMf(k, range) = ...
                (eps_Mf - eps_Ms) ./ (1 + exp(ks * (Tval - T0))) + eps_Ms;
        end

        for k = 1:length(AsAfSig)
            Af_T = AsAfSig(k).Af_T;
            As_T = AsAfSig(k).As_T;
            eps_Af = AsAfSig(k).eps_Af;
            eps_As = AsAfSig(k).eps_As;
            ks = AsAfSig(k).k;
            T0 = AsAfSig(k).T0;
            range = TG(k,:) >= As_T & TG(k,:) <= Af_T;
            Tval = TG(k, range);
            Strain_AsAf(k, range) = ...
                (eps_As - eps_Af) ./ (1 + exp(ks * (Tval - T0))) + eps_Af;
        end

        % === Plot All Surfaces ===
        if Plot == 1
            figure('Name', pn); 
            hold on;
            surf(Tax, Stax, Mgrid_shift, 'FaceColor', 'k', 'EdgeColor', ... 
                'none', 'FaceAlpha', 0.2);       
            surf(Tax, Stax, StrainGrid_A, 'FaceColor', 'k', 'EdgeColor', ... 
                'none', 'FaceAlpha', 0.2);    
            surf(Tax, Stax, Strain_MsMf, 'FaceColor', 'b', 'EdgeColor', ... 
                'none', 'FaceAlpha', 0.3);     
            surf(Tax, Stax, DMgrid_shift, 'FaceColor', 'k', 'EdgeColor', ... 
                'none', 'FaceAlpha', 0.2);
            surf(Tax, Stax, Strain_AsAf, 'FaceColor', 'r', 'EdgeColor', ... 
                'none', 'FaceAlpha', 0.3);         
            xlabel('Temperature (K)');
            ylabel('Stress (MPa)');
            zlabel('Strain');
            view(3);
            box on;
            hold off;
        end

        Paths.(pn).StrainShifted_AM = Strain_MsMf;
        Paths.(pn).StrainShifted_MA = Strain_AsAf;
    end
end