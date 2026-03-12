
%% =============================================================================
% Description:
%{
    TITLE: STRAIN POPULATION (FUNCTION)
    AUTHOR: DERIAN R. MORPHEW / PATRICK J. SHAMBERGER
    VERSION: 1 (2026)      

    SUMMARY:
        [Populates strain values along each path.]

    INPUTS:
        - Paths structure containing all paths.
        - MatProps structure containing relevant material properties.
        - Plotting structure containing grid points.
        - Boolean logic to display plots (1 - Yes, 0 - No).

    OUTPUTS:
        - Updated FullTable in the Paths structure for each path, populated
          with strain values.

    NOTES:
%}
%  =============================================================================
function Paths = strainpop(Paths, MatProps, Plotting, Plot)

    % === Extract Required Inputs ===
    TG = Plotting.TempGrid;
    StG = Plotting.StressGrid;
    Strain_A = Plotting.StrainGrid_A;
    Tax = TG(1,:);
    Stax = StG(:,1);

    % === Loop over each path ===
    pns = fieldnames(Paths);
    pns = pns(contains(pns, 'Path_'));

    for i = 1:length(pns)
        pn = pns{i};

        % === Path Specific Grids ===
        Strain_M = Paths.(pn).StrainShifted_M;
        Strain_DM = Paths.(pn).StrainShifted_DM;
        Strain_AM = Paths.(pn).StrainShifted_AM;
        Strain_MA = Paths.(pn).StrainShifted_MA;

        % === Interpolate Grids ===
        FA = griddedInterpolant({Tax, Stax}, Strain_A.', 'linear', 'none');
        FM = griddedInterpolant({Tax, Stax}, Strain_M.', 'linear', 'none');
        FDM = griddedInterpolant({Tax, Stax}, Strain_DM.', 'linear', 'none');
        F_MA = griddedInterpolant({Tax, Stax}, Strain_MA.', 'linear', 'none');
        F_AM = griddedInterpolant({Tax, Stax}, Strain_AM.', 'linear', 'none');

        % === Add Strain to FullTable ===
        FullTable = Paths.(pn).FullTable;
        T = FullTable.Temperature;
        St = FullTable.Stress;
        PF = FullTable.PhaseFraction;
        N = height(FullTable);
        eps = NaN(N,1);

        % === Extract First Point ===
        T1 = T(1);
        S1 = St(1);
        PF1 = PF(1);

        % === Populate First Strain Point ===
        if PF1 == 1
            eps(1) = FA(T1, S1);
        elseif PF1 > 0 && PF1 < 1
            eps(1) = F_MA(T1, S1);
        end

        % === Conditions ===
        dT = [0; diff(T)];
        dS = [0; diff(St)];
        dPF = [0; diff(PF)];

        %% CASE 1: Fully Austenite (PF = 1)
        i_c1 = PF == 1;
        i_c1(1) = false;

        if any(i_c1)
            Tc1 = T(i_c1);
            Stc1 = St(i_c1);
            epsA = getstrain(FA, Tc1, Stc1, Stax, Tax, Strain_A);
            eps(i_c1) = epsA;
        end

        %% CASE 2: Fully Forward Martensite (PF = 0)
        i_c2 = (PF == 0) & ((dS > 0) | (abs(dS) < 1e-8 & dT < 0));
        i_c2(1) = false;

        if any(i_c2)
            Tc2 = T(i_c2);
            Stc2 = St(i_c2);
            epsM = getstrain(FM, Tc2, Stc2, Stax, Tax, Strain_M);
            eps(i_c2) = epsM;
        end

        %% CASE 3: Fully Detwinned Martensite (PF = 0)
        i_c3 = (PF == 0) & ((dS < 0) | (abs(dS) < 1e-8 & dT > 0));
        i_c3(1) = false;

        if any(i_c3)
            Tc3 = T(i_c3);
            Stc3 = St(i_c3);
            epsDM = getstrain(FDM, Tc3, Stc3, Stax, Tax, Strain_DM);
            eps(i_c3) = epsDM;
        end

        %% CASE 4: Transition A to M (0 < PF < 1 AND dPF < 0)
        i_c4 = (PF > 0) & (PF < 1) & (dPF < 0);
        i_c4(1) = false;

        if any(i_c4)
            Tc4 = T(i_c4);
            Stc4 = St(i_c4);
            epsAM = getstrain(F_AM, Tc4, Stc4, Stax, Tax, Strain_AM);
            eps(i_c4) = epsAM;
        end

        %% CASE 5: Transition M to A (0 < PF < 1 AND dPF > 0)
        i_c5 = (PF > 0) & (PF < 1) & (dPF > 0);
        i_c5(1) = false;

        if any(i_c5)
            Tc5 = T(i_c5);
            Stc5 = St(i_c5);
            epsMA = getstrain(F_MA, Tc5, Stc5, Stax, Tax, Strain_MA);
            eps(i_c5) = epsMA;
        end

        %% CASE 6: Transition Plateau (0 < PF < 1 AND dPF == 0)
        i_c6 = (PF > 0) & (PF < 1) & (dPF == 0);
        i_c6(1) = false;

        if any(i_c6)
            inds = find(i_c6);
            for k = inds.'
                PFk = PF(k);
                Sk = St(k);
                Sk_p = St(k-1);
                dS_loc = Sk - Sk_p;
                Tk = T(k);
                Tk_p = T(k-1);
                dT_loc = Tk - Tk_p;

                % === Nearest Grid Index ===
                [~, iS] = min(abs(Stax - Sk));
                [~, iS_p] = min(abs(Stax - Sk_p));

                % === Stress contribution during constant phase fraction ===
                %{
                    During the constant phase fraction region, there is
                    still an evolution of strain. To account for the stress
                    contribution for the strain evolution, the changes in
                    transformation strain as a function of stress is used.
                %}

                if dS_loc < 0 % Δε_MA

                    Stc_k = Strain_MA(iS, :);
                    Stp_k = Strain_MA(iS_p, :);

                    Tcf_k = find(~isnan(Stc_k),1,'first');
                    Tcl_k = find(~isnan(Stc_k),1,'last');
                    Tpf_k = find(~isnan(Stp_k),1,'first');
                    Tpl_k = find(~isnan(Stp_k),1,'last');

                    As_epsc = Stc_k(Tcf_k);
                    Af_epsc = Stc_k(Tcl_k);
                    As_epsp = Stp_k(Tpf_k);
                    Af_epsp = Stp_k(Tpl_k);

                    % === Current/Previous Δε_MA difference ===
                    Stcon_c = Af_epsc + (1-PFk)*(As_epsc - Af_epsc);
                    Stcon_p = Af_epsp + (1-PFk)*(As_epsp - Af_epsp);

                    % === Thermal offset if T below As ===
                    if Tk <= Tax(Tcf_k)

                        Stc_k = Strain_DM(iS,:);
                        Stp_k = Strain_DM(iS_p,:);

                        Tcl_k = find(~isnan(Stc_k),1,'last');
                        Tpl_k = find(~isnan(Stp_k),1,'last');

                        DM_epsc = getstrain(FDM, Tk, Sk, Stax, Tax, Strain_DM);
                        DM_epsp = getstrain(FDM, Tk_p, Sk_p, Stax, Tax, ...
                                            Strain_DM);
                        As_epsc = Stc_k(Tcl_k);
                        As_epsp = Stp_k(Tpl_k);

                        Stcon_c = Stcon_c + (DM_epsc - As_epsc);
                        Stcon_p = Stcon_p + (DM_epsp - As_epsp);
                    end

                    Stcon = Stcon_c - Stcon_p;

                elseif dS_loc > 0 % Δε_AM
 
                    Stc_k = Strain_AM(iS,:);
                    Stp_k = Strain_AM(iS_p,:);

                    Tcf_k = find(~isnan(Stc_k),1,'first');
                    Tcl_k = find(~isnan(Stc_k),1,'last');
                    Tpf_k = find(~isnan(Stp_k),1,'first');
                    Tpl_k = find(~isnan(Stp_k),1,'last');

                    Mf_epsc = Stc_k(Tcf_k);
                    Ms_epsc = Stc_k(Tcl_k);
                    Mf_epsp = Stp_k(Tpf_k);
                    Ms_epsp = Stp_k(Tpl_k);

                    % === Current/Previous Δε_AM difference ===
                    Stcon_c = Ms_epsc + (1-PFk)*(Mf_epsc - Ms_epsc);
                    Stcon_p = Ms_epsp + (1-PFk)*(Mf_epsp - Ms_epsp);

                    % === Thermal offset if T above Ms ===
                    if Tk >= Tax(Tcf_k)
                        
                        Stc_k = Strain_A(iS,:);
                        Stp_k = Strain_A(iS_p,:);

                        Tcf_k = find(~isnan(Stc_k),1,'first');
                        Tpf_k = find(~isnan(Stp_k),1,'first');

                        A_epsc = getstrain(FA, Tk, Sk, Stax, Tax, Strain_A);
                        A_epsp = getstrain(FA, Tk_p, Sk_p, Stax, Tax, ...
                                            Strain_A);
                        Ms_epsc = Stc_k(Tcf_k);
                        Ms_epsp = Stp_k(Tpf_k);

                        Stcon_c = Stcon_c + (A_epsc - Ms_epsc);
                        Stcon_p = Stcon_p + (A_epsp - Ms_epsp);

                    end

                    Stcon = Stcon_c - Stcon_p;
                else
                    Stcon = 0;
                end

                % === T contribution during constant phase fraction ===
                %{
                    During the constant phase fraction region, there is
                    still an evolution of strain. To account for the
                    temperature contribution for strain evolution, a rule
                    of mixtures between the single phase CTE is used.
                    Ideally, the difference between forward and reverse
                    transformation strains could be used, but that will not
                    work for a simple path prediction here due to the shape
                    memory effect. 
                %}

                cA = MatProps.CTE.Austenite.Coefficients(:).';
                CTEA = polyval(cA, Sk);
                cM = Paths.(pn).DMCTE.Coefficients(:).';
                CTEM = polyval(cM, Sk);

                CTEmix = PFk*CTEA + (1-PFk)*CTEM;
                Tcon = CTEmix*dT_loc;

                % === Strain evolution during constant PF ===
                eps(k) = eps(k-1) + Stcon + Tcon;
            end
        end

        %% Plateau continuity correction
        %{
            Even after accounting for stress/T contributions during the
            phase fraction plateua, the resulting strain at the end of the
            plateua and the strain just after the phase fraction begins
            progressing again can be not smooth. While this effect is
            small, it is corrected for here.
        %}
        if any(i_c6)
            inds = find(i_c6);

            % === Identify Plateaus ===
            breaks = [1; find(diff(inds) > 1) + 1];
            plat_s = inds(breaks);
            plat_e = [inds(breaks(2:end)-1); inds(end)];

            % === Track largest mismatch ===
            max_diff = 0;

            for s = 1:length(plat_s)

                ip1 = plat_s(s); % First plateau point
                ip2 = plat_e(s); % Last plateau point
                i3 = ip2 + 1; % Point after plateau

                if i3 > N % Plateau ends at the end of the path
                    Tsame = abs(T(ip2) - T(1));
                    Stsame = abs(St(ip2) - St(1));
                    PFsame = abs(PF(ip2) - PF(1));

                    if Tsame && Stsame && PFsame
                        eps_after = eps(1);
                    else
                        continue
                    end
                else
                    eps_after = eps(i3);
                end

                eps_ps = eps(ip1);
                eps_pe = eps(ip2);

                diff_end = abs(eps_pe - eps_after);

                if diff_end > max_diff
                    max_diff = diff_end;
                end

                plat = ip1:ip2;
                eps_old = eps(plat);

                old_s = eps_ps;
                old_e = eps_pe;

                new_s = old_s; % Preserve first plateau point
                new_e = eps_after; % Match strain after plateau

                if old_e == old_s % Flat strain plateau
                    if isscalar(plat)
                        eps(plat) = new_e;
                    else
                        eps(plat) = linspace(new_s, new_e, numel(plat)).';
                    end
                else % Preserve relative shape
                    xi = (eps_old - old_s) ./ (old_e - old_s);
                    eps(plat) = new_s + xi .* (new_e - new_s);
                end

                % === Endpoints ===
                eps(ip1) = old_s;
                eps(ip2) = new_e;
            end

            % === Print mismatch ===
            if max_diff > 1e-10
                pnum_str = regexp(pn, '\d+', 'match', 'once');
                pnum = str2double(pnum_str);

                 % === Percent Mismatch ===
                deps = max(eps) - min(eps);
                pd = 100 * max_diff / deps;
                fprintf(['Phase fraction plateau in Path %d results in' ...
                    ' strain mismatch. (Δε mismatch) / (Δε_max of path) =' ...
                         ' %.6f %%\n'], pnum, pd);
            end
        end

        %% NaN Safeguard
        nan_i = isnan(eps);
        if any(nan_i)

            pnum_str = regexp(pn, '\d+', 'match', 'once');
            pnum = str2double(pnum_str);

            n_nans = sum(nan_i);
            % === Print for user awareness ===
            fprintf(['Path %d: %d NaNs found in Strain.' ...
                ' Forward-filling with previous values...\n'], pnum, n_nans);

            % === Fill With Most Recent Non-NaN Value ===
            eps = fillmissing(eps, 'previous');
            if any(isnan(eps))
                eps = fillmissing(eps, 'next');
            end
        end

        % === Save strain column ===
        FullTable.Strain = eps;
        Paths.(pn).FullTable = FullTable;
    end

    %% Optional Plotting
    if Plot == 1

        colormap = [
            161 218 180
            65 182 196
            44 127 184
            37  52 148
            ] / 255;

        nP = length(pns);
        colors = interp1(linspace(0,1,4), colormap, linspace(0,1,nP));

        % === Strain vs Temperature ===
        figure; 
        hold on;
        xlabel('Temperature (K)');
        ylabel('Strain (%)');
        for i = 1:nP
            pn = pns{i};
            FT = Paths.(pn).FullTable;
            T = FT.Temperature;
            eps = FT.Strain;
            plot(T, eps, '-', 'LineWidth', 1, 'Color', colors(i,:));
        end
        hold off;

        % === Strain vs Stress ===
        figure; 
        hold on;
        xlabel('Stress (MPa)');
        ylabel('Strain (%)');
        for i = 1:nP
            pn = pns{i};
            FT = Paths.(pn).FullTable;
            St = FT.Stress;
            eps = FT.Strain;
            plot(St, eps, '-', 'LineWidth', 1.5, 'Color', colors(i,:));
        end
        hold off;
    end
end

%% Helper Function
function eps_o = getstrain(F, Tc, Stc, Stax, Tax, StrainGrid)
    eps_o = F(Tc, Stc);

    % === Fallback ===
    nan_i = find(isnan(eps_o));
    for ii = 1:length(nan_i)
        k = nan_i(ii);
        [~, St_i] = min(abs(Stax - Stc(k)));
        val = find(~isnan(StrainGrid(St_i, :)));
        [~, T_i] = min(abs(Tax(val) - Tc(k)));
        eps_o(k) = StrainGrid(St_i, val(T_i));
    end
end