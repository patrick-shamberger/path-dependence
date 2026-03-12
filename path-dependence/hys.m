
%  =============================================================================
% Description:
%{
    TITLE: FITTING HYSTERESIS ENVELOPES (FUNCTION)
    AUTHOR: DERIAN R. MORPHEW / PATRICK J. SHAMBERGER
    VERSION: 1 (2026)

    SUMMARY:
        [This function takes the raw ISHC experiments and fits a
         two-asymptote logistic function (sigmoid) to model the coexistence
         regions.]

    INPUTS:
        - InitialInputs structure for ISHC data.
        - MatProps structure for phase boundaries and reference points.
        - Plotting structure with sigmoid function properties and grid
          points
        - Plot (boolean logic) to declare if figures associated with this
          fuction should be shown (1 - YES, 0 - NO). 

    OUTPUTS:
        - InitialInputs structure.
        - MatProps structure.
        - Plotting structure to hold the gridpoints of each single phase as
          well as the sigmoid function properties for plotting.
    NOTES:
%}
%  =============================================================================

function [InitialInputs, MatProps, Plotting] = hys(InitialInputs, MatProps, Plotting, Plot)
    
    %% Extract Required Inputs
    % === Experimental Data ===
    ISHCdata = InitialInputs.ISHCdata;
    Stresses = InitialInputs.Stresses;

    % === Sigmoid parameters ===
    FT0_scale = Plotting.Sigmoid.MsMf.T0_scale;
    Fk = Plotting.Sigmoid.MsMf.k;
    RT0_scale = Plotting.Sigmoid.AsAf.T0_scale;
    Rk = Plotting.Sigmoid.AsAf.k;

    % === CTE Coefficients ===
    CTEA_c = MatProps.CTE.Austenite.Coefficients;
    CTEM_c = MatProps.CTE.Martensite.Coefficients;

    % === Compliance ===
    ComA_c = MatProps.Compliance.Austenite.Compliance;
    ComM_c = MatProps.Compliance.Martensite.Compliance;

    % === Reference Points ===
    RA = MatProps.Reference.Austenite;
    RM = MatProps.Reference.Martensite;

    % === Phase Boundaries ===
    PB = MatProps.Phaseboundary;
    val_i = [PB.Stress] > 0;
    PBSt = [PB(val_i).Stress]';
    PBAs = [PB(val_i).As]';
    PBAf = [PB(val_i).Af]';
    PBMs = [PB(val_i).Ms]';
    PBMf = [PB(val_i).Mf]';

    %% Generate Mesh Grid
    N = Plotting.Grid.N;
    St_s = Plotting.Grid.StressRange(1);
    St_f = Plotting.Grid.StressRange(2);
    T_s = Plotting.Grid.TempRange(1);
    T_f = Plotting.Grid.TempRange(2);
    StG = linspace(St_s, St_f, N);
    TG = linspace(T_s, T_f, N);
    [TMesh, StMesh] = meshgrid(TG, StG);

    %% Austenite Grid
    RAT = RA.Temperature;
    RASt = RA.Stress;
    RAeps = RA.Strain;

    dStA = StMesh - RASt;
    epsA = RAeps + ComA_c * dStA; % ε(σ) values at Reference Temp
    CTEA_g = polyval(CTEA_c, StMesh); % α at each σ
    StrainGrid_A = epsA + CTEA_g .* (TMesh - RAT); % ε(T,σ)

    %% Martensite Grid
    RMT = RM.Temperature;
    RMSt = RM.Stress;
    RMeps = RM.Strain;

    a = ComM_c(1); 
    b = ComM_c(2); 
    c = ComM_c(3);
    dStM = StMesh - RMSt;

    epsM = RMeps ...
        + (a/3) * (StMesh.^3 - RMSt^3) ...
        + (b/2) * (StMesh.^2 - RMSt^2) ...
        + c * dStM; % ε(σ) values at Reference Temp

    CTEM_g = polyval(CTEM_c, StMesh); % α at each σ
    StrainGrid_M = epsM + CTEM_g .* (TMesh - RMT); % ε(T,σ)

    % === Store Meshes ===
    Plotting.TempGrid = TMesh;
    Plotting.StressGrid = StMesh;

    %% Clip Single Phase Grids Along Phase Boundaries

    % === Loop over temperature rows (fixed T - varying stress) ===
    for row = 1:size(TMesh, 1)
        St_row = StMesh(row, :);
        T_row = TMesh(row, :);

        % === Interpolate As, Ms, Af ===
        As_row = interp1(PBSt, PBAs, St_row, 'linear', NaN);
        Ms_row = interp1(PBSt, PBMs, St_row, 'linear', NaN);
        Af_row = interp1(PBSt, PBAf, St_row, 'linear', NaN);

        % === For austenite: min(Ms, Af) at each stress ===
        lower_A = min(Ms_row, Af_row);

        % === Clip ===
        StrainGrid_M(row, T_row > As_row) = NaN;
        StrainGrid_A(row, T_row < lower_A) = NaN;
    end

    % === Store Strain Grids ===
    Plotting.StrainGrid_A = StrainGrid_A;
    Plotting.StrainGrid_M = StrainGrid_M;

    %% Sigmoid Visualization: All Datasets
    % === Initialize structures ===
    MsMfSigmoids = struct([]);
    AsAfSigmoids = struct([]);

    % === Nearest grid row for each experimental stress ===
    [~, Stexp_i] = arrayfun(@(s) min(abs(StG - s)), Stresses);

    for i = 1:length(StG)
        stress = StG(i);
        [~, row_i] = min(abs(StG - stress));
        eps_A = StrainGrid_A(row_i, :);
        eps_M = StrainGrid_M(row_i, :);

        % === Interpolate phase boundary temperatures ===
        Mf_T = interp1(PBSt, PBMf, stress, 'linear', NaN);
        Ms_T = interp1(PBSt, PBMs, stress, 'linear', NaN);
        As_T = interp1(PBSt, PBAs, stress, 'linear', NaN);
        Af_T = interp1(PBSt, PBAf, stress, 'linear', NaN);

        % === T0 values for transition regions ===
        FT0 = FT0_scale * Ms_T; 
        RT0 = RT0_scale * Af_T;

        %% Find ε at Mf-Ms-As-Af
        % === ε at Mf ===
        [~, Mf_i] = min(abs(TG - Mf_T));
        eps_Mf = StrainGrid_M(row_i, Mf_i);
        if isnan(eps_Mf) % Fallback
            Mf_ii = find(~isnan(StrainGrid_M(row_i, :)));
            [~, Mf_iii] = min(abs(TG(Mf_ii) - Mf_T));
            eps_Mf = StrainGrid_M(row_i, Mf_ii(Mf_iii));
        end

        % === ε at Ms ===
        [~, Ms_i] = min(abs(TG - Ms_T));
        eps_Ms = StrainGrid_A(row_i, Ms_i);
        if isnan(eps_Ms) % Fallback
            Ms_ii = find(~isnan(StrainGrid_A(row_i, :)));
            [~, Ms_iii] = min(abs(TG(Ms_ii) - Ms_T));
            eps_Ms = StrainGrid_A(row_i, Ms_ii(Ms_iii));
        end

        % === ε at As ===
        [~, As_i] = min(abs(TG - As_T));
        eps_As = StrainGrid_M(row_i, As_i);
        if isnan(eps_As) % Fallback
            As_ii = find(~isnan(StrainGrid_M(row_i, :)));
            [~, As_iii] = min(abs(TG(As_ii) - As_T));
            eps_As = StrainGrid_M(row_i, As_ii(As_iii));
        end

        % === ε at Af ===
        [~, Af_i] = min(abs(TG - Af_T));
        eps_Af = StrainGrid_A(row_i, Af_i);
        if isnan(eps_Af) % Fallback
            Af_ii = find(~isnan(StrainGrid_A(row_i, :)));
            [~, Af_iii] = min(abs(TG(Af_ii) - Af_T));
            eps_Af = StrainGrid_A(row_i, Af_ii(Af_iii));
        end

        %% Plot at experimental stress levels
        if Plot == 1
            match = find(Stexp_i == row_i);
            for m = 1:numel(match)
                Stexp = match(m);

                ISHC = ISHCdata{Stexp};
                T = ISHC(:,2) + 273; % [K]
                eps = ISHC(:,3); % [%]

                PB_i = PB(Stexp + 1); % Ignoring 0 MPa entry

                figure('Name', sprintf('Sigmoid Vis %.1f MPa', ...
                    Stresses(Stexp)));
                hold on;
                plot(T, eps, 'k', 'LineWidth', 1);
                plot(TG, eps_A, 'r--', 'LineWidth', 2);
                plot(TG, eps_M, 'b--', 'LineWidth', 2);

                % === Phase boundary points ===
                for T_bnd = [PB_i.Mf, PB_i.As]
                    [~, idx] = min(abs(TG - T_bnd));
                    if isnan(eps_M(idx)) % Fallback
                        val = find(~isnan(eps_M));
                        [~, ii] = min(abs(TG(val) - T_bnd));
                        idx = val(ii);
                    end
                    scatter(TG(idx), eps_M(idx), 50, 'b', 'filled');
                end
                for T_bnd = [PB_i.Ms, PB_i.Af]
                    [~, idx] = min(abs(TG - T_bnd));
                    if isnan(eps_A(idx)) % Fallback
                        val = find(~isnan(eps_A));
                        [~, ii] = min(abs(TG(val) - T_bnd));
                        idx = val(ii);
                    end
                    scatter(TG(idx), eps_A(idx), 50, 'r', 'filled');
                end

                % === Austentie to Martensite Sigmoids ===
                Ftr = TG >= Mf_T & TG <= Ms_T;
                AtoM = NaN(1, N);
                AtoM(Ftr) = ...
                  (eps_Mf - eps_Ms) ./ (1 + exp(Fk * (TG(Ftr) - FT0))) + eps_Ms;
                plot(TG, AtoM, 'b-', 'LineWidth', 2);

                % === Martensite to Austente Sigmoids ===
                Rtr = TG >= As_T & TG <= Af_T;
                MtoA = NaN(1, N);
                MtoA(Rtr) = ...
                  (eps_As - eps_Af) ./ (1 + exp(Rk * (TG(Rtr) - RT0))) + eps_Af;
                plot(TG, MtoA, 'r-', 'LineWidth', 2);

                xlabel('Temperature (K)');
                ylabel('Strain (%)');
                hold off;
            end
        end
        % === Store parameters ===
        MsMfSigmoids(i).stress = stress;
        MsMfSigmoids(i).Mf_T = Mf_T;
        MsMfSigmoids(i).Ms_T = Ms_T;
        MsMfSigmoids(i).eps_Mf = eps_Mf;
        MsMfSigmoids(i).eps_Ms = eps_Ms;
        MsMfSigmoids(i).T0 = FT0;
        MsMfSigmoids(i).k = Fk;

        AsAfSigmoids(i).stress = stress;
        AsAfSigmoids(i).As_T = As_T;
        AsAfSigmoids(i).Af_T = Af_T;
        AsAfSigmoids(i).eps_As = eps_As;
        AsAfSigmoids(i).eps_Af = eps_Af;
        AsAfSigmoids(i).T0 = RT0;
        AsAfSigmoids(i).k = Rk;
    end

    Plotting.MsMfSigmoids = MsMfSigmoids;
    Plotting.AsAfSigmoids = AsAfSigmoids;
    Plotting.TempGrid = TMesh;
    Plotting.StressGrid = StMesh;
    Plotting.TempAxis = TG;
    Plotting.StressAxis = StG;

    %% All Phase Surfaces
    if Plot == 1
        figure('Name', 'All Strain Grids');
        hold on;

        % Reconstruct using saved parameters
        TempGrid = Plotting.TempGrid;
        StressGrid = Plotting.StressGrid;
        Strain_MsMf = NaN(size(TempGrid));
        Strain_AsAf = NaN(size(TempGrid));
        MsMfSigmoids = Plotting.MsMfSigmoids;
        AsAfSigmoids = Plotting.AsAfSigmoids;
        for i = 1:length(MsMfSigmoids)
            Fpar = Plotting.MsMfSigmoids(i);
            row = i;  % corresponds to stress_grid(i)

            Ftr = TempGrid(row,:) >= Fpar.Mf_T & TempGrid(row,:) <= Fpar.Ms_T;
            T = TempGrid(row, Ftr);

            Strain_MsMf(row, Ftr) = (Fpar.eps_Mf - Fpar.eps_Ms) ./ ...
                (1 + exp(Fpar.k * (T - Fpar.T0))) + Fpar.eps_Ms;
        end
        for i = 1:length(AsAfSigmoids)
            Rpar = Plotting.AsAfSigmoids(i);
            row = i;

            Rtr = TempGrid(row,:) >= Rpar.As_T & TempGrid(row,:) <= Rpar.Af_T;
            T = TempGrid(row, Rtr);

            Strain_AsAf(row, Rtr) = (Rpar.eps_As - Rpar.eps_Af) ./ ...
                (1 + exp(Rpar.k * (T - Rpar.T0))) + Rpar.eps_Af;
        end

        % Surfaces
        surf(TempGrid, StressGrid, StrainGrid_A, ...
            'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        surf(TempGrid, StressGrid, StrainGrid_M, ...
            'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        surf(TempGrid, StressGrid, Strain_MsMf, ...
            'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        surf(TempGrid, StressGrid, Strain_AsAf, ...
            'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

        xlabel('Temperature (K)');
        ylabel('Stress (MPa)');
        zlabel('Strain (%)');
        view(35, 25);
        hold off;
    end
end