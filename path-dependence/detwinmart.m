%% =============================================================================
% Description:
%{
    TITLE: DETWINNED MARTENSITE GRID (FUNCTION)
    AUTHOR: DERIAN R. MORPHEW / PATRICK J. SHAMBERGER
    VERSION: 1 (2026)

    SUMMARY:
        [Creates detwinned martensite strain grid. This grid is
         necessary when unloading or heating in the martensite phase. If
         experimental data under these conditions is unknown, the ε(σ)
         behavior is approximated to be equivalent to the Austenite phase while 
         the α behavior is approximated to be equivalent to the Martensite
         phase at the stress level which the driving force is reversed].

    INPUTS:
        - Paths structure containing all paths.
        - MatProps structure containing relevant material properties
          (including the phase boundary).
        - Plotting structure that contains the strain grid of the
          martensite phase.
        - Boolean logic to display plots (1 - Yes, 0 - No).

    OUTPUTS:
        - Detwinned martensite strain grid for each path. 

    NOTES:
%}
%  =============================================================================
function [Paths, MatProps] = detwinmart(Paths, MatProps, Plotting, Plot)
    %% Extract Required Inputs
    % === Phase boundary ===
    PB = MatProps.Phaseboundary;
    val1 = [PB.Stress] > 0; % Ignore 0 MPa point
    PBSt = [PB(val1).Stress]';
    PBAs = [PB(val1).As]';

    PBtwin = MatProps.TwinPhaseboundary; 
    PBtAs_5 = PBtwin(2).As; % 5 MPa point
    
    % === Plotting grid ===
    StrainGrid_M = Plotting.StrainGrid_M;
    TG = Plotting.TempGrid;
    StG = Plotting.StressGrid;
    Tref_M = MatProps.Reference.Martensite.Temperature;

    % === Axes ===
    S_axis = StG(:,1);
    T_axis = TG(1,:);
    StGmin = min(S_axis);

    %% Detwinned Martensite Properties from Experiments
    % === Import reorientation experiment ===
    Sys = Plotting.Materialsys;
    folder_cur = fileparts(mfilename('fullpath'));
    folder_data = fullfile(folder_cur, 'Experimental Data', 'Detwinning', Sys);
    file = [Sys '_Reversible_Stress.xlsx'];
    fpath = fullfile(folder_data, file);

    if isfile(fpath)

        unload = readmatrix(fpath, 'Sheet', 4); % Unloading leg of data
        heat = readmatrix(fpath, 'Sheet', 5); % Heating leg of data

        eps_u = unload(:,3); % [%]
        St_u = unload(:,4); % [MPa]

        T_h = heat(:,2) + 273; % [K]
        eps_h = heat(:,3); % [%]

        % === DM compliance ===
        DMstrain_p = polyfit(St_u, eps_u, 1);
        DMcomp_p = DMstrain_p(1);
        MatProps.Compliance.DM.Type = 'Constant';
        MatProps.Compliance.DM.Compliance = DMcomp_p;

        % === DM CTE ===
        % Assume CTE of DM is independent of stress
        DMctel_i = (T_h < PBtAs_5);
        DMctel_p = polyfit(T_h(DMctel_i), eps_h(DMctel_i), 1); % Linear ε(T) fit
        DMcte_low = DMctel_p(1);

        % === As temperature for each stress row ===
        As_T = interp1(PBSt, PBAs, S_axis, 'linear', NaN);

        % === Loop over each path ===
        pns = fieldnames(Paths);
        pns = pns(contains(pns, 'Path_'));

        for i = 1:length(pns)
            pn = pns{i};
            
            % === FullTable ===
            FT = Paths.(pn).FullTable;
            S_path = FT.Stress;

            % === Max stress from FullTable ===
            Stmax = max(S_path);
            Paths.(pn).MaxStress = Stmax;

            % === Strain vs T at max stress ===
            [~, Stmax_i] = min(abs(S_axis - Stmax));
            Strain_T = StrainGrid_M(Stmax_i, :);

            % === Fit linear model to strain vs T at max stress ===
            val3 = ~isnan(Strain_T);
            T_val = T_axis(val3);
            eps_val = Strain_T(val3);
            DMcteh_p = polyfit(T_val, eps_val, 1);
            DMcte_high = DMcteh_p(1);

            % === Fit α(σ) ===
            DMcte_vec = [DMcte_low, DMcte_high];
            DMcte_St = [StGmin, Stmax];
            DMcte = polyfit(DMcte_St, DMcte_vec, 1);

            % === Store CTE of detwinned martensite ===
            % Storing for each path to be consistent will the "else" below
            DMCTE.Type = 'Linear';
            DMCTE.Coefficients = DMcte;
            Paths.(pn).DMCTE = DMCTE;

            % === Detwinned martensite grid ===
            StrainGrid_DM = NaN(size(StrainGrid_M));

            % === Reference strain for detwinned martensite ===
            [~, Tref_i] = min(abs(T_axis - Tref_M));
            epsDM_r = StrainGrid_M(Stmax_i, Tref_i);

            for row = 1:length(S_axis)
                St_val = S_axis(row);
                As_row = As_T(row); % As temp at this stress
                mask = (T_axis <= As_row);

                % === Detwinned martensite α(σ) ===
                a_row = polyval(DMcte, St_val);

                % === ε(σ) at reference temperature ===
                dS = St_val - Stmax;
                epsr_row = epsDM_r + DMcomp_p*dS;

                % === Full ε(T,σ) grid for Detwinned martensite ===
                StrainGrid_DM(row, mask) = epsr_row + ...
                    a_row*(T_axis(mask) - T_axis(Tref_i));
            end

            % === Store detwinned martensite grid ===
            Paths.(pn).StrainGrid_DM = StrainGrid_DM;

            % === Optional plotting ===
            if Plot == 1
                figure('Name', sprintf('DetwinnedM & Martensite – %s', pn));
                hold on;
                surf(TG, StG, StrainGrid_M, 'FaceColor', 'b', ...
                    'EdgeColor', 'none', 'FaceAlpha', 0.5);
                surf(TG, StG, StrainGrid_DM, 'FaceColor', 'r', ...
                    'EdgeColor', 'none', 'FaceAlpha', 0.5);
                xlabel('Temperature (K)');
                ylabel('Stress (MPa)');
                zlabel('Strain (%)');
                view(30, 30);
                box on;
            end
        end
    else
        %% Alternative Compliance & CTE for Detwinned Martensite

        fprintf(['Using austenite compliance to approximate detwinned ' ...
                 'martensite compliance for %s because the detwinning ' ...
                 'experiment file was not found\n'], Sys);
        
        % === Austenite compliance ===
        C_A = MatProps.Compliance.Austenite.Compliance; % Scalar

        % === Approximate DM compliance to be equal to austenite compliance ===
        MatProps.Compliance.DM.Type = 'Constant';
        MatProps.Compliance.DM.Compliance = C_A;

        % === As temperature for each stress row ===
        As_T = interp1(PBSt, PBAs, S_axis, 'linear', NaN);

        % === Loop over each path ===
        pns = fieldnames(Paths);
        pns = pns(contains(pns, 'Path_'));

        for i = 1:length(pns)
            pn = pns{i};

            % === FullTable  ===
            FT = Paths.(pn).FullTable;
            S_path = FT.Stress;

            % === Max stress from FullTable ===
            Stmax = max(S_path);
            Paths.(pn).MaxStress = Stmax;

            % === Strain vs T at max stress ===
            [~, Stmax_i] = min(abs(S_axis - Stmax));
            Strain_T = StrainGrid_M(Stmax_i, :);

            % === Fit linear model to strain vs T at max stress ===
            val4 = ~isnan(Strain_T);
            T_val = T_axis(val4);
            eps_val = Strain_T(val4);
            DMcteh_p = polyfit(T_val, eps_val, 1);

            % === Detwinned martensite grid ===
            StrainGrid_DM = NaN(size(StrainGrid_M));

            % === Reference ε(T) at max stress ===
            eps_T = polyval(DMcteh_p, T_axis); % 1 x nT

            for row = 1:length(S_axis)
                St_val = S_axis(row);
                As_row = As_T(row); % As temp at this stress
                mask = T_axis <= As_row; % Only fill below As

                % === Detwinned martensite strain ===
                dS = St_val - Stmax;
                StrainGrid_DM(row, mask) = ...
                    eps_T(mask) + C_A * dS;
            end

            % === Store detwinned martensite grid ===
            Paths.(pn).StrainGrid_DM = StrainGrid_DM;

            % === Store CTE of detwinned martensite ===
            DMCTE.Type = 'Constant';
            DMCTE.Coefficients = DMcteh_p(1);
            Paths.(pn).DMCTE = DMCTE;

            % === Optional plotting ===
            if Plot == 1
                figure('Name', sprintf('DetwinnedM & Martensite – %s', pn));
                hold on;
                surf(TG, StG, StrainGrid_M, 'FaceColor', 'b', ...
                    'EdgeColor', 'none', 'FaceAlpha', 0.5);
                surf(TG, StG, StrainGrid_DM, 'FaceColor', 'r', ...
                    'EdgeColor', 'none', 'FaceAlpha', 0.5);
                xlabel('Temperature (K)');
                ylabel('Stress (MPa)');
                zlabel('Strain (%)');
                view(30, 30);
                box on;
            end
        end
    end
end

