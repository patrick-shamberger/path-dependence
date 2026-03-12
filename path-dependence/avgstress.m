%% =============================================================================
% Description:
%{
    TITLE: AVERAGE STRESS OF TRANSFORMATION (FUNCTION)
    AUTHOR: DERIAN R. MORPHEW / PATRICK J. SHAMBERGER
    VERSION: 1 (2026)    

    SUMMARY:
        [Calculates average stress of transformation from austenite to 
         martensite.]

    INPUTS:
        - Paths structure containing all paths.
        - MatProps structure containing relevant material properties
          (including the phase boundary).
        - Plotting structure that contains grid size.
        - Boolean logic to display plots (1 - Yes, 0 - No).

    OUTPUTS:
        - Average stress of transformation value located within the
          Paths structure for each path.

    NOTES:
%}
%  =============================================================================
function [Paths, PB] = avgstress(Paths, MatProps, Plotting, Plot)
    %% Extract Required Inputs
    % === Phase Boundary ===
    PB = MatProps.Phaseboundary;
    val = [PB.Stress] > 0; % Ignore 0 MPa 
    Mf = [PB(val).Mf]';
    Ms = [PB(val).Ms]';
    Af = [PB(val).Af]';
    As = [PB(val).As]';
    PBSt = [PB(val).Stress]';

    % === Plotting Grid ===
    StG = Plotting.StressGrid;
    S_axis = StG(:,1);

    %% Compute Average Stress
    % === Interpolate Phase Boundaries ===
    Ms_T = interp1(PBSt, Ms, S_axis, 'linear');
    Mf_T = interp1(PBSt, Mf, S_axis, 'linear');
    Af_T = interp1(PBSt, Af, S_axis, 'linear');
    As_T = interp1(PBSt, As, S_axis, 'linear');

    % === Loop over paths ===
    pns = fieldnames(Paths);
    pns = pns(contains(pns, 'Path_'));

    if Plot == 1
        % === Plot Forward Transformation Points ===
        figure('Name', 'A to M Transition Points');
        hold on;
        xlabel('Temperature (K)');
        ylabel('Stress (MPa)');
        fill([Mf_T, fliplr(Ms_T)], ...
            [S_axis, fliplr(S_axis)], ...
            [211 211 211]/255, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        plot(Mf_T, S_axis, 'k-', 'LineWidth', 1);
        plot(Ms_T, S_axis, 'k-', 'LineWidth', 1);
        plot(Af_T, S_axis, 'k--', 'LineWidth', 1);
        plot(As_T, S_axis, 'k--', 'LineWidth', 1);
    end

    for i = 1:length(pns)
        pn = pns{i};

        % === FullTable ===
        FT = Paths.(pn).FullTable;
        S_path = FT.Stress;
        T_path = FT.Temperature;

        % === Initialize ===
        St_sum = 0;
        n_sum = 0;
        valid = false;

        for k = 2:length(T_path)
            T_k = T_path(k);
            T_k0 = T_path(k-1);
            S_k = S_path(k);
            S_k0 = S_path(k-1);
    
            Ms_k = interp1(PBSt, Ms, S_k, 'linear');
            Mf_k = interp1(PBSt, Mf, S_k, 'linear');
            Ms_k0 = interp1(PBSt, Ms, S_k0, 'linear');
            Mf_k0 = interp1(PBSt, Mf, S_k0, 'linear');
    
            % === Forward Transformation Region ===
            tr_k = T_k >= Mf_k && T_k <= Ms_k;
            tr_k0 = T_k0 >= Mf_k0 && T_k0 <= Ms_k0;
    
            if tr_k && tr_k0
                % === Rule of Mixtures Phase Fraction ===
                %{ 
                    Use rule of mixtures for estimating phase fraction here
                    since it is only being used to deterimine if the phase
                    transformation is progressing (and not specifically the
                    magnitude of the progression). 
                %}

                phi_k = (T_k - Ms_k) / (Mf_k - Ms_k);
                phi_k0 = (T_k0 - Ms_k0) / (Mf_k0 - Ms_k0);
    
                if phi_k >= phi_k0 % Transformation is progressing
                    St_sum = St_sum + S_k;
                    n_sum = n_sum + 1;
                    valid = true;
    
                    % === Plot (T,σ) transformation values ===
                    if Plot == 1
                        scatter(T_k, S_k, 25, 'filled', 'MarkerFaceColor', 'k');
                    end
                end
            end
        end

        % === Store ===
        if valid
            avg_stress = St_sum / n_sum;
            Paths.(pn).Avg_Stress = avg_stress;
        else
            fprintf('Path %d: No valid data\n', i);
            Paths.(pn).Avg_Stress = NaN;
        end
    end

    if Plot == 1
        hold off;
    end
end

