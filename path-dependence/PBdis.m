%% =============================================================================
% Description:
%{
    TITLE: PLOTTING PHASE BOUNDARY (FUNCTION)
    AUTHOR: DERIAN R. MORPHEW / PATRICK J. SHAMBERGER
    VERSION: 1 (2026)

    SUMMARY:
        [This function plots the phase boundary points for visualization
         of a proper reference point choice for each single phase.]

    INPUTS:
        - MatProps structure with stored phase boundary points.
        - Plot (boolean logic) to declare if figures associated with this
          fuction should be shown (1 - YES, 0 - NO).

    OUTPUTS:
    NOTES:
%}
%  =============================================================================
function PBdis(MatProps, Plot)
    %% Extract Phase Boundary Points
    PB = MatProps.Phaseboundary;
    n = length(PB);
    
    % === Initialize arrays ===
    St = nan(1, n);
    Mf = nan(1, n);
    Ms = nan(1, n);
    As = nan(1, n);
    Af = nan(1, n);
    
    % === Extract points ===
    for i = 1:n
        St(i) = PB(i).Stress;
        Mf(i) = PB(i).Mf;
        Ms(i) = PB(i).Ms;
        As(i) = PB(i).As;
        Af(i) = PB(i).Af;
    end
    %% Plot Phase Boundaries
    if Plot == 1
        figure('Name', 'Phase Boundaries with Paths');
        hold on;
        plot(Mf, St, 'k-', 'LineWidth', 1, 'DisplayName', 'Mf');
        plot(Ms, St, 'k-', 'LineWidth', 1, 'DisplayName', 'Ms'); 
        plot(As, St, 'k--', 'LineWidth', 1, 'DisplayName', 'As');
        plot(Af, St, 'k--', 'LineWidth', 1, 'DisplayName', 'Af');
        xlabel('Temperature (K)');
        ylabel('Stress (MPa)');
        hold off;
    end
end



