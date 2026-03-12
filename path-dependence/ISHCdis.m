%% =============================================================================
% Description:
%{
    TITLE: DISPLAY EXPERIMENTAL ISHC CURVES (FUNCTION)
    AUTHOR: DERIAN R. MORPHEW / PATRICK J. SHAMBERGER
    VERSION: 1 (2026)

    SUMMARY:
        [This function plots the ε(T) experimental ISHC curves for
         inspection to manually determine the phase boundary points at each
         experimental stress level. Additionally, these curves are used to
         manually determine temperature bounds of consideration for α(σ) in each 
         single phase.]

    INPUTS:
        - InitialInputs structure for ISHC data.
        - Plot (boolean logic) to declare if figures associated with this
          fuction should be shown (1 - YES, 0 - NO).

    OUTPUTS:
    
    NOTES:
%}
%  =============================================================================
function ISHCdis(InitialInputs, Plot)
    %% Extract Inputs
    ISHC = InitialInputs.ISHCdata;
    St = InitialInputs.Stresses;

    % === Combined figure ===
    if Plot == 1
        figAll = figure('Name', 'All ISHC Curves');
        axAll = axes('Parent', figAll);
        hold(axAll, 'on');
        xlabel(axAll, 'Temperature (K)');
        ylabel(axAll, 'Strain (%)');
    end
    %% Loop over each dataset
    for i = 1:length(ISHC)

        Strain = ISHC{i}(:,3); % [%]
        T = ISHC{i}(:,2) + 273; % [K]

        % === Individual figures for each stress ===
        if Plot == 1
            figure('Name', sprintf('%g MPa', St(i)));
            hold on;
            plot(T, Strain, 'LineWidth', 1);
            xlabel('Temperature (K)');
            ylabel('Strain (%)');
            hold off;

            % === Add curve to combined figure ===
            plot(axAll, T, Strain, ...
                 'LineWidth', 1, ...
                 'DisplayName', sprintf('%g MPa', St(i)));
        end
    end

    % === Add legend ===
    if Plot == 1
        legend(axAll, 'show', 'Location', 'best');
        hold(axAll, 'off');
    end
end
