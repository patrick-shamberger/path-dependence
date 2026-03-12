%%  ============================================================================
% Description:
%{
    TITLE: DETWINNED MARTENSITE TO AUSTENITE PHASE BOUNDARY (FUNCTION)
    AUTHOR: DERIAN R. MORPHEW / PATRICK J. SHAMBERGER
    VERSION: 1 (2026)

    SUMMARY:
        [This function replaces the ISHC phase boundary points with the
         detwinned martensite to austenite temperatures where available.
         Particularly, if the appropriate experiment is available, the
         austenite start and finish temperatures are replaced at 5 MPa.]

    INPUTS:
        - InitialInputs structure for material system.
        - MatProps structure for phase boundary.
        - Plot (boolean logic) to declare if figures associated with this
          fuction should be shown (1 - YES, 0 - NO).

    OUTPUTS:
        - MatProps structure after updating phase boundary and storing both
          versions of phase boundaries.
    NOTES:
%}
%  =============================================================================
function MatProps = detPB(InitialInputs, MatProps, Plot)

    % === Extract Inputs ===
    Sys = InitialInputs.Materialsys;
    TwinPB = MatProps.TwinPhaseboundary;
    
    % === Import Reorientation Experiment ===
    folder_cur = fileparts(mfilename('fullpath'));
    folder_data = fullfile(folder_cur, 'Experimental Data', 'Detwinning', Sys);
    file = [Sys '_Reversible_Stress.xlsx'];
    path = fullfile(folder_data, file);

    if isfile(path)
        heat = readmatrix(path, 'Sheet', 5); % Heating leg of data
        T_h = heat(:,2) + 273; % [K]
        eps_h = heat(:,3); % [%]
        if Plot == 1
            figure('Name','Heating at 5 MPa of Detwinned Martensite');
            hold on;
            plot(T_h, eps_h, 'LineWidth', 1);
            xlabel('Temperature (K)');
            ylabel('Strain (%)');
            hold off;
        end
        disp(['Specify the As and Af temperatures from detwinned ' ...
            'martesnite in "ACHANGE_characterize.m". Save then press enter ' ...
            'in the Command Window to continue.']);
    else
        MatProps.Phaseboundary = TwinPB;
        disp(['No Detwinning experiment for', Sys, 'was found. Using ' ...
            'twinned martensite to austenite transition temperatures ' ...
            '(ISHC). Press enter in the Command Window to continue.']);
    end
end