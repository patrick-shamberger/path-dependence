%% =============================================================================
%{
    TITLE: PATH DEPENDENT STRAIN PREDICTION FOR SHAPE MEMORY ALLOYS
    (CHARACTERIZE)
    AUTHOR: DERIAN R. MORPHEW / PATRICK J. SHAMBERGER
    VERSION: 1 (2026)

    SUMMARY:
        [This script is meant to be interactively changed while running the
        "ARUN_characterize.m" script. The variables here are changed for a
        new material based on the manual inspection of the figure outputs
        from ARUN_characterize.m. After a change is made, this script must
        be saved and then hit "enter" in the Command Window to register the
        changes and continue running the other script].

    NOTES:
        - "====CHANGE ME====" regions are meant to be changed for a new
          material
%}
%  =============================================================================

%% MATERIAL SYSTEM
Materialsys = 'Ni50.9Ti49.1'; % ===CHANGE ME===

%% LOAD ISO-STRESS HEATING COOLING (ISHC) DATA
file = 'Ni50.9Ti49.1.xlsx'; % ===CHANGE ME===
Plot = 1; % (1 - Yes, 0 - No)
%  =================================CHANGE ME=================================== 
ISHCtests = {'5 MPa','50 MPa','100 MPa','150 MPa','200 MPa','300 MPa'};
Stresses = [5, 50, 100, 150, 200, 300]; % [MPa]
ISHCselect = [1, 2, 3, 4, 5, 6];
%  =============================================================================

%% CHANGE MODEL PARAMETERs
% ================================ CHANGE ME ===================================
% === DSC Results ===
T_Af_0MPa = NaN; % [K]
T_As_0MPa = NaN; % [K]
T_Mf_0MPa = NaN; % [K]
T_Ms_0MPa = NaN; % [K]
St_0MPa = 0; % [MPa]

% === Transition temperatures (ISHC) ===
Af = [277.7, 288.6, 297.4, 297.7, 303.5, 314.9]; % [K]
As = [240.7, 245.8, 252.3, 255.9, 257.9, 263.5]; % [K]
Mf = [186, 188, 193.6, 201.8, 215.7, 230.2]; % [K]
Ms = [228, 256.0, 265, 272.3, 283.7, 295.3]; % [K]

% === CTE Bounds ===
A_LT = [233.801, 272.401, 289.800, 295.700, 300.600, 318.200]; % [K]
A_HT = [374, 366.900, 365.600, 369.800, 371.200, 376.500]; % [K]
M_LT = [194.201, 237.001, 197.401, 201.001, 200.401, 200.801]; % [K]
M_HT = [239.501, 263.801, 219.601, 243.501, 244.201, 251.901]; % [K]
% ==============================================================================

% === Detwinned transition temperature at 5 MPa ===
% Leave values NaN if unknown
As_D = As;
As_D(1) = 257; % === CHANGE ME ===
Af_D = Af;
Af_D(1) = 283.5; % === CHANGE ME ===

% ================================ CHANGE ME ===================================
% === Reference points ===
Ref_T_A = 310.5; % [K]
Ref_stress_A = 5; % [MPa]
Ref_T_M = 186.5; % [K]
Ref_stress_M = 5; % [MPa]

% === Plotting grid parameters ===
Plotting.Grid.N = 1180; % Number of grid points
Plotting.Grid.TempRange = [150, 350];
Plotting.Grid.StressRange = [5, 300];

% === Transition Region Parameters ===
Plotting.Sigmoid.MsMf.T0_scale = 0.9253; % Scale applied to Ms_T
Plotting.Sigmoid.MsMf.k = 0.22;

Plotting.Sigmoid.AsAf.T0_scale = 0.9560; % Scale applied to Af_T
Plotting.Sigmoid.AsAf.k = 0.6;
% ==============================================================================
fprintf('Edits applied at %s\n', datetime("now"));