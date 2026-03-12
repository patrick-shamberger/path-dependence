%% =============================================================================
%{
    TITLE: PATH DEPENDENT STRAIN PREDICTION FOR SHAPE MEMORY ALLOYS
    AUTHOR: DERIAN R. MORPHEW / PATRICK J. SHAMBERGER
    VERSION: 1 (2026)

    SUMMARY:
        [This script uses the outputs from "ARUN_characterize" to predict
         elastocaloric effect (eCE) type thermomechanical pathways. Strain
         grids are corrected to account for path dependence.]

    NOTES:
        - "====CHANGE ME====" regions are meant to be changed for a new
          material
%}
%  =============================================================================
clear; close all; clc;

%% SPECIFY MATERIAL SYSTEM
Materialsys = 'Ni50.9Ti49.1';  % ===== CHANGE ME =====

%% LOAD STRUCTURES
scriptfolder = fileparts(mfilename('fullpath'));
datafolder = fullfile(scriptfolder, 'Modeled Data', Materialsys);
filename = [Materialsys '_ModelOutputs.mat'];
filepath = fullfile(datafolder, filename);
load(filepath, 'InitialInputs', 'MatProps', 'Plotting');
fprintf('Material system: %s\n', Materialsys);
Plotting.Materialsys = Materialsys;

%% Make Paths
n_leg = 50; % Points per leg

% === Specify Path "Corners" ===
all_paths = {
    [ % Path 1
        295, 5; % Start: 295 K, 5 MPa
        275, 5; % Forward Leg 1 end
        275, 199.6; % Forward Leg 2 end
        230, 199.6; % Forward Leg 3 end
        275, 199.6; % Reverse Leg 3 end
        275, 5; % Reverse Leg 2 end
        295, 5; % End: 295 K, 5 MPa
    ];
    [ % Path 2
        295, 5; % Start: 295 K, 5 MPa
        261.9, 5; % Forward Leg 1 end
        261.9, 138; % Forward Leg 2 end
        239.7, 138; % Forward Leg 3 end
        239.7, 200; % Forward Leg 4 end
        230, 200; % Forward Leg 5 end
        239.9, 200; % Reverse Leg 5 end
        239.9, 138.8; % Reverse Leg 4 end
        262.4, 138.8; % Reverse Leg 3 end
        262.4, 5; % Reverse Leg 2 end
        295, 5; % End: 295 K, 5 MPa
    ];
    [ % Path 3
        295, 5; % Start: 295 K, 5 MPa
        260.8, 5; % Forward Leg 1 end
        260.8, 200; % Forward Leg 2 end
        230, 200; % Forward Leg 3 end
        260.6, 200; % Reverse Leg 3 end
        260.6, 5; % Reverse Leg 2 end
        295, 5; % End: 295 K, 5 MPa
    ];
    [ % Path 4
        295, 5; % Start: 295 K, 5 MPa
        250.3, 5; % Forward Leg 1 end
        250.3, 199.2; % Forward Leg 2 end
        230, 199.2; % Forward Leg 3 end
        251, 199.2; % Reverse Leg 3 end
        251, 5; % Reverse Leg 2 end
        295, 5; % End: 295 K, 5 MPa
    ];
    [ % Path 5
        295, 5; % Start: 295 K, 5 MPa
        241.1, 5; % Forward Leg 1 end
        241.1, 200; % Forward Leg 2 end
        230, 200; % Forward Leg 3 end
        241, 200; % Reverse Leg 3 end
        241, 5; % Reverse Leg 2 end
        295, 5; % End: 295 K, 5 MPa
    ];
    [ % Path 6
        295, 5; % Start: 295 K, 5 MPa
        232.3, 5; % Forward Leg 1 end
        232.3, 200.7; % Forward Leg 2 end
        232.3, 5; % Reverse Leg 2 end
        295, 5; % End: 295 K, 5 MPa
    ];
};

Paths = makepath(all_paths, n_leg);

%% Path Dependence Implementation
[Paths, PB] = avgstress(Paths, MatProps, Plotting, 0);
Type = 1; % Path-Dependence Approach
Paths = pathdep(Paths, MatProps, Plotting, Type, 0);

%% Construct Strain Grids
[Paths, MatProps] = detwinmart(Paths, MatProps, Plotting, 0);
Paths = strainpred(Paths, Plotting, 1, 0);

%% Calculate Phase Fractions
Paths = phasepop(Paths, MatProps, Plotting, 0);

%% Populate Strain Values
Paths = strainpop(Paths, MatProps, Plotting, 1);