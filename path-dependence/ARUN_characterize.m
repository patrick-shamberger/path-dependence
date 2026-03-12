%% =============================================================================
%{
    TITLE: PATH DEPENDENT STRAIN PREDICTION FOR SHAPE MEMORY ALLOYS
    (CHARACTERIZE)
    AUTHOR: DERIAN R. MORPHEW / PATRICK J. SHAMBERGER
    VERSION: 1 (2026)

    SUMMARY:
        [This script characterizes key material properties, and their
         dependencies, from experimental ISHC tests. Linear coefficient of 
         thermal expansion as a function of stress, α(σ), is determined for
         each phase. Strain, ε(σ), is determined at a reference temperature
         for each phase in addition to the compliance. The phase boundary 
         points are manually chosen then linearized for each phase 
         boundary line. The shape of the transformation envelopes are 
         characterized using a two-asymtpote logistic function. Full ε(T,σ) 
         grids are created for each single phase and coexistence boundary. 
         These grids and properties are to be used in "ARUN_predict" to predict
         elastocaloric effect (eCE) type thermomechanical pathways.]

    NOTES:
        - There are pauses in the script to manually inspect key figures to
          determine further inputs to the script.
        - These pauses require you to change values in the
          "ACHANGE_characterize.m" script interactively.
%}
%  =============================================================================
%% MATERIAL SYSTEM
clear; close all; clc;

dirFolder  = fileparts(mfilename('fullpath'));
inputsFile = fullfile(dirFolder, 'ACHANGE_characterize.m');

%% LOAD USER INPUTS
% === Load File ===
dirFolder = fileparts(mfilename('fullpath'));
dataFolder = fullfile(dirFolder, 'Experimental Data', 'ISHC');
run('ACHANGE_characterize.m');
path = fullfile(dataFolder, file);

% === Load Data ===
ISHCdata = cell(numel(ISHCselect),Plot);
for i = 1:length(ISHCselect)
    ISHCdata{i} = readmatrix(path, 'Sheet', ISHCtests{ISHCselect(i)});
end

%% DEFINE INITIAL INPUTS
% == Make Structures ===
InitialInputs.Materialsys = Materialsys;
InitialInputs.ISHCdata = ISHCdata;
InitialInputs.Stresses = Stresses;

% === Store 0 MPa ===
PB0 = struct('Stress', St_0MPa, ...
                'Mf', T_Mf_0MPa, ...
                'Ms', T_Ms_0MPa, ...
                'As', T_As_0MPa, ...
                'Af', T_Af_0MPa);

%% INSPECT DATA
% === Display ISHC curves ===
ISHCdis(InitialInputs, Plot);

% === Inspect the curves and input the PB Point & CTE indices ===
disp(['Specify transformation temps and indices for CTE calculation in ' ...
    '"ACHANGE_characterize.m". Save then press enter in the Command Window ' ...
    'to continue.']);
input('', 's'); rehash; run(inputsFile);

% === Structure for non 0MPA points ===
PB_rest = arrayfun(@(s,mf,ms,as,af) ...
            struct('Stress',s,'Mf',mf,'Ms',ms,'As',as,'Af',af), ...
            Stresses, Mf, Ms, As, Af).';

% === Combine PB ===
MatProps.TwinPhaseboundary = [PB0; PB_rest];

% === Replace ISHC PB with Detwinned values (if available) ===
MatProps = detPB(InitialInputs, MatProps, Plot);
input('', 's'); rehash; run(inputsFile);
if ~isfield(MatProps, 'Phaseboundary') || isempty(MatProps.Phaseboundary)
    PB = MatProps.TwinPhaseboundary;
    for k = 1:numel(PB)-1
        PB(k).As = As_D(k);
        PB(k).Af = Af_D(k);
    end
    MatProps.Phaseboundary = PB;
end

% === Structure for CTE bounds ===
CTE_bounds = arrayfun(@(s,al,au,ml,mu) ...
            struct('Stress',s, ...
                'Austenite', struct('LowerT',al,'UpperT',au), ...
                'Martensite', struct('LowerT',ml,'UpperT',mu)), ...
            Stresses, A_LT, A_HT, M_LT, M_HT).';
InitialInputs.CTEIndices = CTE_bounds;

%% LINEARIZE PHASE BOUNDARIES
PB = MatProps.Phaseboundary(:);
St = [PB.Stress].';
Mf = [PB.Mf].';
Ms = [PB.Ms].';
As = [PB.As].';
Af = [PB.Af].';

% === Fit lines ===
p_Mf = polyfit(St(~isnan(Mf)), Mf(~isnan(Mf)), 1);
p_Ms = polyfit(St(~isnan(Ms)), Ms(~isnan(Ms)), 1);
p_As = polyfit(St(~isnan(As)), As(~isnan(As)), 1);
p_Af = polyfit(St(~isnan(Af)), Af(~isnan(Af)), 1);

% === Evaluate at same stress points ===
Mf_lin = polyval(p_Mf, St);
Ms_lin = polyval(p_Ms, St);
As_lin = polyval(p_As, St);
Af_lin = polyval(p_Af, St);

% === Store linear phase boundaries ===
for k = 1:numel(PB)
    PB(k).Mf = Mf_lin(k);
    PB(k).Ms = Ms_lin(k);
    PB(k).As = As_lin(k);
    PB(k).Af = Af_lin(k);
end
MatProps.Phaseboundary = PB;

% === Display Phase Boundary ===
PBdis(MatProps, Plot);

% === Inspect the curve and input the reference Temp/Stresses ===
disp(['Specify reference temperature and stress for each single phase in ' ...
    '"ACHANGE_characterize.m". Save then press enter in the Command Window ' ...
    'to continue.']);
input('', 's'); rehash; run(inputsFile);

MatProps.Reference.Austenite = struct(...
    'Temperature', Ref_T_A, ...
    'Stress', Ref_stress_A, ...
    'Strain', NaN);
MatProps.Reference.Martensite = struct(...
    'Temperature', Ref_T_M, ...
    'Stress', Ref_stress_M, ...
    'Strain', NaN);

%% DETERMINE SINGLE PHASE PROPERTIES
[InitialInputs, MatProps] = singlephase(InitialInputs, MatProps, Plot);

%% MODEL HYSTERESIS OUTER ENVELOPES
% === Inspect the curve and input the reference Temp/Stresses below ===
disp(['Specify sigmoid parameters T0 and k for both transitions in ' ...
    '"ACHANGE_characterize.m". Save then press enter in the Command Window ' ...
    'to continue.']);
input('', 's'); rehash; run(inputsFile);

[InitialInputs, MatProps, Plotting] = ... 
    hys(InitialInputs, MatProps, Plotting, Plot);

%% SAVE STRUCTURES AS .MAT FILES
parentFolder = fullfile(dirFolder, 'Modeled Data');
materialFolder = fullfile(parentFolder, Materialsys);
if ~exist(materialFolder, 'dir')
    mkdir(materialFolder);
end
saveFile = fullfile(materialFolder, [Materialsys '_ModelOutputs.mat']);
save(saveFile, 'InitialInputs', 'MatProps', 'Plotting');
fprintf('Structures saved to: %s\n', saveFile);