%% prepare workspace before calculation
clear all
clear classes
close all
clc

%% define some variables
lambda = 0.532;     % wavelength of the focused light
NA = 0.95;          % numerical aperture of the lens

rMin = -1;          % size of calculated area (cross section)
rMax = 1;
discretization = 80;% discretization of calculated area

zMin = -0.8;        % size of calculated area (section along propagation axis)
zMax = 0.8;
discretizationZ = discretization;

saveFlag = true;

%% calculation
r = azimuthal(lambda, 0, NA);       % used for calculation of azimuthaly polarized (particularly high-order) beams
r = r.setRBorders(rMin, rMax, discretization);
r = r.setZBorders(zMin, zMax, discretizationZ);

if saveFlag
    temp1 = fix(clock);
    curFolder = strcat(['..\..\data\ver12\', num2str(temp1(1)), '-', num2str(temp1(2)), '-', num2str(temp1(3)), '-', num2str(temp1(4)), '-', num2str(temp1(5)), '-', num2str(temp1(6)) '\']);
    mkdir(curFolder);
    r.storeFlag = false;
end

r.mVortex = 0;      % topological charge of the beam (0 - no vortical phase)
r.mAzimBeam = 2;    % order of the beam (0 - linearly polarized light, 1 - azimuthally polarized light, 2 - second-order CVB with energy backflow on the optical axis)
r.outFolder = curFolder;
r = r.plotAllFocus();     % start calculation of all components in cross section

% r = r.plotAllAlongZ();    % start calculation of all components in the section along propagation axis