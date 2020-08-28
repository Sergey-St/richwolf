%% high-order (2) azimuthal polarization

clear all
clear classes
close all
clc

lambda = 0.532;
NA = 0.95;

rMin = -0.2;
rMax = 0.2;
discretization = 40;

zMin = -0.8;
zMax = 0.8;
discretizationZ = discretization;

saveFlag = true;

r = azimuthal(lambda, 0, NA);
r = r.setRBorders(rMin, rMax, discretization);

r = r.setZBorders(zMin, zMax, discretizationZ);

if saveFlag
    temp1 = fix(clock);
    curFolder = strcat(['..\..\data\ver12\', num2str(temp1(1)), '-', num2str(temp1(2)), '-', num2str(temp1(3)), '-', num2str(temp1(4)), '-', num2str(temp1(5)), '-', num2str(temp1(6)) '\']);
    mkdir(curFolder);
    r.storeFlag = false;
end

r.mVortex = 0;
r.mAzimBeam = 3.5;
r.outFolder = curFolder;
r = r.plotAllFocus();

% r = r.plotAllAlongZ();