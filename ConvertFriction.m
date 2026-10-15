
clear;
close all;

addpath /home/daniel/models/Kori-ULB/subroutines/;
addpath /home/daniel/models/Kori-ULB/;

% Paths.
pathIn      = '/home/daniel/models/Kori-ULB/ice_data/eta1e7/ground_melt_0/';
pathOut     = [pathIn, 'basal_friction/coulomb/u0/m_05/'];

fileIn      = 'INIT_SSA_3.mat';
fileOut     = 'INIT_SSA_3_m_05_u0_150.mat';

filePathIn  = [pathIn, fileIn];
filePathOut = [pathOut, fileOut];

% Settings.
mIn  = 3.0;
mOut = 5.0;
u0   = 150;

mkdir(pathOut)


convPowerLaw(filePathIn, filePathOut, mIn, mOut, u0)
