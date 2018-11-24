clearvars;
clc;

% addpath('the path where you save these functions. Add all related paths if these functions are saved in different paths');

%% =========================== SIMULATION =================================
roiszx = 10240; % nm
roiszy = 10240; % nm
dist_1 = 20; % nm
dist_2 = 50; % nm
NumCoor = 1E4;
[coor_R, coor_G, coor_B] = Simulate_TC(roiszx, roiszy, NumCoor, dist_1, dist_2);
% see details in Simulate_TC function

%% ======================= Triple-Correlation =============================
% The Triple-Correlation will be calculated on a polar system grided every rho_res (nm) and rad_res (rad) in radial and angular directions, respectively.
% The Triple-Correlation will be calculated as far as MaxRhoSize * rho_res (nm).

% Parameter settings
rho_res = 5; % nm
rad_res = pi/45; % rad
MaxRhoSize = 100;
MaxRadSize = round(2*pi/rad_res);

% Triple-Correlation Computing

% 1. If the input sequence of the coordinates is R-G-B as in this demo, The calculated result at this step is a function of distance between R-G (nm, first dim for TC_RGB), distance between R-B (nm, second dim for TC_RGB), and the angle between RG and RB (rad, third dim for TC_RGB) 
[TC_RGB] = smTripleCorrelation(coor_R, coor_G, coor_B, roiszx, roiszy, MaxRhoSize, rho_res, MaxRadSize, 'CPU');
% Changing 'CPU' to 'GPU' will perform calculation on GPU, which requires a GPU suppporting for windows 64-bit, CUDA 8.0, and Matlab 2016 or newer.

% 2. Calculating cross-correlation, which will be subtracted from the TC_RGB
%    The cross-correlation will be calculated along (0 : rho_res : 2*MaxRhoSize*rho_res) 
RhoUpper = (1 : 2*MaxRhoSize).*rho_res; % upper edge of each bin, lower edge of the first bin is 0
PC_RG = smPairCorrelation(coor_R, coor_G, roiszx, roiszy, MaxRadSize, RhoUpper);
PC_RB = smPairCorrelation(coor_R, coor_B, roiszx, roiszy, MaxRadSize, RhoUpper);
PC_GB = smPairCorrelation(coor_G, coor_B, roiszx, roiszy, MaxRadSize, RhoUpper);

% 3. Extracting Cross-Correlation
%    The input of PC_RG, PC_RB, PC_GB must be in a sequence where RG is the 1st dim in TC_RGB, RB is the 2nd dim in TC_RGB
%    The final output triple_trans is the triple cube as a fucntion of the distance RG (1st dim), RB (2nd dim), and GB (3rd dim)
[triple_trans] = TripleTrans(TC_RGB, PC_RG, PC_RB, PC_GB, rho_res);
