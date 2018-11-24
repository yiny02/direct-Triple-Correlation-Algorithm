function [coor_R, coor_G, coor_B] = Simulate_TC(roiszx, roiszy, NumCoor, dist_1, dist_2)

% Given the ROI sizes as roiszx and roiszy, and the number of coordinates for each of the 3 channel, this script generate a list of coordinates for each channel
% The output coodinates have the same unit with the given roisz

% 50% R-G-B and 50% R-B-G linear patterns will be simulated
% for R-G-B pattern, it will have the dist_1 for R-G and dist_2 for R-B
% for R-B-G pattern, it will have the dist_1 for R-B and dist_2 for R-G
% dist_1 and dist_2 are as the same unit of roisz

% randomizing the first 50% of coor_R
x_R_1 = rand(NumCoor/2, 1) .* roiszx;
y_R_1 = rand(NumCoor/2, 1) .* roiszy;

% randomizing the orientation of the patterns
theta_1 = (rand(NumCoor/2, 1) - 0.5) .* 2*pi;

% construct the R-G-B pattern
x_G_1 = x_R_1 + dist_1 .* cos(theta_1);
y_G_1 = y_R_1 + dist_1 .* sin(theta_1);
x_B_1 = x_R_1 + dist_2 .* cos(theta_1);
y_B_1 = y_R_1 + dist_2 .* sin(theta_1);

% randomizing the second 50% of coor_R
x_R_2 = rand(NumCoor/2, 1) .* roiszx;
y_R_2 = rand(NumCoor/2, 1) .* roiszy;

% randomizing the orientation of the patterns
theta_2 = (rand(NumCoor/2, 1) - 0.5) .* 2*pi;

% construct the R-G-B pattern
x_G_2 = x_R_2 + dist_2 .* cos(theta_2);
y_G_2 = y_R_2 + dist_2 .* sin(theta_2);
x_B_2 = x_R_2 + dist_1 .* cos(theta_2);
y_B_2 = y_R_2 + dist_1 .* sin(theta_2);

% generate the output coordinate lists
x_R = cat(1, x_R_1, x_R_2);
y_R = cat(1, y_R_1, y_R_2);
coor_R = [x_R, y_R];

x_G = cat(1, x_G_1, x_G_2);
y_G = cat(1, y_G_1, y_G_2);
coor_G = [x_G, y_G];

x_B = cat(1, x_B_1, x_B_2);
y_B = cat(1, y_B_1, y_B_2);
coor_B = [x_B, y_B];



