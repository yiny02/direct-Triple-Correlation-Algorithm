function [Correlation] = smTripleCorrelation(coor_1, coor_2, coor_3, roiszx, roiszy, MaxRhoSize, rho_res, MaxRadSize, module)
% coors are [x, y] in the unit of nm
% roisz is in the unit of nm
% output is 3D with the indexing of [r_12, r_13, delta_theta] (col major)
% MaxRadSize is predefined as 90 while rad_res is predefined as 2*pi/90
% MaxRhoSize is predefined as 100 while rho_res is predefined as 5 nm

% =================== Switch to single ====================================
coor_1 = single(coor_1);
coor_2 = single(coor_2);
coor_3 = single(coor_3);
roiszx = single(roiszx);
roiszy = single(roiszy);

% ================== Configure MaxRhos and MaxRads ========================
rad_res = 2*pi/MaxRadSize; % rad
rho_center = ((1:MaxRhoSize)' - 0.5).*rho_res; 
rad_center = ((1:MaxRadSize)' - 0.5).*rad_res - pi;

% ================== Construct Normalization Base =========================
BaseX_A = single(bsxfun(@times, rho_center, cos(rad_center')));
BaseY_A = single(bsxfun(@times, rho_center, sin(rad_center')));

TmpX_A = repmat(BaseX_A', MaxRadSize+1, 1);
TmpY_A = repmat(BaseY_A', MaxRadSize+1, 1);

BaseX_B = single(zeros(MaxRadSize+1, MaxRadSize*MaxRhoSize));
BaseX_B(:) = TmpX_A; 
BaseX_B(end, :) = [];

BaseY_B = single(zeros(MaxRadSize+1, MaxRadSize*MaxRhoSize)); 
BaseY_B(:) = TmpY_A; 
BaseY_B(end, :) = [];

clear TmpX_A TmpY_A;

% =========== Converting the matrix into a row major vector for C =========
BaseX_A = reshape(BaseX_A', MaxRhoSize*MaxRadSize, 1);
BaseY_A = reshape(BaseY_A', MaxRhoSize*MaxRadSize, 1);

BaseX_B = reshape(BaseX_B', MaxRhoSize*MaxRadSize*MaxRadSize, 1);
BaseY_B = reshape(BaseY_B', MaxRhoSize*MaxRadSize*MaxRadSize, 1);

% ====================== Calculating spots density ========================
roisz = roiszx * roiszy;
density_1 = size(coor_1, 1)/roisz;
density_2 = size(coor_2, 1)/roisz;
density_3 = size(coor_3, 1)/roisz;

% ===================== Calculating delta theta and r =====================
deltaX_12 = bsxfun(@minus, coor_2(:, 1), coor_1(:, 1)');
deltaY_12 = bsxfun(@minus, coor_2(:, 2), coor_1(:, 2)');
[Theta_12, Rho_12] = cart2pol(deltaX_12, deltaY_12);
clear deltaX_12 deltaY_12;
Theta_12 = reshape(Theta_12, size(coor_2, 1)*size(coor_1, 1), 1);
Rho_12 = reshape(Rho_12, size(coor_2, 1)*size(coor_1, 1), 1);

deltaX_13 = bsxfun(@minus, coor_3(:, 1), coor_1(:, 1)');
deltaY_13 = bsxfun(@minus, coor_3(:, 2), coor_1(:, 2)');
[Theta_13, Rho_13] = cart2pol(deltaX_13, deltaY_13);
clear deltaX_13 deltaY_13;
Theta_13 = reshape(Theta_13, size(coor_3, 1)*size(coor_1, 1), 1);
Rho_13 = reshape(Rho_13, size(coor_3, 1)*size(coor_1, 1), 1);

% ============= Histogram and Normalize at each coor_1_inner ==============
if module == 'GPU'
    Correlation = smTriCorrGPU(coor_1(:, 1), coor_1(:, 2), coor_2(:, 1), coor_2(:, 2), coor_3(:, 1), coor_3(:, 2), roiszx, roiszy, BaseX_A, BaseY_A, BaseX_B, BaseY_B, Rho_12, Theta_12, Rho_13, Theta_13);
else
    Correlation = smTriCorrCPU(coor_1(:, 1), coor_1(:, 2), coor_2(:, 1), coor_2(:, 2), coor_3(:, 1), coor_3(:, 2), roiszx, roiszy, BaseX_A, BaseY_A, BaseX_B, BaseY_B, Rho_12, Theta_12, Rho_13, Theta_13);
end
Correlation = Correlation / (density_1*density_2*density_3) / roisz;

Correlation = reshape(Correlation, 90, 100, 100);
Correlation = fftshift(Correlation, 1);
Correlation = permute(Correlation, [2 3 1]);