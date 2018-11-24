function [PairCorrelation] = smPairCorrelation(coor_1, coor_2, roiszx, roiszy, MaxRadSize, RhoUpper)
% coors are [x, y] in the unit of nm
% roisz is in the unit of nm
% RhoUpper is the upper edge of the radial bins

roisz = roiszx * roiszy;
density_1 = size(coor_1, 1)/roisz;
density_2 = size(coor_2, 1)/roisz;

RadRes = 2*pi / MaxRadSize;

% Constructing bin edges
RhoLower = [0, RhoUpper(1:end-1)];
bins = (RhoUpper+RhoLower)./2;

% Reading coordinates
x_1 = single(coor_1(:, 1));
y_1 = single(coor_1(:, 2));
x_2 = single(coor_2(:, 1));
y_2 = single(coor_2(:, 2));

ind_inner = x_1 >= RhoUpper(end) & x_1 <= roiszx - RhoUpper(end) & y_1 >= RhoUpper(end) & y_1 <= roiszy - RhoUpper(end);
x_1_inner = x_1(ind_inner);
y_1_inner = y_1(ind_inner);
x_1_outer = x_1(~ind_inner);
y_1_outer = y_1(~ind_inner);

% ===================== Preparing normalizations ==========================
NormUnit = bins .* (RhoUpper-RhoLower).*RadRes;

Norm_inner = MaxRadSize .* NormUnit;

if ~isempty(x_1_outer)
    NormMask_outer = zeros(size(x_1_outer, 1), length(RhoUpper));
    for ii = 1 : MaxRadSize
        x_tmp = bsxfun(@plus, x_1_outer, RhoUpper.*cos((ii-0.5)*RadRes));
        y_tmp = bsxfun(@plus, y_1_outer, RhoUpper.*sin((ii-0.5)*RadRes));
        idx = x_tmp >= 0 & x_tmp <= roiszx & y_tmp >= 0 & y_tmp <= roiszy;
        NormMask_outer(idx) = NormMask_outer(idx) + 1;
        clear x_tmp y_tmp idx;
    end
    Norm_outer = NormMask_outer .* NormUnit;
else
    Norm_outer = zeros(1, length(RhoUpper));
end

% ========================= Histogramming =================================
BinEdges = [RhoLower(1), RhoUpper]';
% ----------------------- calculating inners ------------------------------
if ~isempty(x_1_inner)
    PairCorrelation_inner = zeros(size(x_1_inner,1), length(RhoUpper));
    for ii = 1 : size(x_1_inner, 1)
        delta_x_inner = x_2' - x_1_inner(ii);
        delta_y_inner = y_2' - y_1_inner(ii);
        delta_r_inner = sqrt(delta_x_inner.^2 + delta_y_inner.^2);
        clear delta_x_inner delta_y_inner;
        counts_tmp = histcounts(delta_r_inner(delta_r_inner <= RhoUpper(end)), BinEdges);
        PairCorrelation_inner(ii, :) = counts_tmp ./ Norm_inner;
    end
    PairCorrelation_inner = (sum(PairCorrelation_inner, 1))';
    clear delta_r_inner;
else
    PairCorrelation_inner = zeros(length(RhoUpper), 1);
end

% ----------------------- calculating outers ------------------------------
if ~isempty(x_1_outer)
    PairCorrelation_outer = zeros(size(x_1_outer,1), length(RhoUpper));
    for ii = 1 : size(x_1_outer, 1)
        delta_x_outer = x_2' - x_1_outer(ii);
        delta_y_outer = y_2' - y_1_outer(ii);
        delta_r_outer = sqrt(delta_x_outer.^2 + delta_y_outer.^2);
        clear delta_x_outer delta_y_outer;
        counts_tmp = histcounts(delta_r_outer(delta_r_outer <= RhoUpper(end)), BinEdges);
        PairCorrelation_outer(ii, Norm_outer(ii,:)~=0) = counts_tmp(Norm_outer(ii,:)~=0) ./ Norm_outer(ii, Norm_outer(ii,:)~=0);
    end
    PairCorrelation_outer = (sum(PairCorrelation_outer, 1))';
    clear delta_r_outer;
else
    PairCorrelation_outer = zeros(length(RhoUpper), 1);
end
PairCorrelation = (PairCorrelation_outer + PairCorrelation_inner)/(density_1*density_2) / roisz;
bins = bins';

PairCorrelation = [bins, PairCorrelation];

