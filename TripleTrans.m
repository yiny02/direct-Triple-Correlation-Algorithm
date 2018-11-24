function [triple_trans] = TripleTrans(triplemat, PC_1, PC_2, PC_3, rho_res)

SizeofTriple = size(triplemat); % SizofTriple(1) == SizeofTriple(2)
dim = SizeofTriple(1);

rad_res = 2*pi/SizeofTriple(3);
rad_center = ((1:SizeofTriple(3))' - 0.5).*rad_res - pi;

% ======================== Extracting CrossCorrelation ====================

for ii = 1 : SizeofTriple(1)
    edge_1 = PC_1(ii, 1);
    for jj = 1 : SizeofTriple(2)
        edge_2 = PC_2(jj, 1);
        for kk = 1 : SizeofTriple(3)
            theta = rad_center(kk);
            edge_3 = sqrt(edge_1^2 + edge_2^2 - 2*edge_1*edge_2*cos(theta));
            trans_kk = ceil(edge_3 / rho_res);
            triplemat(ii, jj, kk) = triplemat(ii, jj, kk) - (PC_1(ii, 2)-1) - (PC_2(jj, 2)-1) - (PC_3(trans_kk, 2)-1);
        end
    end
end
       
% ==================== building up transforming cube ======================
triple_trans = ones(dim, dim, dim);

trans_cube_1 = repmat((1:dim)', [1, dim, SizeofTriple(3)]);
trans_cube_2 = repmat((1:dim),  [dim, 1, SizeofTriple(3)]);
trans_cube_3 = repmat(permute(rad_center, [2 3 1]), [dim, dim, 1]);
trans_cube = sqrt(trans_cube_1.^2 + trans_cube_2.^2 - 2.*trans_cube_1.*trans_cube_2.*cos(trans_cube_3));

% ====================== transforming and smoothing ======================= 
rho_edge = (0 : dim)'; 
for d1 = 1 : dim
    for d2 = 1 : dim
        tmp_vt = permute(trans_cube(d1, d2, :), [3 2 1]);
        [~, ~, bin] = histcounts(tmp_vt, rho_edge);
        for d3 = 1 : dim
            if ~isempty(triplemat(d1, d2, bin == d3))
                triple_trans(d1, d2, d3) = mean(triplemat(d1, d2, bin == d3));
            end
        end
    end
end
