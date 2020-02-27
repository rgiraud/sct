

function [a_b] = sct_core(a, b_cell, K, filename_a, filename_b)

clc


fprintf('\nProcessing of %s with %s\n', filename_a, filename_b);

%Merging images in B
[b,~,~,lib_img_dims] = merge_img(b_cell);
nb_img = size(b,4);


%% Superpixel decompositions
fprintf('Superpixel decompositions ...');

%Compactness parameter
if (exist('mex/SCALP_mex.cpp'))
    compactness = 0.05;
else %SLIC
    compactness = 10;
end

tic;

if (exist('mex/SCALP_mex.cpp'))
    S = SCALP_mex(uint8(a*255), K, compactness)+1;
else
    [S,~] = SLIC_mex(uint8(a*255),K,compactness); S=S+1;
end
S = scan_order_relabeling_SP(S);


%Superpixel features (normalized cumulative color histogram)
h_bins = 8;
SuperPatch_R = 1; %using only central superpixels
SP_nbra = max(S(:));
[SP_centera, SP_hista] = SP_desc_comp(single(single(a/max(a(:)))), SP_nbra, int32(S), h_bins);

[mat_adja, perima] = SP_adj_graph_building(S);
perima = double((repmat(~perima,1,1,3)));
[mat_adja, ~] = SP_neigh_ordering_fct(int32(mat_adja), single(SP_centera), SP_nbra, SuperPatch_R);
l_vecta = int32(sum(mat_adja > 1,2));
l_maxa = max(l_vecta(:));

%B
SP_nbr_max = 0;
for i=1:nb_img
    if (exist('mex/SCALP_mex.cpp'))
        S_tmp = SCALP_mex(uint8(b_cell{i}*255), K, compactness)+1;
    else
        [S_tmp,~] = SLIC_mex(uint8(b_cell{i}*255),K,compactness); S_tmp=S_tmp+1;
    end
    S_i{i} = scan_order_relabeling_SP(S_tmp);
    SP_nbr_max = max(max(S_i{i}(:)),SP_nbr_max);
end

mat_adjb = zeros(SP_nbr_max,SP_nbr_max,nb_img);
lab_mapb = zeros(size(b,1),size(b,2),nb_img);
SP_centerb = zeros(SP_nbr_max,2,nb_img);
SP_histb = zeros(SP_nbr_max,size(SP_hista,2),nb_img);
l_vectb = int32(zeros(SP_nbr_max, nb_img));

%Superpixel features
for i=1:nb_img
    hb = lib_img_dims(i,1);
    wb = lib_img_dims(i,2);
    S_tmp = S_i{i};
    SP_nbr_i = max(S_tmp(:));
    lab_mapb(1:hb,1:wb,i) = S_tmp;
    [SP_center_i, SP_hist_i] = SP_desc_comp(single(single(b_cell{i}/max(b_cell{i}(:)))), SP_nbr_i, int32(S_tmp), h_bins);
    SP_centerb(1:SP_nbr_i,1:2,i) = SP_center_i;
    SP_histb(1:SP_nbr_i,:,i) = SP_hist_i;
    [mat_adj_i, perim_tmp] = SP_adj_graph_building(S_tmp);
    [mat_adj_i, ~] = SP_neigh_ordering_fct(int32(mat_adj_i), single(SP_center_i), SP_nbr_i, SuperPatch_R);
    mat_adjb(1:SP_nbr_i,1:SP_nbr_i,i) = int32(mat_adj_i);
    l_vectb(1:SP_nbr_i,i) = int32(sum(mat_adj_i > 1,2));
end
l_maxb = max(l_vectb(:));

b_toc = toc;
fprintf('Ok in %f\n', b_toc);


%% Matching

tic;
fprintf('Matching ... ');

epsilon = 3;   %Number of times a superpixel in B can be selected
iterations = 20;

[SP_nnf_ex] = SuperPatchMatch(int32(mat_adja), int32(mat_adjb), single(SP_hista), single(SP_histb), ...
    size(SP_hista,2), single(SP_centera), single(SP_centerb), ...
    SP_nbra, SP_nbr_max, int32(lab_mapb), SuperPatch_R, iterations, 1, ...
    int32(l_vecta), l_maxa, l_maxb, int32(lib_img_dims), epsilon)+1;


fprintf('Ok in %f\n', toc);


% %%Display Matching
% %display_match(S,lab_mapb,SP_nnf_ex,a,b_cell,lib_img_dims,1);
% [avg_map,~,avg_nnf] = avg_match(S,lab_mapb,SP_nnf_ex,b);
%
% figure,
% subplot 121,
% imagesc(avg_map), colormap(gray)
% subplot 122
% imagesc(avg_nnf)
% avg_map_s = (avg_map - min(avg_map(:)))/(max(avg_map(:)) - min(avg_map(:)));
% avg_map_s(1,:) = 1; avg_map_s(end,:) = 1; avg_map_s(:,1) = 1; avg_map_s(:,end) = 1;
% avg_map_s = 1 - avg_map_s;


%% Color transfer

fprintf('Color transfer ... ');

%Parameters
SP_radius = 0.3; % Limits the contributions of neighboring superpixels
sigma_color = 1;
var_color=1e-1;
var_spatiale=1e1;

Epsilon=eye(5)*1e-4;
filtre_decorelation= [1. 1. 0. 0. 0.;  1. 1. 0. 0 0; 0. 0.  1. 1. 1. ;0. 0.  1. 1. 1. ; 0. 0.  1. 1. 1. ];
variance_vector=([var_spatiale var_spatiale var_color var_color var_color]);
filtre_covariance=filtre_decorelation.*(variance_vector'*variance_vector);

%A
[Sj,cov_sp_a] = compute_feat_sp_c(double(a),int32(S),max(S(:)),sigma_color);
max_c = 0;
for iii=1:size(cov_sp_a,3)
    cov_tmp = cov_sp_a(:,:,iii);
    cov_filt = cov_tmp.*filtre_covariance+Epsilon;
    max_c = max(max_c,max(cov_filt(1:2,1:2)));
    cov_tmp = inv(cov_tmp.*filtre_covariance+Epsilon);
    cov_sp_a(:,:,iii) = cov_tmp;
end

%B
Sj_b_cell = cell(size(lab_mapb,3));
Sj_b = zeros(SP_nbr_max,5,size(lab_mapb,3));
for p=1:size(lab_mapb,3)
    [Sk,~] = compute_feat_sp_c(double(b_cell{p}),int32(lab_mapb(1:lib_img_dims(p,1),1:lib_img_dims(p,2),p)),...
        max(max(lab_mapb(1:lib_img_dims(p,1),1:lib_img_dims(p,2),p))),sigma_color);
    Sj_b(1:max(max(lab_mapb(1:lib_img_dims(p,1),1:lib_img_dims(p,2),p))),:,p) = Sk';
    Sj_b_cell{p} = Sk;
end

%Color transfer
tic;
[a_c] = transfer_sp_color(single(a),int32(S),int32(SP_nnf_ex), single(Sj'),single(cov_sp_a),single(Sj_b), sigma_color, SP_radius);
fprintf('Ok in %f\n', toc);

%Regrain
a_c(isnan(a_c)) = a(isnan(a_c));
a_b = color_regrain(a,a_c);


%Display
figure,
kk(1) = subplot(131);
imshow(a)
title('A')
subplot(132)
imshow(b_cell{1})
title('B')
kk(2) = subplot(133);
imshow(a_b)
title('SCT result')
linkaxes(kk)



%% Save results

imwrite(a,sprintf('res/%s_a.png', filename_a(1:end-4)));
imwrite(b_cell{end},sprintf('res/%s_b.png', filename_a(1:end-4)));
imwrite(a_b,sprintf('res/%s_transfer.png', filename_a(1:end-4)));

imwrite(a.*perima,sprintf('res/%s_a_perim.png', filename_a(1:end-4)));
imwrite(b_cell{end}.*repmat(~perim_tmp,[1 1 3]),sprintf('res/%s_b_perim.png', filename_a(1:end-4)));





end