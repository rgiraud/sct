
clear all
close all

addpath(genpath('mex'));
addpath(genpath('utils'));

%Compilation
if(1)
    % Download more accurate SCALP method at   www.remigiraud.fr/research/scalp.php
    % mex  -O CFLAGS="\$CFLAGS -Wall -Wextra -W -std=c99" SCALP_code/SCALP_mex.cpp -outdir mex/
    mex -O CFLAGS="\$CFLAGS -Wall -Wextra -W -std=c99" mex/SLIC_mex.c -outdir mex/
    mex -O CFLAGS="\$CFLAGS -Wall -Wextra -W -std=c99" mex/SuperPatchMatch.c -outdir mex/
    mex -O CFLAGS="\$CFLAGS -Wall -Wextra -W -std=c99" mex/transfer_sp_color.c -outdir mex/
    mex -O CFLAGS="\$CFLAGS -Wall -Wextra -W -std=c99" mex/SP_desc_comp.c -outdir mex/
    mex -O CFLAGS="\$CFLAGS -Wall -Wextra -W -std=c99" mex/SP_neigh_ordering_fct.c -outdir mex/
    mex -O CFLAGS="\$CFLAGS -Wall -Wextra -W -std=c99" mex/compute_feat_sp_c.c -outdir mex/
end


%% MAIN

%Image loading
filename_a = 'scotland_house.jpg';   %creek.jpg     %clown.jpg    %land1.jpg     %flower_red.jpg     %flower_purple.jpg  
a = double(imread(sprintf('img/%s',filename_a)))/255;

filename_b = 'scotland_plain.jpg';   %desert.jpg    %sea.jpg      %stones3.jpg   %flower_purple.jpg  %flower1.jpg  
b_cell{1} = double(imread(sprintf('img/%s',filename_b)))/255;
% filename_b = 'flower2.jpg'; b_cell{2} = double(imread(sprintf('img/%s',filename_b)))/255;
% filename_b = 'flower3.jpg'; b_cell{3} = double(imread(sprintf('img/%s',filename_b)))/255;


%Parameters
K = 300; %Number of superpixels


%Superpixel-based Color Transfer
[a_b] = sct_core(a,b_cell,K,filename_a,filename_b);





