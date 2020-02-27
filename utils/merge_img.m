


function [b,b_g,b_cell_g,lib_img_dims] = merge_img(lib_cell)

nb_img = length(lib_cell);
maxy=0;
maxx=0;
for i=1:nb_img
    if (size(lib_cell{i},1) > maxy)
        maxy = size(lib_cell{1},1);
    end
    if (size(lib_cell{i},2) > maxx)
        maxx = size(lib_cell{1},2);
    end
end

b = zeros(maxy,maxx,3,nb_img);
b_g = zeros(maxy,maxx,3,nb_img);
b_cell_g = cell(1,nb_img);
for i=1:nb_img
    img = lib_cell{i};
    [h,w,z] = size(img);
    if (z == 1)
        b(1:h,1:w,1:3,i) = repmat(img, [1 1 3]);
        b_g(1:h,1:w,1:3,i) = repmat(img, [1 1 3]);
        b_cell_g{i} = repmat(img, [1 1 3]);
    else
        b(1:h,1:w,1:3,i) = img;
        b_g(1:h,1:w,1:3,i) = repmat(rgb2gray(img),[1 1 3]);
        b_cell_g{i} =  repmat(rgb2gray(img),[1 1 3]);
    end
end

lib_img_dims = zeros(nb_img,2);
for p=1:nb_img
    lib_img_dims(p,1) = size(lib_cell{p},1);
    lib_img_dims(p,2) = size(lib_cell{p},2);
end


end