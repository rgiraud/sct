

function [lab_map_new,lab_cor] = scan_order_relabeling_SP(lab_map)

if (min(lab_map(:)) <= 0)
    lab_map = lab_map - min(lab_map(:)) + 1;
end

SP_nbr = max(lab_map(:));
lab_vect = zeros(SP_nbr,1);
[h,w] = size(lab_map);
lab_map_new = zeros(size(lab_map));
lab_cor = zeros(length(unique(lab_map(:))),1);

c = 1;
for i=1:h
    for j=1:w
        label = lab_map(i,j);
        if (lab_vect(label,1) == 0)  
            SP_pos = (lab_map == label);
            lab_map_new(SP_pos) = c;
            lab_cor(c,1) = label;
            c = c + 1;
            lab_vect(label,1) = 1;
        end
    end
end


end