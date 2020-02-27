

function [adj_graph_mat, perim] = SP_adj_graph_building(segments)

if (min(segments(:)) <= 0)
   segments = segments - min(segments(:)) + 1; 
end

[h, w] = size(segments);
nbr_sp = max(segments(:));
adj_graph_mat = zeros(nbr_sp);
perim = zeros(size(segments));

for i=1:h-1
    for j=1:w-1
        
        label = segments(i,j);
        
        if (label ~= segments(max(i-1,1),j))
            adj_graph_mat(label, segments(max(i-1,1),j)) = 1;
            adj_graph_mat(segments(max(i-1,1),j), label) = 1;
            perim(i,j) = 1;
        end
        if (label ~= segments(i,j+1))
            adj_graph_mat(label, segments(i,j+1)) = 1;
            adj_graph_mat(segments(i,j+1), label) = 1;
            perim(i,j) = 1;
        end

    end
end


end

