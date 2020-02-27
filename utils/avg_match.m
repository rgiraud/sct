


function [avg_map,avg_count,avg_nnf] = avg_match(Sa,Sb,SP_nnf,b)

[ha,wa] = size(Sa);
[hb,wb,~,nb] = size(b);
avg_map = zeros(hb,wb,nb);
avg_nnf = zeros(ha,wa);
avg_count = zeros(max(max(Sb(:))),nb);

SP_nbr1 = max(Sa(:));

for i=1:SP_nbr1
    avg_count(SP_nnf(i,1),SP_nnf(i,2)) = avg_count(SP_nnf(i,1),SP_nnf(i,2)) + 1;
end

for n=1:nb
    for p=1:max(max(SP_nnf(:,1)))
        if (avg_count(p,n) > 0)
            sp_pos = Sb(:,:,n) == p;
            tmp = cat(3, repmat(sp_pos*0,[1 1 n-1]), sp_pos, repmat(sp_pos*0,[1 1 nb-n]))>0;
            avg_map(tmp) = avg_count(p,n);
        end
    end
end

for p=1:max(max(Sa(:)))
    sp_pos = Sa == p;
    avg_nnf(sp_pos) = SP_nnf(p,1);
end

end