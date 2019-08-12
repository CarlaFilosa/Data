function [struct_pair] = PairsOrg(pairs_r,NP)
% The function builds a structure that organizes the pairs so that they are
% clustered by binsizes

for k=1:size(pairs_r,2)
pairs_lag_XXX{k}=[];
NP_XXX{k}=[];
end

for k=1:size(pairs_r,2)

pairs_lag_XXX{k}=pairs_r{k};
NP_XXX{k}=NP{k};
% struct for pair to overlap the bin in the figure

struct_pair{k} = struct ('pair',[], 'lag',[], 'pv', [],'bin',[],'ind',[],'new_ind',[],'imean',[]);
j=1; % j is uncremented in the second elseif
   
for i = 1: size(pairs_lag_XXX{k},1)
     if i == 1
         struct_pair{k}.pair{j}(1,1:2) = pairs_lag_XXX{k}(1,1:2);
         struct_pair{k}.lag{j} = pairs_lag_XXX{k}(1,3);
         struct_pair{k}.pv{j} = pairs_lag_XXX{k}(1,4);
         struct_pair{k}.bin{j} = pairs_lag_XXX{k}(1,5);
         struct_pair{k}.ind{j} = pairs_lag_XXX{k}(1,6);
         struct_pair{k}.new_ind{j}(i,:) =j;
         struct_pair{k}.imean{j}(i,:) =i;
    elseif i>1 & ismember(pairs_lag_XXX{k}(i,1:2),pairs_lag_XXX{k}(i-1,1:2)) == [1 ,1],
        struct_pair{k}.pair{j}(i,1:2) = pairs_lag_XXX{k}(i,1:2);
        struct_pair{k}.lag{j}(i,:) = pairs_lag_XXX{k}(i,3);
        struct_pair{k}.pv{j}(i,:) = pairs_lag_XXX{k}(i,4);
        struct_pair{k}.bin{j}(i,:) = pairs_lag_XXX{k}(i,5);
        struct_pair{k}.ind{j}(i,:) = pairs_lag_XXX{k}(i,6);
        struct_pair{k}.new_ind{j}(i,:) =j;
        struct_pair{k}.imean{j}(i,:) =i;
    elseif i>1 & (ismember(pairs_lag_XXX{k}(i,1:2),pairs_lag_XXX{k}(i-1,1:2)) == [0 ,1] | ismember(pairs_lag_XXX{k}(i,1:2),pairs_lag_XXX{k}(i-1,1:2)) == [1,0] | ismember(pairs_lag_XXX{k}(i,1:2),pairs_lag_XXX{k}(i-1,1:2)) == [0,0])
        j=j+1;
       struct_pair{k}.pair{j}(i,1:2) = pairs_lag_XXX{k}(i,1:2);
       struct_pair{k}.lag{j}(i,:) = pairs_lag_XXX{k}(i,3);
       struct_pair{k}.pv{j}(i,:) = pairs_lag_XXX{k}(i,4);
       struct_pair{k}.bin{j}(i,:) = pairs_lag_XXX{k}(i,5);
       struct_pair{k}.ind{j}(i,:) = pairs_lag_XXX{k}(i,6);
       struct_pair{k}.new_ind{j}(i,:) =j;
       struct_pair{k}.imean{j}(i,:) =i;
     end
        
   %j=j+1;
end

end

for k=1:size(struct_pair,2)
    for j=1:size(struct_pair{k}.pair,2)
n= 0; 
for i = size(struct_pair{k}.pair{j},1) - n:-1:1
    if struct_pair{k}.pair{j}(i,1)==0
        struct_pair{k}.pair{j}(i,:)=[]; 
        struct_pair{k}.lag{j}(i,:)=[]; 
        struct_pair{k}.pv{j}(i,:)=[];
        struct_pair{k}.bin{j}(i,:)=[];
        struct_pair{k}.ind{j}(i,:)=[];
        struct_pair{k}.new_ind{j}(i,:)=[]; 
        struct_pair{k}.imean{j}(i,:) =[];
        n=n+1;
    end
end
    end
end


end

