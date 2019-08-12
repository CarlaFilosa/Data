function [struct_pairN] = PairsOrgByPair(pairs_g,NP)
% The function builds a structure that organizes the pairs so that they are
% clustered by binsizes

for k=1:size(pairs_g,2)
pairs_lag_XXX{k}=[];
% NP_XXX{k}=[];
end

 for k=1:size(pairs_r,2)
% k=4;
% pairs_lag_XXX{k}=pairs_r{k};
pairs_lag_XXX{k}=g{k};
% NP_XXX{k}=NP{k};

% struct for pair to overlap the bin in the figure
[ff,pp]=sort(pairs_g{k}(:,1));
sort_g{k}=pairs_g{k}(pp,:);
pairs_lag_XXX{k}=sort_g{k};
struct_pair{k} = struct ('pair',[], 'lag',[], 'pv', [],'bin',[],'ind',[],'new_ind',[],'imean',[]);
struct_pairN{k} = struct ('pair',[], 'lag',[], 'pv', [],'bin',[],'ind',[],'new_ind',[],'imean',[]);
j=1; % j is uncremented in the second elseif
   
for i = 1: size(pairs_lag_XXX{k},1)
     if i ==1 | size(pairs_lag_XXX{k},1)==1
        struct_pair{k}.pair{j}(i,1:2)= pairs_lag_XXX{k}(i,1:2);
        struct_pair{k}.lag{j}(i,:) = pairs_lag_XXX{k}(i,3);
        struct_pair{k}.pv{j}(i,:) = pairs_lag_XXX{k}(i,4);
        struct_pair{k}.bin{j}(i,:) = pairs_lag_XXX{k}(i,5);
        struct_pair{k}.ind{j}(i,:) = pairs_lag_XXX{k}(i,6);
        struct_pair{k}.new_ind{j}(i,:) =j;
        struct_pair{k}.imean{j}(i,:) =i;
     elseif i~=1
    
    if ismember(pairs_lag_XXX{k}(i,1),pairs_lag_XXX{k}(i-1,1))==1
        struct_pair{k}.pair{j}(i,1:2)= pairs_lag_XXX{k}(i,1:2);
        struct_pair{k}.lag{j}(i,:) = pairs_lag_XXX{k}(i,3);
        struct_pair{k}.pv{j}(i,:) = pairs_lag_XXX{k}(i,4);
        struct_pair{k}.bin{j}(i,:) = pairs_lag_XXX{k}(i,5);
        struct_pair{k}.ind{j}(i,:) = pairs_lag_XXX{k}(i,6);
        struct_pair{k}.new_ind{j}(i,:) =j;
        struct_pair{k}.imean{j}(i,:) =i;
    elseif ismember(pairs_lag_XXX{k}(i,1),pairs_lag_XXX{k}(i-1,1))==0
        j=j+1;
        struct_pair{k}.pair{j}(i,1:2)= pairs_lag_XXX{k}(i,1:2);
        struct_pair{k}.lag{j}(i,:) = pairs_lag_XXX{k}(i,3);
        struct_pair{k}.pv{j}(i,:) = pairs_lag_XXX{k}(i,4);
        struct_pair{k}.bin{j}(i,:) = pairs_lag_XXX{k}(i,5);
        struct_pair{k}.ind{j}(i,:) = pairs_lag_XXX{k}(i,6);
        struct_pair{k}.new_ind{j}(i,:) =j;
        struct_pair{k}.imean{j}(i,:) =i;
       
    end
   
     end
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

 
 for k=1:size(struct_pair,2)
 for l=1:length(struct_pair{k}.pair)
    [st,i_st]=sort(struct_pair{k}.pair{l}(:,2));
    struct_pair{k}.pair{l}=struct_pair{k}.pair{l}(i_st,:);
    struct_pair{k}.lag{l}=struct_pair{k}.lag{l}(i_st,:);
    struct_pair{k}.pv{l}=struct_pair{k}.pv{l}(i_st,:);
    struct_pair{k}.bin{l}=struct_pair{k}.bin{l}(i_st,:);
    struct_pair{k}.ind{l}=struct_pair{k}.ind{l}(i_st,:);
 end
 end
 for k=1:size(struct_pair,2)
 j=0; % j is uncremented in the second elseif
 for l=1:length(struct_pair{k}.pair)


for i = 1: size(struct_pair{k}.pair{l},1)
    
    if i ==1 | size(struct_pair{k}.pair{l},1)==1
        struct_pairN{k}.pair{j+l}(i,1:2)= struct_pair{k}.pair{l}(i,1:2);
        struct_pairN{k}.lag{j+l}(i,:) = struct_pair{k}.lag{l}(i,:);
        struct_pairN{k}.pv{j+l}(i,:) = struct_pair{k}.pv{l}(i,:);
        struct_pairN{k}.bin{j+l}(i,:) = struct_pair{k}.bin{l}(i,:);
        struct_pairN{k}.ind{j+l}(i,:) = struct_pair{k}.ind{l}(i,:);
        struct_pairN{k}.new_ind{j+l}(i,:) =j;
        struct_pairN{k}.imean{j+l}(i,:) =i;
    elseif i~=1 & size(struct_pair{k}.pair{l},1)~=1
        
    if ismember(struct_pair{k}.pair{l}(i,2),struct_pair{k}.pair{l}(i-1,2))==1
        
         struct_pairN{k}.pair{j+l}(i,1:2)= struct_pair{k}.pair{l}(i,1:2);
        struct_pairN{k}.lag{j+l}(i,:) = struct_pair{k}.lag{l}(i,:);
        struct_pairN{k}.pv{j+l}(i,:) = struct_pair{k}.pv{l}(i,:);
        struct_pairN{k}.bin{j+l}(i,:) = struct_pair{k}.bin{l}(i,:);
        struct_pairN{k}.ind{j+l}(i,:) = struct_pair{k}.ind{l}(i,:);
        struct_pairN{k}.new_ind{j+l}(i,:) =j;
        struct_pairN{k}.imean{j+l}(i,:) =i;
    elseif ismember(struct_pair{k}.pair{l}(i,2),struct_pair{k}.pair{l}(i-1,2))==0
        
        j=j+1;
         struct_pairN{k}.pair{j+l}(i,1:2)= struct_pair{k}.pair{l}(i,1:2);
         struct_pairN{k}.lag{j+l}(i,:) = struct_pair{k}.lag{l}(i,:);
         struct_pairN{k}.pv{j+l}(i,:) = struct_pair{k}.pv{l}(i,:);
         struct_pairN{k}.bin{j+l}(i,:) = struct_pair{k}.bin{l}(i,:);
         struct_pairN{k}.ind{j+l}(i,:) = struct_pair{k}.ind{l}(i,:);
         struct_pairN{k}.new_ind{j+l}(i,:) =j;
         struct_pairN{k}.imean{j+l}(i,:) =i;

       
    end
    end

 end
 end
 end
 
 for k=1:size(struct_pair,2)
for j=1:size(struct_pairN{k}.pair,2)
n= 0; 
 for i = size(struct_pairN{k}.pair{j},1) - n:-1:1
    if struct_pairN{k}.pair{j}(i,1)==0
        struct_pairN{k}.pair{j}(i,:)=[]; 
        struct_pairN{k}.lag{j}(i,:)=[]; 
        struct_pairN{k}.pv{j}(i,:)=[];
        struct_pairN{k}.bin{j}(i,:)=[];
        struct_pairN{k}.ind{j}(i,:)=[];
        struct_pairN{k}.new_ind{j}(i,:)=[]; 
        struct_pairN{k}.imean{j}(i,:) =[];
        n=n+1;
    end
 end
end
 end
end

