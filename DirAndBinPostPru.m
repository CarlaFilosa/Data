function [new_struct_pair] = DirAndBinPostPru(Dir, struct_pair,k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
count=0;
bb=[];
% struct_pair=copy_struct_pair;
% %  for k=1:TotAn;
% % Dir='vta->vs';
% new_struct_pair=[];
BinDir=[];
 new_struct_pair = struct ('pair',[], 'lag',[], 'pv', [],'bin',[],'ind',[],'new_ind',[],'neuID',[],'labels',[]);
%  struct_pair{k}=copy_struct_pair_small{k};
for i=1:length(struct_pair{k}.lag)
switch Dir 
    
     case 'vs->vta'
        
        cluster_by_lag{i}=find(struct_pair{k}.lag{i}>0);
        
     case 'vta->vs'
        
        cluster_by_lag{i}=find(struct_pair{k}.lag{i}<0);
       
     case 'sync'
        
        cluster_by_lag{i}=find(struct_pair{k}.lag{i}==0);
        
    case 'no'
        cluster_by_lag{i}= i;%find(struct_pair{k}.lag{i}==0) ;
end

% if I want to extend this function to the case of not pruned algorithm I
% have to add a for loop on the size of struct_pair{k}.bin{i}
BinDir{i}=[];
for j=1:length(struct_pair{k}.bin{i})
    
%             if struct_pair{k}.bin{i}(j,1)>=minbin && struct_pair{k}.bin{i}(j,1)<=maxbin
                
        BinDir{i}(j,1)=struct_pair{k}.bin{i}(j,1);
%             end
end



bb=length(cluster_by_lag{i});
if bb~=0
    if ~isempty(BinDir{i})
    count=count+1;
   
    new_struct_pair.pair{count}=struct_pair{k}.pair{i};      
    new_struct_pair.lag{count}=struct_pair{k}.lag{i};
    new_struct_pair.pv{count}=struct_pair{k}.pv{i};
    new_struct_pair.bin{count}=struct_pair{k}.bin{i};
    new_struct_pair.ind{count}=struct_pair{k}.ind{i};
    new_struct_pair.new_ind{count}=count;
%     new_struct_pair.neuID{count}=struct_pair{k}.neuID{i};
%     new_struct_pair.labels{count}=struct_pair{k}.labels{i};
    
    new_struct_pair.neuID(count,:)=struct_pair{k}.neuID(i,:);
    new_struct_pair.labels(count,:)=struct_pair{k}.labels(i,:);
%     new_struct_pair.label_n2{count}=struct_pair{k}.label_n2{i};
    
    
    end
end
end
% %  end
end

