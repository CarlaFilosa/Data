function [new_struct_pair] = DirAndBin(Dir, struct_pair,k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
count=0;
bb=[];
% cluster_by_lag=[];
% struct_pair=copy_struct_pair;
% k=1;
% Dir='vs->vta Dir';
% new_struct_pair=[];
BinDir=[];
new_struct_pair = struct ('pair',[], 'lag',[], 'pv', [],'bin',[],'ind',[],'new_ind',[],'neuID',[],'labels',[]);

for i=1:length(struct_pair{k}.lag)
switch Dir 
    
    case '1->2 small'
        
         minbin=0.01; maxbin=0.25;
        
        cluster_by_lag{i}=find(struct_pair{k}.lag{i}>0);
        
        
    case '1->2 large'
        
         minbin=0.35; maxbin=1.6;
        
        cluster_by_lag{i}=find(struct_pair{k}.lag{i}>0);
        
    case '1->2 AllBins'
        
         minbin=0.01; maxbin=1.6;
        
        cluster_by_lag{i}=find(struct_pair{k}.lag{i}>0);
     case '2->1 large'
        
       minbin=0.35; maxbin=1.6;
        
        cluster_by_lag{i}=find(struct_pair{k}.lag{i}<0);
        
     case '2->1 small'
        
       minbin=0.01; maxbin=0.25;
        
        cluster_by_lag{i}=find(struct_pair{k}.lag{i}<0);
        
        case '2->1 AllBins'
        
       minbin=0.01; maxbin=1.6;
        
        cluster_by_lag{i}=find(struct_pair{k}.lag{i}<0);
        
    case 'sync small'
        
        minbin=0.01; maxbin=0.25;
        
        cluster_by_lag{i}=find(struct_pair{k}.lag{i}==0);
        
   case 'sync large'
        
        minbin=0.35; maxbin=1.6;
        
        cluster_by_lag{i}=find(struct_pair{k}.lag{i}==0);
        
       otherwise
        
        minbin=0.01; maxbin=1.6;
        
        cluster_by_lag{i}=find(struct_pair{k}.lag{i});
  end

% if I want to extend this function to the case of not pruned algorithm I
% have to add a for loop on the size of struct_pair{k}.bin{i}
BinDir{i}=[];
for j=1:length(struct_pair{k}.bin{i})
    
            if struct_pair{k}.bin{i}(j,1)>=minbin && struct_pair{k}.bin{i}(j,1)<=maxbin
                
        BinDir{i}(j,1)=struct_pair{k}.bin{i}(j,1);
            end
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
end

