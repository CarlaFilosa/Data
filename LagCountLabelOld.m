function [LagCountVeryTrue_Sec] = LagCountLabel(struct_pair,c1,c2)
% Labels for the pairs in such a way to took only the couple in wich I am
% interested
% struct_pair = struct_pair_small;
% c1=1;
% c2=4;
%%%%% in my case Data is Data_All.par

if nargin <3
    for k =1 :size(struct_pair,2)
    LagCount_Sec{k} = [];
    jj=0;
    lag{k} = cell2mat(struct_pair{k}.lag);
    bin{k} = cell2mat(struct_pair{k}.bin);
     if ~isempty(struct_pair{k})
     for i = 1:length(struct_pair{k}.lag)
         if struct_pair{k}.labels(i,1) ==c1 && struct_pair{k}.labels(i,2) ==c2 
             jj=jj+1;
             LagCount_Sec{k}(jj,1)= lag{k}(i)*bin{k}(i);
           end
     end
    end
end
LagCountTrue_Sec{1} =[];
kk=0;
if ~isempty(LagCount_Sec)
  for k =1:length(LagCount_Sec)
    if ~isempty(LagCount_Sec{k})
        kk=kk+1;
        LagCountTrue_Sec{kk}=LagCount_Sec{k};
    end
  end
  
    LagCountVeryTrue_Sec = LagCountTrue_Sec{1};
    for i=1:length(LagCountTrue_Sec)-1
        LagCountVeryTrue_Sec = [LagCountVeryTrue_Sec;LagCountTrue_Sec{i+1}];
    end
    else LagCountVeryTrue_Sec{k}= LagCount_Sec{k};
end

elseif nargin==3
for k =1 :size(struct_pair,2)
    LagCount_Sec{k} = [];
    jj=0;
    lag{k} = cell2mat(struct_pair{k}.lag);
    bin{k} = cell2mat(struct_pair{k}.bin);
     if ~isempty(struct_pair{k})
     for i = 1:length(struct_pair{k}.lag)
         if struct_pair{k}.labels(i,1) ==c1 && struct_pair{k}.labels(i,2) ==c2 
             jj=jj+1;
             LagCount_Sec{k}(jj,1)= lag{k}(i)*bin{k}(i);
           end
     end
    end
end
LagCountTrue_Sec{1} =[];
kk=0;
if ~isempty(LagCount_Sec)
  for k =1:length(LagCount_Sec)
    if ~isempty(LagCount_Sec{k})
        kk=kk+1;
        LagCountTrue_Sec{kk}=LagCount_Sec{k};
    end
  end
  
    LagCountVeryTrue_Sec = LagCountTrue_Sec{1};
    for i=1:length(LagCountTrue_Sec)-1
        LagCountVeryTrue_Sec = [LagCountVeryTrue_Sec;LagCountTrue_Sec{i+1}];
    end
    else LagCountVeryTrue_Sec{k}= LagCount_Sec{k};
end

end


end

