function [LagPerBin,LagPerBinNoZero,LagPerBinZero,LagInSecMat,LagInSecNoZero,LagInSecZero] = LagCountTypeMatPairs(pairs,BinSizes,c1,c2)
%% %% this function find lag in secs, for specific neuron typologies or in general for the all typologies;
%%%%% 

 if nargin==2
  for i=1:length(BinSizes) % without the last two bins
    Ind{i}=find(pairs(:,5)==BinSizes(i));
%     LagPerBin{i} = pairs(Ind{i},3);
    if ~isempty(Ind{i})
    LagPerBin{i} = pairs(Ind{i},3);
    LagInSec{i} = pairs(Ind{i},3).*pairs(Ind{i},5);
    else 
       LagPerBin{i} = [];
       LagInSec{i} = [];
    end
  end

elseif nargin==3
  for i=1:length(BinSizes) % without the last two bins
    Ind{i}=find(pairs(:,5)==BinSizes(i) & pairs(:,9)==c1 & pairs(:,10)==c1);
%     LagPerBin{i} = pairs(Ind{i},3);
    if ~isempty(Ind{i})
    LagPerBin{i} = pairs(Ind{i},3);
    LagInSec{i} = pairs(Ind{i},3).*pairs(Ind{i},5);
    else 
        LagPerBin{i} = [];
       LagInSec{i} = [];
    end
  end
  
elseif nargin==4
%     c2=4; pairs=pairs_small;
    for i=1:length(BinSizes) % without the last two bins
    Ind{i}=find(pairs(:,5)==BinSizes(i) & pairs(:,9)==c1 & pairs(:,10)==c2); 
    Ind1{i}=find(pairs(:,5)==BinSizes(i) & pairs(:,9)==c1 & pairs(:,12)==c2);
    if ~isempty(Ind{i})
    LagPerBin{i} = pairs(Ind{i},3);
    LagInSec{i} = pairs(Ind{i},3).*pairs(Ind{i},5);
    elseif ~isempty(Ind1{i})
    LagPerBin{i} = pairs(Ind1{i},3);
    LagInSec{i} = pairs(Ind1{i},3).*pairs(Ind1{i},5);
    else 
        LagPerBin{i} = [];
       LagInSec{i} = [];
% %     lengthLagPerBin(i) =length(LagPerBin{i});
     end
     end
 end


for i =1:length(LagPerBin)
    if ~isempty(LagPerBin{i})
LagPerBinNoZero{i} = LagPerBin{i}(find(LagPerBin{i}~=0));
LagPerBinZero{i} = LagPerBin{i}(find(LagPerBin{i}==0));

    else
        LagPerBinNoZero{i}=[];
        LagPerBinZero{i} =[];
        
    end
end
%%%%%% Lag in seconds 
LagTrue_Sec{1} =[];
kk=0;
if ~isempty(LagInSec)
  for k =1:length(LagInSec)
    if ~isempty(LagInSec{k})
        kk=kk+1;
        LagTrue_Sec{kk}=LagInSec{k};
    end
  end
  
    LagInSecMat = LagTrue_Sec{1};
    for i=1:length(LagTrue_Sec)-1
        LagInSecMat = [LagInSecMat;LagTrue_Sec{i+1}];
    end
%    else LagInSecMat{k}= [];
end

if ~isempty(LagInSecMat)
    for i =1:length(LagInSecMat)
        LagInSecNoZero = LagInSecMat(find(LagInSecMat~=0));
        LagInSecZero = LagInSecMat(find(LagInSecMat==0));
    end
else 
 LagInSecNoZero = [];
 LagInSecZero =[];
end
end

