function [toArrayLagType] = LagCountPerBinType(struct_pair,BinSizes,c1,c2)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
for k= 1: length(struct_pair)

for i =1:length(BinSizes)
    LagPerBin{i}{k} =[];
end
end

%%%%%%% Specificy the neuron of the pair (only one side)

if nargin<4

for k= 1: length(struct_pair)

for i =1:length(BinSizes)
jj=0;
for j = 1:length(struct_pair{k}.lag)
         toArrayLag{k} =cell2mat(struct_pair{k}.lag);
         toArrayBin{k} =cell2mat(struct_pair{k}.bin);
         if struct_pair{k}.labels(j,1) ==c1  || struct_pair{k}.labels(j,2) ==c1
jj=jj+1;
         
         if  toArrayBin{k}(j)==BinSizes(i)
          LagPerBin{i}{k}= toArrayLag{k}(j);
         end
end
end
end
end

for i =1: length(LagPerBin)
 toArrayLagType{i} =cell2mat(LagPerBin{i});
end
%%%%%%% Specify the couple 
elseif nargin ==4
for k= 1: length(struct_pair)

for i =1:length(BinSizes)
jj=0;
for j = 1:length(struct_pair{k}.lag)
         toArrayLag{k} =cell2mat(struct_pair{k}.lag);
         toArrayBin{k} =cell2mat(struct_pair{k}.bin);
         if struct_pair{k}.labels(j,1) ==c1  && struct_pair{k}.labels(j,2) ==c2
jj=jj+1;
         
         if  toArrayBin{k}(j)==BinSizes(i)
          LagPerBin{i}{k}= toArrayLag{k}(j);
         end
end
end
end
end

for k =1: length(LagPerBin)
 toArrayLagType{k} =cell2mat(LagPerBin{k});
end
end
    
end

