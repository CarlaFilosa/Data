function [TabNum,NumPairs] = PairsWSpecNeu_TabNum(TLabels,Region)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

switch Region
case 'vta'
NumPairs = repmat(TLabels{:,5:8},1,4);

TypeI = NumPairs(:,1);
TypeII = NumPairs(:,2);
TypeIII = NumPairs(:,3);
NoIdVta =  NumPairs(:,4);
% time=0:3;
TabNum = table(TypeI,TypeII,TypeIII,NoIdVta,TypeI,TypeII,TypeIII,NoIdVta,TypeI,TypeII,TypeIII,NoIdVta,TypeI,TypeII,TypeIII,NoIdVta);

case 'vs'
NumPairs1 = repmat(TLabels{:,1},1,4);
NumPairs2 = repmat(TLabels{:,2},1,4);
NumPairs3 = repmat(TLabels{:,3},1,4);
NumPairs4 = repmat(TLabels{:,4},1,4);
NumPairs=[NumPairs1, NumPairs2, NumPairs3, NumPairs4];
MSN = NumPairs(:,1);
FSI = NumPairs(:,5);
CIN = NumPairs(:,9);
NoIdVs =  NumPairs(:,13);
TabNum = table( MSN, MSN, MSN, MSN, FSI, FSI,FSI, FSI, CIN, CIN, CIN, CIN, NoIdVs,NoIdVs, NoIdVs, NoIdVs);
end

end

