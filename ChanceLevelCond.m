function [ChanceLevPairs] = ChanceLevelCond(TLabels, given)
% Chance Level Conditional pairs Probability

switch given

case 'vs|vta'
%%% Chance level in case vs|vta is P(Vs)
%%% I have to put this matrix in the rigth order to be the divisor of the
%%% Cond probability
%  TLabels=Tab_PwN_small_vsvta;
Mat_Pairs = TLabels{:,1:4};
Sum_Pairs = sum(Mat_Pairs,2);
% Sum_Pairs(find(Sum_Pairs==0))= 10^-3;
Sum_Pairs = repmat(Sum_Pairs,1,size(Mat_Pairs,2));
% % % ChanceLevPairs = Mat_Pairs./Sum_Pairs;
Chance = Mat_Pairs./Sum_Pairs;
Chance1 = repmat(Chance(:,1),1,size(Mat_Pairs,2));
Chance2 = repmat(Chance(:,2),1,size(Mat_Pairs,2));
Chance3 = repmat(Chance(:,3),1,size(Mat_Pairs,2));
Chance4 = repmat(Chance(:,4),1,size(Mat_Pairs,2));
ChanceLevPairs = [Chance1, Chance2, Chance3, Chance4];
ChanceLevPairs(find(ChanceLevPairs==0)) = nan;

case 'vta|vs'
%%% Chance level in case vs|vta is P(Vta)
%%% I have to put this matrix in the rigth order to be the divisor of the
%%% Cond probability
Mat_Pairs = TLabels{:,5:8};
Sum_Pairs = sum(Mat_Pairs,2);
% Sum_Pairs(find(Sum_Pairs==0))= 10^-3;
Sum_Pairs = repmat(Sum_Pairs,1,size(Mat_Pairs,2));
% % % ChanceLevPairs = Mat_Pairs./Sum_Pairs;

Chance = Mat_Pairs./Sum_Pairs;
% Chance1 = repmat(Chance(:,1),1,size(Mat_Pairs,2));
% Chance2 = repmat(Chance(:,2),1,size(Mat_Pairs,2));
% Chance3 = repmat(Chance(:,3),1,size(Mat_Pairs,2));
% Chance4 = repmat(Chance(:,4),1,size(Mat_Pairs,2));
% ChanceLevPairs = [Chance1, Chance2, Chance3, Chance4];
ChanceLevPairs = repmat(Chance,1, size(Mat_Pairs,2));
ChanceLevPairs(find(ChanceLevPairs==0)) = nan;
end

end

