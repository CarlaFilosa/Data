function [BinZero,BinNonZero] = BinZeroLag(pairs,c1,c2)

Ind =find(pairs(:,3)==0 & pairs(:,9)==c1 & pairs(:,10)==c2)
 if ~isempty(Ind)
    BinZero = pairs(Ind,5);
%     LagZero = pairs(Ind,3);
 else BinZero=[];
 end
    Ind1 = find(pairs(:,3)~=0 & pairs(:,9)==c1 & pairs(:,10)==c2)
 if ~isempty(Ind1)
        BinNonZero = pairs(Ind1,5);
 else BinNonZero=[];
 end
end