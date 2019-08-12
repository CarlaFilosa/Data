function [MeanTotNorm] = NormTotMax(MaxTotAs,MeanTr)
% count=0;
% for k=1:size(MaxA,2);
%     if ~isempty(MaxA{k})
%         count=count+1;
%         MAus(count)=max(MaxA{k});
%     end
% end
% MaxTotAs=max(MAus);

 for k=1:size(MeanTr,2)
    MeanTotNorm{k}=[];
    for i=1:size(MeanTr{k},2)
         if ~isempty(MeanTr{k}{i})
              if MaxTotAs{k}(i)~=0
            MeanTotNorm{k}{i}=MeanTr{k}{i}/MaxTotAs{k}(i);
              elseif MaxTotAs{k}(i)==0
                  MeanTotNorm{k}{i}=MeanTr{k}{i};
         end
    end
    end
 end

end

