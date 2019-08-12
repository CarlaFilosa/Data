function [AsForMean] = MeanNewDim(as_sel)
% mean for those don't have the same dimensions
AsForMean=[];
if ~isempty(as_sel)
    for i=1: length(as_sel.asAct)
        nni=0;
        for j=1:length(as_sel.asAct{i})
            if ~isempty(as_sel.asAct{i}{j})
                nni=nni+1
                for ai=1:length(as_sel.asAct{i}{j})
                 AsForMean{i}(ai,nni)=as_sel.asAct{i}{j}(ai,2)
                end
                end
            end
        end        
  
end

% for count= 1:a
% oo=0;
% nni=0;
% nnt=0;
%    for j=1:length(asTime_Ph{k}{aaa})
%        if (start ~=0) &(new_spM.events{k}(start+j-1).rew_code==rew_code) % hit
%        oo=oo+1;
%    as_sel.asAct{count}{oo}=asTime_Ph{k}{aaa}{j}; 
%    as_sel.asId{count}=aaa;
%        if ~isempty(as_sel.asAct{count}{oo})
%            nni=nni+1;
%    Size(nni)=size(as_sel.asAct{count}{oo},1);
% 
%    for ai=1:Size(nni)
%         AsForMean{count}(ai,nni)= as_sel.asAct{count}{oo}(ai,2);
%    end
%        end
%        end
%    end
%  
%     
%  if ~isempty(AsForMean)
% MeanTr{count}=mean(AsForMean{count},2,'omitnan');
% StETr{count}=std(AsForMean{count},1,2,'omitnan')/sqrt(size(AsForMean{count},2));
%  end
%  
% end
% end
% 
% 
% for i=1:length(MeanTr)
%    if ~isempty(MeanTr{i})
%        MaxA(i)=max(MeanTr{i});
%        if MaxA(i)==0
%         MeanNormTr{i}=MeanTr{i};
%        else
%         MeanNormTr{i}=MeanTr{i}/MaxA(i);
%       end
%    end
% end

end

