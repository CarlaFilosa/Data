function [AsForMean,MeanTr,MeanNormTr,Std,MaxA] = Mean_AsSel_AsAct(as_sel)
% % % % % as_sel=as_sel_MRewNew{4}.asAct;
%Matrix to do the mean of the assembly selected and mean over trials
AsForMean=[];
MeanTr =[];
MeanNormTr=[];
Std=[];
MaxA=[];
for k =1:length(as_sel)
    if ~isempty(as_sel{k}) & ~isempty(as_sel{k}.asAct)
for i=1:length(as_sel{k}.asAct)
    nni=0;
    for j=1:length(as_sel{k}.asAct{i})
        if ~isempty(as_sel{k}.asAct{i}{j})
            nni=nni+1;
        Size(nni)=size(as_sel{k}.asAct{i}{j},1);
        for ai=1:Size(nni)
        AsForMean{k}{i}(ai,nni)=as_sel{k}.asAct{i}{j}(ai,2);
        ML(k)=size(AsForMean{k}{i},1);
        end
        end
    end
end
    end
end
MML=max(ML);
for k=1:length(AsForMean)
    if ~isempty(AsForMean{k})
        for i=1:length(AsForMean{k})
    if isnan(AsForMean{k}{i})==ones(size(AsForMean{k}{i}))
        MeanTr{k}{i} = nan(MML,1);
        Std{k}{i} =nan(MML,1);
    else
        MeanTr{k}{i} = mean(AsForMean{k}{i},2,'omitnan');
        Std{k}{i} = std(AsForMean{k}{i},1,2,'omitnan')/sqrt(size(AsForMean{k}{i},2));
    end 
        end
    end
end
for k=1:length(MeanTr) 
    if~isempty(MeanTr{k})
for i=1:length(MeanTr{k})
   if ~isempty(MeanTr{k}{i})
       MaxA{k}(i)=max(MeanTr{k}{i});
       if MaxA{k}(i)==0
        MeanNormTr{k}{i}=MeanTr{k}{i};
       else
        MeanNormTr{k}{i}=MeanTr{k}{i}/MaxA{k}(i);
      end
   end
end
    end
end

end

