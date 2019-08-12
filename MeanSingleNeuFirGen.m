function [ForMean,Mean_Fir,SEr_Fir] = MeanSingleNeuFirGen(NF)
% Mean of the single firing rate of the neuron
% input is the nueuron firing rate
count=0;
ForMean=[];

for l=1:size(NF,2)
    if EventsMtx(l+start-1).rew_code==H
        count=count+1;
        IndOld(count)=l;
        SizeTr(count)=length(NF{l});
        minSz=min(SizeTr);
        
    
end
 

for l=1:count
    for gg=1:minSz
        ForMean(l,gg)=NF{IndOld(l)}(1,gg);
    end
    
end

% for l=1:countR
%     for gg=1:min(SizeTrR)
%         RejForMean(countR,gg)=NF{IndOldH(countR)}(1,gg);
%     end
%     
% end

Mean_Fir = mean(ForMean,1,'omitnan');
SEr_Fir= std(ForMean,0,1,'omitnan')/sqrt(size(ForMean,1));


end

