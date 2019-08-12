function [ForMean,Mean_Fir,SEr_Fir] = MeanSingleNeuFir(NF,EventsMtx,H,start)
% Mean of the single firing rate of the neuron
count=0;
ForMean=[];
% GG=NF
% clear NF
% NF=GG{5}{4};
% EventsMtx=new_data.events{5};
% H=1
for l=1:size(NF,2)
    if (EventsMtx(l+start(1)-1).rew_code==H) & (~isempty(NF{l}))
        count=count+1;
        IndOld(count)=l;
        SizeTr(count)=length(NF{l});
        minSz=min(SizeTr(~0));
%     elseif  EventsMtx(l+start-1).rew_code==R
%         countR=countR+1;
%         IndOldR(countR)=l;
%         SizeTrR(countR)=length(NF{l});
%         minSzR=min(SizeTrR);
    end
    
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

