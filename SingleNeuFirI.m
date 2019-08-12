function [NeuFirI,BinTrials,NF ] = SingleNeuFirI(new_spM,FindInd,start,minInt)
% Firing rate of single neuron in the phase that i desire

TotAn = size(new_spM.events,2);
 for k= 1:TotAn
    for l = 1:size(FindInd{k},2)
        for i= 1:size(new_spM.spikeT_BegEnd{k},1)
%        Dist_fvP(k,l)=new_spM.events{k}(l+1).fv_on-new_spM.events{k}(l).fv_on;
%        Dist_fv(k,l)=new_spM.events{k}(l).fv_off-new_spM.events{k}(l).fv_on;
         NeuFirI{k}{i}{l}=new_spM.spikeT_BegEnd{k}(i,find(new_spM.spikeT_BegEnd{k}(i,:)>=new_spM.events{k}(l+start(k)-1).fv_on & new_spM.spikeT_BegEnd{k}(i,:)<=new_spM.events{k}(l+start(k)).fv_on));
%        NeuFirI{k}{i}{l}=new_spM.spikeT_BegEnd{k}(i,FindInd{k}{l});
         NeuFirI{k}{i}{l}(2,:)=NeuFirI{k}{i}{l}-new_spM.events{k}(l+start(k)-1).fv_on;
         end
          BinTrials{k}{l} =new_spM.events{k}(l+start(k)-1).fv_on:minInt:new_spM.events{k}(l+start(k)).fv_on;
    end
 end

 for k=1:TotAn
    for i =1: size(new_spM.spiketrain{k},1) 
    for l=1:size(FindInd{k},2)
     
    if size(NeuFirI{k}{i}{l})~=0
NF{k}{i}{l}=histcounts(NeuFirI{k}{i}{l},BinTrials{k}{l});
    end
    end
end
 end

end


% % % minInt=0.1;
% % % for k= 1:size(new_spM.events,2)
% % %     for l = 1: length(new_spM.events{k})-1
% % %          for i= 1:size(new_spM.spikeT_BegEnd{k},1)
% % %        Dist_fvP(k,l)=new_spM.events{k}(l+1).fv_on-new_spM.events{k}(l).fv_on;
% % %        Dist_fv(k,l)=new_spM.events{k}(l).fv_off-new_spM.events{k}(l).fv_on;
% % %        NeuFirI{k}{i}{l}=new_spM.spikeT_BegEnd{k}(i,find(new_spM.spikeT_BegEnd{k}(i,:)>=new_spM.events{k}(l).fv_on & new_spM.spikeT_BegEnd{k}(i,:)<=new_spM.events{k}(l).fv_on+Dist_fvP(k,l)));
% % %        NeuFirI{k}{i}{l}(2,:)=NeuFirI{k}{i}{l}-new_spM.events{k}(l).fv_on;
% % % % a{k}{i}{l}=find(new_spM.spikeT_BegEnd{k}(i,:)>=new_spM.events{k}(l).fv_on);
% % %          end
% % %          BinTrials{k}{l} =new_spM.events{k}(l).fv_onTr:minInt:new_spM.events{k}(l+1).fv_onTr+Dist_fvP(k,l);
% % %     end
% % % end
% % % 
% % % 
% % % clear k i l
