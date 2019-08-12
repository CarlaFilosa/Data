addpath('/zifnas/Carla/CWEM_Project_GenProg');
k=4;
% i=47;
%  i=47;
% cluster_by_lag=[];
count=0;
for i= 1: size(struct_pair{k}.pair,2)
cluster_by_lag{i}=find(struct_pair{k}.lag{i}~=0);
a=size(cluster_by_lag{i},1);
if a~=0
    count=count+1;
    figure(count);
    StrP =strcat('Pair: ',mat2str(struct_pair{k}.pair{i}(end,:)))
    StrB = strcat('Bin: ',mat2str(struct_pair{k}.bin{i}(1:end)))
    firstU=struct_pair{k}.pair{i}(end,1);
    secondU=struct_pair{k}.pair{i}(end,2);
    NeuFirst{l}=histcounts(NeuFirI{k}{firstU}{l},BinTrials{k}{l});
    NeuSecond{l}=histcounts(NeuFirI{k}{secondU}{l},BinTrials{k}{l});
for iii= 1:a
fff{i}=struct_pair{k}.ind{i};
aaa = fff{i}(iii);
in=0;
oo=0;
   for j=1:size(asRes{k}{aaa},2); % j index runs over number of trials
       
       if (new_spM.events{k}(j).rew_code==1) % hit
       

in=in+2;
oo=oo+1;
countH=countH+1;
maxA=max(asRes{k}{fff{i}(iii)}{j}(:,2));
 if maxA==0 
     maxA=1; end


 maxT(j)=max(asRes{k}{aaa}{j}(:,3));

 for gfr=1:length(NF{l})   % Hit for firing rate
        HitForMean(countH,gfr)=NF{l}(1,gfr);
 end
        SZH(countH)= length(NF{l});
       
 
 for ggg=1:size(asRes{k}{aaa}{j},1)  % hit for assembly
   SzA{count}(ggg,1) = size(asRes{k}{aaa}{j}(:,2),1);
   AsForMeanA{count}(ggg,oo)=asRes{k}{aaa}{j}(ggg,2)/maxA;
   AsTimeA{count}(ggg,oo)=asRes{k}{aaa}{j}{ggg,3};
 end

end
    end
 maxTtot=max(maxT);
MeanTrA{count}=mean(AsForMeanA{count}(1:min(SzA{count}),:),2);
StETrA{count}=std(AsForMeanA{count}(1:min(SzA{count}),:),1,2)/sqrt(size(AsForMeanA{count},2));



in=0;
oo=0;
   for j=1:size(asRes{k}{aaa},2);
if (new_spM.events{k}(j).rew_code==3) % correct rejection

     in=in+2;
     oo=oo+1;
     countR=countR+1;
     maxA=max(asRes{k}{fff{i}(iii)}{j}(:,2));
 if maxA==0 
     maxA=1; end

    maxT1(j)=max(asRes{k}{aaa}{j}(:,3));

     
for ggg=1:length(NF{l})  % rej for firing neu
        RejForMean(countR,ggg)=NF{l}(1,ggg);
end
        SZR(countR)= length(NF{l});
    
    
for ggg=1:size(asRes{k}{aaa}{j},1) % rej for assemblies
   SzB{count}(ggg,1) = size(asRes{k}{aaa}{j}(:,2),1);
   AsForMeanB{count}(ggg,oo)=asRes{k}{aaa}{j}(ggg,2)/maxA;
   AsTimeB{count}(ggg,oo)=asRes{k}{aaa}{j}(ggg,3);
end
end
   end
 maxTtot1=max(maxT1);
 MeanTrB{count}=mean(AsForMeanB{count}(1:min(SzB{count}),:),2);
 StETrB{count}=std(AsForMeanB{count}(1:min(SzB{count}),:),1,2)/sqrt(size(AsForMeanB{count},2));
end

minY = min([min(MeanTrA{count}-StETrA{count}),min(MeanTrB{count}-StETrB{count})])-0.01;
maxY = max([max(MeanTrA{count}+StETrA{count}),max(MeanTrB{count}+StETrB{count})])+0.01;
maxX = max([max(asRes{k}{1}{1}(1:min(SzA{count}),3)),max(asRes{k}{1}{1}(1:min(SzB{count}),3))]);

end
end
% pause

%% Single neuron Firing rate %%%%%%%%%%%%%% 
minIntN=0.1;
TotAn=size(new_spM.spiketrain,2);
 [NeuFirI_S,BinTrials_SO,NF_S ] = SingleNeuFirI(new_spM,FindInd_S,startS,minIntN);
 [NeuFirI_U,BinTrials_UO,NF_U ] = SingleNeuFirI(new_spM,FindInd_U,startU,minIntN);
for k=1:TotAn
    for i= 1:size(new_spM.spiketrain{k},1)
       
[ForMean_Hit_S{k}{i},Mean_Fir_Hit_S{k}{i},SEr_Fir_Hit_S{k}{i}] = MeanSingleNeuFir(NF_S{k}{i},new_spM.events{k},hit,startS(k));
[ForMean_Rej_S{k}{i},Mean_Fir_Rej_S{k}{i},SEr_Fir_Rej_S{k}{i}] = MeanSingleNeuFir(NF_S{k}{i},new_spM.events{k},rej,startS(k));
[ForMean_Miss_S{k}{i},Mean_Fir_Miss_S{k}{i},SEr_Fir_Miss_S{k}{i}] = MeanSingleNeuFir(NF_S{k}{i},new_spM.events{k},miss,startS(k));
[ForMean_Fal_S{k}{i},Mean_Fir_Fal_S{k}{i},SEr_Fir_Fal_S{k}{i}] = MeanSingleNeuFir(NF_S{k}{i},new_spM.events{k},fal,startS(k));
[ForMean_Hit_U{k}{i},Mean_Fir_Hit_U{k}{i},SEr_Fir_Hit_U{k}{i}] = MeanSingleNeuFir(NF_U{k}{i},new_spM.events{k},hit,startU(k));
[ForMean_Rej_U{k}{i},Mean_Fir_Rej_U{k}{i},SEr_Fir_Rej_U{k}{i}] = MeanSingleNeuFir(NF_U{k}{i},new_spM.events{k},rej,startU(k));
[ForMean_Miss_U{k}{i},Mean_Fir_Miss_U{k}{i},SEr_Fir_Miss_U{k}{i}] = MeanSingleNeuFir(NF_U{k}{i},new_spM.events{k},miss,startU(k));
[ForMean_Fal_U{k}{i},Mean_Fir_Fal_U{k}{i},SEr_Fir_Fal_U{k}{i}] = MeanSingleNeuFir(NF_U{k}{i},new_spM.events{k},fal,startU(k));
    end
end

%% %%%%%%%% Plot Comparison between stable and unstable


addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla');
% TotAn=size(MeanTrH_S,2); % total animals number
cf=0;
lag='pos';
StrS= strcat('Session: ',mat2str(new_spM.params(1).session_tag),'---','Alg: PruOnePeak');
lick_d=new_spM.events{1}(1).lick_delay/1000; %lick delay (sec)
rew_d=new_spM.events{1}(1).reward_delay/1000; %rew delay (sec)
lick_wind=new_spM.info(1).superflex_parameters.lick_window/1000; %lick window (sec)

cf=5;
% for k=1:TotAn
    k=4;
         nas=size(MeanTrH_S{k},2); 
     if nas~=0
        fc=0;
%         vectN1HS=[];
%         vectN2HS=[];
%         vectN1RS=[];
%         vectN2RS=[];
%         vectN1HU=[];
%         vectN2HU=[];
%         vectN1RU=[];
%         vectN2RU=[];
%     for i=1:nas
        i=9;
        % Single neurons limits and Identity
        
        StrNeu = strcat('Neu1:  ',mat2str(new_struct_pair{k}.pair{i}(end,1)),'---',...
                        'Neu2:  ',mat2str(new_struct_pair{k}.pair{i}(end,2)));
        StrNeu1 = strcat('Neu1: -- ',mat2str(new_struct_pair{k}.pair{i}(end,1)),'--');
        StrNeu2 = strcat('Neu2: -- ',mat2str(new_struct_pair{k}.pair{i}(end,2)),'--');
%         
        neu1{k}{i}=new_struct_pair{k}.pair{i}(end,1);
        neu2{k}{i}=new_struct_pair{k}.pair{i}(end,2);
        vectN1HS=0:minIntN:(size(Mean_Fir_Hit_S{k}{neu1{k}{i}},2)-1)*minIntN;
        vectN2HS=0:minIntN:(size(Mean_Fir_Hit_S{k}{neu2{k}{i}},2)-1)*minIntN;
        vectN1RS=0:minIntN:(size(Mean_Fir_Rej_S{k}{neu1{k}{i}},2)-1)*minIntN;
        vectN2RS=0:minIntN:(size(Mean_Fir_Rej_S{k}{neu2{k}{i}},2)-1)*minIntN;
        vectN1HU=0:minIntN:(size(Mean_Fir_Hit_U{k}{neu1{k}{i}},2)-1)*minIntN;
        vectN2HU=0:minIntN:(size(Mean_Fir_Hit_U{k}{neu2{k}{i}},2)-1)*minIntN;
        vectN1RU=0:minIntN:(size(Mean_Fir_Rej_U{k}{neu1{k}{i}},2)-1)*minIntN;
        vectN2RU=0:minIntN:(size(Mean_Fir_Rej_U{k}{neu2{k}{i}},2)-1)*minIntN;
        
        if(~isempty(vectN1HS) && ~isempty(vectN2HS) && ~isempty(vectN1RS) && ~isempty(vectN2RS)...
                && ~isempty(vectN1HU) && ~isempty(vectN2HU) && ~isempty(vectN1RU) && ~isempty(vectN2RU))
        maxNXS=min([vectN1HS(end),vectN2HS(end),vectN1RS(end),vectN2RS(end)]);
        maxNXU=min([vectN1HU(end),vectN2HU(end),vectN1RU(end),vectN2RU(end)]);
        end
        
        maxNY=max([max(Mean_Fir_Hit_S{k}{neu1{k}{i}}+SEr_Fir_Hit_S{k}{neu1{k}{i}}),...
                   max(Mean_Fir_Hit_S{k}{neu2{k}{i}}+SEr_Fir_Hit_S{k}{neu2{k}{i}}),...
                   max(Mean_Fir_Hit_U{k}{neu1{k}{i}}+SEr_Fir_Hit_U{k}{neu1{k}{i}}),...
                   max(Mean_Fir_Hit_U{k}{neu2{k}{i}}+SEr_Fir_Hit_U{k}{neu2{k}{i}}),...
                   max(Mean_Fir_Rej_S{k}{neu1{k}{i}}+SEr_Fir_Rej_S{k}{neu1{k}{i}}),...
                   max(Mean_Fir_Rej_S{k}{neu2{k}{i}}+SEr_Fir_Rej_S{k}{neu2{k}{i}}),...
                   max(Mean_Fir_Rej_U{k}{neu1{k}{i}}+SEr_Fir_Rej_U{k}{neu1{k}{i}}),...
                   max(Mean_Fir_Rej_U{k}{neu2{k}{i}}+SEr_Fir_Rej_U{k}{neu2{k}{i}})])+0.15;
%                
% %                
        minNY=min([min(Mean_Fir_Hit_S{k}{neu1{k}{i}}-SEr_Fir_Hit_S{k}{neu1{k}{i}}),...
                   min(Mean_Fir_Hit_S{k}{neu2{k}{i}}-SEr_Fir_Hit_S{k}{neu2{k}{i}}),...
                   min(Mean_Fir_Hit_U{k}{neu1{k}{i}}-SEr_Fir_Hit_U{k}{neu1{k}{i}}),...
                   min(Mean_Fir_Hit_U{k}{neu2{k}{i}}-SEr_Fir_Hit_U{k}{neu2{k}{i}}),...
                   min(Mean_Fir_Rej_S{k}{neu1{k}{i}}-SEr_Fir_Rej_S{k}{neu1{k}{i}}),...
                   min(Mean_Fir_Rej_S{k}{neu2{k}{i}}-SEr_Fir_Rej_S{k}{neu2{k}{i}}),...
                   min(Mean_Fir_Rej_U{k}{neu1{k}{i}}-SEr_Fir_Rej_U{k}{neu1{k}{i}}),...
                   min(Mean_Fir_Rej_U{k}{neu2{k}{i}}-SEr_Fir_Rej_U{k}{neu2{k}{i}})])-0.15;
%         
        % Pairs limits and string with identity
        
     StrPBL = strcat('Pair: ',mat2str(new_struct_pair{k}.pair{i}(end,:)),'---',...
                     ' Bin: ',mat2str(new_struct_pair{k}.bin{i}(1:end)),'---',...
                     'Vs->Vta, Lag: ',mat2str(new_struct_pair{k}.lag{i}));
    
    
     vectHS=0:minInt:(size(MeanTrH_S{k}{i},1)-1)*minInt;
     vectRS=0:minInt:(size(MeanTrR_S{k}{i},1)-1)*minInt;
     vectHU=0:minInt:(size(MeanTrH_U{k}{i},1)-1)*minInt;
     vectRU=0:minInt:(size(MeanTrR_U{k}{i},1)-1)*minInt;
     
    
     maxXS=min([vectHS(end),vectRS(end)]);
     maxXU=min([vectHU(end),vectRU(end)]);
     maxY=max([max(MeanTrR_U{k}{i}+StETrR_U{k}{i}),max(MeanTrR_S{k}{i}+StETrR_S{k}{i}),...
          max(MeanTrH_U{k}{i}+StETrH_U{k}{i}),max(MeanTrH_S{k}{i}+StETrH_S{k}{i})])+0.15;
     minY=min([min(MeanTrR_U{k}{i}-StETrR_U{k}{i}),min(MeanTrR_S{k}{i}-StETrR_S{k}{i}),...
          min(MeanTrH_U{k}{i}-StETrH_U{k}{i}),min(MeanTrH_S{k}{i}-StETrH_S{k}{i})])-0.15;
      
      
      % Effective Plots routine
    if new_struct_pair{k}.lag{i}>0
        mylag='pos';
    elseif new_struct_pair{k}.lag{i}<0
        mylag='neg';
    end
    if ismember(mylag,lag)
        fc=fc+1;
    figure(fc+cf);
    subplot(3,2,1);hold on
    box on;
    [l,p]= boundedline(vectHU,MeanTrH_U{k}{i},StETrH_U{k}{i},'-b',...
        vectRU,MeanTrR_U{k}{i},StETrR_U{k}{i},'-m','alpha'); hold on;
   plot([median(Dist_fv_U(k,1:5)) median(Dist_fv_U(k,1:5))], [minY maxY], 'k-.');hold on;
   text(median(Dist_fv_U(k,1:5)+0.1), (minY+maxY/5),{'EndOd'},'Color','k','Fontsize',10)
%     vline(median(Dist_fv_SO(k,:)),'k-.','EndOd');hold on;

    vline(lick_d,'g.-');hold on;
     vline(rew_d,'r.-','RewD'); hold on;
    vline((lick_d+lick_wind),'g.-');hold on;
    text(lick_d+0.1, (maxY-maxY/8),{'LickWindow'},'Color','g','Fontsize',10)
    xlim([0,maxXU])
%     xlabel('Time (sec)','Fontsize',14)
   pst = ylabel('Pair activity','Fontsize',14);
   gpst=get(pst,'Position');
    ylim([minY,maxY])
    title({'Unstable Phase'; 'Pair'})
   
    text(6, (maxY-maxY/10),{'-'},'Color','b','Fontsize',15)
    text(6.35, (maxY-maxY/10),{'Hit'},'Color','k','Fontsize',12)
    text(6, (maxY-1.8*maxY/10),{'-'},'Color','m','Fontsize',15)
    text(6.35, (maxY-1.8*maxY/10),{'Rejection'},'Color','k','Fontsize',12)
    
    
    
    subplot(3,2,2);hold on;
    box on;
    [l,p]= boundedline(vectHS,MeanTrH_S{k}{i},StETrH_S{k}{i},'-b',...
        vectRS,MeanTrR_S{k}{i},StETrR_S{k}{i},'-m','alpha');hold on;
     plot([median(Dist_fv_S(k,:)) median(Dist_fv_S(k,:))], [minY maxY], 'k-.'); hold on;
     text(median(Dist_fv_S(k,:)+0.1), (minY+maxY/5),{'EndOd'},'Color','k','Fontsize',10)
%     vline(median(Dist_fv_SO(k,:)),'k-.','EndOd')
    vline(lick_d,'g.-');hold on;
    vline(rew_d,'r.-','RewD');hold on;
    vline(lick_d+lick_wind,'g.-');hold on;
    text(lick_d+0.1, (maxY-maxY/8),{'LickWindow'},'Color','g','Fontsize',10)
    ylim([minY,maxY])
    xlim([0,maxXS])
%     xlabel('Time (sec)','Fontsize',14)
%     ylabel('Mean Pair activity','Fontsize',14)
    title({'Stable Phase';'Pair'})
    text(maxXS+1.1,(maxY-maxY/15),{StrS, StrPBL},'rotation',270,'Fontsize',12)
    text(6, (maxY-maxY/10),{'-'},'Color','b','Fontsize',15)
    text(6.35, (maxY-maxY/10),{'Hit'},'Color','k','Fontsize',12)
    text(6, (maxY-1.8*maxY/10),{'-'},'Color','m','Fontsize',15)
    text(6.35, (maxY-1.8*maxY/10),{'Rejection'},'Color','k','Fontsize',12)
%     text((maxXS), (maxY+0.2),{'Session:Lastrev1 Alg:Pru Dir: Vs->Vta'},'Color','k','Fontsize',15)
    
subplot(3,2,3);hold on
    box on;
    [l,p]= boundedline(vectN1HU,Mean_Fir_Hit_U{k}{neu1{k}{i}},SEr_Fir_Hit_U{k}{neu1{k}{i}},'-b',...
                       vectN2HU,Mean_Fir_Hit_U{k}{neu2{k}{i}},SEr_Fir_Hit_U{k}{neu2{k}{i}},'-c','alpha');hold on;
% %                vectN1RS,Mean_Fir_Rej_SO{k}{neu1{k}{i}},SEr_Fir_Rej_SO{k}{neu1{k}{i}},'-m',...
% %                vectN2RS,Mean_Fir_Rej_SO{k}{neu2{k}{i}},SEr_Fir_Rej_SO{k}{neu2{k}{i}},'-r','alpha'); hold on;
%    
% %%%%%%%% lines plot of interested episodes
% 
    plot([median(Dist_fv_U(k,1:5)) median(Dist_fv_U(k,1:5))], [minNY maxNY], 'k-.');hold on;
    text(median(Dist_fv_U(k,1:5)+0.1), (minNY+maxNY/5),{'EndOd'},'Color','k','Fontsize',10)
%     vline(median(Dist_fv_SO(k,:)),'k-.','EndOd');hold on;
    vline(lick_d,'g.-');hold on;
    vline((lick_d+lick_wind),'g.-');hold on;
    text(lick_d+0.1, (maxNY-maxNY/8),{'LickWindow'},'Color','g','Fontsize',10)
     
%     xlabel('Time (sec)','Fontsize',14)
%     ylabel('Mean Firing rate','Fontsize',14)
    ylim([minNY,maxNY])
    xlim([0,maxNXU])
    title({'Neurons - Hit'})
%    
      text(6, (maxNY-maxNY/10),{'-'},'Color','b','Fontsize',15)
      text((6.35), (maxNY-maxNY/10),{StrNeu1},'Color','k','Fontsize',12)
      text(6, (maxNY-2*maxNY/10),{'-'},'Color','c','Fontsize',15)
      text((6.35), (maxNY-2*maxNY/10),{StrNeu2},'Color','k','Fontsize',12)


      subplot(3,2,4);hold on
    box on;
    [l,p]= boundedline(vectN1HS,Mean_Fir_Hit_S{k}{neu1{k}{i}},SEr_Fir_Hit_S{k}{neu1{k}{i}},'-b',...
                       vectN2HS,Mean_Fir_Hit_S{k}{neu2{k}{i}},SEr_Fir_Hit_S{k}{neu2{k}{i}},'-c','alpha');hold on;
% % %                vectN1RS,Mean_Fir_Rej_SO{k}{neu1{k}{i}},SEr_Fir_Rej_SO{k}{neu1{k}{i}},'-m',...
% % %                vectN2RS,Mean_Fir_Rej_SO{k}{neu2{k}{i}},SEr_Fir_Rej_SO{k}{neu2{k}{i}},'-r','alpha'); hold on;
% %    
% % %%%%%%%% lines plot of interested episodes
% % 
    plot([median(Dist_fv_U(k,1:5)) median(Dist_fv_U(k,1:5))], [minNY maxNY], 'k-.');hold on;
    text(median(Dist_fv_U(k,1:5)+0.1), (minNY+maxNY/5),{'EndOd'},'Color','k','Fontsize',10)
%     vline(median(Dist_fv_SO(k,:)),'k-.','EndOd');hold on;
    vline(lick_d,'g.-');hold on;
    vline((lick_d+lick_wind),'g.-');hold on;
    text(lick_d+0.1, (maxNY-maxNY/8),{'LickWindow'},'Color','g','Fontsize',10)
%      
% %     xlabel('Time (sec)','Fontsize',14)
% %     ylabel('Mean Firing rate','Fontsize',14)
    ylim([minNY,maxNY])
    xlim([0,maxNXS])
    title({'Neurons - Hit'})
   
      text(6, (maxNY-maxNY/10),{'-'},'Color','b','Fontsize',15)
      text((6.35), (maxNY-maxNY/10),{StrNeu1},'Color','k','Fontsize',12)
      text(6, (maxNY-2*maxNY/10),{'-'},'Color','c','Fontsize',15)
      text((6.35), (maxNY-2*maxNY/10),{StrNeu2},'Color','k','Fontsize',12)
% % 
      subplot(3,2,5);hold on
    box on;
    [l,p]= boundedline(vectN1RU,Mean_Fir_Rej_U{k}{neu1{k}{i}},SEr_Fir_Rej_U{k}{neu1{k}{i}},'-m',...
                       vectN2RU,Mean_Fir_Rej_U{k}{neu2{k}{i}},SEr_Fir_Rej_U{k}{neu2{k}{i}},'-r','alpha');hold on;
% %                vectN1RS,Mean_Fir_Rej_SO{k}{neu1{k}{i}},SEr_Fir_Rej_SO{k}{neu1{k}{i}},'-m',...
% %                vectN2RS,Mean_Fir_Rej_SO{k}{neu2{k}{i}},SEr_Fir_Rej_SO{k}{neu2{k}{i}},'-r','alpha'); hold on;
% %    
% % %%%%%%%% lines plot of interested episodes
% % 
    plot([median(Dist_fv_U(k,1:5)) median(Dist_fv_U(k,1:5))], [minNY maxNY], 'k-.');hold on;
    text(median(Dist_fv_U(k,1:5)+0.1), (minNY+maxNY/5),{'EndOd'},'Color','k','Fontsize',10)
%     vline(median(Dist_fv_SO(k,:)),'k-.','EndOd');hold on;
    vline(lick_d,'g.-');hold on;
    vline((lick_d+lick_wind),'g.-');hold on;
    text(lick_d+0.1, (maxNY-maxNY/8),{'LickWindow'},'Color','g','Fontsize',10)
%      
    xlabel('Time (sec)','Fontsize',14)
    ylim([minNY,maxNY])
    xlim([0,maxNXU])
    title({'Neurons - Rejection'})
    text(gpst(1),minNY+1,'Mean Firing rate','Fontsize',14,'Rotation',90)
   
      text(6, (maxNY-maxNY/10),{'-'},'Color','m','Fontsize',15)
      text((6.35), (maxNY-maxNY/10),{StrNeu1},'Color','k','Fontsize',12)
      text(6, (maxNY-2*maxNY/10),{'-'},'Color','r','Fontsize',15)
      text((6.35), (maxNY-2*maxNY/10),{StrNeu2},'Color','k','Fontsize',12)

    subplot(3,2,6);hold on
    box on;
    [l,p]= boundedline(vectN1RS,Mean_Fir_Rej_S{k}{neu1{k}{i}},SEr_Fir_Rej_S{k}{neu1{k}{i}},'-m',...
                       vectN2RS,Mean_Fir_Rej_S{k}{neu2{k}{i}},SEr_Fir_Rej_S{k}{neu2{k}{i}},'-r','alpha');hold on;
% %                vectN1RS,Mean_Fir_Rej_SO{k}{neu1{k}{i}},SEr_Fir_Rej_SO{k}{neu1{k}{i}},'-m',...
% %                vectN2RS,Mean_Fir_Rej_SO{k}{neu2{k}{i}},SEr_Fir_Rej_SO{k}{neu2{k}{i}},'-r','alpha'); hold on;
% %    
% % %%%%%%%% lines plot of interested episodes
% % 
    plot([median(Dist_fv_U(k,1:5)) median(Dist_fv_U(k,1:5))], [minNY maxNY], 'k-.');hold on;
    text(median(Dist_fv_U(k,1:5)+0.1), (minNY+maxNY/5),{'EndOd'},'Color','k','Fontsize',10)
%     vline(median(Dist_fv_SO(k,:)),'k-.','EndOd');hold on;
    vline(lick_d,'g.-');hold on;
    vline((lick_d+lick_wind),'g.-');hold on;
    text(lick_d+0.1, (maxNY-maxNY/8),{'LickWindow'},'Color','g','Fontsize',10)
     
    xlabel('Time (sec)','Fontsize',14)
%     ylabel('Mean Firing rate','Fontsize',14)
    ylim([minNY,maxNY])
    xlim([0,maxNXS])
    title({'Neurons - Rejection'})
   
      text(6, (maxNY-maxNY/10),{'-'},'Color','m','Fontsize',15)
      text((6.35), (maxNY-maxNY/10),{StrNeu1},'Color','k','Fontsize',12)
      text(6, (maxNY-2*maxNY/10),{'-'},'Color','r','Fontsize',15)
      text((6.35), (maxNY-2*maxNY/10),{StrNeu2},'Color','k','Fontsize',12)
      
      

    end
     end
    cf=cf+fc;
%     end
% end