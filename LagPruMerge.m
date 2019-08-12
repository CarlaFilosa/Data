path = '/zifnas/Carla/CWEM_Project_GenProg'
cd(path)

MergeDataGenPru


%% LAG Two different region...no for the moment

% % % % % % % % % % % % % %%%%%%%%%% Organization of pairs for the Lag distribution
% % % % % % % % % % % % % pairs_lag_XXX = pairs;
% % % % % % % % % % % % % %%%%%%% Bin reduced of the last two bins (post pruned algorithm)
% % % % % % % % % % % % % % BinSizeCopy=BinSizes;
% % % % % % % % % % % % % % clear BinSizes
% % % % % % % % % % % % % % BinSizes = BinSizeCopy(:,1:10); % without the last two bins
% % % % % % % % % % % % % 
% % % % % % % % % % % % % 
% % % % % % % % % % % % % for i=1:length(BinSizes) % without the last two bins
% % % % % % % % % % % % % 
% % % % % % % % % % % % %     LagPerBin{i} = pairs_lag_XXX(find(pairs_lag_XXX(:,5)==BinSizes(i)),3);
% % % % % % % % % % % % %     lengthLagPerBin(i) =length(LagPerBin{i});
% % % % % % % % % % % % % end
% % % % % % % % % % % % % 
% % % % % % % % % % % % %  
% % % % % % % % % % % % % % sum(lengthLagPerBin)
% % % % % % % % % % % % % 
% % % % % % % % % % % % %  %%%%% Different Regions not normalized
% % % % % % % % % % % % % 
% % % % % % % % % % % % %  StrTS=strcat('Regions: VS/VTA--',strcat(nameT));
% % % % % % % % % % % % % 
% % % % % % % % % % % % % %  StrTS=strcat('Regions: VS/VTA Session: ', d_sonew.clust_params(1).session_tag);
% % % % % % % % % % % % % hlagNoNorm=cell(size(BinSizes));
% % % % % % % % % % % % % x_hlag=cell(size(BinSizes));
% % % % % % % % % % % % %  for i=1:size(BinSizes,2)
% % % % % % % % % % % % %    if ~isempty(LagPerBin{i})
% % % % % % % % % % % % % hlagNoNorm{i}=histcounts(LagPerBin{i}(:,1));
% % % % % % % % % % % % % 
% % % % % % % % % % % % % x_hlag{i}=min(LagPerBin{i}(:,1)):max(LagPerBin{i}(:,1))
% % % % % % % % % % % % % maxYY(i)=max(hlagNoNorm{i}(1,:));
% % % % % % % % % % % % %    end
% % % % % % % % % % % % %   end
% % % % % % % % % % % % % 
% % % % % % % % % % % % % %%%% LAG Distribution Plot
% % % % % % % % % % % % % 
% % % % % % % % % % % % % %%%Without last two bins
% % % % % % % % % % % % % % % % % figure(fignum(1));hold on;
% % % % % % % % % % % % % % % % % for i=1:size(BinSizes,2)
% % % % % % % % % % % % % % % % %     if ~isempty(hlagNoNorm{i})
% % % % % % % % % % % % % % % % % subplot(2,5,i);hold on;
% % % % % % % % % % % % % % % % % BinT=strcat('Bin:', mat2str(BinSizes(i)));
% % % % % % % % % % % % % % % % % bar(x_hlag{i},hlagNoNorm{i})
% % % % % % % % % % % % % % % % %  xlabel('Lag','Fontsize',12)
% % % % % % % % % % % % % % % % %  ylabel('# pairs','Fontsize',12)
% % % % % % % % % % % % % % % % %  ylim([0,max(maxYY)+1]);
% % % % % % % % % % % % % % % % %   xlim([-MaxLags(i)-0.5,MaxLags(i)+0.5])
% % % % % % % % % % % % % % % % % %  text(-0.5,3,'...//...','rotation',90,'Fontsize',12)
% % % % % % % % % % % % % % % % %   text(-MaxLags(i)-0.5,max(maxYY)-1,BinT,'Fontsize',12)
% % % % % % % % % % % % % % % % %  plot([0,0],[0,max(maxYY)+1],'k-.','Linewidth',1)
% % % % % % % % % % % % % % % % %   if i==4
% % % % % % % % % % % % % % % % %      text(-MaxLags(i)-0.5,max(maxYY)+0.5,[StrTS],'Fontsize',16) 
% % % % % % % % % % % % % % % % %   end
% % % % % % % % % % % % % % % % %  box on
% % % % % % % % % % % % % % % % %     end
% % % % % % % % % % % % % % % % % end
% % % % % % % % % % % % % 
% % % % % % % % % % % % % %%% With last two bins
% % % % % % % % % % % % % figure(fignum(1));hold on;
% % % % % % % % % % % % % for i=1:size(BinSizes,2)
% % % % % % % % % % % % %     if ~isempty(hlagNoNorm{i})
% % % % % % % % % % % % % subplot(3,4,i);hold on;
% % % % % % % % % % % % % BinT=strcat('Bin:', mat2str(BinSizes(i)));
% % % % % % % % % % % % % bar(x_hlag{i},hlagNoNorm{i})
% % % % % % % % % % % % %  xlabel('Lag','Fontsize',12)
% % % % % % % % % % % % %  ylabel('# pairs','Fontsize',12)
% % % % % % % % % % % % %  ylim([0,max(maxYY)+1]);
% % % % % % % % % % % % %   xlim([-MaxLags(i)-0.5,MaxLags(i)+0.5])
% % % % % % % % % % % % % %  text(-0.5,3,'...//...','rotation',90,'Fontsize',12)
% % % % % % % % % % % % %   text(-MaxLags(i)-0.5,max(maxYY)-1,BinT,'Fontsize',12)
% % % % % % % % % % % % %  plot([0,0],[0,max(maxYY)+1],'k-.','Linewidth',1)
% % % % % % % % % % % % %   if i==4
% % % % % % % % % % % % %      text(-MaxLags(i)-0.5,max(maxYY)+0.5,[StrTS],'Fontsize',16) 
% % % % % % % % % % % % %   end
% % % % % % % % % % % % %  box on
% % % % % % % % % % % % %     end
% % % % % % % % % % % % % end

%% %%%%%%%%%%%%%%%%%%%%%% Bin Distribution

for k = 1: SumTotAn
NP_vsvta{k}= nneuVS_All{k}*nneuVS_All{k};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% %%%%%%% Mean bin distribution
% struct_pair = [struct_pair_small,struct_pair_large];
BinW=cell(size(struct_pair)); 
WW=cell(size(struct_pair));

for k=1:length(struct_pair)
    for i=1:length(struct_pair{k}.bin)
%         for j=1: length(struct_pair{k}.bin{i})
    BinW{k}(i)=struct_pair{k}.bin{i};
%         end
    end
    if size(struct_pair{k}.bin)~=0
    WW{k}=BinW{k}(find(BinW{k}));
    if size(WW{k},2)>1
        WW{k}=WW{k}'
    end
    
    end
end
% % % % % % % % % % % 

NP_XXX=NP_vsvta;

 BinAn=zeros(length(WW),size(BinSizes,2));
for k=1:length(WW)
   
    for i=1:size(BinSizes,2)
        count=find(WW{k}(:)==BinSizes(i));
        if size(count,1)~=0
        BinAn(k,i)=size(count,1);
        end
    end
   BinAnT(k,:)=BinAn(k,:)/sum(BinAn(k,:));
   BinAnP(k,:)=BinAn(k,:)/NP_XXX{k};
end
%%
% % % % % % % % % % % 
MBinAnT = nan(1,length(BinSizes));
STBinAnT = nan(1,length(BinSizes));
MBinAnP = nan(1,length(BinSizes));
STBinAnP = nan(1,length(BinSizes));

MBinAnT(1,1:7)= mean(BinAnT(1:SumTotAn,1:7),'omitnan');
STBinAnT(1,1:7) = std(BinAnT(1:SumTotAn,1:7),0,1,'omitnan')/sqrt(size(BinAnT,1));
MBinAnP(1,1:7) = mean(BinAnP(1:SumTotAn,1:7),'omitnan');
STBinAnP(1,1:7) = std(BinAnP(1:SumTotAn,1:7),0,1,'omitnan')/sqrt(size(BinAnP,1));

MBinAnT(1,8:12)= mean(BinAnT(1:SumTotAn,8:12),'omitnan');
STBinAnT(1,8:12) = std(BinAnT(1:SumTotAn,8:12),0,1,'omitnan')/sqrt(size(BinAnT,1));
MBinAnP(1,8:12) = mean(BinAnP(1:SumTotAn,8:12),'omitnan');
STBinAnP(1,8:12) = std(BinAnP(1:SumTotAn,8:12),0,1,'omitnan')/sqrt(size(BinAnP,1));

% % % % % % % % % % % %%
BinSizesT= mat2str(BinSizes);
% StrTS = strcat(d_sonew.clust_params(1).session_tag);


%%%%%%%%%%%%%%%%%%% PLOT BIN DISTRIBUTION
figure(fignum(2));hold on;
bar(MBinAnP,'r'); hold on;
errorbar(MBinAnP,STBinAnP,'k*');hold on;
% bar(8:10,MBinAnP_Large,'r'); hold on;
% errorbar(8:10,STBinAnP_Large/2,'k*');hold on;
% title({'Mean Bin distribution-VS/VS';'Norm Poss Pairs';StrTS;BinSizesT})
% title({'Mean Bin distribution-VTA/VTA';'Norm Poss Pairs';StrTS;BinSizesT})
 title({'Mean Bin distribution';nameT})
 xticks(1:12)
 xticklabels(BinSizes)
 xlabel('\Delta (sec)')
 ylabel('Averaged pairs number/ # Possible Pairs')
 box on