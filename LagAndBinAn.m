addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM');
% The matrix new_data is the matrix that contains the info, in particular new_data.spikeT_BegEnd has the 
% spiketrain cutted properly in databasePrepro.m
% load('StableOr_rev4.mat')
% load('AsStableOr_rev4.mat')
% load('DispAnStableOrRev4.mat')

% Ex of files
% load('Lastrev1_OrW.mat')
% load('AsStandOrW_lastrev1_rl2.mat')
% load('Displastrev1_Rl2_OrWStand.mat')

% load('Rev2_OrW.mat')
% load('AsPruOrW_rev2_rl2.mat')
% load('DispRev2_Rl2_OrWPru.mat')

% load('Rev2_RevW.mat')
% load('AsPruRevW_rev2_rl2.mat')
% load('DispRev2_Rl2_RevWPru.mat')

% load('Rev2_ExtW.mat')
% load('AsPruExtW_rev2_rl2.mat')
% load('DispRev2_Rl2_ExtWPru.mat')


name1=strcat('Discrimination_W.mat');
load(name1)
names=strcat('An_Disp_Stand_',name1);
load(names)
fignum=[16;17];

%%%%%%% Find number of units of the Striatum and VTA nneuS
clear new_spM
 new_spM=new_data;
 
%%%%%%%% For the pruned algorithm you have to add these two lines
% As_across_bins=As_acr_bins_pru;
% As_order=As_order_pru;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TotAn = size(new_spM.spike_regionNoId,2);
for k=1:TotAn
    nneuVS{k}=0;
    nneuVTA{k}=0;
end

for k=1:TotAn
    for j=1:size(new_spM.spike_regionNoId{k},2)
        [nneuVS{k}] = countNeu(new_spM.spike_regionNoId{k}(1,j),1,nneuVS{k});
        [nneuVTA{k}] = countNeu(new_spM.spike_regionNoId{k}(1,j),1,nneuVTA{k});
    end

end


%%%%% Find Pairs in one region and shared between two differnt regions
% pairs diveded by regions with all info
[pairs_r,pairs_vsvs,pairs_vsvta,pairs_vtavta] = PairsInfo(As_across_bins,As_order,nneuVS, TotAn);

clear a b A B i ii j k
%%%% # of possible pairs in the single regions and intra-region
for k=1:size(Nneu,2)
    NP_vsvs{k}=[];
    NP_vtavta{k}=[];
    NP_vsvta{k}=[];
end

for k=1:size(Nneu,2)
NP_vsvs{k}=(nneuVS{k}*(nneuVS{k}+1))/2;
NP_vtavta{k}=(nneuVTA{k}*(nneuVTA{k}+1))/2;
NP_vsvta{k}=(nneuVS{k}*nneuVTA{k});
end
%%%%%%% Struct pairs classified by bin 

[struct_pair] = PairsOrgByPair(pairs_vsvta);
% % % % % % 
% % % % % % %%%%%%% Lag distribution at fixed bin
% % % % % % 
% % % % % % pairs_lag_XXX=pairs_vsvta;
% % % % % % 
% % % % % % % figure(50);hold on
% % % % % % for i=1:length(BinSizes)
% % % % % %     LagBinFix{i}= cell(size(pairs_lag_XXX));
% % % % % %  for k= 1: size(pairs_lag_XXX,2)
% % % % % %      disp('animal')
% % % % % %      k
% % % % % %       Lag1{k}{i}=pairs_lag_XXX{k}(find(pairs_lag_XXX{k}(:,5)==BinSizes(i)),3);
% % % % % %       LagBinFix{i}{k}=pairs_lag_XXX{k}(find(pairs_lag_XXX{k}(:,5)==BinSizes(i)),3);
% % % % % %       LagBinFixSize{i}{k}=size(LagBinFix{i}{k},1);
% % % % % % 
% % % % % %  end
% % % % % %      
% % % % % %  end
% % % % % % 
% % % % % %  for i=1:length(BinSizes)
% % % % % % 
% % % % % %     LLB{i}=cell2mat(LagBinFix{i}');
% % % % % %     SizeBFix{i}=cell2mat(LagBinFixSize{i});
% % % % % %     SumSize{i} = sum(SizeBFix{i});  % the sum of the size is the total number of pairs detected in that bin
% % % % % %  end




pairs_lag_XXX =pairs_vsvta{1};
for i = 1: length(pairs_vsvta)-1
pairs_lag_XXX=[pairs_lag_XXX;pairs_vsvta{i+1}];
end
% pairs_lag_XXX=pairs_vsvta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % % % % % % for i = 1:length(BinSizes)
% % % % % % % %     LagBinFix{i}= cell(size(pairs_lag_XXX));
% % % % % % % %  for k= 1: size(pairs_lag_XXX,2)
% % % % % % % %      disp('animal')
% % % % % % % %      k
% % % % % % % % %       Lag1{k}{i}=pairs_lag_XXX{k}(find(pairs_lag_XXX{k}(:,5)==BinSizes(i)),3);
% % % % % % % %       LagBinFix{i}{k}=pairs_lag_XXX{k}(find(pairs_lag_XXX{k}(:,5)==BinSizes(i)),3);
% % % % % % % %       LagBinFixSize{i}{k}=size(LagBinFix{i}{k},1);
% % % % % % % % 
% % % % % % % %  end
% % % % % % % %      
% % % % % % % %  end
% % % % % % % % 
% % % % % % % %  for i = 1:length(BinSizes)
% % % % % % % % 
% % % % % % % %     LB{i}=cell2mat(LagBinFix{i}');
% % % % % % % %     SizeBFix{i}=cell2mat(LagBinFixSize{i});
% % % % % % % %     SumSize{i} = sum(SizeBFix{i});  % the sum of the size is the total number of pairs detected in that bin
% % % % % % % %  end


for i=1:length(BinSizes) % without the last two bins

    LagPerBin{i} = pairs_lag_XXX(find(pairs_lag_XXX(:,5)==BinSizes(i)),3);
    lengthLagPerBin(i) =length(LagPerBin{i});
end

 
NumAs_vsvta=sum(lengthLagPerBin)

LLB =LagPerBin;
%%%%% Take the absolute value for the same region distribution
% %  for  i=1:length(BinSizes)
% %      
% %      ABS_Lag{i}(:,1)=abs(LLB{i}(:,1));
% %  end
 
 
 %%%%% Different Regions not normalized


StrTS=strcat('Regions: VS/VTA--', '--Analysis: ',strcat(names));
%  StrTS=strcat('Regions: VS/VTA Session: ', d_sonew.clust_params(1).session_tag);
hlagNoNorm=cell(size(BinSizes));
x_hlag=cell(size(BinSizes));
 for i=1:size(BinSizes,2)
   if ~isempty(LLB{i})
hlagNoNorm{i}=histcounts(LLB{i}(:,1));

x_hlag{i}=min(LLB{i}(:,1)):max(LLB{i}(:,1))
maxY(i)=max(hlagNoNorm{i}(1,:));
% hh=histcounts(LLB{i})
% xx=min(LLB{i}):max(LLB{i})
%  histogram(LLB{i}(:,1))
   end
  end

%%%% LAG Distribution Plot
figure(fignum(1));hold on;
for i=1:size(BinSizes,2)
    if ~isempty(hlagNoNorm{i})
subplot(3,4,i);hold on;
BinT=strcat('Bin:', mat2str(BinSizes(i)));
bar(x_hlag{i},hlagNoNorm{i})
 xlabel('Lag','Fontsize',12)
 ylabel('# pairs','Fontsize',12)
 ylim([0,max(maxY)]);
 xlim([-MaxLags(i)-0.5,MaxLags(i)+0.5])
%  text(-0.5,3,'...//...','rotation',90,'Fontsize',12)
 text(-MaxLags(i)-0.5,max(maxY)-1,BinT,'Fontsize',12)
 plot([0,0],[0,max(maxY)],'k-.','Linewidth',1)
  if i==4
     text(-MaxLags(i)-0.5,max(maxY)+0.5,[StrTS],'Fontsize',16) 
  end
 box on
    end
end


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%% Mean bin distribution
BinW=cell(size(struct_pair)); 
WW=cell(size(struct_pair));

for k=1:size(struct_pair,2)
    
    for i=1:size(struct_pair{k}.bin,2)
        for j=1: size(struct_pair{k}.bin{i},1)
    BinW{k}(i,j)=struct_pair{k}.bin{i}(j);
        end
    end
    if size(struct_pair{k}.bin)~=0
    WW{k}=BinW{k}(find(BinW{k}));
    if size(WW{k},2)>1
        WW{k}=WW{k}'
    end
    
    end
end
% % % % % % % % % % % %%
NP_XXX=NP_vsvta;
 BinAn=zeros(size(WW,2),size(BinSizes,2));
for k=1:size(WW,2)
   
    for i=1:size(BinSizes,2)
        count=find(WW{k}(:)==BinSizes(i));
        if size(count,1)~=0
        BinAn(k,i)=size(count,1);
        end
    end
   BinAnT(k,:)=BinAn(k,:)/sum(BinAn(k,:));
   BinAnP(k,:)=BinAn(k,:)/NP_XXX{k}(:);
end
% % % % % % % % % % % 
MBinAnT = mean(BinAnT,'omitnan');
STBinAnT = std(BinAnT,0,1,'omitnan')/sqrt(size(BinAnT,1));
MBinAnP = mean(BinAnP,'omitnan');
STBinAnP = std(BinAnP,0,1,'omitnan')/sqrt(size(BinAnP,1));
% % % % % % % % % % % %%
BinSizesT= mat2str(BinSizes);
% StrTS = strcat(d_sonew.clust_params(1).session_tag);


%%%%%%%%%%%%%%%%%%% PLOT BIN DISTRIBUTION
% % % % % % figure(figurenum1);hold on;
% % % % % % bar(MAnT);
% % % % % % errorbar(MAnT,STAnT/2,'c*');
% % % % % % % title({'Mean Bin distribution-VS/VS';'Norm Detected Pairs';StrTS;BinSizesT})
% % % % % % title({'Mean Bin distribution-VTA/VTA';'Norm Detected Pairs';StrTS;BinSizesT})
% % % % % % box on
% % % % % % % title({'Mean Bin distribution-VS/VTA';'Norm Detected Pairs';StrTS;BinSizesT})
figure(fignum(2));hold on;
bar(MBinAnP,'r');
errorbar(MBinAnP,STBinAnP,'k*');
% title({'Mean Bin distribution-VS/VS';'Norm Poss Pairs';StrTS;BinSizesT})
% title({'Mean Bin distribution-VTA/VTA';'Norm Poss Pairs';StrTS;BinSizesT})
 title({'Mean Bin distribution';'Norm Poss Pairs';StrTS})
 xticks(1:12)
 xticklabels(BinSizes)
 xlabel('\Delta (sec)')
 ylabel('mean number of pairs per animal/ # possible pairs')
 box on
%%
 
 save BinLag_Stand_Discrimination_W.mat MBinAnP STBinAnP StrTS BinSizesT BinSizes MaxLags x_hlag hlagNoNorm maxY
%  
% % %  load BinLag_Rev2_Rev_Pru.mat
%%%%%%%%%%%%%%%%%%%%%% Trash
% % % % % % % % % % % %%
% % % % % % % % % % % BinSizesT= mat2str(BinSizes);
% % % % % % % % % % % StrTS = strcat(d_sonew.clust_params(1).session_tag);
% % % % % % % % % % % for k=1:size(An,2)
% % % % % % % % % % % figure(k);hold on;
% % % % % % % % % % % StrA = strcat(new_spM.par.animal(k,1));
% % % % % % % % % % % title({'Bin distribution - VS/VS';StrTS;BinSizesT;StrA})
% % % % % % % % % % % bar(An{k})
% % % % % % % % % % % end
% % % % % % % % % % % 
% % % % % % % % % % % %%
% % % % % % % % % % % PWW=cell2mat(WW');
% % % % % % % % % % % %%
% % % % % % % % % % % % for j=1:size(BinB,2)
% % % % % % % % % % %  for i =1: size(BinSizes,2)
% % % % % % % % % % %     
% % % % % % % % % % %     cc=find(PWW(:)==BinSizes(i))
% % % % % % % % % % %     if size(cc,1)~=0
% % % % % % % % % % %     f(i)=size(cc,1)
% % % % % % % % % % %     end
% % % % % % % % % % %  end
% % % % % % % % % % % %  end
% % % % % % % % % % % %%
% % % % % % % % % % % BinSizesT= mat2str(BinSizes);
% % % % % % % % % % % StrTS = strcat(d_sonew.clust_params(1).session_tag);
% % % % % % % % % % % figure(22);hold on;
% % % % % % % % % % % title({'Bin distribution - VS/VS';StrTS;BinSizesT})
% % % % % % % % % % % bar(f)
% % % % % % % % % % % 
% % % % % % % % % % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % 
% % % % % % % % % % % 
% % % % % % % % % % % %for k=1:size(struct_pair,2);
% % % % % % % % % % % k=1;
% % % % % % % % % % % figure(k); hold on;
% % % % % % % % % % % jj{k}=0;
% % % % % % % % % % % aa{k}=[];
% % % % % % % % % % % vni{k}=[];
% % % % % % % % % % % ni{k}=cell(size(struct_pair{k}.pair));
% % % % % % % % % % % for j= 1: (size(struct_pair{k}.pair,2))
% % % % % % % % % % %     ni{k}{j}=find(struct_pair{k}.lag{j}~=0);
% % % % % % % % % % %     if size(ni{k}{j},1)~=0
% % % % % % % % % % %     jj{k}=jj{k}+1;
% % % % % % % % % % %     end
% % % % % % % % % % %  a{k}=size(ni{k}{j},1);
% % % % % % % % % % %  for i= 1:a{k}
% % % % % % % % % % % vni{k}{j}=struct_pair{k}.ind{j};
% % % % % % % % % % % aa{k}(i,:)= vni{k}{j}(i,:);
% % % % % % % % % % % % 
% % % % % % % % % % % % %  assembly_activity{aaa(i)}(:,2);
% % % % % % % % % % % maxA{k}=max(assembly_activity{k}{aa{k}(i)}(:,2));
% % % % % % % % % % % plot(assembly_activity{k}{aa{k}(i)}(:,1), jj{k}+(assembly_activity{k}{aa{k}(i)}(:,2)/maxA{k}),'k');
% % % % % % % % % % % %     ni{j}=struct_pair{k}.ind{j};
% % % % % % % % % % % %      maxA=max(assembly_activity{k}{ni{j}}(:,2));
% % % % % % % % % % % %     
% % % % % % % % % % % %    plot(assembly_activity{k}{ni{j}}(:,1), jj+(assembly_activity{k}{ni{j}}(:,2)/maxA),'b')
% % % % % % % % % % %  end
% % % % % % % % % % % end