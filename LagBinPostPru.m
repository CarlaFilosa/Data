addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM');


nameNew=strcat('Lastrev1_OrRevW.mat');
load(nameNew)
alg = strcat('_Stand_');

nameSmallBig = strcat('PostPru_BinRed_An_Disp',alg,nameNew);
load(nameSmallBig);

fignum=[3,4];

%%%%%%% Find number of units of the Striatum and VTA nneuS
clear c
load('methodbase.mat','c','d')
[c_Stable] = ID_To_StableID(c,d);
d_global = d;
clear c d
% load('classlist.mat');
%  nnn(i,:)=nn(a);
%  aa(i,:)=a;
[new_data] = parId(new_data,a,nn);

TotAn = length(new_data.par);

 for k=1:TotAn
 [new_data] = labels(new_data,c_Stable,k,'stableID');
 end

%%%%%%%% For the pruned algorithm you have to add these two lines
% As_across_bins=As_acr_bins_pru;
% As_order=As_order_pru;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:TotAn
    nneuVS{k}=0;
    nneuVTA{k}=0;
end

for k=1:TotAn
    for j=1:size(new_data.spike_regionNoId{k},2)
        [nneuVS{k}] = countNeu(new_data.spike_regionNoId{k}(1,j),1,nneuVS{k});
        [nneuVTA{k}] = countNeu(new_data.spike_regionNoId{k}(1,j),1,nneuVTA{k});
    end

end


%%%%% Find Pairs in one region and shared between two differnt regions
% pairs diveded by regions with all info

[pairs_r1,pairs_vsvs1,pairs_vsvta1,pairs_vtavta1] = PairsInfo(As_across_smallBins,As_order_smallBins,nneuVS,TotAn);
[pairs_r2,pairs_vsvs2,pairs_vsvta2,pairs_vtavta2] = PairsInfo(As_across_largeBins,As_order_largeBins,nneuVS,TotAn);



%%
clear a b A B i ii j k
%%%% # of possible pairs in the single regions and intra-region
for k=1:TotAn
    NP_vsvs{k}=[];
    NP_vtavta{k}=[];
    NP_vsvta{k}=[];
end

for k=1:TotAn
NP_vsvs{k}=(nneuVS{k}*(nneuVS{k}+1))/2;
NP_vtavta{k}=(nneuVTA{k}*(nneuVTA{k}+1))/2;
NP_vsvta{k}=(nneuVS{k}*nneuVTA{k});
end
%%%%%%% Struct pairs classified by bin 

pairs_vsvta=[pairs_vsvta1,pairs_vsvta2];
 pairs_vsvs=[pairs_vsvs1,pairs_vsvs2];
[struct_pair_Small] = PairsOrgByPair(pairs_vsvta1);
[struct_pair_Large] = PairsOrgByPair(pairs_vsvta2);
%%%%%%% Lag distribution at fixed bin
%%
% pairs_vsvs=[pairs_vsvs_small,pairs_vsvs_large];
% pairs_lag_XXX =pairs_vsvs{1};
% for i = 1: length(pairs_vsvs)-1
% pairs_lag_XXX=[pairs_lag_XXX;pairs_vsvs{i+1}];
% end

% pairs_lag_XXX =pairs_vtavta{1};
% for i = 1: length(pairs_vtavta)-1
% pairs_lag_XXX=[pairs_lag_XXX;pairs_vtavta{i+1}];
% end

% pairs_lag_XXX =pairs_vsvta{1};
% for i = 1: length(pairs_vsvta)-1
% pairs_lag_XXX=[pairs_lag_XXX;pairs_vsvta{i+1}];
% end
pairs_lag_XXX=[pairs_intra_small;pairs_intra_large];
BinSizeCopy=BinSizes;
clear BinSizes
BinSizes = BinSizeCopy(:,1:10); % without the last two bins
%%

for i=1:length(BinSizes) % without the last two bins

    LagPerBin{i} = pairs_lag_XXX(find(pairs_lag_XXX(:,5)==BinSizes(i)),3);
    lengthLagPerBin(i) =length(LagPerBin{i});
end

 
sum(lengthLagPerBin)
 %%
 %%%%% Different Regions not normalized


% StrTS=strcat('Regions: VS/VTA--',strcat(nameSmallBig));

%  StrTS=strcat('Regions: VS/VTA Session: ', d_sonew.clust_params(1).session_tag);
hlagNoNorm=cell(size(BinSizes));
x_hlag=cell(size(BinSizes));
 for i=1:size(BinSizes,2)
   if ~isempty(LagPerBin{i})
hlagNoNorm{i}=histcounts(LagPerBin{i}(:,1));

x_hlag{i}=min(LagPerBin{i}(:,1)):max(LagPerBin{i}(:,1))
maxYY(i)=max(hlagNoNorm{i}(1,:));
   end
  end

%%%% LAG Distribution Plot
%%
BinSizes = BinSizeCopy;
figure(10);hold on;
for i=1:size(BinSizes,2)
    if ~isempty(hlagNoNorm{i})
subplot(2,5,i);hold on;
BinT=strcat('Bin:', mat2str(BinSizes(i)));
bar(x_hlag{i},hlagNoNorm{i})
 xlabel('Lag','Fontsize',12)
 ylabel('# pairs','Fontsize',12)
 set(gca,'FontSize',12)
 ylim([0,max(maxYY)+1]);
  xlim([-MaxLags(i)-0.5,MaxLags(i)+0.5])
%  text(-0.5,3,'...//...','rotation',90,'Fontsize',12)
  text(-MaxLags(i)-0.5,max(maxYY)-1,BinT,'Fontsize',12)
 plot([0,0],[0,max(maxYY)+1],'k-.','Linewidth',1)
%   if i==4
%      text(-MaxLags(i)-0.5,max(maxYY)+0.5,[StrTS],'Fontsize',16) 
%   end
 box on
    end
end
%%

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%% Mean bin distribution
struct_pair = [struct_pair_Small,struct_pair_Large];
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

NP_XXX=repmat(NP_vsvs,1,2);

 BinAn=zeros(size(WW,2),size(BinSizes,2));
for k=1:size(WW,2)
   
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

MBinAnT(1,1:7)= mean(BinAnT(1:TotAn,1:7),'omitnan');
STBinAnT(1,1:7) = std(BinAnT(1:TotAn,1:7),0,1,'omitnan')/sqrt(size(BinAnT,1));
MBinAnP(1,1:7) = mean(BinAnP(1:TotAn,1:7),'omitnan');
STBinAnP(1,1:7) = std(BinAnP(1:TotAn,1:7),0,1,'omitnan')/sqrt(size(BinAnP,1));

MBinAnT(1,8:10)= mean(BinAnT(TotAn+1:2*TotAn,8:10),'omitnan');
STBinAnT(1,8:10) = std(BinAnT(TotAn+1:2*TotAn,8:10),0,1,'omitnan')/sqrt(size(BinAnT,1));
MBinAnP(1,8:10) = mean(BinAnP(TotAn+1:2*TotAn,8:10),'omitnan');
STBinAnP(1,8:10) = std(BinAnP(TotAn+1:2*TotAn,8:10),0,1,'omitnan')/sqrt(size(BinAnP,1));

% % % % % % % % % % % %%
BinSizesT= mat2str(BinSizes);
% StrTS = strcat(d_sonew.clust_params(1).session_tag);


%%%%%%%%%%%%%%%%%%% PLOT BIN DISTRIBUTION



figure(3);hold on;
bar(MBinAnP,'r'); hold on;
errorbar(MBinAnP,STBinAnP/2,'k*');hold on;
% title({'Mean Bin distribution';nameSmallBig})
xticks(1:10);
xticklabels(BinSizes)
 box on