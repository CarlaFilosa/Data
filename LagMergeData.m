%%%% Programs that does the lag distribution for each specific couple of
%%%% units

nameMerge{1} = strcat('Lastrev1_OrRevW.mat');
nameMerge{2} = strcat('Revconc_OrRevW.mat');
nameMerge{3} = strcat('Succession_2PhW.mat');
nameMerge{4} = strcat('Rev2_OrRevW.mat');
nameMerge{5} = strcat('Rev3_OrRevW.mat');

nameT = strcat('Lastrev1; Revconc; Succession; Rev2; Rev3');


alg= strcat('Stand_');
% algPru= strcat('Pru_');
fignum =[3,4,5,6,7,8,9,10];

for i =1: length(nameMerge)
%     nameSmallBigMerge{i} = strcat('SmallBig_An_Disp',alg,nameMerge{i})     %postpruned algorithm
    nameSmallBigMerge{i} = strcat('PostPru_BinRed_An_Disp',alg,nameMerge{i})
    load(nameSmallBigMerge{i})


clear c
load('classlist.mat');
%  nnn(i,:)=nn(a);
%  aa(i,:)=a;
[new_data] = parId(new_data,a,nn);


 for k=1:length(new_data.par)
 [new_data] = labels(new_data,c,k);
 end
 
 very_data{i}=new_data;
%%%% Postpruned
As_across_smallBinsT{i} = As_across_smallBins;
As_order_smallBinsT{i} = As_order_smallBins;
As_across_largeBinsT{i} = As_across_largeBins;
As_order_largeBinsT{i} = As_order_largeBins; 


clear new_data a nn As_across_smallBins As_order_smallBins As_across_largeBins As_order_largeBins 
end
%%
 clearvars -except very_data As_across_smallBinsT As_order_smallBinsT As_across_largeBinsT As_order_largeBinsT c nameMerge nameT fignum BinSizes MaxLags nameT
 
 %%
 for i =1:length(very_data)
    TotAn(i)=length(very_data{i}.par);
for k=1:TotAn(i)
    nneuVS{i}{k}=0;
    nneuVTA{i}{k}=0;
end
end

SumTotAn=sum(TotAn);

for i =1:length(very_data)
for k=1:TotAn(i)
    for j=1:size(very_data{i}.spike_regionNoId{k},2)
        [nneuVS{i}{k}] = countNeu(very_data{i}.spike_regionNoId{k}(1,j),1,nneuVS{i}{k});
        [nneuVTA{i}{k}] = countNeu(very_data{i}.spike_regionNoId{k}(1,j),2,nneuVTA{i}{k});
    end

end
end

As_across_smallBins_All=As_across_smallBinsT{1};
As_across_largeBins_All=As_across_largeBinsT{1};
As_order_smallBins_All=As_order_smallBinsT{1};
As_order_largeBins_All=As_order_largeBinsT{1};
nneuVS_All = nneuVS{1};
nneuVTA_All = nneuVTA{1};
Data_All.par=very_data{1}.par;
for i=1:length(very_data)-1
    [As_across_smallBins_All]=[As_across_smallBins_All,As_across_smallBinsT{i+1}];
    [As_across_largeBins_All]=[As_across_largeBins_All,As_across_largeBinsT{i+1}];
    [As_order_smallBins_All]=[As_order_smallBins_All,As_order_smallBinsT{i+1}];
    [As_order_largeBins_All]=[As_order_largeBins_All,As_order_largeBinsT{i+1}];
    [nneuVS_All]=[nneuVS_All,nneuVS{i+1}];
    [nneuVTA_All]=[nneuVTA_All,nneuVTA{i+1}];
    [Data_All.par]=[Data_All.par,very_data{i+1}.par];
end



[pairs_r_small,pairs_vsvs_small,pairs_vsvta_small,pairs_vtavta_small] = PairsInfo(As_across_smallBins_All,As_order_smallBins_All,nneuVS_All,SumTotAn);
[pairs_r_large,pairs_vsvs_large,pairs_vsvta_large,pairs_vtavta_large] = PairsInfo(As_across_largeBins_All,As_order_largeBins_All,nneuVS_All,SumTotAn);

clear a b A B i ii j k

[struct_pair_small] = PairsOrgByPair(pairs_vsvta_small);
[struct_pair_large] = PairsOrgByPair(pairs_vsvta_large);


[struct_pair_small] = PairNeuLabelsID(struct_pair_small,Data_All);
[struct_pair_large] = PairNeuLabelsID(struct_pair_large,Data_All);
copy_struct_pair_small=struct_pair_small;
copy_struct_pair_large=struct_pair_large;
clear struct_pair

%%


%%%%%%%%%% Organization of pairs for the Lag distribution

pairs_vsvta = [pairs_vsvta_small,pairs_vsvta_large];

pairs_lag_XXX =pairs_vsvta{1};
for i = 1: length(pairs_vsvta)-1
pairs_lag_XXX=[pairs_lag_XXX;pairs_vsvta{i+1}];
end


%%%%%%% Bin reduced of the last two bins (post pruned algorithm)
BinSizeCopy=BinSizes;
clear BinSizes
BinSizes = BinSizeCopy(:,1:10); % without the last two bins


for i=1:length(BinSizes) % without the last two bins

    LagPerBin{i} = pairs_lag_XXX(find(pairs_lag_XXX(:,5)==BinSizes(i)),3);
    lengthLagPerBin(i) =length(LagPerBin{i});
end

 
% sum(lengthLagPerBin)

 %%%%% Different Regions not normalized

 StrTS=strcat('Regions: VS/VTA--',strcat(nameT));

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

figure(fignum(1));hold on;
for i=1:size(BinSizes,2)
    if ~isempty(hlagNoNorm{i})
subplot(2,5,i);hold on;
BinT=strcat('Bin:', mat2str(BinSizes(i)));
bar(x_hlag{i},hlagNoNorm{i})
 xlabel('Lag','Fontsize',12)
 ylabel('# pairs','Fontsize',12)
 ylim([0,max(maxYY)]);
  xlim([-MaxLags(i)-0.5,MaxLags(i)+0.5])
%  text(-0.5,3,'...//...','rotation',90,'Fontsize',12)
  text(-MaxLags(i)-0.5,max(maxYY)-1,BinT,'Fontsize',12)
 plot([0,0],[0,max(maxYY)],'k-.','Linewidth',1)
  if i==4
     text(-MaxLags(i)-0.5,max(maxYY)+0.5,[StrTS],'Fontsize',16) 
  end
 box on
    end
end

%% %%%%%%%%%%%%%%%%%%%%%% Bin Distribution

for k = 1: SumTotAn
NP_vsvta{k}= nneuVS_All{k}*nneuVTA_All{k};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% %%%%%%% Mean bin distribution
struct_pair = [struct_pair_small,struct_pair_large];
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

NP_XXX=repmat(NP_vsvta,1,2);

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
% % % % % % % % % % % 
MBinAnT = nan(1,length(BinSizes));
STBinAnT = nan(1,length(BinSizes));
MBinAnP = nan(1,length(BinSizes));
STBinAnP = nan(1,length(BinSizes));

MBinAnT(1,1:7)= mean(BinAnT(1:SumTotAn,1:7),'omitnan');
STBinAnT(1,1:7) = std(BinAnT(1:SumTotAn,1:7),0,1,'omitnan')/sqrt(size(BinAnT,1));
MBinAnP(1,1:7) = mean(BinAnP(1:SumTotAn,1:7),'omitnan');
STBinAnP(1,1:7) = std(BinAnP(1:SumTotAn,1:7),0,1,'omitnan')/sqrt(size(BinAnP,1));

MBinAnT(1,8:10)= mean(BinAnT(SumTotAn+1:end,8:10),'omitnan');
STBinAnT(1,8:10) = std(BinAnT(SumTotAn+1:end,8:10),0,1,'omitnan')/sqrt(size(BinAnT,1));
MBinAnP(1,8:10) = mean(BinAnP(SumTotAn+1:end,8:10),'omitnan');
STBinAnP(1,8:10) = std(BinAnP(SumTotAn+1:end,8:10),0,1,'omitnan')/sqrt(size(BinAnP,1));

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
 xticks(1:10)
 xticklabels(BinSizes)
 
 box on