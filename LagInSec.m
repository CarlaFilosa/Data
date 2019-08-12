% addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM');
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


% name1=strcat('Lastrev1_OrRevW.mat');
% load(name1)
% names=strcat('An_Disp_Pru_',name1);
% load(names)


%%
for i=1:length(BinSizes)
LagInSec{i} = BinSizes(i)*(1:MaxLags(i));
end
clear new_spM
new_spM=new_data;
TotAn=size(new_spM.spike_regionNoId,2);


fignum=33;
%%%%%%% Find number of units of the Striatum and VTA nneuS

 
%%%%%%%% For the pruned algorithm you have to add these two lines
As_across_bins=As_acr_bins_pru;
As_order=As_order_pru;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
[pairs_r,pairs_vsvs,pairs_vsvta,pairs_vtavta] = PairsInfo(As_across_bins,As_order,nneuVS,TotAn);

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

%%%%%%% Lag distribution at fixed bin
pairs_lag_XXX =pairs_vsvta{1};
for i = 1: length(pairs_vsvta)-1
pairs_lag_XXX=[pairs_lag_XXX;pairs_vsvta{i+1}];
end
% pairs_lag_XXX=pairs_vsvta;

%%

% for i = 1:length(BinSizes)
%     LagBinFix{i}= cell(size(pairs_lag_XXX));
%  for k= 1: size(pairs_lag_XXX,2)
%      disp('animal')
%      k
% %       Lag1{k}{i}=pairs_lag_XXX{k}(find(pairs_lag_XXX{k}(:,5)==BinSizes(i)),3);
%       LagBinFix{i}{k}=pairs_lag_XXX{k}(find(pairs_lag_XXX{k}(:,5)==BinSizes(i)),3);
%       LagBinFixSize{i}{k}=size(LagBinFix{i}{k},1);
% 
%  end
%      
%  end
% 
%  for i = 1:length(BinSizes)
% 
%     LB{i}=cell2mat(LagBinFix{i}');
%     SizeBFix{i}=cell2mat(LagBinFixSize{i});
%     SumSize{i} = sum(SizeBFix{i});  % the sum of the size is the total number of pairs detected in that bin
%  end
% %%%% Trasform Lag in Seconds

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


LB = LagPerBin;

for i = 1: length(LB)
    LagSec{i}=LB{i}*BinSizes(i);
end


%%% Division in two time scales (from 0.01 to 0.25  and from 0.35 and 1.6)
%%% %0.25 is the 7th BinSizes

LagSecTwoScales{1} = LagSec{1};
for i = 1:find(BinSizes==0.25)-1
    if ~isempty (LagSec{i+1})
    LagSecTwoScales{1} =[LagSecTwoScales{1};LagSec{i+1}];
    end
end

LagSecTwoScales{2} = LagSec{find(BinSizes==0.35)};
for i = find(BinSizes==0.35):find(BinSizes==1.6)-1
    if ~isempty (LagSec{i+1})
        LagSecTwoScales{2}= [LagSecTwoScales{2};LagSec{i+1}];
    end
end



StrTS=strcat('Regions: VS/VTA--', '--Analysis: ',strcat(names));
%  StrTS=strcat('Regions: VS/VTA Session: ', d_sonew.clust_params(1).session_tag);
hLagNoNorm=cell(1,2);
x_hlag=cell(1,2);
 
[hlagNoNorm{1},edges{1}]=histcounts(LagSecTwoScales{1}(:,1));
[hlagNoNorm{2},edges{2}]=histcounts(LagSecTwoScales{2}(:,1));

x_hlag{1} = -0.6:0.12:0.6;
x_hlag_zero{1} = -0.06:0.12:0.06;
x_hlag{2} = -3:0.6:3;
x_hlag_zero{2} = -0.3:0.6:0.3;
% maxY(i)=max(hhhNoNorm{i}(1,:));

hhlagNoNorm{1} = histcounts(LagSecTwoScales{1}(:,1),x_hlag{1});
LagNoZeroInd{1}= find(LagSecTwoScales{1}(:,1)~=0);
LagZeroInd{1}= find(LagSecTwoScales{1}(:,1)==0);
LagNoZeroInd{2}= find(LagSecTwoScales{2}(:,1)~=0);
LagZeroInd{2}= find(LagSecTwoScales{2}(:,1)==0);
% % figure(20);hold on;
% % histogram(LagSecTwoScales{2}(:,1),x_hlag{2})

figure(36);hold on;
subplot(1,2,1);hold on;
 h_lag_smallNozero = histogram(LagSecTwoScales{1}(LagNoZeroInd{1},1),x_hlag{1},'FaceAlpha',0.5); hold on;
 h_lag_smallZero = histogram(LagSecTwoScales{1}(LagZeroInd{1},1),x_hlag_zero{1},'FaceAlpha',0.5); hold on;
 h_lag_smallNozero.FaceColor='b';
 h_lag_smallZero.FaceColor = 'g';
 xlabel('Lag (sec)','Fontsize',12);hold on; 
 ylabel('# pairs','Fontsize',12); hold on;
 ylim([0 12]);hold on;
 title('Small Bins')
 box on;
% % % histogram(LagSecTwoScales{1}(:,1),edges{1})
% % %  xlabel('Lag','Fontsize',12)
% % %  ylabel('# pairs','Fontsize',12)
 subplot(1,2,2);hold on;
 
 h_lag_smallNozero = histogram(LagSecTwoScales{2}(LagNoZeroInd{2},1),x_hlag{2},'FaceAlpha',0.5); hold on;
 h_lag_smallZero = histogram(LagSecTwoScales{2}(LagZeroInd{2},1),x_hlag_zero{2},'FaceAlpha',0.5); hold on;
 h_lag_smallNozero.FaceColor='b';
 h_lag_smallZero.FaceColor = 'g';
 xlabel('Lag (sec)','Fontsize',12);hold on; 
 ylabel('# pairs','Fontsize',12); hold on;
 ylim([0 12]);hold on;
 title('Large Bins')
 box on;
% BinT=strcat('Bin:', mat2str(BinSizes(i)));
% bar(x_hlag{2},hlagNoNorm{2})
%  xlabel('Lag','Fontsize',12)
%  ylabel('# pairs','Fontsize',12)
%  ylim([0,max(maxY)]);
%  xlim([-MaxLags(i)-0.5,MaxLags(i)+0.5])
% %  text(-0.5,3,'...//...','rotation',90,'Fontsize',12)
%  text(-MaxLags(i)-0.5,max(maxY)-1,BinT,'Fontsize',12)
%  plot([0,0],[0,max(maxY)],'k-.','Linewidth',1)
%   if i==3
%      text(-MaxLags(i)-0.5,max(maxY)+0.5,[StrTS],'Fontsize',16) 
%   end
 box on
    
  
  