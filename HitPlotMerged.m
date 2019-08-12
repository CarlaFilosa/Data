MergeDataGenPostPru
%% Indices for start and end

ChapOr=1; % these are the phases Original
ChapRev=2;

%%%% GivInEnd gives me the indices in a particular chapter and differently
%%%% for the phases 'stable', 'unstable' and 'all'
[cinOr, cendOr, ccrOr, inOr, finOr] = GiveInEnd( Data_All,ChapOr, 'all' );  % indices for the whole phases 
[cinRev,cendRev, ccrRev, inRev, finRev] = GiveInEnd( Data_All,ChapRev, 'all' );

start = cinOr;
stop = cendRev;


startOr = cinOr;
stopOr = cendOr;

startRev = cinRev;
stopRev = cendRev;

%% Recale assembly activity and find indices for stable and non-stable phase

%%%%% For small bins
As_activity_All = As_activity_smallBins_All;
struct_pair =struct_pair_small;

%%%%% For large Bins
% As_activity_All = As_activity_largeBins_All;
% struct_pair =struct_pair_large;

minInt=0.01;
for k=1:SumTotAn
  if ~isempty(As_activity_All{k}) && ~isempty(Data_All.spikeT_BegEnd{k})
  binref.B{k}=min(min(Data_All.spikeT_BegEnd{k},[],2)):minInt:max(max(Data_All.spikeT_BegEnd{k},[],2));

   for i=1:size(As_activity_All{k},1)
   
      [AsTime{k}{i}(:,1)]=TimeRescaleMin(As_activity_All{k}{i}(:,1),As_activity_All{k}{i}(:,1), binref.B{k});
      [AsTime{k}{i}(:,2)]=TimeRescaleMin(As_activity_All{k}{i}(:,2),As_activity_All{k}{i}(:,1), binref.B{k});
    
   end
 end
end

%% Whole trials for rewarded and not
% % for k= 1:SumTotAn
% % [FindIndT{k}] = FindIndicesForAssem(Data_All,binref.B,start,stop,k,'trial');
% % end
% % 
% % %%% the whole trial
% % [asTime_US] = AsOnInd(AsTime,FindIndT,Data_All,start);

%% the whole trial
% % Region='vsvta';
% % % Region='vsvs';
% % % Region='vtavta';
% % 
% % %Dir stays for "Directionality", can assume these values
% % %'vs->vta Dir','vta->vs Dir','vs->vta Inv', 'vta->vs Inv'0
% % %  Dir='vs->vta Dir';
% % %Dir='vs->vta Inv';
% % % Dir='vs->vta AllBins';
% % 
% % %  Dir='vta->vs Dir';
% % %  Dir='vta->vs Inv';
% % %  Dir='vta->vs AllBins';
% % %    Dir=0;
% % % Dir='small';
% % Dir='large';
% % 
% % hit = 1;
% % rej = 3;
% % fal = 2;
% % miss = 4;
% % 
% % 
% % for k=1:SumTotAn
% %     if ~isempty(struct_pair{k}.pair)  % & (length(new_spM.eventsOld{k})~=length(new_spM.events{k}))
% % [as_sel_H{k},MeanTr_H{k},MNormTr_H{k},StETr_H{k},MaxA_H{k},AsForMean_H{k},new_struct_pairAll{k},BinDirAll{k}] = As_TaskCodePru(Dir,struct_pair,Data_All,asTime_US,hit,start(k),k);
% % [as_sel_R{k},MeanTr_R{k},MNormTr_R{k},StETr_R{k},MaxA_R{k},AsForMean_R{k}] = As_TaskCodePru(Dir,struct_pair,Data_All,asTime_US,rej,start(k),k);
% % [as_sel_F{k},MeanTr_F{k},MNormTr_F{k},StETr_F{k},MaxA_F{k},AsForMean_F{k}] = As_TaskCodePru(Dir,struct_pair,Data_All,asTime_US,fal,start(k),k);
% % [as_sel_M{k},MeanTr_M{k},MNormTr_M{k},StETr_M{k},MaxA_M{k},AsForMean_M{k}] = As_TaskCodePru(Dir,struct_pair,Data_All,asTime_US,miss,start(k),k);
% %     end 
% % end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% All Phases Together

% % for k=1:size(MaxA_H,2)
% %     [MaxATot{k}] = MaxMaxTr(MaxA_H{k},MaxA_R{k},MaxA_F{k},MaxA_M{k});
% % end
% %  
% % 
% % [MNorm_H] = NormTotMax(MaxATot,MeanTr_H);
% % [MNorm_R] = NormTotMax(MaxATot,MeanTr_R);
% % [MNorm_F] = NormTotMax(MaxATot,MeanTr_F);
% % [MNorm_M] = NormTotMax(MaxATot,MeanTr_M);
% % % [MNT_HU] = NormTotMax(MaxATot,MeanTr_HU);
% % % [MNT_RU] = NormTotMax(MaxATot,MeanTr_RU);
% % % [MNT_FU] = NormTotMax(MaxATot,MeanTr_FU);
% % % [MNT_MU] = NormTotMax(MaxATot,MeanTr_MU);
% % 
% % 
% % [Mean_P] = PrDim(MNorm_H,MNorm_R,MNorm_M,MNorm_F,MNorm_H,MNorm_R,MNorm_M,MNorm_F);  % two times the input arguments simply because the function was made to have 8 input
% % 
% % [ActH,MaxInd1,SortMaxR1] = ActAsAllNew(length(MNorm_H), MNorm_H,Mean_P);
% % [ActSortH] = SortAct(ActH,SortMaxR1);
% % 
% % [ActR,MaxInd,SortMaxR] = ActAsAllNew(length(MNorm_H), MNorm_R,Mean_P);
% % [ActSortR] = SortAct(ActR,SortMaxR);
% % 
% % 
% % [ActM,MaxInd,SortMaxR] = ActAsAllNew(length(MNorm_H), MNorm_M,Mean_P);
% % [ActSortM] = SortAct(ActM,SortMaxR);
% % 
% % 
% % [ActF,MaxInd,SortMaxR] = ActAsAllNew(length(MNorm_H), MNorm_F,Mean_P);
% % [ActSortF] = SortAct(ActF,SortMaxR);
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%% Each Panel sorted o the maxima of the first panel
% % [ActSortH1] = SortAct(ActH,SortMaxR1);
% % 
% % [ActSortR1] = SortAct(ActR,SortMaxR1);
% % 
% % [ActSortM1] = SortAct(ActM,SortMaxR1);
% % 
% % [ActSortF1] = SortAct(ActF,SortMaxR1);
%%
%%%%%%%%%%%%%%%%%%%%%%% Hit Plot
% % ActAllH=[];  ActAllR=[]; ActAllM=[];  ActAllF=[]; 
% % ActAllH=ActSortH;  ActAllR=ActSortR; ActAllM=ActSortM;  ActAllF=ActSortF; 
% % %  ActAllH=ActSortH1; ActAllR=ActSortR1; ActAllM=ActSortM1; ActAllF=ActSortF1;
% % 
% % ActAllH(:,end+1)=nan; ActAllR(:,end+1)=nan; ActAllM(:,end+1)=nan; ActAllF(:,end+1)=nan;
% % % h = pcolor(X,Y,C);
% % % set(h, 'EdgeColor', 'none');
% % 
% % switch Dir
% %     case 'small'
% % BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
% %     case 'large'
% % BinDirText=strcat(Dir,'--BinSize Range: [0.35, 1.6]');
% %     case 'vs->vta Dir'
% % BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
% %     case 'vs->vta Inv'
% % BinDirText=strcat(Dir,'--BinSize Range: [0.35, 1.6]');
% %     case 'vta->vs Dir'
% % BinDirText=strcat(Dir,'--BinSize Range: [0.35, 1.6]');
% %     case 'vta->vs Inv'
% % BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
% % % BinDirText=strcat(Dir,'--BinSize Range: [0.01, 1.6]');
% % case 'vta->vs AllBins' 
% %     BinDirText=strcat(Dir);
% % case 'vs->vta AllBins'
% % BinDirText=strcat(Dir);
% %     otherwise
% %         BinDirText=strcat('Region--',Region,'--All Bin');
% % end
% % addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla');
% % figure(fignum(6));
% % subplot(2,2,1)
% % pp=pcolor(ActAllH);hold on;
% % set(pp, 'EdgeColor', 'none');
% % grid on;
% % vline(1/minInt,'-.r');hold on;
% % vline(1.5/minInt,'-.r');hold on;
% % ylabel('Pairs Actvity')
% % title('Hit')
% % 
% % 
% % box on;
% % 
% % subplot(2,2,2)
% % pp=pcolor(ActAllR);hold on;
% % set(pp, 'EdgeColor', 'none');hold on;
% % grid on;hold on;
% % vline(1/minInt,'-.r');hold on;
% % vline(1.5/minInt,'-.r');hold on;
% % title('Rejection')
% % text(1000,size(ActAllM,1)/2,{nameT;BinDirText},'rotation',270,'Fontsize',12);
% % box on;
% % 
% % subplot(2,2,3)
% % pp=pcolor(ActAllM);hold on;
% % set(pp, 'EdgeColor', 'none');
% % grid on; hold on;
% % vline(1/minInt,'-.r');hold on;
% % vline(1.5/minInt,'-.r');hold on;
% % xlabel('Time (10^{-2} sec)')
% % ylabel('Pairs activity')
% % title('Miss')
% % box on;
% % % text(-250,size(ActAllM,1)/2,{'Sorted on the 1st Panel'},'rotation',90,'Fontsize',12);
% %  text(-250,size(ActAllM,1)/2,{'Each panel sorted'},'rotation',90,'Fontsize',12);
% % 
% % 
% % 
% % subplot(2,2,4)
% % pp=pcolor(ActAllF);hold on;
% % set(pp, 'EdgeColor', 'none');
% % grid on;
% % vline(1/minInt,'-.r');hold on;
% % vline(1.5/minInt,'-.r');hold on;
% % xlabel('Time (10^{-2} sec)')
% % title('False Alarm')
% % box on;
% %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Around reward and odour
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k= 1:SumTotAn
[FindIndOd{k}] = FindIndicesForAssem(Data_All,binref.B,start,stop,k,'odor',0.705,0.705);  %around odour
end
for k= 1:SumTotAn
[FindIndRew{k}] = FindIndicesForAssem(Data_All,binref.B,start,stop,k,'reward',0.705,0.705);  %around reward
end

for k= 1:SumTotAn
[FindIndNoRew{k}] = FindIndicesForAssem(Data_All,binref.B,start,stop,k,'HitNoRew',0.705,2.255); % supposed rewarded, not rewarded
end

for k= 1:SumTotAn
[FindIndRewT{k}] = FindIndicesForAssem(Data_All,binref.B,start,stop,k,'HitRew',0.705,2.255);  %rewarded trials one sec before odor onset, one sec after end odour
end
for k= 1:SumTotAn
[FindIndFLick{k}] = FindIndicesForAssem(Data_All,binref.B,start,stop,k,'FirstLick',0.705,0.705);  %around first lick
end

%%
%%% Odour 
[asTime_US_Od] = AsOnInd(AsTime,FindIndOd);
%%% Reward
[asTime_US_Rew] = AsOnInd(AsTime,FindIndRew);
%%%%%% Hit Not rewarded in the interval odour plus
[asTime_US_NoRew] = AsOnInd(AsTime,FindIndNoRew);
%%%% Hit rewarded in the interval odour plus
[asTime_US_RewT] = AsOnInd(AsTime,FindIndRewT);
%%%% Around the first Lick
[asTime_US_FLick] = AsOnInd(AsTime,FindIndFLick);


%%%%%%%%%%%%%%%%%%%%%%%%%
%%
Region='vsvta';
% Region='vsvs';
% Region='vtavta';

%Dir stays for "Directionality", can assume these values
%'vs->vta Dir','vta->vs Dir','vs->vta Inv', 'vta->vs Inv'0
%  Dir='vs->vta Dir';
%Dir='vs->vta Inv';
% Dir='vs->vta AllBins';

%  Dir='vta->vs Dir';
%  Dir='vta->vs Inv';
%  Dir='vta->vs AllBins';
%    Dir=0;
Dir='small';
% Dir='large';
hit = 1;
rej = 3;
fal = 2;
miss = 4;

%% %%%%%%% Odour

for k=1:SumTotAn
    if ~isempty(struct_pair{k}.pair)  % & (length(new_spM.eventsOld{k})~=length(new_spM.events{k}))
[as_sel_HOd{k},MeanTr_HOd{k},MNormTr_HOd{k},StETr_HOd{k},MaxA_HOd{k},AsForMean_HOd{k},~,~] = As_TaskCodePru(Dir,struct_pair,Data_All,asTime_US_Od,hit,start(k),k);
[as_sel_ROd{k},MeanTr_ROd{k},MNormTr_ROd{k},StETr_ROd{k},MaxA_ROd{k},AsForMean_ROd{k},~,~] = As_TaskCodePru(Dir,struct_pair,Data_All,asTime_US_Od,rej,start(k),k);  
[as_sel_FOd{k},MeanTr_FOd{k},MNormTr_FOd{k},StETr_FOd{k},MaxA_FOd{k},AsForMean_FOd{k},~,~] = As_TaskCodePru(Dir,struct_pair,Data_All,asTime_US_Od,fal,start(k),k);    
[as_sel_MOd{k},MeanTr_MOd{k},MNormTr_MOd{k},StETr_MOd{k},MaxA_MOd{k},AsForMean_MOd{k},~,~] = As_TaskCodePru(Dir,struct_pair,Data_All,asTime_US_Od,miss,start(k),k); 
    end 
end
%% First Lick

for k=1:SumTotAn
    if ~isempty(struct_pair{k}.pair)  % & (length(new_spM.eventsOld{k})~=length(new_spM.events{k}))
[as_sel_HLick{k},MeanTr_HLick{k},MNormTr_HLick{k},StETr_HLick{k},MaxA_HLick{k},AsForMean_HLick{k},~,~] = As_TaskCodePru(Dir,struct_pair,Data_All,asTime_US_FLick,hit,start(k),k);
% [as_sel_RLick{k},MeanTr_RLick{k},MNormTr_RLick{k},StETr_RLick{k},MaxA_RLick{k},AsForMean_RLick{k},~,~] = As_TaskCodePru(Dir,struct_pair,Data_All,asTime_US_FLick,rej,start(k),k);  
% [as_sel_FLick{k},MeanTr_FLick{k},MNormTr_FLick{k},StETr_FLick{k},MaxA_FLick{k},AsForMean_FLick{k},~,~] = As_TaskCodePru(Dir,struct_pair,Data_All,asTime_US_FLick,fal,start(k),k);    
% [as_sel_MLick{k},MeanTr_MLick{k},MNormTr_MLick{k},StETr_MLick{k},MaxA_MLick{k},AsForMean_MLick{k},~,~] = As_TaskCodePru(Dir,struct_pair,Data_All,asTime_US_FLick,miss,start(k),k); 
    end 
end
%% %% Reward

for k=1:SumTotAn
    if ~isempty(struct_pair{k}.pair)  % & (length(new_spM.eventsOld{k})~=length(new_spM.events{k}))
[as_sel_Rew{k},MeanTr_HRew{k},MNormTr_HRew{k},StETr_HRew{k},MaxA_HRew{k},AsForMean_HRew{k},~,~] = As_TaskCodePru(Dir,struct_pair,Data_All,asTime_US_Rew,hit,start(k),k);

    end 
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% All Phases Together

for k=1:length(MaxA_HOd)
    [MaxATot{k}] = MaxMaxTr(MaxA_HOd{k},MaxA_FOd{k},MaxA_HLick{k},MaxA_FLick{k});
end
 

[MNorm_HOd] = NormTotMax(MaxATot,MeanTr_HOd);
[MNorm_ROd] = NormTotMax(MaxATot,MeanTr_ROd);
[MNorm_FOd] = NormTotMax(MaxATot,MeanTr_FOd);
[MNorm_MOd] = NormTotMax(MaxATot,MeanTr_MOd);

[MNorm_HLick] = NormTotMax(MaxATot,MeanTr_HLick);
[MNorm_RLick] = NormTotMax(MaxATot,MeanTr_RLick);
[MNorm_FLick] = NormTotMax(MaxATot,MeanTr_FLick);
[MNorm_MLick] = NormTotMax(MaxATot,MeanTr_MLick);
%%%%%%%% Hit/ Different sorting
[ActH_Od,MaxI_HO1,SortH_O1] = ActAsAllSameLengthInt(length(MNorm_HOd), MNorm_HOd);
[ActSortH_Od] = SortAct(ActH_Od,SortH_O1);

[ActH_Lick,MaxI_L1,SortH_L1] = ActAsAllSameLengthInt(length(MNorm_HLick), MNorm_HLick);
[ActSortH_Lick] = SortAct(ActH_Lick,SortH_O1);
[ActSortH_Lick1] = SortAct(ActH_Lick,SortH_L1);

[ActSortH_Od1] = SortAct(ActH_Od,SortH_L1);

%%%%%%% False Alarm
[ActF_Od,MaxI_FO1,SortF_O1] = ActAsAllSameLengthInt(length(MNorm_FOd), MNorm_FOd);
[ActSortF_Od] = SortAct(ActF_Od,SortF_O1);

[ActF_Lick,MaxFI_L1,SortF_L1] = ActAsAllSameLengthInt(length(MNorm_FLick), MNorm_FLick);
[ActSortF_Lick] = SortAct(ActF_Lick,SortF_O1);
[ActSortF_Lick1] = SortAct(ActF_Lick,SortF_L1);
[ActSortF_Lick2] = SortAct(ActF_Lick,SortH_O1);  % on the first panel odor
[ActSortF_Lick3] = SortAct(ActF_Lick,SortH_L1); % on the first panel lick

[ActSortF_Od1] = SortAct(ActF_Od,SortF_L1);
[ActSortF_Od2] = SortAct(ActF_Od,SortH_O1);   % on the first panel odor
[ActSortF_Od3] = SortAct(ActF_Od,SortH_L1);  % on the first panel lick

[ActSortH_Od22] = SortAct(ActH_Od,SortF_O1);
[ActSortH_Lick22] = SortAct(ActH_Lick,SortF_O1);


gvect=-70:70;

NoVect=nan(size(ActSortH_Od,1),2);
%%%%%% Hit
[ActH_LickOdS]=[ActSortH_Od,NoVect,ActSortH_Lick]; % ordered on the Odour
ActH_LickOdS(:,end+1)=nan;

[ActH_LickOdS1]=[ActSortH_Od1,NoVect,ActSortH_Lick1]; % Ordered on Lick
 ActH_LickOdS1(:,end+1)=nan;

%%%%%% False Alarm

[ActF_LickOdS]=[ActSortF_Od,NoVect,ActSortF_Lick]; % ordered on the Odour
 ActF_LickOdS(:,end+1)=nan;

[ActF_LickOdS1]=[ActSortF_Od1,NoVect,ActSortF_Lick1]; % Ordered on Lick
 ActF_LickOdS1(:,end+1)=nan;

[ActF_LickOdS2]=[ActSortF_Od2,NoVect,ActSortF_Lick2]; % Ordered on First Odor (Hit)
 ActF_LickOdS2(:,end+1)=nan;
 
 [ActF_LickOdS3]=[ActSortF_Od3,NoVect,ActSortF_Lick3]; % Ordered on First lick (Hit)
 ActF_LickOdS3(:,end+1)=nan;

 [ActH_LickOdS2]=[ActSortH_Od22,NoVect,ActSortH_Lick22]; % Ordered on second odour (Hit)
 ActH_LickOdS2(:,end+1)=nan;

tick=[find(gvect==-50),find(gvect==0),find(gvect==50),...
      size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-50),...
      size(ActH_Lick,2)+size(NoVect,2)+find(gvect==0),...
      size(ActH_Lick,2)+size(NoVect,2)+find(gvect==50)];
ticklabel =[-0.5,0,0.5,-0.5,0,0.5];
%%
ActH1 = ActH_LickOdS2; ActH2 = ActH_LickOdS1; ActF1 = ActF_LickOdS; ActF2 = ActF_LickOdS3; 
addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla')
figure(3);hold on;
subplot(2,2,1)
pp=pcolor(ActH1);hold on;
xlim([1,size(ActH1,2)])
ylim([1,size(ActH1,1)])
xticks(tick)
xticklabels(ticklabel)
text(find(gvect==-10),-18,'Time(sec)')
text(find(gvect==-10),size(ActH_Lick,1)+10,'Odor onset')
text(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-10),-18,'Time(sec)')
text(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-10),size(ActH_Lick,1)+10,'First Lick')
% text(find(gvect==-70)-10,size(ActH_Lick,1)/2,'Sorted on Odour','Rotation',90)
vline(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==0),'-.r');
vline(find(gvect==0),'-.r');
set(pp, 'EdgeColor', 'none');
grid on;
ylabel('Pairs Actvity')
title('Hit')

subplot(2,2,2)
pp=pcolor(ActH2);hold on;
xlim([1,size(ActH2,2)])
ylim([1,size(ActH2,1)])
xticks(tick)
xticklabels(ticklabel)
text(find(gvect==-10),-18,'Time(sec)')
text(find(gvect==-10),size(ActH_Lick,1)+10,'Odor onset')
text(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-10),-18,'Time(sec)')
text(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-10),size(ActH_Lick,1)+10,'First Lick')
% text(find(gvect==-70)-10,size(ActH_Lick,1)/2,'Sorted on First Lick','Rotation',90)
vline(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==0),'-.r');
vline(find(gvect==0),'-.r');
set(pp, 'EdgeColor', 'none');
grid on;
ylabel('Pairs Actvity')
title('Hit')


subplot(2,2,3)
pp=pcolor(ActF1);hold on;
xlim([1,size(ActF1,2)])
ylim([1,size(ActF1,1)])
xticks(tick)
xticklabels(ticklabel)
text(find(gvect==-10),-18,'Time(sec)')
text(find(gvect==-10),size(ActH_Lick,1)+10,'Odor onset')
text(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-10),-18,'Time(sec)')
text(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-10),size(ActH_Lick,1)+10,'First Lick')
text(find(gvect==-70)-20,size(ActH_Lick,1)/2,'Sorted on Odour','Rotation',90)
vline(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==0),'-.r');
vline(find(gvect==0),'-.r');
set(pp, 'EdgeColor', 'none');
grid on;
ylabel('Pairs Actvity')
title('False Alarm')

subplot(2,2,4)
pp=pcolor(ActF2);hold on;
xlim([1,size(ActF2,2)])
ylim([1,size(ActF2,1)])
xticks(tick)
xticklabels(ticklabel)
text(find(gvect==-10),-18,'Time(sec)')
text(find(gvect==-10),size(ActH_Lick,1)+10,'Odor onset')
text(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-10),-18,'Time(sec)')
text(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-10),size(ActH_Lick,1)+10,'First Lick')
text(find(gvect==-70)-20,size(ActH_Lick,1)/2,'Sorted on First Lick','Rotation',90)
vline(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==0),'-.r');
vline(find(gvect==0),'-.r');
set(pp, 'EdgeColor', 'none');
grid on;
ylabel('Pairs Actvity')
title('False Alarm')


% % 
% % 
% % [Mean_P] = PrDim(MNorm_HOd,MNorm_ROd,MNorm_MOd,MNorm_FOd,MNorm_HLick,MNorm_RLick,MNorm_MLick,MNorm_FLick);  % two times the input arguments simply because the function was made to have 8 input
% % 
% % [ActH,MaxInd1,SortMaxR1] = ActAsAllNew(length(MNorm_H), MNorm_H,Mean_P);
% % [ActSortH] = SortAct(ActH,SortMaxR1);
% % 
% % [ActR,MaxInd,SortMaxR] = ActAsAllNew(length(MNorm_H), MNorm_R,Mean_P);
% % [ActSortR] = SortAct(ActR,SortMaxR);
% % 
% % 
% % [ActM,MaxInd,SortMaxR] = ActAsAllNew(length(MNorm_H), MNorm_M,Mean_P);
% % [ActSortM] = SortAct(ActM,SortMaxR);
% % 
% % 
% % [ActF,MaxInd,SortMaxR] = ActAsAllNew(length(MNorm_H), MNorm_F,Mean_P);
% % [ActSortF] = SortAct(ActF,SortMaxR);
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%% Each Panel sorted o the maxima of the first panel
% % [ActSortH1] = SortAct(ActH,SortMaxR1);
% % 
% % [ActSortR1] = SortAct(ActR,SortMaxR1);
% % 
% % [ActSortM1] = SortAct(ActM,SortMaxR1);
% % 
% % [ActSortF1] = SortAct(ActF,SortMaxR1);
%%
%%%%%%%%%%%%%%%%%%%%%%% Hit Plot
% % ActAllH=[];  ActAllR=[]; ActAllM=[];  ActAllF=[]; 
% % ActAllH=ActSortH;  ActAllR=ActSortR; ActAllM=ActSortM;  ActAllF=ActSortF; 
% % %  ActAllH=ActSortH1; ActAllR=ActSortR1; ActAllM=ActSortM1; ActAllF=ActSortF1;
% % 
% % ActAllH(:,end+1)=nan; ActAllR(:,end+1)=nan; ActAllM(:,end+1)=nan; ActAllF(:,end+1)=nan;
% % % h = pcolor(X,Y,C);
% % % set(h, 'EdgeColor', 'none');
% % 
% % switch Dir
% %     case 'small'
% % BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
% %     case 'large'
% % BinDirText=strcat(Dir,'--BinSize Range: [0.35, 1.6]');
% %     case 'vs->vta Dir'
% % BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
% %     case 'vs->vta Inv'
% % BinDirText=strcat(Dir,'--BinSize Range: [0.35, 1.6]');
% %     case 'vta->vs Dir'
% % BinDirText=strcat(Dir,'--BinSize Range: [0.35, 1.6]');
% %     case 'vta->vs Inv'
% % BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
% % % BinDirText=strcat(Dir,'--BinSize Range: [0.01, 1.6]');
% % case 'vta->vs AllBins' 
% %     BinDirText=strcat(Dir);
% % case 'vs->vta AllBins'
% % BinDirText=strcat(Dir);
% %     otherwise
% %         BinDirText=strcat('Region--',Region,'--All Bin');
% % end
% % addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla');
% % figure(fignum(6));
% % subplot(2,2,1)
% % pp=pcolor(ActAllH);hold on;
% % set(pp, 'EdgeColor', 'none');
% % grid on;
% % vline(1/minInt,'-.r');hold on;
% % vline(1.5/minInt,'-.r');hold on;
% % ylabel('Pairs Actvity')
% % title('Hit')
% % 
% % 
% % box on;
% % 
% % subplot(2,2,2)
% % pp=pcolor(ActAllR);hold on;
% % set(pp, 'EdgeColor', 'none');hold on;
% % grid on;hold on;
% % vline(1/minInt,'-.r');hold on;
% % vline(1.5/minInt,'-.r');hold on;
% % title('Rejection')
% % text(1000,size(ActAllM,1)/2,{nameT;BinDirText},'rotation',270,'Fontsize',12);
% % box on;
% % 
% % subplot(2,2,3)
% % pp=pcolor(ActAllM);hold on;
% % set(pp, 'EdgeColor', 'none');
% % grid on; hold on;
% % vline(1/minInt,'-.r');hold on;
% % vline(1.5/minInt,'-.r');hold on;
% % xlabel('Time (10^{-2} sec)')
% % ylabel('Pairs activity')
% % title('Miss')
% % box on;
% % % text(-250,size(ActAllM,1)/2,{'Sorted on the 1st Panel'},'rotation',90,'Fontsize',12);
% %  text(-250,size(ActAllM,1)/2,{'Each panel sorted'},'rotation',90,'Fontsize',12);
% % 
% % 
% % 
% % subplot(2,2,4)
% % pp=pcolor(ActAllF);hold on;
% % set(pp, 'EdgeColor', 'none');
% % grid on;
% % vline(1/minInt,'-.r');hold on;
% % vline(1.5/minInt,'-.r');hold on;
% % xlabel('Time (10^{-2} sec)')
% % title('False Alarm')
% % box on;
% %  

%%
%%%%%%%%% Rewarded around reward 
for k=1:SumTotAn
    if ~isempty(struct_pair{k}.pair)  % & (length(new_spM.eventsOld{k})~=length(new_spM.events{k}))
[as_sel_HRew{k},MeanTr_HRew{k},MNormTr_HRew{k},StETr_HRew{k},MaxA_HRew{k},AsForMean_HRew{k},struct_pairRewAll{k},~] = As_TaskCodePru(Dir,struct_pair,Data_All,asTime_US_Rew,hit,start(k),k);
    end 
end


%%%%%%%%% Not Rewarded in the interval (fv_on-0.5; fv_on+2.25) (the maximum lick delay is 1250 and after 1sec the reward is supposed to be there)
for k=1:SumTotAn
    if ~isempty(struct_pair{k}.pair)  % & (length(new_spM.eventsOld{k})~=length(new_spM.events{k}))
[as_sel_HNoRew{k},MeanTr_HNoRew{k},MNormTr_HNoRew{k},StETr_HNoRew{k},MaxA_HNoRew{k},AsForMean_HNoRew{k},struct_pairNoRewAll{k},~] = As_TaskCodePru(Dir,struct_pair,Data_All,asTime_US_NoRew,hit,start(k),k);
    end 
end


%%%%%%%%% Rewarded in the interval (fv_on-0.5; fv_on+2.25) (the maximum lick delay is 1250 and after 1sec the reward is supposed to be there)
for k=1:SumTotAn
    if ~isempty(struct_pair{k}.pair)  % & (length(new_spM.eventsOld{k})~=length(new_spM.events{k}))
[as_sel_HRewT{k},MeanTr_HRewT{k},MNormTr_HRewT{k},StETr_HRewT{k},MaxA_HRewT{k},AsForMean_HRewT{k},struct_pairRewTAll{k},~] = As_TaskCodePru(Dir,struct_pair,Data_All,asTime_US_RewT,hit,start(k),k);
    end 
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% All Phases Together


for k=1:size(MaxA_HOd,2)
    [MaxAOd{k}] = MaxMaxTr(MaxA_HOd{k},MaxA_ROd{k},MaxA_FOd{k},MaxA_MOd{k});
end
 

[MNorm_HOd] = NormTotMax(MaxAOd,MeanTr_HOd);
[MNorm_ROd] = NormTotMax(MaxAOd,MeanTr_ROd);
[MNorm_FOd] = NormTotMax(MaxAOd,MeanTr_FOd);
[MNorm_MOd] = NormTotMax(MaxAOd,MeanTr_MOd);

%%
[Mean_POd] = PrDim(MNorm_HOd,MNorm_ROd,MNorm_MOd,MNorm_FOd,MNorm_HOd,MNorm_ROd,MNorm_MOd,MNorm_FOd);  % two times the input arguments simply because the function was made to have 8 input

[ActHOd,MaxIndOd1,SortMaxROd1] = ActAsAllNew(length(MNorm_HOd), MNorm_HOd,Mean_POd);
[ActSortHOd] = SortAct(ActHOd,SortMaxROd1);

[ActROd,MaxIndOd,SortMaxROd] = ActAsAllNew(length(MNorm_HOd), MNorm_ROd,Mean_POd);
[ActSortROd] = SortAct(ActROd,SortMaxROd);


[ActMOd,MaxIndOd,SortMaxROd] = ActAsAllNew(length(MNorm_HOd), MNorm_MOd,Mean_POd);
[ActSortMOd] = SortAct(ActMOd,SortMaxROd);


[ActFOd,MaxIndOd,SortMaxROd] = ActAsAllNew(length(MNorm_HOd), MNorm_FOd,Mean_POd);
[ActSortFOd] = SortAct(ActFOd,SortMaxROd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Each Panel sorted o the maxima of the first panel
[ActSortHOd1] = SortAct(ActHOd,SortMaxROd1);

[ActSortROd1] = SortAct(ActROd,SortMaxROd1);

[ActSortMOd1] = SortAct(ActMOd,SortMaxROd1);

[ActSortFOd1] = SortAct(ActFOd,SortMaxROd1);

%%

% tick = -0.5:minInt:0.5;
%%%%%%%%%%%%%%%%%%%%%%% Hit Plot
ActAllH=ActSortHOd;  ActAllR=ActSortROd; ActAllM=ActSortMOd;  ActAllF=ActSortFOd; 
%  ActAllH=ActSortHOd1; ActAllR=ActSortROd1; ActAllM=ActSortMOd1; ActAllF=ActSortFOd1;

ActAllH(:,end+1)=nan; ActAllR(:,end+1)=nan; ActAllM(:,end+1)=nan; ActAllF(:,end+1)=nan;
% h = pcolor(X,Y,C);
% set(h, 'EdgeColor', 'none');

switch Dir
    case 'small'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
    case 'large'
BinDirText=strcat(Dir,'--BinSize Range: [0.35, 1.6]');
    case 'vs->vta Dir'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
    case 'vs->vta Inv'
BinDirText=strcat(Dir,'--BinSize Range: [0.35, 1.6]');
    case 'vta->vs Dir'
BinDirText=strcat(Dir,'--BinSize Range: [0.35, 1.6]');
    case 'vta->vs Inv'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
% BinDirText=strcat(Dir,'--BinSize Range: [0.01, 1.6]');
case 'vta->vs AllBins' 
    BinDirText=strcat(Dir);
case 'vs->vta AllBins'
BinDirText=strcat(Dir);
    otherwise
        BinDirText=strcat('Region--',Region,'--All Bin');
end
addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla');
figure(fignum(2));
subplot(2,2,1)
pp=pcolor(ActAllH);hold on;
set(pp, 'EdgeColor', 'none');
grid on;
xticks([1,50,101])
xticklabels([-0.5,0,0.5])
% vline(1/minInt,'-.r');hold on;
vline(0.5/minInt,'-.r');hold on;
ylabel('Pairs Actvity')
title('Hit')


box on;

subplot(2,2,2)
pp=pcolor(ActAllR);hold on;
set(pp, 'EdgeColor', 'none');hold on;
grid on;hold on;
% xticklabels(tick)
% vline(1/minInt,'-.r');hold on;
vline(0.5/minInt,'-.r');hold on;
xticks([1,50,101])
xticklabels([-0.5,0,0.5])
title('Rejection')
text(110,size(ActAllM,1)/2,{nameT;BinDirText},'rotation',270,'Fontsize',12);
box on;

subplot(2,2,3)
pp=pcolor(ActAllM);hold on;
set(pp, 'EdgeColor', 'none');
grid on; hold on;
% xticklabels(tick)
% vline(1/minInt,'-.r');hold on;
vline(0.5/minInt,'-.r');hold on;
xticks([1,50,101])
xticklabels([-0.5,0,0.5])
xlabel('Time around odour onset(sec)')
ylabel('Pairs activity')
title('Miss')
box on;
% text(-110,size(ActAllM,1)/2,{'Sorted on the 1st Panel'},'rotation',90,'Fontsize',12);
 text(-20,size(ActAllM,1)/2,{'Each panel sorted'},'rotation',90,'Fontsize',12);



subplot(2,2,4)
pp=pcolor(ActAllF);hold on;
set(pp, 'EdgeColor', 'none');
grid on;
% xticklabels(tick)
% vline(1/minInt,'-.r');hold on;
vline(0.5/minInt,'-.r');hold on;
xticks([1,50,101])
xticklabels([-0.5,0,0.5])
xlabel('Time around odour onset(sec)')
title('False Alarm')
box on;
 %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:size(MaxA_HRew,2)
    [MaxARewOd{k}] = MaxMaxTr(MaxA_HRew{k},MaxA_HOd{k});
end


[MNorm_HOd_2] = NormTotMax(MaxARewOd,MeanTr_HOd);
[MNorm_HRew_2] = NormTotMax(MaxARewOd,MeanTr_HRew);


[Mean_POdRew] = PrDim(MNorm_HRew_2,MNorm_HOd_2);  % two times the input arguments simply because the function was made to have 8 input

[ActHOd_2,MaxIndOd1_2,SortMaxOd1_2] = ActAsAllNew(length(MNorm_HOd_2), MNorm_HOd_2,Mean_POdRew);
[ActSortHOd_2] = SortAct(ActHOd_2,SortMaxOd1_2);


[ActHRew_2,MaxIndRew1_2,SortMaxRew1_2] = ActAsAllNew(length(MNorm_HRew_2), MNorm_HRew_2,Mean_POdRew);
[ActSortHRew_2] = SortAct(ActHRew_2,SortMaxRew1_2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Each Panel sorted o the maxima of the first panel
[ActSortHOd1_2] = SortAct(ActHOd_2,SortMaxOd1_2);

[ActSortHRew1_2] = SortAct(ActHRew_2,SortMaxOd1_2);

[ActSortHOd2_2] = SortAct(ActHOd_2,SortMaxRew1_2);

% tick = -0.5:minInt:0.5;
%% %%%%%%%%%%%%%%%%%%%%% Hit Plot

ActAllH=[]; ActAllR=[];
 ActAllH=ActSortHOd_2;  ActAllR=ActSortHRew_2; 
% ActAllH=ActSortHOd1_2;  ActAllR=ActSortHRew1_2; 
% ActAllH=ActSortHOd2_2;  ActAllR=ActSortHRew_2;

ActAllH(:,end+1)=nan; ActAllR(:,end+1)=nan; 
% h = pcolor(X,Y,C);
% set(h, 'EdgeColor', 'none');

switch Dir
    case 'small'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
    case 'large'
BinDirText=strcat(Dir,'--BinSize Range: [0.35, 1.6]');
    case 'vs->vta Dir'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
    case 'vs->vta Inv'
BinDirText=strcat(Dir,'--BinSize Range: [0.35, 1.6]');
    case 'vta->vs Dir'
BinDirText=strcat(Dir,'--BinSize Range: [0.35, 1.6]');
    case 'vta->vs Inv'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
% BinDirText=strcat(Dir,'--BinSize Range: [0.01, 1.6]');
case 'vta->vs AllBins' 
    BinDirText=strcat(Dir);
case 'vs->vta AllBins'
BinDirText=strcat(Dir);
    otherwise
        BinDirText=strcat('Region--',Region,'--All Bin');
end
addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla');
figure(fignum(4));
subplot(1,2,1)
pp=pcolor(ActAllH);hold on;
set(pp, 'EdgeColor', 'none');
grid on;
xticks([1,50,101])
xticklabels([-0.5,0,0.5])
% vline(1/minInt,'-.r');hold on;
vline(0.5/minInt,'-.r');hold on;
ylabel('Pairs Actvity')
title('Odour')
text(-20,size(ActAllR,1)/2,{'Each panel sorted'},'rotation',90,'Fontsize',12);
xlabel('Time(sec)')
box on;

subplot(1,2,2)
pp=pcolor(ActAllR);hold on;
set(pp, 'EdgeColor', 'none');hold on;
grid on;hold on;
% xticklabels(tick)
% vline(1/minInt,'-.r');hold on;
xticks([1,50,101])
xticklabels([-0.5,0,0.5])
vline(0.5/minInt,'-.r');hold on;
title('Reward')
xlabel('Time(sec)')

text(110,size(ActAllH,1)/2,{nameT, BinDirText},'rotation',270,'Fontsize',12);
box on;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



