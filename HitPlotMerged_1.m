MergeDataGenPostPru

ChapOr=1; % these are the phases Original
ChapRev=2;

%%%% GivInEnd gives me the indices in a particular chapter and differently
%%%% for the phases 'stable', 'unstable' and 'all'
[cinOr, cendOr, ccrOr, inOr, finOr] = GiveInEnd( Data_All,ChapOr, 'all' );  % indices for the whole phases 
[cinRev,cendRev, ccrRev, inRev, finRev] = GiveInEnd( Data_All,ChapRev, 'all' );

startOr = cinOr;
stopOr = cendOr;

startRev = cinRev;
stopRev = cendRev;

%%%%%%%%%%%%%%%%%%
start = startOr;
stop = stopOr;
% start = startRev;
% stop = stopRev;

%%%%% For small bins
As_activity_All = As_activity_smallBins_All;
struct_pair =struct_pair_small;

FindIndicesScript
pairs=pairs_vsvta_small;
As_TaskCodePruScript

%% %%%%%%
gvect=-70:70;
gvect1= -70:470;
NoVect=nan(size(ActSortH_OdOd,1),2);

%%
%%%%%% Hit
[ActH_LickOdRew_Od]=[ActSortH_OdOd,NoVect,ActSortH_LickOd,NoVect,ActSortH_RewOd]; % ordered on the Odour
[ActH_LickOdRew_Lick]=[ActSortH_OdLick,NoVect,ActSortH_LickLick,NoVect,ActSortH_RewLick]; % ordered on the first lick
[ActH_LickOdRew_Rew]=[ActSortH_OdRew,NoVect,ActSortH_LickRew,NoVect,ActSortH_RewRew]; % ordered on the Odour

[ActH_LickOdRew_ES]=[ActSortH_OdOd,NoVect,ActSortH_LickLick,NoVect,ActSortH_RewRew]; % each panel sorted


% ActH_LickOdRewS(:,end+1)=nan;

%%

% [ActH_LickOdS1]=[ActSortH_Od1,NoVect,ActSortH_Lick1]; % Ordered on Lick
% ActH_LickOdS1(:,end+1)=nan;
% 
% [ActH_LickOdS2]=[ActSortH_Od22,NoVect,ActSortH_Lick22]; % Ordered on second odour (Hit)
% ActH_LickOdS2(:,end+1)=nan;

tick=[find(gvect==-50),find(gvect==0),find(gvect==50),...
      size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-50),...
      size(ActH_Lick,2)+size(NoVect,2)+find(gvect==0),...
      size(ActH_Lick,2)+size(NoVect,2)+find(gvect==50)...
      2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect1==-50),...
      2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect1==0),...
      2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect1==50),...
      2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect1==100),...
      2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect1==150),...
      2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect1==200),...
      2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect1==250),...
      2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect1==300),...
      2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect1==350),...
      2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect1==400)];
ticklabel =[-0.5,0,0.5,-0.5,0,0.5,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4];
%%

ActH1=ActH_LickOdRew_Od;
ActH2=ActH_LickOdRew_Lick;
ActH3=ActH_LickOdRew_Rew;
ActH4=ActH_LickOdRew_ES;
% ActH1(:,end+1)=nan;
% ActH2(:,end+1)=nan;
% ActH3(:,end+1)=nan;
% ActH4(:,end+1)=nan;


%% Three plot
%  addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla')

switch Dir
    case 'small'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
 case 'smallSync'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
    case 'large'
BinDirText=strcat(Dir,'--BinSize Range: [0.35, 0.6]');
  case 'largeSync'
BinDirText=strcat(Dir,'--BinSize Range: [0.35, 0.6]');
    case 'vs->vta Dir'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
    case 'vs->vta Inv'
BinDirText=strcat(Dir,'--BinSize Range: [0.35, 0.6]');
    case 'vta->vs Dir'
BinDirText=strcat(Dir,'--BinSize Range: [0.35, 0.6]');
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

if (ismember(start,startOr) & ismember(stop,stopOr))
    PhaseStr = strcat('Original');
elseif (ismember(start,startRev) & ismember(stop,stopRev))
    PhaseStr = strcat('Reversal');
elseif (ismember(start,startOr) & ismember(stop,stopRev))
    PhaseStr = strcat('Or and Rev');   
end

figure(fignum(5));hold on;
subplot(3,1,1)
pp=pcolor(ActH1);hold on;
xlim([1,size(ActH1,2)])
ylim([1,size(ActH1,1)])
xticks(tick)
xticklabels(ticklabel)
% text(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-10),-18,'Time(sec)')
% text(find(gvect==-10),-18,'Time(sec)')
% text(2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==-10),-18,'Time(sec)')
text(find(gvect==-10),size(ActH_Lick,1)+10,'Odor onset','Fontsize',12)
text(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-10),size(ActH_Lick,1)+10,'First Lick','Fontsize',12)
text(2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==-10),size(ActH_Lick,1)+10,'Reward','Fontsize',12)
plot([size(ActH_Lick,2)+size(NoVect,2)+find(gvect==0),size(ActH_Lick,2)+size(NoVect,2)+find(gvect==0)],[0,size(ActH2,1)],'-.r','Linewidth',1);
plot([2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==0),2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==0)],[0,size(ActH2,1)],'-.r','Linewidth',1);
plot([find(gvect==0),find(gvect==0)],[0,size(ActH2,1)],'-.r','Linewidth',1);
% text(find(gvect==-70)-10,size(ActH_Lick,1)/2,'Sorted on Odour','Rotation',90)
% vline(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==0),'-.r');
% vline(2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==0),'-.r');
% vline(find(gvect==0),'-.r');
set(pp, 'EdgeColor', 'none');
grid on;
ylabel('Pairs Actvity')
title('Hit - Odour Sorting -')


subplot(3,1,2)
pp=pcolor(ActH2);hold on;
xlim([1,size(ActH2,2)])
ylim([1,size(ActH2,1)])
xticks(tick)
xticklabels(ticklabel)
% text(find(gvect==-10),-18,'Time(sec)')
% text(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-10),-18,'Time(sec)')
% text(2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==-10),-18,'Time(sec)')
text(find(gvect==-10),size(ActH_Lick,1)+10,'Odor onset','Fontsize',12)
text(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-10),size(ActH_Lick,1)+10,'First Lick','Fontsize',12)
text(2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==-10),size(ActH_Lick,1)+10,'Reward','Fontsize',12)
plot([size(ActH_Lick,2)+size(NoVect,2)+find(gvect==0),size(ActH_Lick,2)+size(NoVect,2)+find(gvect==0)],[0,size(ActH2,1)],'-.r','Linewidth',1);
plot([2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==0),2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==0)],[0,size(ActH2,1)],'-.r','Linewidth',1);
plot([find(gvect==0),find(gvect==0)],[0,size(ActH2,1)],'-.r','Linewidth',1);
% vline(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==0),'-.r');
% vline(2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==0),'-.r');
% vline(find(gvect==0),'-.r');
set(pp, 'EdgeColor', 'none');
grid on;
ylabel('Pairs Actvity')
title('Hit - First Lick Sorting -')
text(size(ActH1,2)+20,size(ActH_Lick,1),{nameT, BinDirText,PhaseStr},'rotation',270,'Fontsize',12);

subplot(3,1,3)
pp=pcolor(ActH3);hold on;
xlim([1,size(ActH3,2)])
ylim([1,size(ActH3,1)])
xticks(tick)
xticklabels(ticklabel)
text(find(gvect==-10),-50,'Time(sec)','Fontsize',12)
text(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-10),-50,'Time(sec)','Fontsize',12)
text(2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==-10),-50,'Time(sec)','Fontsize',12)
text(find(gvect==-10),size(ActH_Lick,1)+10,'Odor onset','Fontsize',12)
text(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-10),size(ActH_Lick,1)+10,'First Lick','Fontsize',12)
text(2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==-10),size(ActH_Lick,1)+10,'Reward','Fontsize',12)
plot([size(ActH_Lick,2)+size(NoVect,2)+find(gvect==0),size(ActH_Lick,2)+size(NoVect,2)+find(gvect==0)],[0,size(ActH2,1)],'-.r','Linewidth',1);
plot([2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==0),2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==0)],[0,size(ActH2,1)],'-.r','Linewidth',1);
plot([find(gvect==0),find(gvect==0)],[0,size(ActH2,1)],'-.r','Linewidth',1);
% text(find(gvect==-70)-10,size(ActH_Lick,1)/2,'Sorted on Odour','Rotation',90)
% vline(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==0),'-.r');
% vline(2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==0),'-.r');
% vline(find(gvect==0),'-.r');
set(pp, 'EdgeColor', 'none');
grid on;
ylabel('Pairs Actvity')
title('Hit - Reward Sorting -')
%%   
addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla')
figure(4);hold on;
pp=pcolor(ActH1);hold on;
xlim([1,size(ActH1,2)])
ylim([1,size(ActH1,1)])
xticks(tick)
xticklabels(ticklabel)
text(find(gvect==-10),-18,'Time(sec)')
text(find(gvect==-10),size(ActH_Lick,1)+10,'Odor onset')
text(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-10),-18,'Time(sec)')
text(2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==-10),-18,'Time(sec)')
text(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-10),size(ActH_Lick,1)+10,'First Lick')
text(2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==-10),size(ActH_Lick,1)+10,'Reward')

% text(find(gvect==-70)-10,size(ActH_Lick,1)/2,'Sorted on Odour','Rotation',90)
vline(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==0),'-.r');
vline(2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==0),'-.r');
vline(find(gvect==0),'-.r');
set(pp, 'EdgeColor', 'none');
% grid on;
ylabel('Pairs Actvity')
title('Hit')

%%
%%%%%%%%%%%%%%%%
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



