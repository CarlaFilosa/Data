clear all;
%%% Non Parametric version of a repeated measures Anova
MergeDataGenPostPru_1

MergeDataGenPostPru_2

%%%%%%%%%%% Chapters division %%%%%%%%%%%%%%%%%%%%%%%%%5
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
start = startOr+1; % because of the baseline I can't take the first trial of the original phase
stop = stopOr;
% start = startRev;
% stop = stopRev;

%%%%% For small bins
As_activity_All = As_activity_smallBins_All;
%%%% to pass to the stable id, that I chooce as universal ID
Data_All.clust_params = Data_Par;
[new_rates.id_Stable] = ID_To_StableID_Gen(new_rates.id,Data_All);
%% Cross REGION    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
struct_pair =struct_pair_small;
Region='vsvta';
Dir='small';
hit=1;
cvta = [4,5,6,4.1,4.2,4.3,4.4,0,20];
cvs=[1,2.1,2.2,3];  

col1=9;
col2= [10,10,10,12,12,12,12,12,10];
scriptDivisionPairs;

pairs = pairs_vsvta_small;

SelectCouple=[];
SelectCouple = pairs_MT{8};
TitleText = {'MSN No Repetions','NonLaser No Repetions','MSN with Repetions','NonLaser with Repetitions',};


SelectCouple_Matx=[];
SelectCouple_Matx = SelectCouple{1};
for k=1:length(SelectCouple)-1
    SelectCouple_Matx =[SelectCouple_Matx; SelectCouple{k+1}];
end


VTA_un_NoRepVec=[];
VS_un_NoRepVec=[];
VTA_un_NoRep=[];
VS_un_NoRep=[];
VTA_un=[];
VS_un=[];

VTA_un_NoRepVec = unique(SelectCouple_Matx(:,8));
VS_un_NoRepVec = unique(SelectCouple_Matx(:,7));

for i = 1: size(VTA_un_NoRepVec,1)
    VTA_un_NoRep(i) = new_rates.value(new_rates.id_Stable == VTA_un_NoRepVec(i));
end

for i = 1: size(VS_un_NoRepVec,1)
    VS_un_NoRep(i) = new_rates.value(new_rates.id_Stable == VS_un_NoRepVec(i));
end


for i = 1: size(SelectCouple_Matx,1)
    VTA_un(i) = new_rates.value(new_rates.id_Stable == SelectCouple_Matx(i,8));
end

for i = 1: size(SelectCouple_Matx,1)
    VS_un(i) = new_rates.value(new_rates.id_Stable == SelectCouple_Matx(i,7));
end
clear i
step=2;
xticksVec = 0:5:round(max(VTA_un_NoRep))+step;
XLim = round(max([max(VTA_un_NoRep),max(VS_un_NoRep),max(VTA_un),max(VS_un)]))+step;
figure(17)
subplot(2,2,1)
histogram(VS_un_NoRep,min(VS_un_NoRep):step:max(VS_un_NoRep)+step,'FaceColor',[0 0.5 0.5]);
xlabel('Baseline Firing Rate (Hz)')
xlim([0 XLim])
ylim([0 15])
xticks(xticksVec)
set(gca,'FontSize',12)
% xticklabels(xticksVec)
title(TitleText{1})
subplot(2,2,2)
histogram(VTA_un_NoRep,min(VTA_un_NoRep):step:max(VTA_un_NoRep)+step,'FaceColor',[0.7 0 0.7]);
xlim([0 XLim])
ylim([0 15])
xticks(xticksVec)
xlabel('Baseline Firing Rate (Hz)')
set(gca,'FontSize',12)
title(TitleText{2})
subplot(2,2,3)
histogram(VS_un,min(VS_un):step:max(VS_un)+step,'FaceColor',[0 0.5 0.5]);
xlim([0 XLim])
ylim([0 15])
xticks(xticksVec)
title(TitleText{3})
xlabel('Baseline Firing Rate (Hz)')
set(gca,'FontSize',12)
subplot(2,2,4)
histogram(VTA_un,min(VTA_un):step:max(VTA_un)+step,'FaceColor',[0.7 0 0.7]);
xlim([0 XLim])
xticks(xticksVec)
xlabel('Baseline Firing Rate (Hz)')
set(gca,'FontSize',12)
ylim([0 15])
title(TitleText{4})

DiffFr_VTAVS=[];
DiffFr_VTAVS=abs(VTA_un-VS_un);
xticksL=0:2:round(max(DiffFr_VTAVS))+step;
figure(19)
histogram(DiffFr_VTAVS,round(min(DiffFr_VTAVS)):step:round(max(DiffFr_VTAVS)+step),'FaceColor',[0 0 0.5]);
xticks([xticksL])
title('Absolute Firing Rate Difference between VS and VTA units forming a couple')
xlabel('Absolute Firing Rate Difference (Hz)')
text(24, 10, 'MSN - NonLasert.VTA','FontSize',12)
set(gca,'FontSize',12)



DiffFr_VTAVS_Pos=[];
DiffFr_VTAVS_Neg=[];
DiffFr_VTAVS_Null=[];
DiffFr_VTAVS_Pos = VS_un(find(SelectCouple_Matx(:,3)>0))-VTA_un(find(SelectCouple_Matx(:,3)>0));
DiffFr_VTAVS_Neg = VS_un(find(SelectCouple_Matx(:,3)<0))-VTA_un(find(SelectCouple_Matx(:,3)<0));
DiffFr_VTAVS_Null = VS_un(find(SelectCouple_Matx(:,3)==0))-VTA_un(find(SelectCouple_Matx(:,3)==0));
figure(8)
histogram(abs(DiffFr_VTAVS_Pos),20)
figure(9)
histogram(abs(DiffFr_VTAVS_Neg),20)
figure(10)
histogram(abs(DiffFr_VTAVS_Null),20)
%%