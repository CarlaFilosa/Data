% This script intend to put toghether all the means
 addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM');
% The matrix new_spM is built in the program MaxDataPrepro.m
%  load('basedatabese.mat')


load('Lastrev1_OrRevW.mat')
load('AsPruOrRevW_lastrev1_rl2.mat')
load('DispLastrev1_Rl2_OrRevWPru.mat');
%%
[new_data] = parId(new_data,newInd,nn);
%% Chose the region and the directionality
% Region of interest, it could be: 'vsvta','vsvs','vtavta';
% Region  ='all';
Region = 'vsvta';

%Dir stays for "Directionality", can assume these values
%'vs->vta Dir','vta->vs Dir','vs->vta Inv', 'vta->vs Inv'0
%  Dir='vs->vta Dir';
% Dir='vs->vta Inv';
Dir='vs->vta AllBins';


% as_act=assembly_activity;
% as_acr_bins=As_across_bins;
% as_order=As_order;

as_act=ass_act_pru;
as_acr_bins=As_acr_bins_pru;
as_order=As_order_pru;

ChapOr=1; % these are the phases Original
ChapRev=2; % reversal phase

[cin, cend, ccr, inS, finS] = GiveInEnd( new_data,ChapOr, 'stable' );

[~, ~, ~, inU, finU] = GiveInEnd( new_data,ChapOr, 'unstable' );

[cin_rev, cend_rev, ccr_rev, inrevSS, finrevS] = GiveInEnd( new_data,ChapRev, 'stable' );

[~, ~, ~, inrevU, finrevU] = GiveInEnd( new_data,ChapRev, 'unstable' );

new_spM=new_data;
new_spM.eventsOld=new_spM.events;
new_spM.events=[];
new_spM.events=new_spM.eventsCut;

TotAn=size(new_spM.events,2);
%%%%%%% Unstable Original or reversal
startU=cin;  % original
stopU=ccr-1;
% startU=cin_rev;  % reversal
% stopU=ccr_rev-1;
%%%%%%%% Stable- original or reversal
startS=ccr;
stopS=cend-1;  %original
% startS=ccr_rev;
% stopS=cend_rev-1;  %reversal


new_spM=new_data;
new_spM.eventsOld=new_spM.events;
new_spM.events=[];
new_spM.events=new_spM.eventsCut;

minInt=0.01;
for k=1:size(new_spM.events,2)
    
binref.B{k}=min(min(new_spM.spikeT_BegEnd{k},[],2)):minInt:max(max(new_spM.spikeT_BegEnd{k},[],2));
if ~isempty(as_act{k})
 for i=1:size(as_act{k},1)
   
      [AsTime{k}{i}(:,1)]=TimeRescaleMin(as_act{k}{i}(:,1),as_act{k}{i}(:,1), binref.B{k});
      [AsTime{k}{i}(:,2)]=TimeRescaleMin(as_act{k}{i}(:,2),as_act{k}{i}(:,1), binref.B{k});
    
 end
end
end

 %% Big matrix all together
 
 prepost = 1; % one second befor of after the reward / or, if I need, before/after the first/last lick of the trial
[ new_spM.events] = LicksRewInt( new_spM.events,TotAn,prepost );

Base=0.7; % 0.7 sec after starting odor-onset is my baseline
 
  for k=1:TotAn
     for l=startU(k):stopS(k)
 FindInd.Odor{k}{l}=[];
 FindInd.Base{k}{l}=[];
 FindInd.LickPre{k}{l}=[]; 
 FindInd.LickPost{k}{l}=[];
 FindInd.RewPre{k}{l}=[]; 
 FindInd.RewPost{k}{l}=[];
 FindInd.Code{k}{l}=[];
     end
  end
 lu=[];

 
 for k=1:TotAn
     for l=startU(k)+1:stopS(k)
     lu(l)=length(new_spM.events{k}(l).licks);
     lr=length(new_spM.events{k}(l).reward_time);
             
             FindInd.Odor{k}{l}= find(binref.B{k}>=new_spM.events{k}(l).fv_on & binref.B{k}<=new_spM.events{k}(l).fv_off);
             FindInd.Base{k}{l}= find(binref.B{k}>=new_spM.events{k}(l).fv_on-Base & binref.B{k}<=new_spM.events{k}(l).fv_on);
             FindInd.Code{k}{l}=new_spM.events{k}(l).rew_code;
             
             if lr==1  % rewards 
             FindInd.RewPre{k}{l}= find(binref.B{k}>=new_spM.events{k}(l).rewardPre & binref.B{k}<=new_spM.events{k}(l).reward_time);
             FindInd.RewPost{k}{l}= find(binref.B{k}>=new_spM.events{k}(l).reward_time & binref.B{k}<=new_spM.events{k}(l).rewardPost);
             end
             
             if lu(l)~=0 %licks
             for i=1:lu(l)
             FindInd.LickPre{k}{l}{i}= find(binref.B{k}>=new_spM.events{k}(l).lickPre(i) & binref.B{k}<=new_spM.events{k}(l).licks(i));
             FindInd.LickPost{k}{l}{i}= find(binref.B{k}>=new_spM.events{k}(l).licks(i) & binref.B{k}<=new_spM.events{k}(l).lickPost(i));
             end
             end
            
 
     end
 end
 
 %% Indices for the licks of the same trial...Put all togheter this indices to use the same function
 for k=1:TotAn
    for l=1:length(FindInd.LickPre{k})
        FindInd.LickPreVect{k}{l}=[];
    end
 end
 
 
 for k=1:TotAn
    for l=1:length(FindInd.LickPre{k})
in=0;
for i=1:length(FindInd.LickPre{k}{l})
for j=1:length(FindInd.LickPre{k}{l}{i})
    FindInd.LickPreVect{k}{l}(j+in)=FindInd.LickPre{k}{l}{i}(j);
end
in=in+length(FindInd.LickPre{k}{l}{i});
end
    end
 end
 
 for k=1:TotAn
    for l=1:length(FindInd.LickPost{k})
        FindInd.LickPostVect{k}{l}=[];
    end
 end
 
 
 for k=1:TotAn
    for l=1:length(FindInd.LickPost{k})
in=0;
for i=1:length(FindInd.LickPost{k}{l})
for j=1:length(FindInd.LickPost{k}{l}{i})
    FindInd.LickPostVect{k}{l}(j+in)=FindInd.LickPost{k}{l}{i}(j);
end
in=in+length(FindInd.LickPost{k}{l}{i});
end
    end
 end
 

%% Find assemblies activity on indices
% IN FIND_IND I have 
[asTimeInt.US] = AsOnInd(AsTime,FindInd.US,new_spM,startU);
[asTimeInt.USBase] = AsOnInd(AsTime,FindInd.USBase,new_spM,startU);

[asTimeInt.RewPre] = AsOnInd(AsTime,FindInd.RewPre,new_spM,startU);
[asTimeInt.RewPost] = AsOnInd(AsTime,FindInd.RewPost,new_spM,startU);

[asTimeInt.LickPre] = AsOnInd(AsTime,FindInd.LickPreVect,new_spM,startU);
[asTimeInt.LickPost] = AsOnInd(AsTime,FindInd.LickPostVect,new_spM,startU);


 %% At the end before the plots
%%%%% I build my usual struct_pair

% Find number of units of the Striatum and VTA nneuS




for k=1:size(new_spM.spike_regionNoId,2)
    nneuVS{k}=0;
    nneuVTA{k}=0;
end

for k=1:size(new_spM.spike_regionNoId,2)
    for j=1:size(new_spM.spike_regionNoId{k},2)
        [nneuVS{k}] = countNeu(new_spM.spike_regionNoId{k}(1,j),1,nneuVS{k});
        [nneuVTA{k}] = countNeu(new_spM.spike_regionNoId{k}(1,j),1,nneuVTA{k});
    end

end

% pairs diveded by regions with all info
[pairs_r,pairs_vsvs,pairs_vsvta,pairs_vtavta] = PairsInfo(as_acr_bins,as_order,nneuVS);

clear a b A B i ii j k
%%%%%%%%%%%%%%%%%%%%%%% # of possible pairs in the single regions and intra-region
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
switch Region
    case 'vsvta'
[struct_pair] = PairsOrg(pairs_vsvta,NP_vsvta);
    case 'vsvs'
        [struct_pair] = PairsOrg(pairs_vsvs,NP_vsvs);
    case 'vtavta'
        [struct_pair] = PairsOrg(pairs_vtavta,NP_vtavta);
    case 'all'
         [struct_pair] = PairsOrg(pairs_r,Nneu);
end

%% Mean-activitity on each trials



% lick='no';

for k=1:length(asTimeInt.US)
% [as_sel.Trial{k},MeanWind.Trial{k},NewStructPair.T{k},BinDir.T{k}] = AsActCode_WMeanOld(new_spM,minbin,maxbin,struct_pair,asTimeInt.US,k);
[as_sel.Trial{k},MeanWind.Trial{k},NewStructPair.T{k},BinDir.T{k}] = AsActCode_WMean(Dir,new_spM,struct_pair,asTimeInt.US,k);

[as_sel.Base{k},MeanWind.Base{k},~,~] = AsActCode_WMean(Dir,new_spM,struct_pair,asTimeInt.USBase,k);
end
%% mean activity reward period in each trial
% lick='no'

for k=1:length(asTimeInt.RewPre)
[as_sel.RewPre{k},MeanWind.RewPre{k},NewStructPair.RewPre{k},BinDir.RewPre{k}] = AsActCode_WMean(Dir,new_spM,struct_pair,asTimeInt.RewPre,k);
end

for k=1:length(asTimeInt.RewPost)
% [as_selRewPost{k},MeanWind_RewPost{k},~,~] = AsActCode_WMeanWLi(new_spM,minbin,maxbin,struct_pair,asTime_RewPost,k,lick);
[as_sel.RewPost{k},MeanWind.RewPost{k},~,~] = AsActCode_WMean(Dir,new_spM,struct_pair,asTimeInt.RewPost,k);
end
%% Mean of each Lick period within the trial
% lick='yes';
minbin=0.01;
maxbin=1.6;
for k=1:length(asTimeInt.LickPre)
[as_sel.LickPre{k},MeanWind.LickPre{k},NewStructPair.LickPre{k},BinDir.LickPre{k}] = AsActCode_WMean(Dir,new_spM,struct_pair,asTimeInt.LickPre,k);
end

for k=1:length(asTimeInt.LickPost)
[as_sel.LickPost{k},MeanWind.LickPost{k},~,~] = AsActCode_WMean(Dir,new_spM,struct_pair,asTimeInt.LickPost,k);
end

%% Composition of the huge matrix
% for k=1:length(MeanWind.Trial)
%     for i=1:length(MeanWind.Trial{k}.act)
%         Mtx{k}{i}=[];
%     end
% end


for k=1:length(MeanWind.Trial)
    if ~isempty(MeanWind.Trial{k})
        for i=1:length(MeanWind.Trial{k}.act)
        
Mtx{k}{i}=[MeanWind.Trial{k}.code{i},MeanWind.Trial{k}.act{i},MeanWind.Base{k}.act{i},...
           MeanWind.LickPre{k}.act{i},MeanWind.LickPost{k}.act{i},...
           MeanWind.RewPre{k}.act{i},MeanWind.RewPost{k}.act{i}];
       [S,R]=sort(Mtx{k}{i}(:,1));
       Mtx{k}{i}=Mtx{k}{i}(R,:);
        end
    end
end
%%


save Mtx_Anova_VsVtaAllB.mat Mtx NewStructPair new_spM struct_pair
load Mtx_Anova_VsVtaAllB.mat

 %% From Matlab to excel/ SPSS
 temp=0;
for k=1:length(Mtx)
   if ~isempty(Mtx{k})
     for i=1:length(Mtx{k})
    
Code=(Mtx{k}{i}(:,1)); Trial=(Mtx{k}{i}(:,2)); Base=(Mtx{k}{i}(:,3));
LickPre=(Mtx{k}{i}(:,4)); LickPost=(Mtx{k}{i}(:,5));
RewPre=(Mtx{k}{i}(:,6)); RewPost=(Mtx{k}{i}(:,7));
% Pair = NewStructPair.T{k}.pair{i};
% Bin = NewStructPair.T{k}.bin{i};
% Lag = NewStructPair.T{k}.lag{i};
% Ind = NewStructPair.T{k}.ind{i};
% column = 2;
% columnLetters = char(ExcelCol(column));
% % Convert to Excel A1 format
% cellReference = sprintf('%s1', columnLetters);
% writetable(t1,'myfile.xlsx','Sheet',1,'Range', cellReference)



TableAn{k}{i} = table(Code,Trial,Base,LickPre,LickPost,RewPre,RewPost);
TableAn{k}{i}(:,:);

% B=[1 1 1];
filename = 'Mtx_Anova_VsVtaAllB.xlsx';
% Write the table
writetable(TableAn{k}{i},filename,'Sheet',i+temp);
      end
   end
temp=temp+length(Mtx{k});
end
