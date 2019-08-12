% This script intend to put toghether all the means
 addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM');
% The matrix new_spM is built in the program MaxDataPrepro.m
%  load('basedatabese.mat')


load('Lastrev1_OrRevW.mat')
load('AsPruOrRevW_lastrev1_rl2.mat')
load('DispLastrev1_Rl2_OrRevWPru.mat');
%%

% Region of interest, it could be: 'vsvta','vsvs','vtavta';
Region  ='all';

%Dir stays for "Directionality", can assume these values: 'vs->vta','vta->vs',0
Dir=0;

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




% Find indices in stable and unstable phase for the chapter that I took

% FindInd_US=[];   
% for k=1:TotAn
%     ll=0;
%     for l=startU(k):stopS(k)
%        
%           ll=ll+1;
%         FindInd_US{k}{ll}= find(binref.B{k}>=new_spM.events{k}(l).fv_on & binref.B{k}<=new_spM.events{k}(l).fv_off);
%         FindInd_USBase{k}{ll}= find(binref.B{k}>=new_spM.events{k}(l).fv_on & binref.B{k}<=new_spM.events{k}(l).fv_on+Base);
%         
%         
%     end
% end



 %% Big matrix all together
 
 prepost = 1; % one second befor of after the reward / or, if I need, before/after the first/last lick of the trial
[ new_spM.events] = LicksRewInt( new_spM.events,TotAn,prepost );

Base=0.7; % 0.7 sec after starting odor-onset is my baseline
 
  for k=1:TotAn
     for l=startU(k):stopS(k)
 FindInd.US{k}{l}=[];
 FindInd.USBase{k}{l}=[];
 FindInd.LickPre{k}{l}=[]; 
 FindInd.LickPost{k}{l}=[];
 FindInd.RewPre{k}{l}=[]; 
 FindInd.RewPost{k}{l}=[];
 FindInd.Code{k}{l}=[];
     end
  end
 lu=[];

 
 for k=1:TotAn
     for l=startU(k):stopS(k)
     lu(l)=length(new_spM.events{k}(l).licks);
     lr=length(new_spM.events{k}(l).reward_time);
             
             FindInd.US{k}{l}= find(binref.B{k}>=new_spM.events{k}(l).fv_on & binref.B{k}<=new_spM.events{k}(l).fv_off);
             FindInd.USBase{k}{l}= find(binref.B{k}>=new_spM.events{k}(l).fv_on & binref.B{k}<=new_spM.events{k}(l).fv_on+Base);
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
 
 %% Indices for the licks of the same 
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
 
%% 
%  end
 

% for k=1:TotAn
%     ll=0;
%         new_spM.events{k}=rmfield(new_spM.events{k}, 'lickDiff');
%           new_spM.events{k}=rmfield(new_spM.events{k}, 'lickPre');
%          new_spM.events{k}=rmfield(new_spM.events{k}, 'lickPost');
%         new_spM.events{k}= rmfield(new_spM.events{k}, 'rewardPre');
%         new_spM.events{k}= rmfield(new_spM.events{k}, 'rewardPost');
%          new_spM.events{k}=rmfield(new_spM.events{k}, 'lickFirstPre');
%          new_spM.events{k}=rmfield(new_spM.events{k}, 'lickLastPost');
% end

%% Find assemblies activity on indices
% IN FIND_IND I have 
[asTimeInt.US] = AsOnInd(AsTime,FindInd.US,new_spM,startU);
[asTimeInt.USBase] = AsOnInd(AsTime,FindInd.USBase,new_spM,startU);

[asTimeInt.RewPre] = AsOnInd(AsTime,FindInd.RewPre,new_spM,startU);
[asTimeInt.RewPost] = AsOnInd(AsTime,FindInd.RewPost,new_spM,startU);

[asTimeInt.LickPre] = AsOnInd(AsTime,FindInd.LickPreVect,new_spM,startU);
[asTimeInt.LickPost] = AsOnInd(AsTime,FindInd.LickPostVect,new_spM,startU);

% FindInd=FindInd_LickPre;
% for k=1:length(AsTime) % number of animals
%     for n=1:length(AsTime{k}) % number of pairs
%         for t=1: length(FindInd{k}) % number of trials
%             for l=1:length(FindInd{k}{t}) % number of licks in the trials
%  [asTime_LickPre{k}{n}{t}{l}] = AsOnIndwoFor(AsTime{k}{n},FindInd{k}{t}{l});
%  
%             end
%         end
%     end
% end
% 
% FindInd=FindInd_LickPost;
% for k=1:length(AsTime)
%     for n=1:length(AsTime{k})
%         for t=1: length(FindInd{k})
%             for l=1:length(FindInd{k}{t})
%  [asTime_LickPost{k}{n}{t}{l}] = AsOnIndwoFor(AsTime{k}{n},FindInd{k}{t}{l});
%             end
%         end
%     end
% end
% I divide the pairs in a struct coerently with the subset that I choose

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
minbin=0.01;
maxbin=1.6;
for k=1:length(asTimeInt.US)
[as_sel.Trial{k},MeanWind.Trial{k},NewStructPair.T{k},BinDir.T{k}] = AsActCode_WMean(new_spM,minbin,maxbin,struct_pair,asTimeInt.US,k);
[as_sel.Base{k},MeanWind.Base{k},~,~] = AsActCode_WMean(new_spM,minbin,maxbin,struct_pair,asTimeInt.USBase,k);
end
%% mean activity reward period in each trial
% lick='no'
minbin=0.01;
maxbin=1.6;
for k=1:length(asTimeInt.RewPre)
[as_sel.RewPre{k},MeanWind.RewPre{k},NewStructPair.RewPre{k},BinDir.RewPre{k}] = AsActCode_WMean(new_spM,minbin,maxbin,struct_pair,asTimeInt.RewPre,k);
end

for k=1:length(asTimeInt.RewPost)
% [as_selRewPost{k},MeanWind_RewPost{k},~,~] = AsActCode_WMeanWLi(new_spM,minbin,maxbin,struct_pair,asTime_RewPost,k,lick);
[as_sel.RewPost{k},MeanWind.RewPost{k},~,~] = AsActCode_WMean(new_spM,minbin,maxbin,struct_pair,asTimeInt.RewPost,k);
end
%% Mean of each Lick period within the trial
% lick='yes';
minbin=0.01;
maxbin=1.6;
for k=1:length(asTimeInt.LickPre)
[as_sel.LickPre{k},MeanWind.LickPre{k},NewStructPair.LickPre{k},BinDir.LickPre{k}] = AsActCode_WMean(new_spM,minbin,maxbin,struct_pair,asTimeInt.LickPre,k);
end

for k=1:length(asTimeInt.LickPost)
[as_sel.LickPost{k},MeanWind.LickPost{k},~,~] = AsActCode_WMean(new_spM,minbin,maxbin,struct_pair,asTimeInt.LickPost,k);
end
% %% Mean of mean licks period in the same trial
% MeanWind=MeanWind_LickPre;
% 
% for k=1:length(MeanWind) % k runs over the mice
% for a=1:length(MeanWind{k})% a runs over the assemblies
%  for j= 1:length(MeanWind{k}.act{a}) % j runs over the trials
%     
% 
%         MeanWindTrial_LickPre{k}.act{a}(j)=mean(MeanWind{k}.act{a}{j});
%           MeanWindTrial_LickPre{k}.code{a}(j)=MeanWind{k}.code{a}{j}(1);
%     
%    
%   end
% end
% end
% 
% MeanWind=MeanWind_LickPost;
% 
% for k=1:length(MeanWind) % k runs over the mice
% for a=1:length(MeanWind{k})% a runs over the assemblies
%  for j= 1:length(MeanWind{k}.act{a}) % j runs over the trials
%     
% 
%         MeanWindTrial_LickPost{k}.act{a}(j)=mean(MeanWind{k}.act{a}{j});
%         MeanWindTrial_LickPost{k}.code{a}(j)=MeanWind{k}.code{a}{j}(1);
%     
%    
%   end
% end
% end
%%

MTX=[MeanWind_SH{1}{1};MeanWind_SH{1}{2}];
MM=MTX';
% size(MTX)
% export2spss(vars,MTX', titl);
save MTX.mat MM
load MTX.mat
MeanA=MM(:,1);
MeanB=MM(:,2);
T = table(MeanA,MeanB);
T(:,1:2)

% B=[1 1 1];
filename = 'MTX.xlsx';
writetable(T,filename,'Sheet',1)


%%
% k=1;
% as_sel=as_selLickPre;
% for k=1:length(as_sel)
% for m=1:length(as_sel{k}.act) % m runs over the mice
% if ~isempty(as_sel{k}.act{m})
% %     bb=0;
%     for j= 1:length(as_sel{k}.act{m})
%         for l=1:length(as_sel{k}.act{m}{j})
% %         bb=bb+1;
%         MeanWind{k}.act{m}{j}(l)=mean(as_sel{k}.act{m}{j}{l}(:,2));
% %         MeanWind.code{m}{j}(l)=as_sel.code{m}{j}{l}(:);
%         end
%     end
% end
% end
% end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Trash but useful in case the function 
%  for k=1:TotAn
% %     ll=0;
% %     Int{k}=[];
%     for l=1:length(new_spM.events{k})
%        lu= length(new_spM.events{k}(l).licks);
%        
% %        new_spM.events{k}(l).lickIn=nan(size(new_spM.events{k}(l).licks));
%        new_spM.events{k}(l).lickIn=[];
%        if lu~=0
%        if lu>1
% %           ll=ll+1;
%           new_spM.events{k}(l).lickDiff=diff(new_spM.events{k}(l).licks);
%            elseif lu==1
%             new_spM.events{k}(l).lickDiff=0.3;
%        end
%           new_spM.events{k}(l).lickIn(1)=new_spM.events{k}(l).licks(1)-new_spM.events{k}(l).lickDiff(1)/2;
%           new_spM.events{k}(l).lickIn(2:lu)=new_spM.events{k}(l).licks(2:lu)-new_spM.events{k}(l).lickDiff(1:end)/2;
%        
%           new_spM.events{k}(l).lickEnd(1:lu-1)=new_spM.events{k}(l).licks(1:lu-1)+new_spM.events{k}(l).lickDiff(:)/2;
%           new_spM.events{k}(l).lickEnd(lu)=new_spM.events{k}(l).licks(end)+new_spM.events{k}(l).lickDiff(end)/2;
% %         FindInd_Lick{k}{ll}= find(binref.B{k}>=new_spM.events{k}(l).fv_on & binref.B{k}<=new_spM.events{k}(l).fv_off);
% %         FindInd_LickBase{k}{ll}= find(binref.B{k}>=new_spM.events{k}(l).fv_on & binref.B{k}<=new_spM.events{k}(l).fv_on+Base);
%        end
%     end
%     
%  end



% % %  %% Only rewarded
% % %  for k=1:TotAn
% % %      ll=0;
% % %     
% % %      for l=startU(k):stopS(k)
% % %      lu(l)=length(new_spM.events{k}(l).licks);
% % % %         if length(new_spM.events{k}(l).licks)~=0
% % %      lr=length(new_spM.events{k}(l).reward_time);
% % %           if length(new_spM.events{k}(l).reward_time)==1
% % %               ll=ll+1;
% % %              FindInd.Rew.RewPre{k}{l}= find(binref.B{k}>=new_spM.events{k}(l).rewardPre & binref.B{k}<=new_spM.events{k}(l).reward_time);
% % %              FindInd.Rew.RewPost{k}{l}= find(binref.B{k}>=new_spM.events{k}(l).reward_time & binref.B{k}<=new_spM.events{k}(l).rewardPost);
% % % 
% % %              FindInd.Rew.IndOld{k}{ll}=l;
% % %              for i=1:lu(l)
% % %              FindInd.Rew.LickPre{k}{l}{i}= find(binref.B{k}>=new_spM.events{k}(l).lickPre(i) & binref.B{k}<=new_spM.events{k}(l).licks(i));
% % %              FindInd.Rew.LickPost{k}{l}{i}= find(binref.B{k}>=new_spM.events{k}(l).licks(i) & binref.B{k}<=new_spM.events{k}(l).lickPost(i));
% % % %           FindInd_Lick.Ind{k}{ll}=l;
% % %              end
% % % %         end
% % % %         ri=0;
% % %         
% % % %         if length(new_spM.events{k}(l).reward_time)==1
% % % %             ri=ri+1;
% % %            
% % % %         end
% % %            end
% % %      end
% % %  end



%%%%%%%%%%%%%%%%%%%%%%%%%

%% Rewarded and not-rewarded separeted
 
% % %  FindInd.Rew.LickPre=[]; % rewarded
% % %  FindInd.Rew.LickPost=[];
% % %  FindInd.Rew.RewPre=[]; 
% % %  FindInd.Rew.RewPost=[];
% % %  FindInd.Rew.IndOld=[]; 
% % %  
% % %  
% % %  FindInd.Unrew_L.LickPre=[]; % unrewarded with lick
% % %  FindInd.Unrew_L.LickPost=[];
% % %  FindInd.Unrew_L.RewPre=[]; 
% % %  FindInd.Unrew_L.RewPost=[];
% % %  FindInd.Unrew_L.IndOld=[]; 
% % %  
% % %  FindInd.Unrew_N.LickPre=[]; % unrewarded No lick
% % %  FindInd.Unrew_N.LickPost=[];
% % %  FindInd.Unrew_N.RewPre=[]; 
% % %  FindInd.Unrew_N.RewPost=[];
% % %  FindInd.Unrew_N.IndOld=[]; 
% % %  
% % %  
% % %  lu=[];
% % %  
% % %  
% % %  
% % %  
% % %  
% % %  for k=1:TotAn
% % %      h=0;
% % %      f=0;
% % %      c=0;
% % %      for l=startU(k):stopS(k)
% % %      lu(l)=length(new_spM.events{k}(l).licks);
% % %      lr=length(new_spM.events{k}(l).reward_time);
% % %      
% % %           if new_spM.events{k}(l).rew_code==1 % rewarded case % HIT
% % %               
% % %               h=h+1;
% % %              FindInd.Hit.RewPre{k}{h}= find(binref.B{k}>=new_spM.events{k}(l).rewardPre & binref.B{k}<=new_spM.events{k}(l).reward_time);
% % %              FindInd.Hit.RewPost{k}{h}= find(binref.B{k}>=new_spM.events{k}(l).reward_time & binref.B{k}<=new_spM.events{k}(l).rewardPost);
% % %              FindInd.Hit.US{k}{h}= find(binref.B{k}>=new_spM.events{k}(l).fv_on & binref.B{k}<=new_spM.events{k}(l).fv_off);
% % %              FindInd.Hit.USBase{k}{h}= find(binref.B{k}>=new_spM.events{k}(l).fv_on & binref.B{k}<=new_spM.events{k}(l).fv_on+Base);
% % %              FindInd.Hit.IndOld{k}{h}=l;
% % %              for i=1:lu(l)
% % %              FindInd.Hit.LickPre{k}{h}{i}= find(binref.B{k}>=new_spM.events{k}(l).lickPre(i) & binref.B{k}<=new_spM.events{k}(l).licks(i));
% % %              FindInd.Hit.LickPost{k}{h}{i}= find(binref.B{k}>=new_spM.events{k}(l).licks(i) & binref.B{k}<=new_spM.events{k}(l).lickPost(i));
% % %              end
% % %             
% % %           elseif  new_spM.events{k}(l).rew_code==2  %nonrewarded mouse went, False Alarm
% % %               if length(lu(l))~=0  %non rewardedLick
% % %               f=f+1;
% % %               FindInd.FA.IndOld{k}{f}=l;
% % %               FindInd.FA.US{k}{f}= find(binref.B{k}>=new_spM.events{k}(l).fv_on & binref.B{k}<=new_spM.events{k}(l).fv_off);
% % %               FindInd.FA.USBase{k}{f}= find(binref.B{k}>=new_spM.events{k}(l).fv_on & binref.B{k}<=new_spM.events{k}(l).fv_on+Base);
% % %               for i=1:lu(l)
% % %              FindInd.FA.LickPre{k}{f}{i}= find(binref.B{k}>=new_spM.events{k}(l).lickPre(i) & binref.B{k}<=new_spM.events{k}(l).licks(i));
% % %              FindInd.FA.LickPost{k}{f}{i}= find(binref.B{k}>=new_spM.events{k}(l).licks(i) & binref.B{k}<=new_spM.events{k}(l).lickPost(i));
% % %               end
% % %               end
% % %               
% % %               
% % %            elseif  new_spM.events{k}(l).rew_code==3  %nonrewarded mouse sat quiet, Correct Rejection
% % %               
% % %               c=c+1;
% % %               FindInd.FA.IndOld{k}{c}=l;
% % %               FindInd.FA.US{k}{c}= find(binref.B{k}>=new_spM.events{k}(l).fv_on & binref.B{k}<=new_spM.events{k}(l).fv_off);
% % %               FindInd.FA.USBase{k}{c}= find(binref.B{k}>=new_spM.events{k}(l).fv_on & binref.B{k}<=new_spM.events{k}(l).fv_on+Base);
% % %               if length(lu(l))~=0  %non rewardedLick
% % %               for i=1:lu(l)
% % %              FindInd.FA.LickPre{k}{c}{i}= find(binref.B{k}>=new_spM.events{k}(l).lickPre(i) & binref.B{k}<=new_spM.events{k}(l).licks(i));
% % %              FindInd.FA.LickPost{k}{c}{i}= find(binref.B{k}>=new_spM.events{k}(l).licks(i) & binref.B{k}<=new_spM.events{k}(l).lickPost(i));
% % %               end
% % %               end
% % %               
% % %               elseif  new_spM.events{k}(l).rew_code==4  %rewarded mouse didn't go for it, MISS
% % %               m=m+1;
% % %               FindInd.Miss.IndOld{k}{m}=l;
% % %               FindInd.Miss.US{k}{m}= find(binref.B{k}>=new_spM.events{k}(l).fv_on & binref.B{k}<=new_spM.events{k}(l).fv_off);
% % %               FindInd.Miss.USBase{k}{m}= find(binref.B{k}>=new_spM.events{k}(l).fv_on & binref.B{k}<=new_spM.events{k}(l).fv_on+Base);
% % %               if length(lu(l))~=0  %in case there are licks
% % %               for i=1:lu(l)
% % %              FindInd.Miss.LickPre{k}{m}{i}= find(binref.B{k}>=new_spM.events{k}(l).lickPre(i) & binref.B{k}<=new_spM.events{k}(l).licks(i));
% % %              FindInd.Miss.LickPost{k}{m}{i}= find(binref.B{k}>=new_spM.events{k}(l).licks(i) & binref.B{k}<=new_spM.events{k}(l).lickPost(i));
% % %               end
% % %               end
% % % 
% % %            end
% % %      end
% % %  end
