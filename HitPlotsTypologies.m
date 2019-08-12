
clear all;
MergeDataGenPostPru_1  %% Script that merges data 

MergeDataGenPostPru_2  %%% pairs small takes in account the directionality vs->vta (be careful if you want to use the othere directionality)
%%
%%%%%% Indices for start and end

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

%%%%% For large Bins
% As_activity_All = As_activity_largeBins_All;
% struct_pair =struct_pair_large;


 FindIndicesScript   %% Script that find indices in assemblies activity after the time rescaling


Region='vsvta';

%Dir stays for "Directionality", can assume these values
%'vs->vta','vta->vs','small', 'large','smallSync','largeSync'0

% Dir='vs->vta';
% Dir='vta->vs';
%    Dir=0;
Dir='small';
% Dir='large';
% Dir = 'largeSync';
% Dir = 'smallSync';

hit = 1;
rej = 3;
fal = 2;
miss = 4;


cvs1=1;  cvs2=2; cvs3=3;
cvta = [4, 5,6, 7.1, 0, 7.3,20];
col1 = 9;
col2= [10, 10, 10, 12, 12, 12, 10];
for k=1:length(As_activity_All)
    if ~isempty(pairs_vsvta_small{k})
for i =1:length(cvta)
    
[pairs_MT{i}{k}] = pairs_NeuSelection(pairs_vsvta_small{k},cvs1,cvta(i),col1,col2(i));
[pairs_FT{i}{k}] = pairs_NeuSelection(pairs_vsvta_small{k},cvs2,cvta(i),col1,col2(i));
[pairs_CT{i}{k}] = pairs_NeuSelection(pairs_vsvta_small{k},cvs3,cvta(i),col1,col2(i));
end
% [pairs_MT{4}{k}] = pairs_NeuSelection(pairs_vsvta_small{k},cvs1,20);
% [pairs_FT{4}{k}] = pairs_NeuSelection(pairs_vsvta_small{k},cvs2,20);
% [pairs_CT{4}{k}] = pairs_NeuSelection(pairs_vsvta_small{k},cvs3,20);
end
end

pairs = pairs_MT{1};
%    pairs =pairs_vsvta_small;
% pairs =pairs_vsvta_large;
%%
% for i=1:length(pairs_MT{4})
% pairs{i} =[pairs_FT{4}{i};pairs_MT{4}{i}]
% end

%%
  
  pairs_Mat=pairs{1};
for i = 1: length(pairs)-1
    if ~isempty(pairs{i+1})
pairs_Mat=[pairs_Mat;pairs{i+1}];
    end
end
clear i

As_TaskCodePruScript   % Assembly activity divided by reward code and mean over trials
%%
%  [as_sel_HLick1] = Comp_RewCode(as_sel_HLick,as_sel_HOd);
%  for k=1:length(as_sel_HLick1)
%  AsForMean{k} = MeanNewDim(as_sel_HLick1{k});
%  end
% [as_sel_MRew1] = Comp_RewCode(as_sel_MRew,as_sel_MOd);
% [as_sel_MLick1] = Comp_RewCode(as_sel_MLick,as_sel_MOd);
% [as_sel_FLick1] = Comp_RewCode(as_sel_FLick,as_sel_FOd);
% [as_sel_RLick1] = Comp_RewCode(as_sel_RLick,as_sel_ROd);
% 
% 
% [as_sel_MLick2] = Comp_RewCode(as_sel_Sec_MLick,as_sel_MOd);
% [as_sel_FLick2] = Comp_RewCode(as_sel_Sec_FLick,as_sel_FOd);
% [as_sel_RLick2] = Comp_RewCode(as_sel_Sec_RLick,as_sel_ROd);
%  for k=1:length(as_sel_HRew)
%     if ~isempty(as_sel_MOd{k})
%         if ~isempty(as_sel_MOd{k}.asAct)
% [as_sel_MRewNew{k}.asAct] = Comp_RewCode(as_sel_MRew{k}.asAct,as_sel_MOd{k}.asAct);
% [as_sel_MLickNew{k}.asAct] = Comp_RewCode(as_sel_MLick{k}.asAct,as_sel_MOd{k}.asAct);
%         end
%     end
%  end


 %%%%%%%%%% Hit Reward code
 
 
         
 MaxA_Od = MaxA_HOd; MaxA_Lick=MaxA_HLick; MaxA_Rew=MaxA_HRew; 
 MeanTr_Od = MeanTr_HOd; MeanTr_Lick= MeanTr_HLick; MeanTr_Rew = MeanTr_HRew;
% % 
 %%%%%%%%%%%%%%%%%%%%% For Original
  SortingCodeScript  %%%%%% Sort the Things
%   save Sort_MSNNoId_Or_Firstrev1_short Sort_O1 Sort_L1 Sort_R1
%   save Sort_MSNNoId_Or_Lastrev1_short Sort_O1 Sort_L1 Sort_R1
% save Sort_MSNNoId_Or_Rev2_short Sort_O1 Sort_L1 Sort_R1
% save Sort_FSINoId_Or_Rev2_short Sort_O1 Sort_L1 Sort_R1
% save Sort_FSITypeII_Or_Rev2_short Sort_O1 Sort_L1 Sort_R1
% save Sort_MSNTypeII_Or_Rev3_short Sort_O1 Sort_L1 Sort_R1
%   save Sort_MSNNoId_NewOr_3Data Sort_O1 Sort_L1 Sort_R1
% %   save Sort_MSNNoId_DisFirst_Data Sort_O1 Sort_L1 Sort_R1
% % save Sort_MSNTypeI_NewOr_3Data Sort_O1 Sort_L1 Sort_R1
%  % save Sort_MSNTypeI_NewOr_4Data Sort_O1 Sort_L1 Sort_R1
% save Sort_MSNTypeI_BinRed_Or_3Data Sort_O1 Sort_L1 Sort_R1




%%%%%%%%%%%%%%%%% For Reversal
% % % load Sort_MSNTypeI_Or_5Data
% % % load Sort_MSNTypeI_NewOr_4Data
% % % load Sort_MSNTypeII_Or_4Data
% % % load Sort_MSNTypeI_NewOr_3Data
% % %  load Sort_MSNTypeII_NewOr_3Data
% % % load Sort_FSITypeI_NewOr_3Data
% % % load Sort_FSITypeII_NewOr_3Data
% % % load Sort_MSNNoId_NewOr_3Data
% % % load Sort_FSINoId_NewOr_3Data
% % load Sort_FSINoId_Or_Firstrev1_short Sort_O1 Sort_L1 Sort_R1
% % load Sort_FSITypeI_Or_Firstrev1_short Sort_O1 Sort_L1 Sort_R1
% % for k=1:length(MaxA_Od)
% %     [MaxA_Tot{k}] = MaxMaxTr(MaxA_Od{k},MaxA_Lick{k},MaxA_Rew{k},MaxA_Od{k});
% % end
% % 
% % [MNorm_Od] = NormTotMax(MaxA_Tot,MeanTr_Od);
% % [MNorm_Lick] = NormTotMax(MaxA_Tot,MeanTr_Lick);
% % [MNorm_Rew] = NormTotMax(MaxA_Tot,MeanTr_Rew);
% % 
% % %%%%%%%% Hit/ Different sorting
% % [Act_Od,~,~] = ActAsAllSameLengthInt(length(MNorm_Od), MNorm_Od);
% % [Act_Lick,~,~] = ActAsAllSameLengthInt(length(MNorm_Lick), MNorm_Lick);
% % [Act_Rew,~,~] = ActAsAllSameLengthInt(length(MNorm_Rew), MNorm_Rew);
% % 
% % [ActSort_OdOd] = SortAct(Act_Od,Sort_O1);         % Sorted on the odour
% % [ActSort_LickOd] = SortAct(Act_Lick,Sort_O1);
% % [ActSort_RewOd] = SortAct(Act_Rew,Sort_O1);
% % 
% % [ActSort_OdLick] = SortAct(Act_Od,Sort_L1);        % sorted on the first lick
% % [ActSort_LickLick] = SortAct(Act_Lick,Sort_L1);
% % [ActSort_RewLick] = SortAct(Act_Rew,Sort_L1);
% % 
% % [ActSort_OdRew] = SortAct(Act_Od,Sort_R1);       % sorted on the reward
% % [ActSort_LickRew] = SortAct(Act_Lick,Sort_R1);   % Not active in false alarm and rejection
% % [ActSort_RewRew] = SortAct(Act_Rew,Sort_R1);
% 
%%%%%%%%%%%%%%%%%%% end of reversal paradigma %%%%%%
 Act_HLick = Act_Lick; Act_HOd = Act_Od; Act_HRew = Act_Rew;
 ActSort_HOdOd =  ActSort_OdOd;  ActSort_HOdRew =  ActSort_OdRew;  ActSort_HOdLick = ActSort_OdLick;
 ActSort_HRewOd =  ActSort_RewOd;  ActSort_HRewRew =  ActSort_RewRew;  ActSort_HRewLick = ActSort_RewLick;
 ActSort_HLickOd =  ActSort_LickOd;  ActSort_HLickRew =  ActSort_LickRew;  ActSort_HLickLick = ActSort_LickLick;
 
 
NoVect=nan(size(ActSort_HOdOd,1),2);
 %%%% Hit
[ActH_LickOdRew_Od]=[ActSort_HOdOd,NoVect,ActSort_HLickOd,NoVect,ActSort_HRewOd]; % ordered on the Odour
[ActH_LickOdRew_Lick]=[ActSort_HOdLick,NoVect,ActSort_HLickLick,NoVect,ActSort_HRewLick]; % ordered on the first lick
[ActH_LickOdRew_Rew]=[ActSort_HOdRew,NoVect,ActSort_HLickRew,NoVect,ActSort_HRewRew]; % ordered on the Odour

[ActH_LickOdRew_ES]=[ActSort_HOdOd,NoVect,ActSort_HLickLick,NoVect,ActSort_HRewRew]; % each panel sorted

Act1=ActH_LickOdRew_Od;
Act2=ActH_LickOdRew_Lick;
Act3=ActH_LickOdRew_Rew;
Act4=ActH_LickOdRew_ES;

% 
% 
% 
% num=1;
binwidthType = 'small';
TrialsT ='Hit';
HitPlotVisualization
 
 
 %% %%%%%%%%%%%% Miss Reward Code
 
%  [AsForMean_MRewNew,Mean_MRewNew,MNorm_MRewNew,NewStd_MRewNew,MaxA_MRewNew] = Mean_AsSel(as_sel_MRew1);
%  [AsForMean_MLickNew,Mean_MLickNew,MNorm_MLickNew,Std_MLickNew,MaxA_MLickNew] = Mean_AsSel(as_sel_MLick1);
%  MaxA_Od = MaxA_MOd; MaxA_Lick=MaxA_MLickNew; MaxA_Rew=MaxA_MRewNew; 
%  MeanTr_Od = MeanTr_MOd; MeanTr_Lick= Mean_MLickNew; MeanTr_Rew = Mean_MRewNew;
%  
%  
% %%%%%% If I sort the miss on the itselves indices  
% % % SortingCodeScript  %%%%%% Sort the Things
% % %  
% % %  Act_MLick = Act_Lick; Act_MOd = Act_Od; Act_MRew = Act_Rew;
% % %  ActSort_MOdOd =  ActSort_OdOd;  ActSort_MOdRew =  ActSort_OdRew;  ActSort_MOdLick = ActSort_OdLick;
% % %  ActSort_MRewOd =  ActSort_RewOd;  ActSort_MRewRew =  ActSort_RewRew;  ActSort_MRewLick = ActSort_RewLick;
% % %  ActSort_MLickOd =  ActSort_LickOd;  ActSort_MLickRew =  ActSort_LickRew;  ActSort_MLickLick = ActSort_LickLick;
% % %%%%%% end of sorting on itselves indeces 
%  
% 
%  load Sort_MSNTypeI_Or_5Data
%  for k=1:length(MaxA_Od)
%     [MaxA_Tot{k}] = MaxMaxTr(MaxA_Od{k},MaxA_Lick{k},MaxA_Rew{k},MaxA_Od{k});
% end
% 
% [MNorm_Od] = NormTotMax(MaxA_Tot,MeanTr_Od);
% [MNorm_Lick] = NormTotMax(MaxA_Tot,MeanTr_Lick);
% [MNorm_Rew] = NormTotMax(MaxA_Tot,MeanTr_Rew);
% 
% %%%%%%%% Miss/ Different sorting
% [Act_Od,~,~] = ActAsAllSameLengthInt(length(MNorm_Od), MNorm_Od);
% [Act_Lick,~,~] = ActAsAllSameLengthInt(length(MNorm_Lick), MNorm_Lick);
% [Act_Rew,~,~] = ActAsAllSameLengthInt(length(MNorm_Rew), MNorm_Rew);
% 
% 
% 
% [ActSort_OdOd] = SortAct(Act_Od,Sort_O1);         % Sorted on the odour
% [ActSort_LickOd] = SortAct(Act_Lick,Sort_O1);
% [ActSort_RewOd] = SortAct(Act_Rew,Sort_O1);
% 
% [ActSort_OdLick] = SortAct(Act_Od,Sort_L1);        % sorted on the first lick
% [ActSort_LickLick] = SortAct(Act_Lick,Sort_L1);
% [ActSort_RewLick] = SortAct(Act_Rew,Sort_L1);
% 
% [ActSort_OdRew] = SortAct(Act_Od,Sort_R1);       % sorted on the reward
% [ActSort_LickRew] = SortAct(Act_Lick,Sort_R1);   % Not active in false alarm and rejection
% [ActSort_RewRew] = SortAct(Act_Rew,Sort_R1);
% 
%  Act_MLick = Act_Lick; Act_MOd = Act_Od; Act_MRew = Act_Rew;
%  ActSort_MOdOd =  ActSort_OdOd;  ActSort_MOdRew =  ActSort_OdRew;  ActSort_MOdLick = ActSort_OdLick;
%  ActSort_MRewOd =  ActSort_RewOd;  ActSort_MRewRew =  ActSort_RewRew;  ActSort_MRewLick = ActSort_RewLick;
%  ActSort_MLickOd =  ActSort_LickOd;  ActSort_MLickRew =  ActSort_LickRew;  ActSort_MLickLick = ActSort_LickLick;
%  
% %%%%% Miss
% NoVect=nan(size(ActSort_MOdOd,1),2);
% 
% 
% [ActM_LickOdRew_Od]=[ActSort_MOdOd,NoVect,ActSort_MLickOd,NoVect,ActSort_MRewOd]; % ordered on the Odour
% [ActM_LickOdRew_Lick]=[ActSort_MOdLick,NoVect,ActSort_MLickLick,NoVect,ActSort_MRewLick]; % ordered on the first lick
% [ActM_LickOdRew_Rew]=[ActSort_MOdRew,NoVect,ActSort_MLickRew,NoVect,ActSort_MRewRew]; % ordered on the Odour
% 
% [ActM_LickOdRew_ES]=[ActSort_MOdOd,NoVect,ActSort_MLickLick,NoVect,ActSort_MRewRew]; % each panel sorted
% 
% % [ActM_LickOdS1]=[ActSortM_Od1,NoVect,ActSortM_Lick1]; % Ordered on Lick
% % ActM_LickOdS1(:,end+1)=nan;
% % 
% % [ActM_LickOdS2]=[ActSortM_Od22,NoVect,ActSortM_Lick22]; % Ordered on second odour (Hit)
% % ActM_LickOdS2(:,end+1)=nan;
% 
% Act1=ActM_LickOdRew_Od;
% Act2=ActM_LickOdRew_Lick;
% Act3=ActM_LickOdRew_Rew;
% Act4=ActM_LickOdRew_ES;
% 
% num=5;
% binwidthType = 'small';
% TrialsT ='Miss';
% HitPlotVisualization
%% %%%%%%%%% False Alarm Reward Code
% % 
 [AsForMean_FLickNew,Mean_FLickNew,MNorm_FLickNew,Std_FLickNew,MaxA_FLickNew] = Mean_AsSel(as_sel_FLick1);
[AsForMean_FSecLickNew,Mean_FSecLickNew,MNorm_FSecLickNew,Std_FSecLickNew,MaxA_FSecLickNew] = Mean_AsSel(as_sel_FLick2);
 MaxA_Od = MaxA_FOd; MaxA_Lick=MaxA_FLickNew; MaxA_SecLick=MaxA_FSecLickNew;
 MeanTr_Od = MeanTr_FOd; MeanTr_Lick= Mean_FLickNew; MeanTr_SecLick= Mean_FSecLickNew; 
% % % % % %  
SortingCodeScriptFalseRej  %%%%%% Sort the Things
 
% % %%%%%%%%%%% reward code
% % for k=1:length(MaxA_Od)
% %     [MaxA_Tot{k}] = MaxMaxTr(MaxA_Od{k},MaxA_Lick{k},MaxA_SecLick{k},MaxA_Od{k});
% % end
% % 
% % [MNorm_Od] = NormTotMax(MaxA_Tot,MeanTr_Od);
% % [MNorm_Lick] = NormTotMax(MaxA_Tot,MeanTr_Lick);
% % [MNorm_SecLick] = NormTotMax(MaxA_Tot,MeanTr_SecLick);


%%%%%%%% Hit/ Different sorting
% % [Act_Od,~,~] = ActAsAllSameLengthInt(length(MNorm_Od), MNorm_Od);
% % [Act_Lick,~,~] = ActAsAllSameLengthInt(length(MNorm_Lick), MNorm_Lick);
% % [Act_SecLick,~,~] = ActAsAllSameLengthInt(length(MNorm_SecLick), MNorm_SecLick);

% % load Sort_MSNTypeI_DisFirst_Data
% % load Sort_FSITypeI_DisFirst_Data
% % load Sort_FSITypeII_DisFirst_Data

% % load Sort_MSNTypeI_Or_3Data
% % %load Sort_MSNTypeI_NewOr_3Data
% % % load Sort_MSNTypeII_NewOr_3Data
% % %load Sort_FSITypeI_NewOr_3Data
% % % load Sort_FSITypeII_NewOr_3Data
% % % load Sort_MSNNoId_Or_3Data
% % %load Sort_FSINoId_NewOr_3Data
% load Sort_MSNNoId_NewOr_3Data

% % [ActSort_OdOd] = SortAct(Act_Od,Sort_O1);         % Sorted on the odour
% % [ActSort_LickOd] = SortAct(Act_Lick,Sort_O1);
% % [ActSort_SecLickOd] = SortAct(Act_SecLick,Sort_O1);
% % 
% % 
% % [ActSort_OdLick] = SortAct(Act_Od,Sort_L1);        % sorted on the first lick
% % [ActSort_LickLick] = SortAct(Act_Lick,Sort_L1);
% % [ActSort_SecLickLick] = SortAct(Act_SecLick,Sort_L1);
% % 
% % [ActSort_OdSecLick] = SortAct(Act_Od,Sort_L2);        % sorted on the second lick
% % [ActSort_LickSecLick] = SortAct(Act_Lick,Sort_L2);
% % [ActSort_SecLickSecLick] = SortAct(Act_SecLick,Sort_L2);


 Act_FLick = Act_Lick; Act_FOd = Act_Od; 
% %  ActSort_FOdOd =  ActSort_OdOd; ActSort_FOdLick = ActSort_OdLick; ActSort_FOdLick = ActSort_OdLick;
% %  ActSort_FLickOd =  ActSort_LickOd; ActSort_FLickLick = ActSort_LickLick;
% %  ActSort_FSecLickOd =  ActSort_SecLickOd; ActSort_FsecLickLick = ActSort_SecLickLick;
 


NoVect=nan(size(ActSort_OdOd,1),2);


%%%%% False Alarm
[ActF_LickOd_Od]=[ActSort_OdOd,NoVect,ActSort_LickOd,NoVect,ActSort_SecLickOd]; % ordered on the Odour
[ActF_LickOd_Lick]=[ActSort_OdLick,NoVect,ActSort_LickLick,NoVect,ActSort_SecLickLick]; % ordered on the first lick
[ActF_LickOd_SecLick]=[ActSort_OdSecLick,NoVect,ActSort_LickSecLick,NoVect,ActSort_SecLickSecLick]; % ordered on the Second lick

% [ActF_LickOd_ES]=[ActSort_FOdOd,NoVect,ActSort_FLickLick]; % each panel sorted

Act1=ActF_LickOd_Od;
Act2=ActF_LickOd_Lick;
Act3=ActF_LickOd_SecLick;


num=5;
binwidthType = 'small';
TrialsT ='False Alarm';
HitPlotVisualization
% HeatPlotVisu_FalseRej
%% %%%%%%%%%%%%%%%%%% Correct rejection reward code


% [AsForMean_RLickNew,Mean_RLickNew,MNorm_RLickNew,Std_RLickNew,MaxA_RLickNew] = Mean_AsSel(as_sel_RLick1);
%  MaxA_Od = MaxA_ROd; MaxA_Lick=MaxA_RLickNew; 
%  MeanTr_Od = MeanTr_ROd; MeanTr_Lick= Mean_RLickNew; 
%  
%  
%  
% %%% Original phase
% SortingCodeScriptFalseRej  %%%%%% Sort the Things
% Act_RLick = Act_Lick; Act_ROd = Act_Od; 
% ActSort_ROdOd =  ActSort_OdOd; ActSort_ROdLick = ActSort_OdLick;
% ActSort_RLickOd =  ActSort_LickOd; ActSort_RLickLick = ActSort_LickLick;
 


%%%%%%%%%%%%%%%%%%%%%reversal
%%% when i use reversal phase and i want to sort the pairs on the original
%%% phase of the hit


% % load Sort_MSNTypeI_Or_5Data
% % %%%%%%%% reward code
% % for k=1:length(MaxA_Od)
% %     [MaxA_Tot{k}] = MaxMaxTr(MaxA_Od{k},MaxA_Lick{k},MaxA_Od{k},MaxA_Od{k});
% % end
% % 
% % [MNorm_Od] = NormTotMax(MaxA_Tot,MeanTr_Od);
% % [MNorm_Lick] = NormTotMax(MaxA_Tot,MeanTr_Lick);


%%%%%%%% Hit/ Different sorting
% % [Act_Od,~,~] = ActAsAllSameLengthInt(length(MNorm_Od), MNorm_Od);
% % [Act_Lick,~,~] = ActAsAllSameLengthInt(length(MNorm_Lick), MNorm_Lick);
% % 
% % 
% % [ActSort_OdOd] = SortAct(Act_Od,Sort_O1);         % Sorted on the odour
% % [ActSort_LickOd] = SortAct(Act_Lick,Sort_O1);
% % 
% % 
% % [ActSort_OdLick] = SortAct(Act_Od,Sort_L1);        % sorted on the first lick
% % [ActSort_LickLick] = SortAct(Act_Lick,Sort_L1);
% % 
% %  
% %  Act_RLick = Act_Lick; Act_ROd = Act_Od; 
% %  ActSort_ROdOd =  ActSort_OdOd; ActSort_ROdLick = ActSort_OdLick;
% %  ActSort_RLickOd =  ActSort_LickOd; ActSort_RLickLick = ActSort_LickLick;
% %  
%%%%%%%%%%%%%%%% end of reversal case

% NoVect=nan(size(ActSort_ROdOd,1),2);
% 
% 
% %%%%% Correct rejection
% [ActR_LickOd_Od]=[ActSort_ROdOd,NoVect,ActSort_RLickOd]; % ordered on the Odour
% [ActR_LickOd_Lick]=[ActSort_ROdLick,NoVect,ActSort_RLickLick]; % ordered on the first lick
% 
% [ActR_LickOd_ES]=[ActSort_ROdOd,NoVect,ActSort_RLickLick]; % each panel sorted
% 
% Act1=ActR_LickOd_Od;
% Act2=ActR_LickOd_Lick;
% Act3=ActR_LickOd_ES;
% 
% % ActH1(:,end+1)=nan;
% % ActH2(:,end+1)=nan;
% % ActH3(:,end+1)=nan;
% % ActH4(:,end+1)=nan;
% num=4;
% binwidthType = 'small';
% TrialsT ='Correct Rejection';
% HeatPlotVisu_FalseRej
%%
% Act1=ActSortM_RewRew;
% num=2;
% figure(fignum(num));hold on;
% pp=pcolor(Act1);hold on;
% xlim([1,size(Act1,2)])
% ylim([1,size(Act1,1)])
% xticks(tick)
% xticklabels(ticklabel)
% text(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-10),-18,'Time(sec)')
% text(find(gvect==-10),-18,'Time(sec)')
% text(2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==-10),-18,'Time(sec)')
% text(find(gvect==-10),size(Act_Lick,1)+10,'Odor onset','Fontsize',12)
% text(size(Act_Lick,2)+size(NoVect,2)+find(gvect==-10),size(Act_Lick,1)+10,'First Lick','Fontsize',12)
% text(2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect==-10),size(Act_Lick,1)+10,'Reward','Fontsize',12)
% plot([size(Act_Lick,2)+size(NoVect,2)+find(gvect==0),size(Act_Lick,2)+size(NoVect,2)+find(gvect==0)],[0,size(Act2,1)],'-.r','Linewidth',1);
% plot([2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect==0),2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect==0)],[0,size(Act2,1)],'-.r','Linewidth',1);
% plot([find(gvect==0),find(gvect==0)],[0,size(Act2,1)],'-.r','Linewidth',1);
% text(find(gvect==-70)-10,size(ActH_Lick,1)/2,'Sorted on Odour','Rotation',90)
% vline(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==0),'-.r');
% vline(2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==0),'-.r');
% vline(find(gvect==0),'-.r');
% set(pp, 'EdgeColor', 'none');
% grid on;
% ylabel('Pairs Actvity')
% title('Miss-reward -')