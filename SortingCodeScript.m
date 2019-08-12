%% reward code
for k=1:length(MaxA_Od)
    [MaxA_Tot{k}] = MaxMaxTr(MaxA_Od{k},MaxA_Lick{k},MaxA_Rew{k},MaxA_Od{k});
end

[MNorm_Od] = NormTotMax(MaxA_Tot,MeanTr_Od);
[MNorm_Lick] = NormTotMax(MaxA_Tot,MeanTr_Lick);
[MNorm_Rew] = NormTotMax(MaxA_Tot,MeanTr_Rew);

%%%%%%%% Hit/ Different sorting
[Act_Od,MaxI_O1,Sort_O1] = ActAsAllSameLengthInt(length(MNorm_Od), MNorm_Od);
[Act_Lick,MaxI_L1,Sort_L1] = ActAsAllSameLengthInt(length(MNorm_Lick), MNorm_Lick);
[Act_Rew,MaxI_R1,Sort_R1] = ActAsAllSameLengthInt(length(MNorm_Rew), MNorm_Rew); % Not active in false alarm and rejection


[ActSort_OdOd] = SortAct(Act_Od,Sort_O1);         % Sorted on the odour
[ActSort_LickOd] = SortAct(Act_Lick,Sort_O1);
[ActSort_RewOd] = SortAct(Act_Rew,Sort_O1);

[ActSort_OdLick] = SortAct(Act_Od,Sort_L1);        % sorted on the first lick
[ActSort_LickLick] = SortAct(Act_Lick,Sort_L1);
[ActSort_RewLick] = SortAct(Act_Rew,Sort_L1);

[ActSort_OdRew] = SortAct(Act_Od,Sort_R1);       % sorted on the reward
[ActSort_LickRew] = SortAct(Act_Lick,Sort_R1);   % Not active in false alarm and rejection
[ActSort_RewRew] = SortAct(Act_Rew,Sort_R1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% Miss reward code
% % % % % for k=1:length(MaxA_HOd)
% % % % %     [MaxAM_Tot{k}] = MaxMaxTr(MaxA_MOd{k},MaxA_MLick{k},MaxA_MRew{k},MaxA_MOd{k});
% % % % % end
% % % % % 
% % % % % [MNorm_MOd] = NormTotMax(MaxAM_Tot,MeanTr_MOd);
% % % % % [MNorm_MLick] = NormTotMax(MaxAM_Tot,MeanTr_MLick);
% % % % % [MNorm_MRew] = NormTotMax(MaxAM_Tot,MeanTr_MRew);
% % % % % 
% % % % % %%%%%%%% Hit/ Different sorting
% % % % % [ActM_Od,MaxI_MO1,SortM_O1] = ActAsAllSameLengthInt(length(MNorm_MOd), MNorm_MOd);
% % % % % [ActM_Lick,MaxI_ML1,SortM_L1] = ActAsAllSameLengthInt(length(MNorm_MLick), MNorm_MLick);
% % % % % [ActM_Rew,MaxI_MR1,SortM_R1] = ActAsAllSameLengthInt(length(MNorm_MRew), MNorm_MRew);
% % % % % 
% % % % % 
% % % % % [ActSortM_OdOd] = SortAct(ActM_Od,SortM_O1);         % Sorted on the odour
% % % % % % [ActSortM_LickOd] = SortAct(ActM_Lick,SortM_O1);
% % % % % % [ActSortM_RewOd] = SortAct(ActM_Rew,SortM_O1);
% % % % % % 
% % % % % % [ActSortM_OdLick] = SortAct(ActM_Od,SortM_L1);        % sorted on the first lick
% % % % % [ActSortM_LickLick] = SortAct(ActM_Lick,SortM_L1);
% % % % % % [ActSortM_RewLick] = SortAct(ActM_Rew,SortM_L1);
% % % % % % 
% % % % % % [ActSortM_OdRew] = SortAct(ActM_Od,SortM_R1);       % sorted on the reward
% % % % % % [ActSortM_LickRew] = SortAct(ActM_Lick,SortM_R1);
% % % % % [ActSortM_RewRew] = SortAct(ActM_Rew,SortM_R1);
