%% reward code
for k=1:length(MaxA_Od)
    [MaxA_Tot{k}] = MaxMaxTr(MaxA_Od{k},MaxA_Lick{k},MaxA_SecLick{k},MaxA_Od{k});
end

[MNorm_Od] = NormTotMax(MaxA_Tot,MeanTr_Od);
[MNorm_Lick] = NormTotMax(MaxA_Tot,MeanTr_Lick);
[MNorm_SecLick] = NormTotMax(MaxA_Tot,MeanTr_SecLick);

%%%%%%%% Hit/ Different sorting
[Act_Od,MaxI_O1,Sort_O1] = ActAsAllSameLengthInt(length(MNorm_Od), MNorm_Od);
[Act_Lick,MaxI_L1,Sort_L1] = ActAsAllSameLengthInt(length(MNorm_Lick), MNorm_Lick);
[Act_SecLick,MaxI_L2,Sort_L2] = ActAsAllSameLengthInt(length(MNorm_SecLick), MNorm_SecLick);


[ActSort_OdOd] = SortAct(Act_Od,Sort_O1);         % Sorted on the odour
[ActSort_LickOd] = SortAct(Act_Lick,Sort_O1);
[ActSort_SecLickOd] = SortAct(Act_SecLick,Sort_O1);

[ActSort_OdLick] = SortAct(Act_Od,Sort_L1);        % sorted on the first lick
[ActSort_LickLick] = SortAct(Act_Lick,Sort_L1);
[ActSort_SecLickLick] = SortAct(Act_SecLick,Sort_L1);


[ActSort_OdSecLick] = SortAct(Act_Od,Sort_L2);        % sorted on the second lick
[ActSort_LickSecLick] = SortAct(Act_Lick,Sort_L2);
[ActSort_SecLickSecLick] = SortAct(Act_SecLick,Sort_L2);
