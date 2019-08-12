

BinSizesPostPru = BinSizes(:,1:7);
cvs=1;
for i=1:3
    cvta=i+3;
[BinZeroMT{i}] = BinZeroLag(new_pairs_small,cvs,cvta)
end

cvs=2;
for i=1:3
    cvta=i+3;
[BinZeroFT{i}] = BinZeroLag(new_pairs_small,cvs,cvta)
end

for z =1 :length(BinZeroMT)
for i=1:length(BinSizesPostPru)
%     for j=1:length(BinZero)
   countBinZeroMT{z}{i}= find(BinZeroMT{z}(:)==BinSizesPostPru(i))
   CountZeroBinMT{z}(i)=length(countBinZeroMT{z}{i})
   
%     end
end
end


for z =1 :length(BinZeroFT)
for i=1:length(BinSizesPostPru)
%     for j=1:length(BinZero)
   countBinZeroFT{z}{i}= find(BinZeroFT{z}(:)==BinSizesPostPru(i))
   CountZeroBinFT{z}(i)=length(countBinZeroFT{z}{i})
   
%     end
end
end



xx=1:10;

figure(7); 
subplot(2,2,1);
 b=bar(xx,CountZeroBinMT{1},'r'); hold on;
 xticks(1:7);
xticklabels(BinSizesPostPru)
xlim([0,8])
ylim([0,6])
xlabel('\Delta (sec)')
ylabel('# Pairs')
%  b.FaceAlpha=0.5;
subplot(2,2,3);
 b1=bar(xx,CountZeroBinMT{2},'g'); hold on;
 xticks(1:10);
xticklabels(BinSizesPostPru)
xlim([0,8])
ylim([0,6])
xlabel('\Delta (sec)')
ylabel('# Pairs')
%  b1.FaceAlpha=0.5;
subplot(2,2,2);
 b2=bar(xx,CountZeroBinFT{1},'b'); hold on;
 xticks(1:10);
xticklabels(BinSizesPostPru)
xlim([0,8])
ylim([0,6])
xlabel('\Delta (sec)')
ylabel('# Pairs')
% ylim([0,4])
%  b2.FaceAlpha=0.5;
subplot(2,2,4);
 b3=bar(xx,CountZeroBinFT{2},'y'); hold on;
%  b3.FaceAlpha=0.5;
% bar([CountZeroBinMT{1};CountZeroBinMT{2};CountZeroBinFT{1};CountZeroBinFT{2}]); hold on;
xticks(1:10);
xticklabels(BinSizesPostPru)
xlim([0,8])
ylim([0,6])
xlabel('\Delta (sec)')
ylabel('# Pairs')
% ylim([0,1])
%% Intra region bin distribution
MergeDataGenPostPru_1;
MergeDataGenPostPru_Intra2;

pairs_intra_small=pairs_small;
pairs_intra_large=pairs_large;
pairs_intra = [pairs_intra_small;pairs_intra_large];
BinSizes = BinSizeCopy(:,1:10);
vs1=1; vs2=2;

[BinZeroMF{k},BinNonZeroMF{k}] = BinZeroLag(pairs_intra,vs1,vs2);

for k =1:length(Data_All.par)
    nneuMSN(k)=0;
    nneuFSI(k)=0;
end
for k =1:length(Data_All.par)
    for j=1:size(Data_All.par{k},2)
        [nneuMSN(k)] = countNeu(Data_All.par{k}(1,j).labels,1,nneuMSN(k));
      [nneuFSI(k)] = countNeu(Data_All.par{k}(1,j).labels,2,nneuFSI(k));
    end

end


for i=1:length(BinSizes)
   countBinZeroMF{i}= find(BinZeroMF(:)==BinSizes(i))
   CountZeroBinMF(i)=length(countBinZeroMF{i})/(sum(nneuMSN)*sum(nneuFSI));
    
end

for i=1:length(BinSizes)
   countBinNonZeroMF{i}= find(BinNonZeroMF(:)==BinSizes(i))
   CountNonZeroBinMF(i)=length(countBinNonZeroMF{i})/(sum(nneuMSN)*sum(nneuFSI));
end
%%%%%%% Zero and Non Zero Lag Toghether
CountAllMF=CountZeroBinMF+CountNonZeroBinMF;
xx=1:12;

figure(15); 

 b=bar(xx,CountAllMF,'r'); hold on;
 xticks(1:12);
xticklabels(BinSizes)
% xlim([0,8])
% ylim([0,6])
set(gca,'FontSize',12)
xlabel('\Delta (sec)')
ylabel('# Pairs')
%  b.FaceAlpha=0.5;


%%%%% Figures zero non zero separeted


figure(11); 

 b=bar(xx,CountZeroBinMF,'r'); hold on;
 xticks(1:10);
xticklabels(BinSizes)
% xlim([0,8])
% ylim([0,6])
set(gca,'FontSize',12)
xlabel('\Delta (sec)')
ylabel('# Pairs')
%  b.FaceAlpha=0.5;

figure(12); 

 b=bar(xx,CountNonZeroBinMF,'r'); hold on;
 xticks(1:10);
xticklabels(BinSizes)
% xlim([0,8])
% ylim([0,6])
set(gca,'FontSize',12)
xlabel('\Delta (sec)')
ylabel('# Pairs')
%  b.FaceAlpha=0.5;
%% Mean across animals
MergeDataGenPostPru_1;
MergeDataGenPostPru_Intra2;

pairs_intra=pairs_vsvs_small;
% pairs_intra_large=pairs_large;
% pairs_intra = [pairs_intra_small;pairs_intra_large];
BinSizeCopy=BinSizes;
clear BinSizes
BinSizes = BinSizeCopy(:,1:7);
vs1=1; vs2=2;

for k=1:SumTotAn
    if ~isempty(pairs_intra{k})
[BinZeroMF{k},BinNonZeroMF{k}] = BinZeroLag(pairs_intra{k},vs1,vs2);
    end
end


for k =1:length(Data_All.par)
    nneuMSN(k)=0;
    nneuFSI(k)=0;
end
for k =1:length(Data_All.par)
    for j=1:size(Data_All.par{k},2)
        [nneuMSN(k)] = countNeu(Data_All.par{k}(1,j).labels,1,nneuMSN(k));
      [nneuFSI(k)] = countNeu(Data_All.par{k}(1,j).labels,2,nneuFSI(k));
    end

end

for k=1:SumTotAn
for i=1:length(BinSizes)
   countBinZeroMF{k}{i}= find(BinZeroMF{k}(:)==BinSizes(i))
   CountZeroBinMF{k}(i)=length(countBinZeroMF{k}{i})./(nneuMSN(k)*nneuFSI(k));
    
end
end

for k=1:SumTotAn
for i=1:length(BinSizes)
   countBinNonZeroMF{k}{i}= find(BinNonZeroMF{k}(:)==BinSizes(i))
   CountNonZeroBinMF{k}(i)=length(countBinNonZeroMF{k}{i})./(nneuMSN(k)*nneuFSI(k));
end
end

for k=1:SumTotAn
CountAllMF{k}=CountZeroBinMF{k}
end
BinDistr=CountAllMF{1};
for k=1:SumTotAn-1
    BinDistr=[BinDistr;CountAllMF{k+1}];
end

BinDistrM = mean(BinDistr,'omitnan');
BinDistrSt = std(BinDistr,'omitnan')/sqrt(size(BinDistr,1));

figure(17);hold on;
bar(BinDistrM,'r'); hold on;
errorbar(BinDistrM,BinDistrSt/2,'k*');hold on;
 title({'Mean Bin Zero Lag distribution';nameT})
 ylabel('Assembly pairs/All possible pairs')
 xlabel('\Delta (sec)')
xticks(1:7);
xticklabels(BinSizes);
 box on

 clear CountAllMF BinDistr
 for k=1:SumTotAn
CountAllMF{k}=CountNonZeroBinMF{k}+CountZeroBinMF{k}
end
BinDistr=CountAllMF{1};
for k=1:SumTotAn-1
    BinDistr=[BinDistr;CountAllMF{k+1}];
end

BinDistrM = mean(BinDistr,'omitnan');
BinDistrSt = std(BinDistr,'omitnan')/sqrt(size(BinDistr,1));

figure(20);hold on;
bar(BinDistrM,'r'); hold on;
errorbar(BinDistrM,BinDistrSt/2,'k*');hold on;
 title({'Mean Bin distribution';nameT})
 ylabel('Assembly pairs/All possible pairs')
 xlabel('\Delta (sec)')
xticks(1:7);
xticklabels(BinSizes);
 box on