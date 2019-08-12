%%%%%%  %%%%%%%%%%%%

%%%%% Bin intra-regional
path = '/zifnas/Carla/CWEM_Project_GenProg'
cd(path)
%  Merge_Pru_IntraR;
% MergeDataGenPostPru_1;
% MergeDataGenPostPru_Intra2;
clear i j l
%% Sessions togheter
for k =1:length(Data_All.par)
    nneuMSN(k)=0;
    nneuFSIlow(k)=0;
    nneuFSIhigh(k)=0;
end
for k =1:length(Data_All.par)
    for j=1:size(Data_All.par{k},2)
        [nneuMSN(k)] = countNeu(Data_All.par{k}(1,j).labels,1,nneuMSN(k));
      [nneuFSIlow(k)] = countNeu(Data_All.par{k}(1,j).labels,2.1,nneuFSIlow(k));
      [nneuFSIhigh(k)] = countNeu(Data_All.par{k}(1,j).labels,2.2,nneuFSIhigh(k));
    end

end

% pairs_vsvs=pairs_vsvs_small
col1=9;col2=10;
for k =1:SumTotAn
    pairs_smallMFlow{k}=[];
    pairs_smallMFhigh{k}=[];
end
for k =1:SumTotAn
    if~isempty(pairs_vsvs{k})
    
    [pairs_smallMFlow{k}] = pairs_NeuSelection(pairs_vsvs{k},1,2.1,col1,col2);
    [pairs_smallMFhigh{k}] = pairs_NeuSelection(pairs_vsvs{k},1,2.2,col1,col2);
    end
end
%%
for k=1:SumTotAn
    count{k}=[];
end

for k=1:SumTotAn
    if~isempty(pairs_smallMFlow{k})
for i =1:length(BinSizes)
Ind{k}{i}=find(pairs_smallMFlow{k}(:,5)==BinSizes(i) & pairs_smallMFlow{k}(:,3)~=0)
% count{k}(i) = length(Ind{k}{i})./(sum(nneuFSIlow)*sum(nneuMSN));
 count{k}(i) = length(Ind{k}{i})./(nneuFSIlow(k)*nneuMSN(k)); %depending on the couplings that you are considering
% count{k}(i) = length(Ind{k}{i})./(nneuFSIlow(k)*nneuFSIlow(k)); %fsi-fsi low
end
    end
end


for k=1:SumTotAn
    count_high{k}=[];
end
for k=1:SumTotAn
    if~isempty(pairs_smallMFhigh{k})
for i =1:length(BinSizes)
Ind_h{k}{i}=find(pairs_smallMFhigh{k}(:,5)==BinSizes(i) & pairs_smallMFhigh{k}(:,3)~=0)
% count_high{k}(i) = length(Ind_h{k}{i})./(sum(nneuFSIhigh)*sum(nneuMSN));
 count_high{k}(i) = length(Ind_h{k}{i})./(nneuFSIhigh(k)*nneuMSN(k));
%  count_high{k}(i) = length(Ind_h{k}{i})./(nneuFSIhigh(k)*nneuFSIhigh(k)); % fsi fsi high


end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% for counting: important change right values when I pass from low to
%%%%%%% high firing and vice versa
% ccc=count_high;
% nneuFSIc=nneuFSIhigh;
ccc=count;
nneuFSIc=nneuFSIlow;
for k=1:SumTotAn
ccc_n{k}=[];
end
for k=1:SumTotAn
    if~isempty(ccc{k})
for i =1:length(BinSizes)
ccc_n{k}(i) = ccc{k}(i).*(nneuFSIc(k)*nneuMSN(k));
%  ccc_n{k}(i) = ccc{k}(i).*(sum(nneuFSIc)*sum(nneuMSN));
% ccc_n{k}(i) = ccc{k}(i).*(nneuFSIc(k)*nneuFSIc(k));  % fsi-fsi
% ccc_n{k}(i) = ccc{k}(i).*(nneuMSN(k)*nneuMSN(k));   %msn-msn
end
    end
end
ccc_num=ccc_n{1};
for i=1:SumTotAn-1
    ccc_num = [ccc_num; ccc_n{i+1}];
end
   ccc_num= sum(ccc_num);
count_low_num=ccc_num;
% count_high_num=ccc_num;
%%%%%%%%%%%%%%%%%%%%


BinDistr = count{1};
for i=1:SumTotAn-1
    BinDistr = [BinDistr;count{i+1}];
    
end

% BinDistr = count_high{1};
% for i=1:SumTotAn-1
%     BinDistr = [BinDistr;count_high{i+1}];
%     
% end

BinDistrM = mean(BinDistr,'omitnan');
 BinDistrSt = std(BinDistr,0,1,'omitnan')/sqrt(SumTotAn);
 
 %%
figure(3);hold on;
bar(BinDistrM,'r'); hold on;
 errorbar(BinDistrM,BinDistrSt/2,'k*');hold on;
%  title({'Mean Bin distribution';nameT})
 ylabel('Assembly pairs/All possible pairs')
 xlabel('\Delta (sec)')
xticks(1:12);
xticklabels(BinSizes);
 box on
ax1 = gca; % current axes
ax1.XColor = 'k';
ax1.YColor = 'k';
 
% ax1.YColor = 'none';
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');hold on;
% set the same Limits and Ticks on ax2 as on ax1;
set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
set(ax2, 'XTick', get(ax1, 'XTick'), 'YTick', get(ax1, 'YTick'));
OppTickLabels = ccc_num;
% Set the x-tick and y-tick  labels for the second axes
set(ax2, 'XTickLabel', OppTickLabels,'YTickLabel',OppTickLabels);
ax2.YColor = 'none';
% set(ax1, 'XTickLabel', BinSizes,'YTickLabel',gca);



%% Session pulled apart
% for k=1:TotAn(1)
%     NP_vsvs{k}=[];
%     NP_vtavta{k}=[];
%     NP_vsvta{k}=[];
% end
% 
% for k=1:TotAn(1)
% NP_vsvs{k}=(nneuVS_All{k}*(nneuVS_All{k}+1))/2;
% NP_vtavta{k}=(nneuVTA_All{k}*(nneuVTA_All{k}+1))/2;
% NP_vsvta{k}=(nneuVS_All{k}*nneuVTA_All{k});
% end
% 
% for k=1:TotAn(1)
% count{k}=[];
% Ind{k} =[];
% end

% % for k=1:TotAn(1)
% %     if~isempty(pairs_vsvs{k})
% % for i =1:length(BinSizes)
% % Ind{k}{i}=find(pairs_vsvs{k}(:,5)==BinSizes(i))
% %  count{k}(i) = length(Ind{k}{i})./NP_vtavta{k};
% % end
% %     end
% % end

BinDistr = count{1};
for i=1:SumTotAn-1
    BinDistr = [BinDistr;count{i+1}];
    
end

BinDistrM = mean(BinDistr,'omitnan');
BinDistrSt = std(BinDistr,0,1,'omitnan')/sqrt(size(BinDistr,1));








figure(fignum(5));hold on;
bar(BinDistrM,'r'); hold on;
errorbar(BinDistrM,BinDistrSt/2,'k*');hold on;
 title({'Mean Bin distribution';nameT})
 ylabel('Assembly pairs/All possible pairs')
 xlabel('\Delta (sec)')
xticks(1:12);
xticklabels(BinSizes);
 box on