
% %%% Pie chart program for pruned and post pruned algorithm...
% %
% %% %%%% Postpruned algorithm
% 
% 
% nameMerge{1} = strcat('Lastrev1_OrRevW.mat');
% nameMerge{2} = strcat('Revconc_OrRevW.mat');
% nameMerge{3} = strcat('Succession_2PhW.mat');
% nameMerge{4} = strcat('Rev2_OrRevW.mat');
% nameMerge{5} = strcat('Rev3_OrRevW.mat');
% nameT = strcat('Lastrev1; Revconc; Succession_2Ph');
% 
% alg= strcat('Stand_');
% 
% fignum =[15,16,17,18,19,20,21,22];
% 
% for i =1: length(nameMerge)
% %     nameSmallBigMerge{i} = strcat('SmallBig_An_Disp',alg,nameMerge{i})
% nameSmallBigMerge{i} = strcat('PostPru_BinRed_An_Disp',alg,nameMerge{i})
% load(nameSmallBigMerge{i})
% clear c
% load('classlist.mat');
% %  nnn(i,:)=nn(a);
% %  aa(i,:)=a;
% [new_data] = parId(new_data,a,nn);
% 
% 
%  for k=1:length(new_data.par)
%  [new_data] = labels(new_data,c,k);
%  end
%  
%  very_data{i}=new_data;
% 
% As_across_smallBinsT{i} = As_across_smallBins;
% As_order_smallBinsT{i} = As_order_smallBins;
% As_across_largeBinsT{i} = As_across_largeBins;
% As_order_largeBinsT{i} = As_order_largeBins; 
% clear new_data a nn As_across_smallBins As_order_smallBins As_across_largeBins As_order_largeBins
% end
% %%
% clearvars -except very_data As_across_smallBinsT As_order_smallBinsT As_across_largeBinsT As_order_largeBinsT c nameMerge nameT fignum
% %%
% for i =1:length(very_data)
%     TotAn(i)=length(very_data{i}.par);
% for k=1:TotAn(i)
%     nneuVS{i}{k}=0;
%     nneuVTA{i}{k}=0;
% end
% end
% 
% SumTotAn=sum(TotAn);
% 
% for i =1:length(very_data)
% for k=1:TotAn(i)
%     for j=1:size(very_data{i}.spike_regionNoId{k},2)
%         [nneuVS{i}{k}] = countNeu(very_data{i}.spike_regionNoId{k}(1,j),1,nneuVS{i}{k});
%         [nneuVTA{i}{k}] = countNeu(very_data{i}.spike_regionNoId{k}(1,j),2,nneuVTA{i}{k});
%     end
% 
% end
% end
% 
% 
% 
% As_across_smallBins_All=As_across_smallBinsT{1};
% As_across_largeBins_All=As_across_largeBinsT{1};
% As_order_smallBins_All=As_order_smallBinsT{1};
% As_order_largeBins_All=As_order_largeBinsT{1};
% nneuVS_All = nneuVS{1};
% nneuVTA_All = nneuVTA{1};
% Data_All.par=very_data{1}.par;
% for i=1:length(very_data)-1
%     [As_across_smallBins_All]=[As_across_smallBins_All,As_across_smallBinsT{i+1}];
%     [As_across_largeBins_All]=[As_across_largeBins_All,As_across_largeBinsT{i+1}];
%     [As_order_smallBins_All]=[As_order_smallBins_All,As_order_smallBinsT{i+1}];
%     [As_order_largeBins_All]=[As_order_largeBins_All,As_order_largeBinsT{i+1}];
%     [nneuVS_All]=[nneuVS_All,nneuVS{i+1}];
%     [nneuVTA_All]=[nneuVTA_All,nneuVTA{i+1}];
%     [Data_All.par]=[Data_All.par,very_data{i+1}.par];
% end
% % 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Postpruned algorithm
%%%%%%% Find Pairs in one region and shared between two differnt regions
% %%%pairs divided by regions with all info (then counting of type of neurons in pairs)
% [pairs_r_small,pairs_vsvs_small,pairs_vsvta_small,pairs_vtavta_small] = PairsInfo(As_across_smallBins_All,As_order_smallBins_All,nneuVS_All,SumTotAn);
% [pairs_r_large,pairs_vsvs_large,pairs_vsvta_large,pairs_vtavta_large] = PairsInfo(As_across_largeBins_All,As_order_largeBins_All,nneuVS_All,SumTotAn);
% 
% clear a b A B i ii j k
% 
% [struct_pair_small] = PairsOrgByPair(pairs_vsvta_small);
% [struct_pair_large] = PairsOrgByPair(pairs_vsvta_large);
% 
% 
% [struct_pair_small] = PairNeuLabelsID(struct_pair_small,Data_All);
% [struct_pair_large] = PairNeuLabelsID(struct_pair_large,Data_All);
% copy_struct_pair_small=struct_pair_small;
% copy_struct_pair_large=struct_pair_large;
% clear struct_pair

%%
addpath('/zifnas/Carla/CWEM_Project_GenProg');
% cd(path_pn)
MergeDataGenPostPru_1;
MergeDataGenPostPru_2;
%%% Create the struct coerently with the directionality
%%%%%%% Here the directionality doesn't cointan the info about the bins
%%%%%%% because the pairs are already divided for bins

% %Dir stays for "Directionality", can assume these values

 Dir1 = 'vs->vta';
 Dir2 = 'vta->vs';
 Dir3 = 'sync';
 Dir4 = 'no';

addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM');

for k=1:SumTotAn
    if ~isempty(copy_struct_pair_small)
[struct_pair_small_vsvta{k}] = DirAndBinPostPru(Dir1, copy_struct_pair_small,k);
 [struct_pair_small_vtavs{k}] = DirAndBinPostPru(Dir2, copy_struct_pair_small,k);
 [struct_pair_small_sync{k}] = DirAndBinPostPru(Dir3, copy_struct_pair_small,k);
    end
    if ~isempty(copy_struct_pair_large)
    [struct_pair_large_vsvta{k}] = DirAndBinPostPru(Dir1, copy_struct_pair_large,k);
    [struct_pair_large_vtavs{k}] = DirAndBinPostPru(Dir2, copy_struct_pair_large,k);
    [struct_pair_large_sync{k}] = DirAndBinPostPru(Dir3, copy_struct_pair_large,k);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%% Pruned Algorithm
%%%%%%% PieChart for pruned data set

% % % % % nameMerge{1} = strcat('Lastrev1_OrRevW.mat');
% % % % % nameMerge{2} = strcat('Revconc_OrRevW.mat');
% % % % % nameMerge{3} = strcat('Succession_2PhW.mat');
% % % % % nameMerge{4} = strcat('Rev2_OrRevW.mat');
% % % % % nameMerge{5} = strcat('Rev3_OrRevW.mat');
% % % % % nameT = strcat('Lastrev1; Revconc; Succession_2Ph');
% % % % % 
% % % % % 
% % % % % algPru= strcat('Pru_');
% % % % % 
% % % % % fignum =[1,2,3,4,5,6,7,8];
% % % % % 
% % % % % for i =1: length(nameMerge)
% % % % % %     nameSmallBigMerge{i} = strcat('SmallBig_An_Disp',alg,nameMerge{i})
% % % % % nameSmallBigMerge{i} = strcat('An_Disp_',algPru,nameMerge{i})
% % % % %     load(nameSmallBigMerge{i})
% % % % % clear c
% % % % % load('classlist.mat');
% % % % % %  nnn(i,:)=nn(a);
% % % % % %  aa(i,:)=a;
% % % % % [new_data] = parId(new_data,a,nn);
% % % % % 
% % % % % 
% % % % %  for k=1:length(new_data.par)
% % % % %  [new_data] = labels(new_data,c,k);
% % % % %  end
% % % % %  
% % % % %  very_data{i}=new_data;
% % % % % 
% % % % % As_across_binsT{i} = As_acr_bins_pru;
% % % % % As_order_T{i} = As_order_pru;
% % % % %  
% % % % % clear new_data a nn As_acr_bins_pru As_order_pru
% % % % % end
% % % % % 
% % % % %  clearvars -except very_data As_across_binsT As_order_T c nameMerge nameT fignum
% % % % % 
% % % % % for i =1:length(very_data)
% % % % %     TotAn(i)=length(very_data{i}.par);
% % % % % for k=1:TotAn(i)
% % % % %     nneuVS{i}{k}=0;
% % % % %     nneuVTA{i}{k}=0;
% % % % % end
% % % % % end
% % % % % 
% % % % % SumTotAn=sum(TotAn);
% % % % % 
% % % % % for i =1:length(very_data)
% % % % % for k=1:TotAn(i)
% % % % %     for j=1:size(very_data{i}.spike_regionNoId{k},2)
% % % % %         [nneuVS{i}{k}] = countNeu(very_data{i}.spike_regionNoId{k}(1,j),1,nneuVS{i}{k});
% % % % %         [nneuVTA{i}{k}] = countNeu(very_data{i}.spike_regionNoId{k}(1,j),2,nneuVTA{i}{k});
% % % % %     end
% % % % % 
% % % % % end
% % % % % end
% % % % % 
% % % % % 
% % % % % 
% % % % % As_across_bins_All=As_across_binsT{1};
% % % % % As_order_All=As_order_T{1};
% % % % % nneuVS_All = nneuVS{1};
% % % % % nneuVTA_All = nneuVTA{1};
% % % % % Data_All.par=very_data{1}.par;
% % % % % for i=1:length(very_data)-1
% % % % %     [As_across_bins_All]=[As_across_bins_All,As_across_binsT{i+1}];
% % % % %     [As_order_All]=[As_order_All,As_order_T{i+1}];
% % % % %     [nneuVS_All]=[nneuVS_All,nneuVS{i+1}];
% % % % %     [nneuVTA_All]=[nneuVTA_All,nneuVTA{i+1}];
% % % % %     [Data_All.par]=[Data_All.par,very_data{i+1}.par];
% % % % % end
% % % % % 
% % % % % 
% % % % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pruned algorithm Pairs
% % % % % %%%%% Find Pairs in one region and shared between two differnt regions
% % % % % % pairs diveded by regions with all info
% % % % % [pairs_r,pairs_vsvs,pairs_vsvta,pairs_vtavta] = PairsInfo(As_across_bins_All,As_order_All,nneuVS_All,SumTotAn);
% % % % % 
% % % % % clear a b A B i ii j k
% % % % % 
% % % % % [struct_pair] = PairsOrgByPair(pairs_vsvta);
% % % % % [struct_pair] = PairNeuLabelsID(struct_pair,Data_All);
% % % % % 
% % % % % % copy_struct_pair_small=struct_pair_small;
% % % % % % copy_struct_pair_large=struct_pair_large;
% % % % % % clear struct_pair
% % % % % 
% % % % % 
% % % % % 
% % % % % %%%%%% Create the struct coerently with the directionality
% % % % % %%%%%%%%%% Here the directionality doesn't cointan the info about the bins
% % % % % %%%%%%%%%% because the pairs are already divided for bins
% % % % % 
% % % % % % % % % %Dir stays for "Directionality", can assume these values
% % % % % 
% % % % %  Dir1 = 'vs->vta small';
% % % % %  Dir2 = 'vta->vs small';
% % % % %  Dir3 = 'sync small';
% % % % %  Dir4 = 'no';
% % % % %  Dir5 = 'vs->vta large';
% % % % %  Dir6 = 'vta->vs large';
% % % % %  Dir7 = 'sync large';
% % % % %  
% % % % % 
% % % % % 
% % % % % for k=1:SumTotAn
% % % % %     if ~isempty(struct_pair{k})
% % % % %  [struct_pair_small_vsvta{k}] = DirAndBin(Dir1, struct_pair,k);
% % % % %  [struct_pair_small_vtavs{k}] = DirAndBin(Dir2, struct_pair,k);
% % % % %  [struct_pair_small_sync{k}] = DirAndBin(Dir3, struct_pair,k);
% % % % %     end
% % % % %     if ~isempty(struct_pair{k})
% % % % %     [struct_pair_large_vsvta{k}] = DirAndBin(Dir5, struct_pair,k);
% % % % %     [struct_pair_large_vtavs{k}] = DirAndBin(Dir6, struct_pair,k);
% % % % %     [struct_pair_large_sync{k}] = DirAndBin(Dir7, struct_pair,k);
% % % % %     end
% % % % % end


%% %%%%%%%%%%%% Counting of the different quantities
CountType_All = zeros(SumTotAn,8);
for k =1: SumTotAn
    for i= 1:length(Data_All.par{k})
[TLabNeu,CountType_All] = CountTypeNeuGen(Data_All.par{k}(i).labels,k,CountType_All);
    end
end

CountTypeVS = sum(TLabNeu{:,1:4});
CountTypeVTA = sum(TLabNeu{:,5:8});
CountTypeAll = sum(TLabNeu{:,:});
%% Plot Global Pie

TextVs = {'SPN ';'FSN ';'CIN'; 'NoID'}; % strings
TextVta = {'DAN ';'GABA ';'GLU'; 'NoID'}; % strings
figure(fignum(1));
subplot(1,4,[1,2]); hold on;
Pie_Vs_Tot = pie(CountTypeVS(find(CountTypeVS~=0)));hold on;
Pie_Vs_TotText = findobj(Pie_Vs_Tot,'Type','text'); % text object handles
percentValues = get(Pie_Vs_TotText,'String'); % percent values
% txt = TextVs(find(CountTypeVS));
% combinedtxt = strcat(txt,percentValues); % strings and percent values
% oldExtents_cell = get(Pie_Vs_TotText,'Extent'); % cell array
% oldExtents = cell2mat(oldExtents_cell);
% for i =1:length(txt)
% Pie_Vs_TotText(i).String = combinedtxt(i);
% endLEG = findobj(AX,'type','text');
lgd=legend(TextVs,'Location','southoutside','Orientation','horizontal');
lgd.FontSize=18;
axis off

subplot(1,4,[3,4]); hold on;
Pie_Vta_Tot = pie(CountTypeVTA(find(CountTypeVTA~=0)));hold on;
Pie_Vta_TotText = findobj(Pie_Vta_Tot,'Type','text'); % text object handles
percentValues = get(Pie_Vta_TotText,'String'); % percent values
% txt = TextVta(find(CountTypeVTA))
% combinedtxt = strcat(txt,percentValues); % strings and percent values
% oldExtents_cell = get(Pie_Vta_TotText,'Extent'); % cell array
% oldExtents = cell2mat(oldExtents_cell);
% for i=1:length(txt)
% Pie_Vta_TotText(i).String = combinedtxt(i);
% end
lgd1=legend(TextVta,'Location','southoutside','Orientation','horizontal')
lgd1.FontSize=18;
axis off; 

%% Counting for neurons in pairs without repetions



for k=1:SumTotAn
    Un_small_vsvta{k}=struct('neuID',[],'ia',[],'ic',[],'labels',[]);
    [Un_small_vsvta{k}.neuID,Un_small_vsvta{k}(:,1).ia(:,1),Un_small_vsvta{k}.ic(:,1)]=unique(struct_pair_small_vsvta{k}.neuID);
    Un_small_vsvta{k}.labels(:,1) = struct_pair_small_vsvta{k}.labels(Un_small_vsvta{k}.ia); % I create here the structure with the index unique
       
    
    Un_small_vtavs{k}=struct('neuID',[],'ia',[],'ic',[],'labels',[]);
    [Un_small_vtavs{k}.neuID,Un_small_vtavs{k}(:,1).ia(:,1),Un_small_vtavs{k}.ic(:,1)]=unique(struct_pair_small_vtavs{k}.neuID);
    Un_small_vtavs{k}.labels(:,1) = struct_pair_small_vtavs{k}.labels(Un_small_vtavs{k}.ia); % I create here the structure with the index unique
    
       
    Un_small_sync{k}=struct('neuID',[],'ia',[],'ic',[],'labels',[]);
    [Un_small_sync{k}.neuID,Un_small_sync{k}(:,1).ia(:,1),Un_small_sync{k}.ic(:,1)]=unique(struct_pair_small_sync{k}.neuID);
    Un_small_sync{k}.labels(:,1) = struct_pair_small_sync{k}.labels(Un_small_sync{k}.ia); % I create here the structure with the index unique
       
    
    
    
    Un_large_vsvta{k}=struct('neuID',[],'ia',[],'ic',[],'labels',[]);
    [Un_large_vsvta{k}.neuID,Un_large_vsvta{k}(:,1).ia(:,1),Un_large_vsvta{k}.ic(:,1)]=unique(struct_pair_large_vsvta{k}.neuID);
    Un_large_vsvta{k}.labels(:,1) = struct_pair_large_vsvta{k}.labels(Un_large_vsvta{k}.ia); % I create here the structure with the index unique
    
    Un_large_vtavs{k}=struct('neuID',[],'ia',[],'ic',[],'labels',[]);
    [Un_large_vtavs{k}.neuID,Un_large_vtavs{k}(:,1).ia(:,1),Un_large_vtavs{k}.ic(:,1)]=unique(struct_pair_large_vtavs{k}.neuID);
    Un_large_vtavs{k}.labels(:,1) = struct_pair_large_vtavs{k}.labels(Un_large_vtavs{k}.ia); % I create here the structure with the index unique
    
    Un_large_sync{k}=struct('neuID',[],'ia',[],'ic',[],'labels',[]);
    [Un_large_sync{k}.neuID,Un_large_sync{k}(:,1).ia(:,1),Un_large_sync{k}.ic(:,1)]=unique(struct_pair_large_sync{k}.neuID);
    Un_large_sync{k}.labels(:,1) = struct_pair_large_sync{k}.labels(Un_large_sync{k}.ia); % I create here the structure with the index unique
end


%%%%% Count Small Bins


CountNeuP_small_vsvta=zeros(SumTotAn,8);
for k=1:SumTotAn
    if~isempty(Un_small_vsvta{k}.neuID)
        for i=1:length(Un_small_vsvta{k}.neuID)
[TabNeuP_small_vsvta,CountNeuP_small_vsvta] = CountTypeNeuGen(Un_small_vsvta{k}.labels(i),k,CountNeuP_small_vsvta); 
        end
    end
end

SumCountNeuP_small_vsvta_VS = sum(TabNeuP_small_vsvta{:,1:4});
SumCountNeuP_small_vsvta_VTA = sum(TabNeuP_small_vsvta{:,5:8});
SumCountNeuP_small_vsvta = sum(TabNeuP_small_vsvta{:,:});


CountNeuP_small_vtavs=zeros(SumTotAn,8);
for k=1:SumTotAn
    if~isempty(Un_small_vtavs{k}.neuID)
        for i=1:length(Un_small_vtavs{k}.neuID)
[TabNeuP_small_vtavs,CountNeuP_small_vtavs] = CountTypeNeuGen(Un_small_vtavs{k}.labels(i),k,CountNeuP_small_vtavs); 
        end
    end
end

SumCountNeuP_small_vtavs_VS = sum(TabNeuP_small_vtavs{:,1:4});
SumCountNeuP_small_vtavs_VTA = sum(TabNeuP_small_vtavs{:,5:8});
SumCountNeuP_small_vtavs = sum(TabNeuP_small_vtavs{:,:});


CountNeuP_small_sync=zeros(SumTotAn,8);
for k=1:SumTotAn
    if~isempty(Un_small_sync{k}.neuID)
        for i=1:length(Un_small_sync{k}.neuID)
[TabNeuP_small_sync,CountNeuP_small_sync] = CountTypeNeuGen(Un_small_sync{k}.labels(i),k,CountNeuP_small_sync); 
        end
    end
end

SumCountNeuP_small_sync_VS = sum(TabNeuP_small_sync{:,1:4});
SumCountNeuP_small_sync_VTA = sum(TabNeuP_small_sync{:,5:8});
SumCountNeuP_small_sync = sum(TabNeuP_small_sync{:,:});


%%%%%% Count Neurons in Pairs Large Bins


CountNeuP_large_vsvta=zeros(SumTotAn,8);
for k=1:SumTotAn
    if~isempty(Un_large_vsvta{k}.neuID)
        for i=1:length(Un_large_vsvta{k}.neuID)
[TabNeuP_large_vsvta,CountNeuP_large_vsvta] = CountTypeNeuGen(Un_large_vsvta{k}.labels(i),k,CountNeuP_large_vsvta); 
        end
    end
end


SumCountNeuP_large_vsvta_VS = sum(TabNeuP_large_vsvta{:,1:4});
SumCountNeuP_large_vsvta_VTA = sum(TabNeuP_large_vsvta{:,5:8});
SumCountNeuP_large_vsvta = sum(TabNeuP_large_vsvta{:,:});

% SumCountNeuP_large_vsvta_VS = mean(TabNeuP_large_vsvta{:,1:4},'omitnan');
% SumCountNeuP_large_vsvta_VTA = mean(TabNeuP_large_vsvta{:,5:8},'omitnan');
% SumCountNeuP_large_vsvta = mean(TabNeuP_large_vsvta{:,:},'omitnan');


CountNeuP_large_vtavs=zeros(SumTotAn,8);
for k=1:SumTotAn
    if~isempty(Un_large_vtavs{k}.neuID)
        for i=1:length(Un_large_vtavs{k}.neuID)
[TabNeuP_large_vtavs,CountNeuP_large_vtavs] = CountTypeNeuGen(Un_large_vtavs{k}.labels(i),k,CountNeuP_large_vtavs); % general...for each type of stucture,then I have to pass the .labels
        end
    end
end


SumCountNeuP_large_vtavs_VS = sum(TabNeuP_large_vtavs{:,1:4});
SumCountNeuP_large_vtavs_VTA = sum(TabNeuP_large_vtavs{:,5:8});
SumCountNeuP_large_vtavs = sum(TabNeuP_large_vtavs{:,:});


% SumCountNeuP_large_vtavs_VS = mean(TabNeuP_large_vtavs{:,1:4},'omitnan');
% SumCountNeuP_large_vtavs_VTA = mean(TabNeuP_large_vtavs{:,5:8},'omitnan');
% SumCountNeuP_large_vtavs = mean(TabNeuP_large_vtavs{:,:},'omitnan');

CountNeuP_large_sync=zeros(SumTotAn,8);
for k=1:SumTotAn
    if~isempty(Un_large_sync{k}.neuID)
        for i=1:length(Un_large_sync{k}.neuID)
[TabNeuP_large_sync,CountNeuP_large_sync] = CountTypeNeuGen(Un_large_sync{k}.labels(i),k,CountNeuP_large_sync); 
        end
    end
end

% SumCountNeuP_large_sync_VS = mean(TabNeuP_large_sync{:,1:4},'omitnan');
% SumCountNeuP_large_sync_VTA = mean(TabNeuP_large_sync{:,5:8},'omitnan');
% SumCountNeuP_large_sync = mean(TabNeuP_large_sync{:,:},'omitnan');

SumCountNeuP_large_sync_VS = sum(TabNeuP_large_sync{:,1:4});
SumCountNeuP_large_sync_VTA = sum(TabNeuP_large_sync{:,5:8});
SumCountNeuP_large_sync = sum(TabNeuP_large_sync{:,:});
%% Plot w/orep only small bin size
figure(68);hold on;
subplot(1,4,[1,2]); hold on;
PVS_small_vsvta = pie(SumCountNeuP_small_vsvta_VS(find(SumCountNeuP_small_vsvta_VS~=0)));hold on;
PVS_small_vsvtaText = findobj(PVS_small_vsvta,'Type','text'); % text object handles
percentValues = get(PVS_small_vsvtaText,'String'); % percent values
% txt = TextVs(find(SumCountNeuP_small_vsvta_VS)); % strings
% combinedtxt = strcat(txt,percentValues); % strings and percent values
% oldExtents_cell = get(PVS_small_vsvtaText,'Extent'); % cell array
% oldExtents = cell2mat(oldExtents_cell);
% for i=1:length(txt)
% PVS_small_vsvtaText(i).String = combinedtxt(i);hold on; 
% end
lgd1=legend(TextVs,'Location','southoutside','Orientation','horizontal')
lgd1.FontSize=18;
axis off;hold on; 


subplot(1,4,[3,4]); hold on;
PVS_small_vtavs = pie(SumCountNeuP_small_vtavs_VS(find(SumCountNeuP_small_vtavs_VS~=0)));hold on;
PVS_small_vtavsText = findobj(PVS_small_vtavs,'Type','text'); % text object handles
percentValues = get(PVS_small_vtavsText,'String'); % percent values
% txt = TextVs(find(SumCountNeuP_small_vtavs_VS)); % strings
% combinedtxt = strcat(txt,percentValues); % strings and percent values
% oldExtents_cell = get(PVS_small_vtavsText,'Extent'); % cell array
% oldExtents = cell2mat(oldExtents_cell);
% for i=1:length(txt)
% PVS_small_vtavsText(i).String = combinedtxt(i);
% end
lgd1=legend(TextVs,'Location','southoutside','Orientation','horizontal')
lgd1.FontSize=18;
axis off;hold on; 
axis off


%% Plot Neu in Pairs without repetitions


%%%%%%%%%%%%% VS
figure(fignum(2));hold on;
subplot(2,3,1); hold on;
PVS_small_vsvta = pie(SumCountNeuP_small_vsvta_VS(find(SumCountNeuP_small_vsvta_VS~=0)));hold on;
PVS_small_vsvtaText = findobj(PVS_small_vsvta,'Type','text'); % text object handles
percentValues = get(PVS_small_vsvtaText,'String'); % percent values
txt = TextVs(find(SumCountNeuP_small_vsvta_VS)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PVS_small_vsvtaText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
PVS_small_vsvtaText(i).String = combinedtxt(i);hold on; 
end
lgd1=legend(TextVs,'Location','southoutside','Orientation','horizontal')
lgd1.FontSize=12;
axis off;hold on; 


subplot(2,3,2); hold on;
PVS_small_vtavs = pie(SumCountNeuP_small_vtavs_VS(find(SumCountNeuP_small_vtavs_VS~=0)));hold on;
PVS_small_vtavsText = findobj(PVS_small_vtavs,'Type','text'); % text object handles
percentValues = get(PVS_small_vtavsText,'String'); % percent values
txt = TextVs(find(SumCountNeuP_small_vtavs_VS)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PVS_small_vtavsText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
PVS_small_vtavsText(i).String = combinedtxt(i);
end
axis off



subplot(2,3,3); hold on;
PVS_small_sync = pie(SumCountNeuP_small_sync_VS(find(SumCountNeuP_small_sync_VS~=0)));hold on;
PVS_small_syncText = findobj(PVS_small_sync,'Type','text'); % text object handles
percentValues = get(PVS_small_syncText,'String'); % percent values
txt = TextVs(find(SumCountNeuP_small_sync_VS)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PVS_small_syncText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
PVS_small_syncText(i).String = combinedtxt(i);
end
% PVS_small_syncText(4).String = combinedtxt(4);

axis off

subplot(2,3,4); hold on;
PVS_large_vsvta = pie(SumCountNeuP_large_vsvta_VS(find(SumCountNeuP_large_vsvta_VS~=0)));hold on;
PVS_large_vsvtaText = findobj(PVS_large_vsvta,'Type','text'); % text object handles
percentValues = get(PVS_large_vsvtaText,'String'); % percent values
txt = TextVs(find(SumCountNeuP_large_vsvta_VS)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PVS_large_vsvtaText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i =1:length(txt)
PVS_large_vsvtaText(i).String = combinedtxt(i);
end

axis off


subplot(2,3,5); hold on;
PVS_large_vtavs = pie(SumCountNeuP_large_vtavs_VS(find(SumCountNeuP_large_vtavs_VS~=0)));hold on;
PVS_large_vtavsText = findobj(PVS_large_vtavs,'Type','text'); % text object handles
percentValues = get(PVS_large_vtavsText,'String'); % percent values
txt = TextVs(find(SumCountNeuP_large_vtavs_VS)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PVS_large_vtavsText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
PVS_large_vtavsText(i).String = combinedtxt(i);
end

axis off


subplot(2,3,6); hold on;
PVS_large_sync = pie(SumCountNeuP_large_sync_VS(find(SumCountNeuP_large_sync_VS~=0)));hold on;
PVS_large_syncText = findobj(PVS_large_sync,'Type','text'); % text object handles
percentValues = get(PVS_large_syncText,'String'); % percent values
txt = TextVs(find(SumCountNeuP_large_sync_VS)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PVS_large_syncText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
PVS_large_syncText(i).String = combinedtxt(i);
end

axis off

%% %%%%%%%%%%%%%%%%%%%%%% VTA
TextVTA_ = {'DAN: ';'GABA: ';'GLU:'; 'NoID:'};
figure(fignum(3));hold on;
 subplot(2,3,1); hold on;
PVTA_small_vsvta = pie(SumCountNeuP_small_vsvta_VTA(find(SumCountNeuP_small_vsvta_VTA~=0)));hold on;
PVTA_small_vsvtaText = findobj(PVTA_small_vsvta,'Type','text'); % text object handles
percentValues = get(PVTA_small_vsvtaText,'String'); % percent values
txt =TextVTA_(find(SumCountNeuP_small_vsvta_VTA)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PVTA_small_vsvtaText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
PVTA_small_vsvtaText(i).String = combinedtxt(i);
end
% PVTA_small_vsvtaText(1).String = combinedtxt(1);hold on; 
% PVTA_small_vsvtaText(2).String = combinedtxt(2);hold on; 
% PVTA_small_vsvtaText(3).String = combinedtxt(3);hold on; 
% PVTA_small_vsvtaText(4).String = combinedtxt(4);hold on; 
%  title('vs->vta')
axis off;hold on; 


subplot(2,3,2); hold on;
PVTA_small_vtavs = pie(SumCountNeuP_small_vtavs_VTA(find(SumCountNeuP_small_vtavs_VTA~=0)));hold on;
PVTA_small_vtavsText = findobj(PVTA_small_vtavs,'Type','text'); % text object handles
percentValues = get(PVTA_small_vtavsText,'String'); % percent values
txt = TextVTA_(find(SumCountNeuP_small_vtavs_VTA)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PVTA_small_vtavsText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
PVTA_small_vtavsText(i).String = combinedtxt(i);
end
% PVTA_small_vtavsText(1).String = combinedtxt(1);hold on; 
% PVTA_small_vtavsText(2).String = combinedtxt(2);hold on; 
% PVTA_small_vtavsText(3).String = combinedtxt(3);hold on; 
% PVTA_small_vtavsText(4).String = combinedtxt(4);hold on; 
%  title('vs->vta')
axis off;hold on; 


subplot(2,3,3); hold on;
PVTA_small_sync = pie(SumCountNeuP_small_sync_VTA(find(SumCountNeuP_small_sync_VTA~=0)));hold on;
PVTA_small_syncText = findobj(PVTA_small_sync,'Type','text'); % text object handles
percentValues = get(PVTA_small_syncText,'String'); % percent values
txt = TextVTA_(find(SumCountNeuP_small_sync_VTA)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PVTA_small_syncText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
PVTA_small_syncText(i).String = combinedtxt(i);
end
% PVTA_small_syncText(1).String = combinedtxt(1);hold on; 
% PVTA_small_syncText(2).String = combinedtxt(2);hold on; 
% PVTA_small_syncText(3).String = combinedtxt(3);hold on; 
% PVTA_small_syncText(4).String = combinedtxt(4);hold on; 
%  title('vs->vta')
axis off;hold on; 


subplot(2,3,4); hold on;
PVTA_large_vsvta = pie(SumCountNeuP_large_vsvta_VTA(find(SumCountNeuP_large_vsvta_VTA~=0)));hold on;
PVTA_large_vsvtaText = findobj(PVTA_large_vsvta,'Type','text'); % text object handles
percentValues = get(PVTA_large_vsvtaText,'String'); % percent values
txt = TextVTA_(find(SumCountNeuP_large_vsvta_VTA)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PVTA_large_vsvtaText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
PVTA_large_vsvtaText(i).String = combinedtxt(i);
end
% PVTA_large_vsvtaText(1).String = combinedtxt(1);hold on; 
% PVTA_large_vsvtaText(2).String = combinedtxt(2);hold on; 
% PVTA_large_vsvtaText(3).String = combinedtxt(3);hold on; 
% PVTA_large_vsvtaText(4).String = combinedtxt(4);hold on; 
%  title('vs->vta')
axis off;hold on; 


subplot(2,3,5); hold on;
PVTA_large_vtavs = pie(SumCountNeuP_large_vtavs_VTA(find(SumCountNeuP_large_vtavs_VTA~=0)));hold on;
PVTA_large_vtavsText = findobj(PVTA_large_vtavs,'Type','text'); % text object handles
percentValues = get(PVTA_large_vtavsText,'String'); % percent values
txt = TextVTA_(find(SumCountNeuP_large_vtavs_VTA)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PVTA_large_vtavsText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
PVTA_large_vtavsText(i).String = combinedtxt(i);
end
% PVTA_large_vtavsText(1).String = combinedtxt(1);hold on; 
% PVTA_large_vtavsText(2).String = combinedtxt(2);hold on; 
% PVTA_large_vtavsText(3).String = combinedtxt(3);hold on; 
% PVTA_large_vtavsText(4).String = combinedtxt(4);hold on; 
%  title('vs->vta')
axis off;hold on; 

subplot(2,3,6); hold on;
PVTA_large_sync = pie(SumCountNeuP_large_sync_VTA(find(SumCountNeuP_large_sync_VTA~=0)));hold on;
PVTA_large_syncText = findobj(PVTA_large_sync,'Type','text'); % text object handles
percentValues = get(PVTA_large_syncText,'String'); % percent values
txt = TextVTA_(find(SumCountNeuP_large_sync_VTA)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PVTA_large_syncText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
PVTA_large_syncText(i).String = combinedtxt(i);
end
% PVTA_large_syncText(1).String = combinedtxt(1);hold on; 
% PVTA_large_syncText(2).String = combinedtxt(2);hold on; 
% PVTA_large_syncText(3).String = combinedtxt(3);hold on; 
% PVTA_large_syncText(4).String = combinedtxt(4);hold on; 
%  title('vs->vta')
axis off;hold on; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Do the assemblies prefer specific units?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Total Number of assemblies

for k=1:SumTotAn
   % if ~isempty(struct_pair{k}.pair)
        count_Pairs_small_vsvta(k)=length(struct_pair_small_vsvta{k}.pair);
        count_Pairs_small_vtavs(k)=length(struct_pair_small_vtavs{k}.pair);
        count_Pairs_small_sync(k)=length(struct_pair_small_sync{k}.pair);
        
        count_Pairs_large_vsvta(k)=length(struct_pair_large_vsvta{k}.pair);
        count_Pairs_large_vtavs(k)=length(struct_pair_large_vtavs{k}.pair);
        count_Pairs_large_sync(k)=length(struct_pair_large_sync{k}.pair);
   % end
end

TotAs_small_vsvta=sum(count_Pairs_small_vsvta);
TotAs_small_vtavs=sum(count_Pairs_small_vtavs);
TotAs_small_sync=sum(count_Pairs_small_sync);

TotAs_large_vsvta=sum(count_Pairs_large_vsvta);
TotAs_large_vtavs=sum(count_Pairs_large_vtavs);
TotAs_large_sync=sum(count_Pairs_large_sync);
%%%%%%%%%%%%%%
%%%%%%%%%% Number of assembly with a specific class of neurons %%%%%%% !!!!! 

%%% Small bins different directionality
count_PwN_small_vsvta=zeros(SumTotAn,8);
 for k=1:SumTotAn
[Tab_PwN_small_vsvta,count_PwN_small_vsvta] = CountTypeNeuPair(struct_pair_small_vsvta{k},k,count_PwN_small_vsvta); % pairs with specific units, animal per animal
 end
 
 count_PwN_small_vtavs=zeros(SumTotAn,8);
 for k=1:SumTotAn
[Tab_PwN_small_vtavs,count_PwN_small_vtavs] = CountTypeNeuPair(struct_pair_small_vtavs{k},k,count_PwN_small_vtavs); % pairs with specific units, animal per animal
 end
 
 count_PwN_small_sync=zeros(SumTotAn,8);
 for k=1:SumTotAn
[Tab_PwN_small_sync,count_PwN_small_sync] = CountTypeNeuPair(struct_pair_small_sync{k},k,count_PwN_small_sync); % pairs with specific units, animal per animal
 end
 
 %%% Large bins different directionality
count_PwN_large_vsvta=zeros(SumTotAn,8);
 for k=1:SumTotAn
[Tab_PwN_large_vsvta,count_PwN_large_vsvta] = CountTypeNeuPair(struct_pair_large_vsvta{k},k,count_PwN_large_vsvta); % pairs with specific units, animal per animal
 end
 
 count_PwN_large_vtavs=zeros(SumTotAn,8);
 for k=1:SumTotAn
[Tab_PwN_large_vtavs,count_PwN_large_vtavs] = CountTypeNeuPair(struct_pair_large_vtavs{k},k,count_PwN_large_vtavs); % pairs with specific units, animal per animal
 end
 
 count_PwN_large_sync=zeros(SumTotAn,8);
 for k=1:SumTotAn
[Tab_PwN_large_sync,count_PwN_large_sync] = CountTypeNeuPair(struct_pair_large_sync{k},k,count_PwN_large_sync); % pairs with specific units, animal per animal
 end
 
 %% %%%%%%
 SumCount_PwN_small_vsvta = sum(Tab_PwN_small_vsvta{:,:});
 SumCount_PwN_small_vsvtaVS = sum(Tab_PwN_small_vsvta{:,1:4});
 SumCount_PwN_small_vsvtaVTA = sum(Tab_PwN_small_vsvta{:,5:8});
 
 SumCount_PwN_small_vtavs = sum(Tab_PwN_small_vtavs{:,:});
 SumCount_PwN_small_vtavsVS = sum(Tab_PwN_small_vtavs{:,1:4});
 SumCount_PwN_small_vtavsVTA = sum(Tab_PwN_small_vtavs{:,5:8});
 
 SumCount_PwN_small_sync = sum(Tab_PwN_small_sync{:,:});
 SumCount_PwN_small_syncVS = sum(Tab_PwN_small_sync{:,1:4});
 SumCount_PwN_small_syncVTA = sum(Tab_PwN_small_sync{:,5:8});
 
 SumCount_PwN_large_vsvta = sum(Tab_PwN_large_vsvta{:,:});
 SumCount_PwN_large_vsvtaVS = sum(Tab_PwN_large_vsvta{:,1:4});
 SumCount_PwN_large_vsvtaVTA = sum(Tab_PwN_large_vsvta{:,5:8});
 
 SumCount_PwN_large_vtavs = sum(Tab_PwN_large_vtavs{:,:});
 SumCount_PwN_large_vtavsVS = sum(Tab_PwN_large_vtavs{:,1:4});
 SumCount_PwN_large_vtavsVTA = sum(Tab_PwN_large_vtavs{:,5:8});
 
 SumCount_PwN_large_sync = sum(Tab_PwN_large_sync{:,:});
 SumCount_PwN_large_syncVS = sum(Tab_PwN_large_sync{:,1:4});
 SumCount_PwN_large_syncVTA = sum(Tab_PwN_large_sync{:,5:8});
 
 %% %%%%% Plot Pie Chart with repetitions
 TextAsRep = {'SPN: ';'FSN: ';'CIN:';'NoID:' }; % strings
figure(fignum(4));hold on;
subplot(2,3,1); hold on;
PwN_small_vsvtaVS = pie(SumCount_PwN_small_vsvtaVS(find(SumCount_PwN_small_vsvtaVS~=0)));hold on;
PwN_small_vsvtaVSText = findobj(PwN_small_vsvtaVS,'Type','text'); % text object handles
percentValues = get(PwN_small_vsvtaVSText,'String'); % percent values
txt = TextAsRep(find(SumCount_PwN_small_vsvtaVS));
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PwN_small_vsvtaVSText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i =1:length(txt)
PwN_small_vsvtaVSText(i).String = combinedtxt(i);hold on; 
end 

% PVS_small_vsvtaText(4).String = combinedtxt(4);hold on; 
%  title('vs->vta')
axis off;hold on; 

subplot(2,3,2); hold on;
PwN_small_vtavsVS = pie(SumCount_PwN_small_vtavsVS(find(SumCount_PwN_small_vtavsVS~=0)));hold on;
PwN_small_vtavsVSText = findobj(PwN_small_vtavsVS,'Type','text'); % text object handles
percentValues = get(PwN_small_vtavsVSText,'String'); % percent values
txt = TextAsRep(find(SumCount_PwN_small_vtavsVS)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PwN_small_vtavsVSText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i =1:length(txt)
PwN_small_vtavsVSText(i).String = combinedtxt(i);hold on; 
end
axis off;hold on; 


subplot(2,3,3); hold on;
PwN_small_syncVS = pie(SumCount_PwN_small_syncVS(find(SumCount_PwN_small_syncVS~=0)));hold on;
PwN_small_syncVSText = findobj(PwN_small_syncVS,'Type','text'); % text object handles
percentValues = get(PwN_small_syncVSText,'String'); % percent values
txt = TextAsRep(find(SumCount_PwN_small_syncVS)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PwN_small_syncVSText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
PwN_small_syncVSText(i).String = combinedtxt(i);hold on;  
end
axis off;hold on; 

subplot(2,3,4); hold on;
PwN_large_vsvtaVS = pie(SumCount_PwN_large_vsvtaVS(find(SumCount_PwN_large_vsvtaVS~=0)));hold on;
PwN_large_vsvtaVSText = findobj(PwN_large_vsvtaVS,'Type','text'); % text object handles
percentValues = get(PwN_large_vsvtaVSText,'String'); % percent values
txt = TextAsRep(find(SumCount_PwN_large_vsvtaVS)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PwN_large_vsvtaVSText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
PwN_large_vsvtaVSText(i).String = combinedtxt(i);hold on; 
end
axis off;hold on; 


subplot(2,3,5); hold on;
PwN_large_vtavsVS = pie(SumCount_PwN_large_vtavsVS(find(SumCount_PwN_large_vtavsVS~=0)));hold on;
PwN_large_vtavsVSText = findobj(PwN_large_vtavsVS,'Type','text'); % text object handles
percentValues = get(PwN_large_vtavsVSText,'String'); % percent values
txt = TextAsRep(find(SumCount_PwN_large_vtavsVS)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PwN_large_vtavsVSText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i =1:length(txt)
PwN_large_vtavsVSText(i).String = combinedtxt(i);hold on; 
end 
axis off;hold on; 


subplot(2,3,6); hold on;
PwN_large_syncVS = pie(SumCount_PwN_large_syncVS(find(SumCount_PwN_large_syncVS~=0)));hold on;
PwN_large_syncVSText = findobj(PwN_large_syncVS,'Type','text'); % text object handles
percentValues = get(PwN_large_syncVSText,'String'); % percent values
txt = TextAsRep(find(SumCount_PwN_large_syncVS)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PwN_large_syncVSText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
PwN_large_syncVSText(i).String = combinedtxt(i);hold on; 
end
axis off;hold on; 

%% %%%%%%%%%%%% VTA Region with repetions
TextVTA= {'DAN: ';'GABA: ';'GLU:'; 'NoID:'};
figure(fignum(5));hold on;
subplot(2,3,1); hold on;
PwN_small_vsvtaVTA = pie(SumCount_PwN_small_vsvtaVTA(find(SumCount_PwN_small_vsvtaVTA~=0)));hold on;
PwN_small_vsvtaVTAText = findobj(PwN_small_vsvtaVTA,'Type','text'); % text object handles
percentValues = get(PwN_small_vsvtaVTAText,'String'); % percent values
txt = TextVTA(find(SumCount_PwN_small_vsvtaVTA)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PwN_small_vsvtaVTAText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
PwN_small_vsvtaVTAText(i).String = combinedtxt(i);hold on; 
end
% PwN_small_vsvtaVTAText(1).String = combinedtxt(1);hold on; 
% PwN_small_vsvtaVTAText(2).String = combinedtxt(2);hold on; 
% PwN_small_vsvtaVTAText(3).String = combinedtxt(3);hold on; 
% PwN_small_vsvtaVTAText(4).String = combinedtxt(4);hold on; 
% PVS_small_vsvtaText(4).String = combinedtxt(4);hold on; 
%  title('vs->vta')
axis off;hold on; 

subplot(2,3,2); hold on;
PwN_small_vtavsVTA = pie(SumCount_PwN_small_vtavsVTA(find(SumCount_PwN_small_vtavsVTA~=0)));hold on;
PwN_small_vtavsVTAText = findobj(PwN_small_vtavsVTA,'Type','text'); % text object handles
percentValues = get(PwN_small_vtavsVTAText,'String'); % percent values
txt = TextVTA(find(SumCount_PwN_small_vtavsVTA)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PwN_small_vtavsVTAText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
PwN_small_vtavsVTAText(i).String = combinedtxt(i);hold on; 
end
% PwN_small_vtavsVTAText(1).String = combinedtxt(1);hold on; 
% PwN_small_vtavsVTAText(2).String = combinedtxt(2);hold on; 
% PwN_small_vtavsVTAText(3).String = combinedtxt(3);hold on; 
% PwN_small_vtavsVTAText(4).String = combinedtxt(4);hold on; 
axis off;hold on; 


subplot(2,3,3); hold on;
PwN_small_syncVTA = pie(SumCount_PwN_small_syncVTA(find(SumCount_PwN_small_syncVTA~=0)));hold on;
PwN_small_syncVTAText = findobj(PwN_small_syncVTA,'Type','text'); % text object handles
percentValues = get(PwN_small_syncVTAText,'String'); % percent values
txt = TextVTA(find(SumCount_PwN_small_syncVTA)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PwN_small_syncVTAText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
PwN_small_syncVTAText(i).String = combinedtxt(i);hold on; 
end
% PwN_small_syncVTAText(1).String = combinedtxt(1);hold on; 
% PwN_small_syncVTAText(2).String = combinedtxt(2);hold on; 
% PwN_small_syncVTAText(3).String = combinedtxt(3);hold on; 
% PwN_small_syncVTAText(4).String = combinedtxt(4);hold on; 
axis off;hold on; 

subplot(2,3,4); hold on;
PwN_large_vsvtaVTA = pie(SumCount_PwN_large_vsvtaVTA(find(SumCount_PwN_large_vsvtaVTA~=0)));hold on;
PwN_large_vsvtaVTAText = findobj(PwN_large_vsvtaVTA,'Type','text'); % text object handles
percentValues = get(PwN_large_vsvtaVTAText,'String'); % percent values
txt = TextVTA(find(SumCount_PwN_large_vsvtaVTA)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PwN_large_vsvtaVTAText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
PwN_large_vsvtaVTAText(i).String = combinedtxt(i);hold on; 
end
% PwN_large_vsvtaVTAText(1).String = combinedtxt(1);hold on; 
% PwN_large_vsvtaVTAText(2).String = combinedtxt(2);hold on; 
% PwN_large_vsvtaVTAText(3).String = combinedtxt(3);hold on; 
% PwN_large_vsvtaVTAText(4).String = combinedtxt(4);hold on
axis off;hold on; 


subplot(2,3,5); hold on;
PwN_large_vtavsVTA = pie(SumCount_PwN_large_vtavsVTA(find(SumCount_PwN_large_vtavsVTA~=0)));hold on;
PwN_large_vtavsVTAText = findobj(PwN_large_vtavsVTA,'Type','text'); % text object handles
percentValues = get(PwN_large_vtavsVTAText,'String'); % percent values
txt = TextVTA(find(SumCount_PwN_large_vtavsVTA)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PwN_large_vtavsVTAText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
PwN_large_vtavsVTAText(i).String = combinedtxt(i);hold on; 
end
% PwN_large_vtavsVTAText(1).String = combinedtxt(1);hold on; 
% PwN_large_vtavsVTAText(2).String = combinedtxt(2);hold on; 
% PwN_large_vtavsVTAText(3).String = combinedtxt(3);hold on; 
% PwN_large_vtavsVTAText(4).String = combinedtxt(4);hold on; 
axis off;hold on; 


subplot(2,3,6); hold on;
PwN_large_syncVTA = pie(SumCount_PwN_large_syncVTA(find(SumCount_PwN_large_syncVTA~=0)));hold on;
PwN_large_syncVTAText = findobj(PwN_large_syncVTA,'Type','text'); % text object handles
percentValues = get(PwN_large_syncVTAText,'String'); % percent values
txt = TextVTA(find(SumCount_PwN_large_syncVTA)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(PwN_large_syncVTAText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
PwN_large_syncVTAText(i).String = combinedtxt(i);hold on; 
end
% PwN_large_syncVTAText(1).String = combinedtxt(1);hold on; 
% PwN_large_syncVTAText(2).String = combinedtxt(2);hold on; 
% PwN_large_syncVTAText(3).String = combinedtxt(3);hold on; 
% PwN_large_syncVTAText(4).String = combinedtxt(4);hold on; 
axis off;hold on; 

%%

% % % % %%%%%% VTA %%%%%%%%%
% % % % figure(fignum(6));hold on;
% % % % subplot(2,3,1); hold on;
% % % % PwN_small_vsvtaVTA = pie(SumCount_PwN_small_vsvtaVTA);hold on;
% % % % PwN_small_vsvtaVTAText = findobj(PwN_small_vsvtaVTA,'Type','text'); % text object handles
% % % % percentValues = get(PwN_small_vsvtaVTAText,'String'); % percent values
% % % % txt = {'Type I: ';'Type II: ';'Type III:'; 'No-Id:'}; % strings
% % % % combinedtxt = strcat(txt,percentValues); % strings and percent values
% % % % oldExtents_cell = get(PwN_small_vsvtaVTAText,'Extent'); % cell array
% % % % oldExtents = cell2mat(oldExtents_cell);
% % % % PwN_small_vsvtaVTAText(1).String = combinedtxt(1);hold on; 
% % % % PwN_small_vsvtaVTAText(2).String = combinedtxt(2);hold on; 
% % % % PwN_small_vsvtaVTAText(3).String = combinedtxt(3);hold on; 
% % % % PwN_small_vsvtaVTAText(4).String = combinedtxt(4);hold on; 
% % % % % PVS_small_vsvtaText(4).String = combinedtxt(4);hold on; 
% % % % %  title('vs->vta')
% % % % axis off;hold on; 
% % % % 
% % % % subplot(2,3,2); hold on;
% % % % PwN_small_vtavsVTA = pie(SumCount_PwN_small_vtavsVTA);hold on;
% % % % PwN_small_vtavsVTAText = findobj(PwN_small_vtavsVTA,'Type','text'); % text object handles
% % % % percentValues = get(PwN_small_vtavsVTAText,'String'); % percent values
% % % % txt = {'Type I: ';'Type II: ';'Type III:'; 'No-Id:'}; % strings
% % % % combinedtxt = strcat(txt,percentValues); % strings and percent values
% % % % oldExtents_cell = get(PwN_small_vtavsVTAText,'Extent'); % cell array
% % % % oldExtents = cell2mat(oldExtents_cell);
% % % % PwN_small_vtavsVTAText(1).String = combinedtxt(1);hold on; 
% % % % PwN_small_vtavsVTAText(2).String = combinedtxt(2);hold on; 
% % % % PwN_small_vtavsVTAText(3).String = combinedtxt(3);hold on; 
% % % % PwN_small_vtavsVTAText(4).String = combinedtxt(4);hold on; 
% % % % axis off;hold on; 
% % % % 
% % % % 
% % % % subplot(2,3,3); hold on;
% % % % PwN_small_syncVTA = pie(SumCount_PwN_small_syncVTA);hold on;
% % % % PwN_small_syncVTAText = findobj(PwN_small_syncVTA,'Type','text'); % text object handles
% % % % percentValues = get(PwN_small_syncVTAText,'String'); % percent values
% % % % txt = {'Type I: ';'Type II: ';'Type III:'; 'No-Id:'}; % strings
% % % % combinedtxt = strcat(txt,percentValues); % strings and percent values
% % % % oldExtents_cell = get(PwN_small_syncVTAText,'Extent'); % cell array
% % % % oldExtents = cell2mat(oldExtents_cell);
% % % % PwN_small_syncVTAText(1).String = combinedtxt(1);hold on; 
% % % % PwN_small_syncVTAText(2).String = combinedtxt(2);hold on; 
% % % % PwN_small_syncVTAText(3).String = combinedtxt(3);hold on; 
% % % % PwN_small_syncVTAText(4).String = combinedtxt(4);hold on; 
% % % % axis off;hold on; 
% % % % 
% % % % subplot(2,3,4); hold on;
% % % % PwN_large_vsvtaVTA = pie(SumCount_PwN_large_vsvtaVTA);hold on;
% % % % PwN_large_vsvtaVTAText = findobj(PwN_large_vsvtaVTA,'Type','text'); % text object handles
% % % % percentValues = get(PwN_large_vsvtaVTAText,'String'); % percent values
% % % % txt = {'Type I: ';'Type II: ';'Type III:'; 'No-Id:'}; % strings
% % % % combinedtxt = strcat(txt,percentValues); % strings and percent values
% % % % oldExtents_cell = get(PwN_large_vsvtaVTAText,'Extent'); % cell array
% % % % oldExtents = cell2mat(oldExtents_cell);
% % % % PwN_large_vsvtaVTAText(1).String = combinedtxt(1);hold on; 
% % % % PwN_large_vsvtaVTAText(2).String = combinedtxt(2);hold on; 
% % % % PwN_large_vsvtaVTAText(3).String = combinedtxt(3);hold on; 
% % % % PwN_large_vsvtaVTAText(4).String = combinedtxt(4);hold on
% % % % axis off;hold on; 
% % % % 
% % % % 
% % % % subplot(2,3,5); hold on;
% % % % PwN_large_vtavsVTA = pie(SumCount_PwN_large_vtavsVTA);hold on;
% % % % PwN_large_vtavsVTAText = findobj(PwN_large_vtavsVTA,'Type','text'); % text object handles
% % % % percentValues = get(PwN_large_vtavsVTAText,'String'); % percent values
% % % % txt = {'Type I: ';'Type II: ';'Type III:'; 'No-Id:'}; % strings
% % % % combinedtxt = strcat(txt,percentValues); % strings and percent values
% % % % oldExtents_cell = get(PwN_large_vtavsVTAText,'Extent'); % cell array
% % % % oldExtents = cell2mat(oldExtents_cell);
% % % % PwN_large_vtavsVTAText(1).String = combinedtxt(1);hold on; 
% % % % PwN_large_vtavsVTAText(2).String = combinedtxt(2);hold on; 
% % % % PwN_large_vtavsVTAText(3).String = combinedtxt(3);hold on; 
% % % % PwN_large_vtavsVTAText(4).String = combinedtxt(4);hold on; 
% % % % axis off;hold on; 
% % % % 
% % % % 
% % % % subplot(2,3,6); hold on;
% % % % PwN_large_syncVTA = pie(SumCount_PwN_large_syncVTA);hold on;
% % % % PwN_large_syncVTAText = findobj(PwN_large_syncVTA,'Type','text'); % text object handles
% % % % percentValues = get(PwN_large_syncVTAText,'String'); % percent values
% % % % txt = {'Type I: ';'Type II: ';'Type III:'; 'No-Id:'}; % strings
% % % % combinedtxt = strcat(txt,percentValues); % strings and percent values
% % % % oldExtents_cell = get(PwN_large_syncVTAText,'Extent'); % cell array
% % % % oldExtents = cell2mat(oldExtents_cell);
% % % % PwN_large_syncVTAText(1).String = combinedtxt(1);hold on; 
% % % % PwN_large_syncVTAText(2).String = combinedtxt(2);hold on; 
% % % % PwN_large_syncVTAText(3).String = combinedtxt(3);hold on; 
% % % % PwN_large_syncVTAText(4).String = combinedtxt(4);hold on; 
% % % % axis off;hold on; 
% % % % 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Conditional Probability = Joint Probability/ Probability to have assembly with a specific neuron-type
count_Joint_small_vsvta = zeros(SumTotAn, 16);
count_Joint_small_vtavs = zeros(SumTotAn, 16);
count_Joint_small_sync = zeros(SumTotAn, 16);
count_Joint_large_vsvta = zeros(SumTotAn, 16);
count_Joint_large_vtavs = zeros(SumTotAn, 16);
count_Joint_large_sync = zeros(SumTotAn, 16);


% Probability to have assembly with a specific neuron-type (already calculated but here repmat)
Region = 'vs';
[TabNumPwN_VS_small_vsvta,NumPwN_VS_small_vsvta] = PairsWSpecNeu_TabNum(Tab_PwN_small_vsvta,Region);  % Number of pairs with specific neurons repeated in such a way to do the ratio with the joint
[TabNumPwN_VS_small_vtavs,NumPwN_VS_small_vtavs] = PairsWSpecNeu_TabNum(Tab_PwN_small_vtavs,Region);
[TabNumPwN_VS_small_sync,NumPwN_VS_small_sync] = PairsWSpecNeu_TabNum(Tab_PwN_small_sync,Region);

[TabNumPwN_VS_large_vsvta,NumPwN_VS_large_vsvta] = PairsWSpecNeu_TabNum(Tab_PwN_large_vsvta,Region);  % Number of pairs with specific neurons repeated in such a way to do the ratio with the joint
[TabNumPwN_VS_large_vtavs,NumPwN_VS_large_vtavs] = PairsWSpecNeu_TabNum(Tab_PwN_large_vtavs,Region);
[TabNumPwN_VS_large_sync,NumPwN_VS_large_sync] = PairsWSpecNeu_TabNum(Tab_PwN_large_sync,Region);

% [TabNumPwN_VS_large,NumPwN_VS_large] = PairsWSpecNeu_TabNum(Tab_PwN_small_vtavs,Region);
 Region = 'vta';
[TabNumPwN_VTA_small_vsvta,NumPwN_VTA_small_vsvta] = PairsWSpecNeu_TabNum(Tab_PwN_small_vsvta,Region);  % Number of pairs with specific neurons repeated in such a way to do the ratio with the joint
[TabNumPwN_VTA_small_vtavs,NumPwN_VTA_small_vtavs] = PairsWSpecNeu_TabNum(Tab_PwN_small_vtavs,Region);
[TabNumPwN_VTA_small_sync,NumPwN_VTA_small_sync] = PairsWSpecNeu_TabNum(Tab_PwN_small_sync,Region);

[TabNumPwN_VTA_large_vsvta,NumPwN_VTA_large_vsvta] = PairsWSpecNeu_TabNum(Tab_PwN_large_vsvta,Region);  % Number of pairs with specific neurons repeated in such a way to do the ratio with the joint
[TabNumPwN_VTA_large_vtavs,NumPwN_VTA_large_vtavs] = PairsWSpecNeu_TabNum(Tab_PwN_large_vtavs,Region);
[TabNumPwN_VTA_large_sync,NumPwN_VTA_large_sync] = PairsWSpecNeu_TabNum(Tab_PwN_large_sync,Region);


% Joint Probability
for k=1:SumTotAn
[TLabelsJoint_small_vsvta,count_Joint_small_vsvta]=CountTypeNeuPairTwoSides(struct_pair_small_vsvta{k},k,count_Joint_small_vsvta); %joint probability to have vs_k and vta_k
[TLabelsJoint_small_vtavs,count_Joint_small_vtavs]=CountTypeNeuPairTwoSides(struct_pair_small_vtavs{k},k,count_Joint_small_vtavs);
[TLabelsJoint_small_sync,count_Joint_small_sync]=CountTypeNeuPairTwoSides(struct_pair_small_sync{k},k,count_Joint_small_sync);


[TLabelsJoint_large_vsvta,count_Joint_large_vsvta]=CountTypeNeuPairTwoSides(struct_pair_large_vsvta{k},k,count_Joint_large_vsvta);
[TLabelsJoint_large_vtavs,count_Joint_large_vtavs]=CountTypeNeuPairTwoSides(struct_pair_large_vtavs{k},k,count_Joint_large_vtavs);
[TLabelsJoint_large_sync,count_Joint_large_sync]=CountTypeNeuPairTwoSides(struct_pair_large_sync{k},k,count_Joint_large_sync);
end
%%
SumCount_Joint_small_vsvta = sum(count_Joint_small_vsvta);
SumCount_Joint_small_vtavs = sum(count_Joint_small_vtavs);
SumCount_Joint_small_sync = sum(count_Joint_small_sync);

SumCount_Joint_small=SumCount_Joint_small_vsvta + SumCount_Joint_small_vtavs + SumCount_Joint_small_sync


SumCount_Joint_large_vsvta = sum(count_Joint_large_vsvta);
SumCount_Joint_large_vtavs = sum(count_Joint_large_vtavs);
SumCount_Joint_large_sync = sum(count_Joint_large_sync);

SumCount_Joint_large=SumCount_Joint_large_vsvta+SumCount_Joint_large_vtavs+SumCount_Joint_large_sync;

[TabSumCount_Joint_small_vsvta] = MakeCoupleTable(SumCount_Joint_small_vsvta);
[TabSumCount_Joint_small_vtavs] = MakeCoupleTable(SumCount_Joint_small_vtavs);
[TabSumCount_Joint_small_sync] = MakeCoupleTable(SumCount_Joint_small_sync);

TabSumCount_Joint_small=TabSumCount_Joint_small_vsvta{:,:}+TabSumCount_Joint_small_vtavs{:,:}+TabSumCount_Joint_small_sync{:,:};


[TabSumCount_Joint_large_vsvta] = MakeCoupleTable(SumCount_Joint_large_vsvta);
[TabSumCount_Joint_large_vtavs] = MakeCoupleTable(SumCount_Joint_large_vtavs);
[TabSumCount_Joint_large_sync] = MakeCoupleTable(SumCount_Joint_large_sync);

TabSumCount_Joint_large=TabSumCount_Joint_large_vsvta{:,:}+TabSumCount_Joint_large_vtavs{:,:}+TabSumCount_Joint_large_sync{:,:};


pathfile='/zifnas/Carla/CWEM_Project_DATA/AssemblyOcc/';
filename='AssemblyOcc.mat';
save(fullfile(pathfile,filename),'SumCount_Joint_small_vsvta','SumCount_Joint_small_vtavs','SumCount_Joint_small_sync','SumCount_Joint_small',...
    'SumCount_Joint_large_vsvta','SumCount_Joint_large_vtavs','SumCount_Joint_large_sync','SumCount_Joint_large',...
    'TabSumCount_Joint_small_vsvta','TabSumCount_Joint_small_vtavs','TabSumCount_Joint_small_sync','TabSumCount_Joint_small',...
    'TabSumCount_Joint_large_vsvta','TabSumCount_Joint_large_vtavs','TabSumCount_Joint_large_sync','TabSumCount_Joint_large');



textCouple={'SPN-DAN: '; 'SPN-GABA: '; 'SPN-GLU: '; 'SPN-NoID: '; 'FSN-DAN: '; 'FSN-GABA: '; 'FSN-GLU: '; 'FSN-NoID: ';...
                 'CIN-DAN: '; 'CIN-GABA: '; 'CIN-GLU: ';'CIN-NoID: '; 'NoID-DAN: '; 'NoID-GABA: '; 'NoID-GLU: '; 'NoID-NoID: '};
% ZerosFind = find(SumCount_Joint_small_vsvta==0)


%%   No distinction directionality (only small bins)
SumAll_Count_Joint_small=sum([SumCount_Joint_small_vsvta;SumCount_Joint_small_vtavs; SumCount_Joint_small_sync]);
SumAll_Count_Joint_small_NoSync=sum([SumCount_Joint_small_vsvta;SumCount_Joint_small_vtavs]);
textCouple={'SPN-DAN'; 'SPN-GABA'; 'SPN-GLU'; 'SPN-NoID'; 'FSN-DAN'; 'FSN-GABA '; 'FSN-GLU'; 'FSN-NoID';...
                 'CIN-DAN'; 'CIN-GABA'; 'CIN-GLU';'CIN-NoID'; 'NoID-DAN'; 'NoID-GABA'; 'NoID-GLU'; 'NoID-NoID'};
figure(91);hold on;
Pie_All_PerType = pie(SumAll_Count_Joint_small(find(SumAll_Count_Joint_small~=0)));hold on;
% Pie_All_PerTypeText =findobj(Pie_All_PerType,'Type','text');
% percentValues = get(Pie_All_PerTypeText,'String'); % percent values
% txt = textCouple(find(SumAll_Count_Joint_small)); % strings
% combinedtxt = strcat(txt,percentValues); % strings and percent values
% oldExtents_cell = get(Pie_All_PerTypeText,'Extent'); % cell array
% oldExtents = cell2mat(oldExtents_cell);
% for i=1:length(txt)
% Pie_All_PerTypeText(i).String = combinedtxt(i);hold on; 
%  
% end
lgd1=legend(textCouple,'Location','eastoutside','Orientation','vertical')
lgd1.FontSize=14;
axis off;hold on;


%%

figure(90);hold on;
Pie_All_PerType1 = pie(SumAll_Count_Joint_small_NoSync(find(SumAll_Count_Joint_small_NoSync~=0)));hold on;
Pie_All_PerTypeText1 =findobj(Pie_All_PerType1,'Type','text');
percentValues = get(Pie_All_PerTypeText1,'String'); % percent values
txt = textCouple(find(SumAll_Count_Joint_small_NoSync)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(Pie_All_PerTypeText1,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
Pie_All_PerTypeText1(i).String = combinedtxt(i);hold on; 
 
end
axis off;hold on; 
%%
figure(fignum(7));hold on;
subplot(2,3,1); hold on;
P2S_small_vsvta = pie(SumCount_Joint_small_vsvta(find(SumCount_Joint_small_vsvta~=0)));hold on;
P2S_small_vsvtaText = findobj(P2S_small_vsvta,'Type','text'); % text object handles
percentValues = get(P2S_small_vsvtaText,'String'); % percent values
txt = textCouple(find(SumCount_Joint_small_vsvta)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(P2S_small_vsvtaText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
P2S_small_vsvtaText(i).String = combinedtxt(i);hold on; 
 
end
axis off;hold on; 

subplot(2,3,2); hold on;
P2S_small_vtavs = pie(SumCount_Joint_small_vtavs(find(SumCount_Joint_small_vtavs~=0)));hold on;
P2S_small_vtavsText = findobj(P2S_small_vtavs,'Type','text'); % text object handles
percentValues = get(P2S_small_vtavsText,'String'); % percent values
txt = textCouple(find(SumCount_Joint_small_vtavs)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(P2S_small_vtavsText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
P2S_small_vtavsText(i).String = combinedtxt(i);hold on; 
 
end
axis off;hold on; 

subplot(2,3,3); hold on;
P2S_small_sync = pie(SumCount_Joint_small_sync(find(SumCount_Joint_small_sync~=0)));hold on;
P2S_small_syncText = findobj(P2S_small_sync,'Type','text'); % text object handles
percentValues = get(P2S_small_syncText,'String'); % percent values
txt = textCouple(find(SumCount_Joint_small_sync)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(P2S_small_syncText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
P2S_small_syncText(i).String = combinedtxt(i);hold on; 
 
end
axis off;hold on; 

subplot(2,3,4); hold on;
P2S_large_vsvta = pie(SumCount_Joint_large_vsvta(find(SumCount_Joint_large_vsvta~=0)));hold on;
P2S_large_vsvtaText = findobj(P2S_large_vsvta,'Type','text'); % text object handles
percentValues = get(P2S_large_vsvtaText,'String'); % percent values
txt = textCouple(find(SumCount_Joint_large_vsvta)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(P2S_large_vsvtaText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
P2S_large_vsvtaText(i).String = combinedtxt(i);hold on; 
 
end
axis off;hold on; 

subplot(2,3,5); hold on;
P2S_large_vtavs = pie(SumCount_Joint_large_vtavs(find(SumCount_Joint_large_vtavs~=0)));hold on;
P2S_large_vtavsText = findobj(P2S_large_vtavs,'Type','text'); % text object handles
percentValues = get(P2S_large_vtavsText,'String'); % percent values
txt = textCouple(find(SumCount_Joint_large_vtavs)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(P2S_large_vtavsText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
P2S_large_vtavsText(i).String = combinedtxt(i);hold on; 
 
end
axis off;hold on; 


subplot(2,3,6); hold on;
P2S_large_sync = pie(SumCount_Joint_large_sync(find(SumCount_Joint_large_sync~=0)));hold on;
P2S_large_syncText = findobj(P2S_large_sync,'Type','text'); % text object handles
percentValues = get(P2S_large_syncText,'String'); % percent values
txt = textCouple(find(SumCount_Joint_large_sync)); % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(P2S_large_syncText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell);
for i=1:length(txt)
P2S_large_syncText(i).String = combinedtxt(i);hold on; 
 
end
axis off;hold on; 

%% %%%%%%%%%%




 %%%%%%%%%%%%%%%%% Distict Animals
% countMat_Pairs_small_vsvta=repmat(count_Pairs_small_vsvta,8,1);
% countMat_Pairs_small_vtavs=repmat(count_Pairs_small_vtavs,8,1);
% countMat_Pairs_small_sync=repmat(count_Pairs_small_sync,8,1);
% 
% countMat_Pairs_large_vsvta=repmat(count_Pairs_large_vsvta,8,1);
% countMat_Pairs_large_vtavs=repmat(count_Pairs_large_vtavs,8,1);
% countMat_Pairs_large_sync=repmat(count_Pairs_large_sync,8,1);
% 
% % Percentage of assemblies with a specific neuron
% 
% Perc_PwN_pAn_small_vsvta=count_PwN_small_vsvta./countMat_Pairs_small_vsvta'; % # assemblies with V_{ri}/ # assemblies Animal per animal
% Perc_PwN_pAn_small_vtavs=count_PwN_small_vtavs./countMat_Pairs_small_vtavs';
% Perc_PwN_pAn_small_sync=count_PwN_small_sync./countMat_Pairs_small_sync';
% 
% Perc_PwN_pAn_large_vsvta=count_PwN_large_vsvta./countMat_Pairs_large_vsvta'; % # assemblies with V_{ri}/ # assemblies Animal per animal
% Perc_PwN_pAn_large_vtavs=count_PwN_large_vtavs./countMat_Pairs_large_vtavs';
% Perc_PwN_pAn_large_sync=count_PwN_large_sync./countMat_Pairs_large_sync';
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % % figure(fignum(4));hold on;
% % % % subplot(2,3,1); hold on;
% % % % PwN_small_vsvta = pie(SumCount_PwN_small_vsvta(:,[1:3,5:8]));hold on;
% % % % PwN_small_vsvtaText = findobj(PwN_small_vsvta,'Type','text'); % text object handles
% % % % percentValues = get(PwN_small_vsvtaText,'String'); % percent values
% % % % txt = {'MSN: ';'FSI: ';'CNI:';'Type I:'; 'Type II:', 'Type III:', 'No-Id:'}; % strings
% % % % combinedtxt = strcat(txt,percentValues); % strings and percent values
% % % % oldExtents_cell = get(PwN_small_vsvtaText,'Extent'); % cell array
% % % % oldExtents = cell2mat(oldExtents_cell);
% % % % PwN_small_vsvtaText(1).String = combinedtxt(1);hold on; 
% % % % PwN_small_vsvtaText(2).String = combinedtxt(2);hold on; 
% % % % PwN_small_vsvtaText(3).String = combinedtxt(3);hold on; 
% % % % PwN_small_vsvtaText(4).String = combinedtxt(4);hold on; 
% % % % PwN_small_vsvtaText(5).String = combinedtxt(5);hold on; 
% % % % PwN_small_vsvtaText(6).String = combinedtxt(6);hold on; 
% % % % PwN_small_vsvtaText(7).String = combinedtxt(7);hold on; 
% % % % % PVS_small_vsvtaText(4).String = combinedtxt(4);hold on; 
% % % % %  title('vs->vta')
% % % % axis off;hold on; 
% % % % 
% % % % subplot(2,3,2); hold on;
% % % % PwN_small_vtavs = pie(SumCount_PwN_small_vtavs(:,[1:3,5:8]));hold on;
% % % % PwN_small_vtavsText = findobj(PwN_small_vtavs,'Type','text'); % text object handles
% % % % percentValues = get(PwN_small_vtavsText,'String'); % percent values
% % % % txt = {'MSN: ';'FSI: ';'CNI:';'Type I:'; 'Type II:', 'Type III:', 'No-Id:'}; % strings
% % % % combinedtxt = strcat(txt,percentValues); % strings and percent values
% % % % oldExtents_cell = get(PwN_small_vtavsText,'Extent'); % cell array
% % % % oldExtents = cell2mat(oldExtents_cell);
% % % % PwN_small_vtavsText(1).String = combinedtxt(1);hold on; 
% % % % PwN_small_vtavsText(2).String = combinedtxt(2);hold on; 
% % % % PwN_small_vtavsText(3).String = combinedtxt(3);hold on; 
% % % % PwN_small_vtavsText(4).String = combinedtxt(4);hold on; 
% % % % PwN_small_vtavsText(5).String = combinedtxt(5);hold on; 
% % % % PwN_small_vtavsText(6).String = combinedtxt(6);hold on; 
% % % % PwN_small_vtavsText(7).String = combinedtxt(7);hold on; 
% % % % % PVS_small_vsvtaText(4).String = combinedtxt(4);hold on; 
% % % % %  title('vs->vta')
% % % % axis off;hold on; 
% % % % 
% % % % 
% % % % subplot(2,3,3); hold on;
% % % % PwN_small_vtavs = pie(SumCount_PwN_small_vtavs(:,[1:3,5:8]));hold on;
% % % % PwN_small_vtavsText = findobj(PwN_small_vtavs,'Type','text'); % text object handles
% % % % percentValues = get(PwN_small_vtavsText,'String'); % percent values
% % % % txt = {'MSN: ';'FSI: ';'CNI:';'Type I:'; 'Type II:', 'Type III:', 'No-Id:'}; % strings
% % % % combinedtxt = strcat(txt,percentValues); % strings and percent values
% % % % oldExtents_cell = get(PwN_small_vtavsText,'Extent'); % cell array
% % % % oldExtents = cell2mat(oldExtents_cell);
% % % % PwN_small_vtavsText(1).String = combinedtxt(1);hold on; 
% % % % PwN_small_vtavsText(2).String = combinedtxt(2);hold on; 
% % % % PwN_small_vtavsText(3).String = combinedtxt(3);hold on; 
% % % % PwN_small_vtavsText(4).String = combinedtxt(4);hold on; 
% % % % PwN_small_vtavsText(5).String = combinedtxt(5);hold on; 
% % % % PwN_small_vtavsText(6).String = combinedtxt(6);hold on; 
% % % % PwN_small_vtavsText(7).String = combinedtxt(7);hold on; 
% % % % % PVS_small_vsvtaText(4).String = combinedtxt(4);hold on; 
% % % % %  title('vs->vta')
% % % % axis off;hold on; 
% % % % 
% % % % subplot(2,3,4); hold on;
% % % % PwN_large_vsvta = pie(SumCount_PwN_large_vsvta(:,[1:3,5:8]));hold on;
% % % % PwN_small_vtavsText = findobj(PwN_small_vtavs,'Type','text'); % text object handles
% % % % percentValues = get(PwN_small_vtavsText,'String'); % percent values
% % % % txt = {'MSN: ';'FSI: ';'CNI:';'Type I:'; 'Type II:', 'Type III:', 'No-Id:'}; % strings
% % % % combinedtxt = strcat(txt,percentValues); % strings and percent values
% % % % oldExtents_cell = get(PwN_small_vtavsText,'Extent'); % cell array
% % % % oldExtents = cell2mat(oldExtents_cell);
% % % % PwN_small_vtavsText(1).String = combinedtxt(1);hold on; 
% % % % PwN_small_vtavsText(2).String = combinedtxt(2);hold on; 
% % % % PwN_small_vtavsText(3).String = combinedtxt(3);hold on; 
% % % % PwN_small_vtavsText(4).String = combinedtxt(4);hold on; 
% % % % PwN_small_vtavsText(5).String = combinedtxt(5);hold on; 
% % % % PwN_small_vtavsText(6).String = combinedtxt(6);hold on; 
% % % % PwN_small_vtavsText(7).String = combinedtxt(7);hold on; 
% % % % % PVS_small_vsvtaText(4).String = combinedtxt(4);hold on; 
% % % % %  title('vs->vta')
% % % % axis off;hold on; 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






% % % % % % % figure(fignum(2));hold on;
% % % % % % % subplot(3,2,1); hold on;
% % % % % % % PVS_small_vsvta = pie(SumCountNeuP_small_vsvta_VS(:,1:3));hold on;
% % % % % % % PVS_small_vsvtaText = findobj(PVS_small_vsvta,'Type','text'); % text object handles
% % % % % % % percentValues = get(PVS_small_vsvtaText,'String'); % percent values
% % % % % % % txt = {'MSN: ';'FSI: ';'CNI:'}; % strings
% % % % % % % combinedtxt = strcat(txt,percentValues); % strings and percent values
% % % % % % % oldExtents_cell = get(PVS_small_vsvtaText,'Extent'); % cell array
% % % % % % % oldExtents = cell2mat(oldExtents_cell);
% % % % % % % PVS_small_vsvtaText(1).String = combinedtxt(1);hold on; 
% % % % % % % PVS_small_vsvtaText(2).String = combinedtxt(2);hold on; 
% % % % % % % PVS_small_vsvtaText(3).String = combinedtxt(3);hold on; 
% % % % % % % % PVS_small_vsvtaText(4).String = combinedtxt(4);hold on; 
% % % % % % % % title('Neurons Types in pairs; Small Bins v->vta')
% % % % % % % axis off;hold on; 
% % % % % % % 
% % % % % % % 
% % % % % % % subplot(3,2,3); hold on;
% % % % % % % PVS_small_vtavs = pie(SumCountNeuP_small_vtavs_VS(:,1:3));hold on;
% % % % % % % PVS_small_vtavsText = findobj(PVS_small_vtavs,'Type','text'); % text object handles
% % % % % % % percentValues = get(PVS_small_vtavsText,'String'); % percent values
% % % % % % % txt = {'MSN: ';'FSI: ';'CNI:'}; % strings
% % % % % % % combinedtxt = strcat(txt,percentValues); % strings and percent values
% % % % % % % oldExtents_cell = get(PVS_small_vtavsText,'Extent'); % cell array
% % % % % % % oldExtents = cell2mat(oldExtents_cell);
% % % % % % % PVS_small_vtavsText(1).String = combinedtxt(1);
% % % % % % % PVS_small_vtavsText(2).String = combinedtxt(2);
% % % % % % % PVS_small_vtavsText(3).String = combinedtxt(3);
% % % % % % % % PVS_small_vtavsText(4).String = combinedtxt(4);
% % % % % % % 
% % % % % % % axis off
% % % % % % % 
% % % % % % % 
% % % % % % % 
% % % % % % % subplot(3,2,5); hold on;
% % % % % % % PVS_small_sync = pie(SumCountNeuP_small_sync_VS(:,1:3));hold on;
% % % % % % % PVS_small_syncText = findobj(PVS_small_sync,'Type','text'); % text object handles
% % % % % % % percentValues = get(PVS_small_syncText,'String'); % percent values
% % % % % % % txt = {'MSN: ';'FSI: ';'CNI:'}; % strings
% % % % % % % combinedtxt = strcat(txt,percentValues); % strings and percent values
% % % % % % % oldExtents_cell = get(PVS_small_syncText,'Extent'); % cell array
% % % % % % % oldExtents = cell2mat(oldExtents_cell);
% % % % % % % PVS_small_syncText(1).String = combinedtxt(1);
% % % % % % % PVS_small_syncText(2).String = combinedtxt(2);
% % % % % % % PVS_small_syncText(3).String = combinedtxt(3);
% % % % % % % % PVS_small_syncText(4).String = combinedtxt(4);
% % % % % % % 
% % % % % % % axis off
% % % % % % % 
% % % % % % % subplot(3,2,2); hold on;
% % % % % % % PVS_large_vsvta = pie(SumCountNeuP_large_vsvta_VS);hold on;
% % % % % % % PVS_large_vsvtaText = findobj(PVS_large_vsvta,'Type','text'); % text object handles
% % % % % % % percentValues = get(PVS_large_vsvtaText,'String'); % percent values
% % % % % % % txt = {'MSN: ';'FSI: ';'CNI:';'NoId'}; % strings
% % % % % % % combinedtxt = strcat(txt,percentValues); % strings and percent values
% % % % % % % oldExtents_cell = get(PVS_large_vsvtaText,'Extent'); % cell array
% % % % % % % oldExtents = cell2mat(oldExtents_cell);
% % % % % % % PVS_large_vsvtaText(1).String = combinedtxt(1);
% % % % % % % PVS_large_vsvtaText(2).String = combinedtxt(2);
% % % % % % % PVS_large_vsvtaText(3).String = combinedtxt(3);
% % % % % % % PVS_large_vsvtaText(4).String = combinedtxt(4);
% % % % % % % 
% % % % % % % axis off
% % % % % % % 
% % % % % % % 
% % % % % % % subplot(3,2,4); hold on;
% % % % % % % PVS_large_vtavs = pie(SumCountNeuP_large_vtavs_VS(:,1:3));hold on;
% % % % % % % PVS_large_vtavsText = findobj(PVS_large_vtavs,'Type','text'); % text object handles
% % % % % % % percentValues = get(PVS_large_vtavsText,'String'); % percent values
% % % % % % % txt = {'MSN: ';'FSI: ';'CNI:'}; % strings
% % % % % % % combinedtxt = strcat(txt,percentValues); % strings and percent values
% % % % % % % oldExtents_cell = get(PVS_large_vtavsText,'Extent'); % cell array
% % % % % % % oldExtents = cell2mat(oldExtents_cell);
% % % % % % % PVS_large_vtavsText(1).String = combinedtxt(1);
% % % % % % % PVS_large_vtavsText(2).String = combinedtxt(2);
% % % % % % % PVS_large_vtavsText(3).String = combinedtxt(3);
% % % % % % % % PVS_large_vtavsText(4).String = combinedtxt(4);
% % % % % % % 
% % % % % % % axis off
% % % % % % % 
% % % % % % % 
% % % % % % % subplot(3,2,6); hold on;
% % % % % % % PVS_large_sync = pie(SumCountNeuP_large_sync_VS(:,1:3));hold on;
% % % % % % % PVS_large_syncText = findobj(PVS_large_sync,'Type','text'); % text object handles
% % % % % % % percentValues = get(PVS_large_syncText,'String'); % percent values
% % % % % % % txt = {'MSN: ';'FSI: ';'CNI:'}; % strings
% % % % % % % combinedtxt = strcat(txt,percentValues); % strings and percent values
% % % % % % % oldExtents_cell = get(PVS_large_syncText,'Extent'); % cell array
% % % % % % % oldExtents = cell2mat(oldExtents_cell);
% % % % % % % PVS_large_syncText(1).String = combinedtxt(1);
% % % % % % % PVS_large_syncText(2).String = combinedtxt(2);
% % % % % % % PVS_large_syncText(3).String = combinedtxt(3);
% % % % % % % % PVS_large_syncText(4).String = combinedtxt(4);
% % % % % % % 
% % % % % % % axis off



% % % % % figure(fignum(1));
% % % % % subplot(4,4,[1,2,5,6]); hold on;
% % % % % Pie_Vs_Tot = pie(CountTypeVS);hold on;
% % % % % Pie_Vs_TotText = findobj(Pie_Vs_Tot,'Type','text'); % text object handles
% % % % % percentValues = get(Pie_Vs_TotText,'String'); % percent values
% % % % % txt = {'MSN: ';'FSI: ';'CNI:'; 'NoId:'}; % strings
% % % % % combinedtxt = strcat(txt,percentValues); % strings and percent values
% % % % % oldExtents_cell = get(Pie_Vs_TotText,'Extent'); % cell array
% % % % % oldExtents = cell2mat(oldExtents_cell);
% % % % % Pie_Vs_TotText(1).String = combinedtxt(1);
% % % % % Pie_Vs_TotText(2).String = combinedtxt(2);
% % % % % Pie_Vs_TotText(3).String = combinedtxt(3);
% % % % % Pie_Vs_TotText(4).String = combinedtxt(4);
% % % % % % newExtents_cell = get(Pie_Vs_TotText,'Extent'); % cell array
% % % % % % newExtents = cell2mat(newExtents_cell); % numeric array 
% % % % % % width_change = newExtents(:,4)-oldExtents(:,4);
% % % % % % signValues = sign(oldExtents(:,1));
% % % % % % offset = signValues.*(width_change/2);
% % % % % % textPositions_cell = get(Pie_Vs_TotText,{'Position'}); % cell array
% % % % % % textPositions = cell2mat(textPositions_cell); % numeric array
% % % % % % textPositions(:,1) = textPositions(:,1) + offset; % add offset 
% % % % % % 
% % % % % % PieVsTotText(1).Position = textPositions(1,:);
% % % % % % PieVsTotText(2).Position = textPositions(2,:);
% % % % % % PieVsTotText(3).Position = textPositions(3,:);
% % % % % % PieVsTotText(4).Position = textPositions(4,:);
% % % % % % title('Pie Total VS Neurons')
% % % % % axis off
% % % % % 
% % % % % subplot(4,4,[3,4,7,8]); hold on;
% % % % % Pie_Vta_Tot = pie(CountTypeVTA);hold on;
% % % % % Pie_Vta_TotText = findobj(Pie_Vta_Tot,'Type','text'); % text object handles
% % % % % percentValues = get(Pie_Vta_TotText,'String'); % percent values
% % % % % txt = {'Type I: ';'Type II: ';'Type III:'; 'NoId:'}; % strings
% % % % % combinedtxt = strcat(txt,percentValues); % strings and percent values
% % % % % oldExtents_cell = get(Pie_Vta_TotText,'Extent'); % cell array
% % % % % oldExtents = cell2mat(oldExtents_cell);
% % % % % Pie_Vta_TotText(1).String = combinedtxt(1);
% % % % % Pie_Vta_TotText(2).String = combinedtxt(2);
% % % % % Pie_Vta_TotText(3).String = combinedtxt(3);
% % % % % Pie_Vta_TotText(4).String = combinedtxt(4);
% % % % % % title('Pie Total VTA Neurons')
% % % % % axis off; 
% % % % % 
% % % % % subplot(4,4,[10,11,14,15]); hold on;
% % % % % Pie_Neu_Tot = pie(CountTypeAll);hold on;
% % % % % Pie_Neu_TotText = findobj(Pie_Neu_Tot,'Type','text'); % text object handles
% % % % % percentValues = get(Pie_Neu_TotText,'String'); % percent values
% % % % % txt = {'MSN: ';'FSI: ';'CNI:'; 'NoIdVs:';'Type I: ';'Type II: ';'Type III:'; 'NoIdVta:'}; % strings
% % % % % combinedtxt = strcat(txt,percentValues); % strings and percent values
% % % % % oldExtents_cell = get(Pie_Neu_TotText,'Extent'); % cell array
% % % % % oldExtents = cell2mat(oldExtents_cell);
% % % % % Pie_Neu_TotText(1).String = combinedtxt(1);
% % % % % Pie_Neu_TotText(2).String = combinedtxt(2);
% % % % % Pie_Neu_TotText(3).String = combinedtxt(3);
% % % % % Pie_Neu_TotText(4).String = combinedtxt(4);
% % % % % Pie_Neu_TotText(5).String = combinedtxt(5);
% % % % % Pie_Neu_TotText(6).String = combinedtxt(6);
% % % % % Pie_Neu_TotText(7).String = combinedtxt(7);
% % % % % Pie_Neu_TotText(8).String = combinedtxt(8);
% % % % % % title('Pie Total VTA Neurons')
% % % % % axis off; 
