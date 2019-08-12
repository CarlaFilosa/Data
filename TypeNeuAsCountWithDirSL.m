%% Script to count whether ensembles do prefer specific cathegories of units
addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM');
% % % % % % % % The matrix new_data is the matrix that contains the info, in particular new_data.spikeT_BegEnd has the 
% % % % % % % % spiketrain cutted properly

% name1=strcat('Admixture_OrRevW.mat');
% % load(name1)
% alg= strcat('Stand_');
% % alg=strcat('Pru_');
% % name=strcat('An_Disp_',alg,name1);
% % % name=strcat('An_Disp_',name);
% % load(name);
% % name_bins=strcat('BinLag_',alg,name1);
% % load(name_bins)
% % MatBin=-MBinAnP;
% nameSmallBig=strcat('SmallBig_An_Disp',alg,name1)
% load(nameSmallBig)
% 
% 
% fignum1=4;
% fignum2=5;
% fignum3=6;
% 
% clear c
% load('classlist.mat');
% 
% 
% Region='vsvta';
% % % % % % Region='vsvs';
% % % % % % Region='vtavta';
% % % % % 
% 
% 
% nnr=nn(a);
% [new_data] = parId(new_data,a,nnr);
% 
% 
% for k=1:length(new_data.par)
%  [new_data] = labels(new_data,c,k);
%  end
%  clear tt i events licks laser info spikes map params spM spikes spike_region j k new_spM
% 
% %%%%%%%%%%%%%%%% Find number of units of the Striatum and VTA nneuS
% 
% new_spM=new_data;
% 
% % TotAn=size(new_spM.spike_regionNoId,2);
% TotAn=length(new_data.par);
% for k=1:TotAn
%     nneuVS{k}=0;
%     nneuVTA{k}=0;
% end
% 
% for k=1:TotAn
%     for j=1:size(new_spM.spike_regionNoId{k},2)
%         [nneuVS{k}] = countNeu(new_spM.spike_regionNoId{k}(1,j),1,nneuVS{k});
%         [nneuVTA{k}] = countNeu(new_spM.spike_regionNoId{k}(1,j),2,nneuVTA{k});
%     end
% 
% end
% %%
% %%%%% Find Pairs in one region and shared between two differnt regions
% % pairs diveded by regions with all info
% [pairs_r_small,pairs_vsvs_small,pairs_vsvta_small,pairs_vtavta_small] = PairsInfo(As_across_smallBins,As_order_smallBins,nneuVS,TotAn);
% [pairs_r_large,pairs_vsvs_large,pairs_vsvta_large,pairs_vtavta_large] = PairsInfo(As_across_largeBins,As_order_largeBins,nneuVS,TotAn);
% 
% clear a b A B i ii j k
% 
% [struct_pair_small] = PairsOrgByPair(pairs_vsvta_small);
% [struct_pair_large] = PairsOrgByPair(pairs_vsvta_large);
% 
% 
% [struct_pair_small] = PairNeuLabelsID(struct_pair_small,new_data);
% [struct_pair_large] = PairNeuLabelsID(struct_pair_large,new_data);
% copy_struct_pair_small=struct_pair_small;
% copy_struct_pair_large=struct_pair_large;
% clear struct_pair
%%%%%% Create the struct coerently with the directionality
%%%%%%%%%% Here the directionality doesn't cointan the info about the bins
%%%%%%%%%% because the pairs are already divided for bins

% % % % %Dir stays for "Directionality", can assume these values

 

path_p= '/zifnas/Carla/CWEM_Project_GenProg';
cd(path_p)
MergeDataGenPostPru_1;
MergeDataGenPostPru_2;

path_po='/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM';
cd(path_po)
copy_struct_pair_small=struct_pair_small;
copy_struct_pair_large=struct_pair_large;

Dir1 = 'vs->vta';
 Dir2 = 'vta->vs';
 Dir3 = 'sync';
 Dir4 = 'no';
TotAnPerS = TotAn;
clear TotAn;
TotAn = SumTotAn;
for k=1:TotAn
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


%%%%% Do the assemblies prefer specific units?
% Total Number of assemblies

for k=1:TotAn
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
count_PwN_small_vsvta=zeros(TotAn,8);
 for k=1:TotAn
[Tab_PwN_small_vsvta,count_PwN_small_vsvta] = CountTypeNeuPair(struct_pair_small_vsvta{k},k,count_PwN_small_vsvta); % pairs with specific units, animal per animal
 end
 
 count_PwN_small_vtavs=zeros(TotAn,8);
 for k=1:TotAn
[Tab_PwN_small_vtavs,count_PwN_small_vtavs] = CountTypeNeuPair(struct_pair_small_vtavs{k},k,count_PwN_small_vtavs); % pairs with specific units, animal per animal
 end
 
 count_PwN_small_sync=zeros(TotAn,8);
 for k=1:TotAn
[Tab_PwN_small_sync,count_PwN_small_sync] = CountTypeNeuPair(struct_pair_small_sync{k},k,count_PwN_small_sync); % pairs with specific units, animal per animal
 end
 
 %%% Large bins different directionality
count_PwN_large_vsvta=zeros(TotAn,8);
 for k=1:TotAn
[Tab_PwN_large_vsvta,count_PwN_large_vsvta] = CountTypeNeuPair(struct_pair_large_vsvta{k},k,count_PwN_large_vsvta); % pairs with specific units, animal per animal
 end
 
 count_PwN_large_vtavs=zeros(TotAn,8);
 for k=1:TotAn
[Tab_PwN_large_vtavs,count_PwN_large_vtavs] = CountTypeNeuPair(struct_pair_large_vtavs{k},k,count_PwN_large_vtavs); % pairs with specific units, animal per animal
 end
 
 count_PwN_large_sync=zeros(TotAn,8);
 for k=1:TotAn
[Tab_PwN_large_sync,count_PwN_large_sync] = CountTypeNeuPair(struct_pair_large_sync{k},k,count_PwN_large_sync); % pairs with specific units, animal per animal
 end
 
 
 %%%%%%%%%%%%%%%%% Distict Animals
countMat_Pairs_small_vsvta=repmat(count_Pairs_small_vsvta,8,1);
countMat_Pairs_small_vtavs=repmat(count_Pairs_small_vtavs,8,1);
countMat_Pairs_small_sync=repmat(count_Pairs_small_sync,8,1);

countMat_Pairs_large_vsvta=repmat(count_Pairs_large_vsvta,8,1);
countMat_Pairs_large_vtavs=repmat(count_Pairs_large_vtavs,8,1);
countMat_Pairs_large_sync=repmat(count_Pairs_large_sync,8,1);

% Percentage of assemblies with a specific neuron

Perc_PwN_pAn_small_vsvta=count_PwN_small_vsvta./countMat_Pairs_small_vsvta'; % # assemblies with V_{ri}/ # assemblies Animal per animal
Perc_PwN_pAn_small_vtavs=count_PwN_small_vtavs./countMat_Pairs_small_vtavs';
Perc_PwN_pAn_small_sync=count_PwN_small_sync./countMat_Pairs_small_sync';

Perc_PwN_pAn_large_vsvta=count_PwN_large_vsvta./countMat_Pairs_large_vsvta'; % # assemblies with V_{ri}/ # assemblies Animal per animal
Perc_PwN_pAn_large_vtavs=count_PwN_large_vtavs./countMat_Pairs_large_vtavs';
Perc_PwN_pAn_large_sync=count_PwN_large_sync./countMat_Pairs_large_sync';


% Function that count how many neurons of specific cathegories are present
new_data=Data_All;
count_Neu=zeros(TotAn,8);
for k=1:TotAn
[TabCount_Neu,count_Neu] = CountTypeNeuAll(Data_All.par{k},k,count_Neu);
end


ChLevNeu_pAn= NaN(TotAn,8);
copy_nneuVS=nneuVS_All;
copy_nneuVTA=nneuVTA_All;
for k=1:TotAn
    if copy_nneuVS{k}==0, copy_nneuVS{k}=nan; end
    if copy_nneuVTA{k}==0, copy_nneuVTA{k}=nan; end
    ChLevNeu_pAn(k,1:4) = count_Neu(k,1:4)./copy_nneuVS{k};
    ChLevNeu_pAn(k,5:8) = count_Neu(k,5:8)./copy_nneuVTA{k};
end

% Ratio with chance level 
RatioAs_small_vsvta=Perc_PwN_pAn_small_vsvta./ChLevNeu_pAn;
RatioAs_small_vtavs=Perc_PwN_pAn_small_vtavs./ChLevNeu_pAn;
RatioAs_small_sync=Perc_PwN_pAn_small_sync./ChLevNeu_pAn;


RatioAs_large_vsvta=Perc_PwN_pAn_large_vsvta./ChLevNeu_pAn;
RatioAs_large_vtavs=Perc_PwN_pAn_large_vtavs./ChLevNeu_pAn;
RatioAs_large_sync=Perc_PwN_pAn_large_sync./ChLevNeu_pAn;



%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot AS with Directionality
%%
Text=strcat('Post Pruned---',nameT);
type{1}=strcat('MSN');type{2}=strcat('FSI');type{3}=strcat('CIN');type{4}= strcat('NoId-Vs');type{5}= strcat('type1');type{6}=strcat('type2');
type{7}=strcat('type3');type{8}=strcat('NoId-Vta');

ysup=max([max(max(Perc_PwN_pAn_small_vsvta./ChLevNeu_pAn)),max(max(Perc_PwN_pAn_small_vtavs./ChLevNeu_pAn)),max(max(Perc_PwN_pAn_small_sync./ChLevNeu_pAn)),...
    max(max(Perc_PwN_pAn_large_vsvta./ChLevNeu_pAn)),max(max(Perc_PwN_pAn_large_vtavs./ChLevNeu_pAn)),max(max(Perc_PwN_pAn_large_sync./ChLevNeu_pAn))]);
yinf =min([min(min(Perc_PwN_pAn_small_vsvta./ChLevNeu_pAn)),min(min(Perc_PwN_pAn_small_vtavs./ChLevNeu_pAn)),min(min(Perc_PwN_pAn_small_sync./ChLevNeu_pAn)),...
    min(min(Perc_PwN_pAn_large_vsvta./ChLevNeu_pAn)),min(min(Perc_PwN_pAn_large_vtavs./ChLevNeu_pAn)),min(min(Perc_PwN_pAn_large_sync./ChLevNeu_pAn))]);
figure(10);hold on;
% boxplot(Perc_PwN_pAn)
subplot(3,2,1)
for i=1:size(ChLevNeu_pAn,2)
    text(i-0.2,ysup,type{i});hold on;
end
plot([0 size(ChLevNeu_pAn,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((RatioAs_small_vsvta)','b*');hold on;
boxplot(RatioAs_small_vsvta);hold on;
ylim([(yinf-0.5) (ysup+0.5)]); hold on;
title('Small Bins vs->vta');

subplot(3,2,2)

for i=1:size(ChLevNeu_pAn,2)
    text(i-0.2,ysup,type{i});hold on;
end
plot([0 size(ChLevNeu_pAn,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((RatioAs_large_vsvta)','b*');hold on;
boxplot(RatioAs_large_vsvta);hold on;
ylim([(yinf-0.5) (ysup+0.5)]);hold on;
title('Large Bins vs->vta')
% ylim([(yinf-0.5) (ysup+0.5)]);

subplot(3,2,3)
for i=1:size(ChLevNeu_pAn,2)
    text(i-0.2,ysup,type{i});hold on;
end
plot([0 size(ChLevNeu_pAn,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((RatioAs_small_vtavs)','b*');hold on;
boxplot(RatioAs_small_vtavs);hold on;
ylim([(yinf-0.5) (ysup+0.5)]);hold on;
text(-0.5,0,Text,'Rotation',90,'Fontsize',14);hold on;
title('Small Bins vta->vs')
% ylim([(yinf-0.5) (ysup+0.5)]);

subplot(3,2,4)
for i=1:size(ChLevNeu_pAn,2)
    text(i-0.2,ysup,type{i});hold on;
end
plot([0 size(ChLevNeu_pAn,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((RatioAs_large_vtavs)','b*');hold on;
boxplot(RatioAs_large_vtavs);hold on;
ylim([(yinf-0.5) (ysup+0.5)]);hold on;
title('Large Bins vta->vs')


subplot(3,2,5)
for i=1:size(ChLevNeu_pAn,2)
    text(i-0.2,ysup,type{i});hold on;
end
plot([0 size(ChLevNeu_pAn,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((RatioAs_small_sync)','b*');hold on;
boxplot(RatioAs_small_sync);hold on;
ylim([(yinf-0.5) (ysup+0.5)]);hold on;
title('Small Bins sync')


subplot(3,2,6)
for i=1:size(ChLevNeu_pAn,2)
    text(i-0.2,ysup,type{i});hold on;
end
plot([0 size(ChLevNeu_pAn,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((RatioAs_large_sync)','b*');hold on;
boxplot(RatioAs_large_sync);hold on; 
ylim([(yinf-0.5) (ysup+0.5)]);hold on;
title('Large Bins sync')

box on;
%%
%%%%% Are there typologies of neurons most probably deputed than others to form assemblies?

% functions to count how many neurons of a cathegories are part of one
% assembly (the neurons are taken only once), the first funtion eliminates
% double elements, the second one counts

for k=1:TotAn
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


CountNeuP_small_vsvta=zeros(TotAn,8);
for k=1:TotAn
    if~isempty(Un_small_vsvta{k}.neuID)
        for i=1:length(Un_small_vsvta{k}.neuID)
[TabNeuP_small_vsvta,CountNeuP_small_vsvta] = CountTypeNeuGen(Un_small_vsvta{k}.labels(i),k,CountNeuP_small_vsvta); 
        end
    end
end


CountNeuP_small_vtavs=zeros(TotAn,8);
for k=1:TotAn
    if~isempty(Un_small_vtavs{k}.neuID)
        for i=1:length(Un_small_vtavs{k}.neuID)
[TabNeuP_small_vtavs,CountNeuP_small_vtavs] = CountTypeNeuGen(Un_small_vtavs{k}.labels(i),k,CountNeuP_small_vtavs); 
        end
    end
end



CountNeuP_small_sync=zeros(TotAn,8);
for k=1:TotAn
    if~isempty(Un_small_sync{k}.neuID)
        for i=1:length(Un_small_sync{k}.neuID)
[TabNeuP_small_sync,CountNeuP_small_sync] = CountTypeNeuGen(Un_small_sync{k}.labels(i),k,CountNeuP_small_sync); 
        end
    end
end


%%%%%% Large Bins


CountNeuP_large_vsvta=zeros(TotAn,8);
for k=1:TotAn
    if~isempty(Un_large_vsvta{k}.neuID)
        for i=1:length(Un_large_vsvta{k}.neuID)
[TabNeuP_large_vsvta,CountNeuP_large_vsvta] = CountTypeNeuGen(Un_large_vsvta{k}.labels(i),k,CountNeuP_large_vsvta); 
        end
    end
end


CountNeuP_large_vtavs=zeros(TotAn,8);
for k=1:TotAn
    if~isempty(Un_large_vtavs{k}.neuID)
        for i=1:length(Un_large_vtavs{k}.neuID)
[TabNeuP_large_vtavs,CountNeuP_large_vtavs] = CountTypeNeuGen(Un_large_vtavs{k}.labels(i),k,CountNeuP_large_vtavs); 
        end
    end
end



CountNeuP_large_sync=zeros(TotAn,8);
for k=1:TotAn
    if~isempty(Un_large_sync{k}.neuID)
        for i=1:length(Un_large_sync{k}.neuID)
[TabNeuP_large_sync,CountNeuP_large_sync] = CountTypeNeuGen(Un_large_sync{k}.labels(i),k,CountNeuP_large_sync); 
        end
    end
end


% Function that count how many neurons of specific cathegories are present

count_Neu=zeros(TotAn,8);
for k=1:TotAn
[TabCount_Neu,count_Neu] = CountTypeNeuAll(new_data.par{k},k,count_Neu);
end


%%%%%% Make sum and ratio for each categories of neurons
%%%%%%% Goal:  #neu_c in any assembly/ #neu_c Tot 
%%%%%%% To see whether there exist typologies of neurons more deputed to
%%%%%%% take part in any assembly

SumNeuP_small_vsvta = sum(TabNeuP_small_vsvta{:,:});
SumNeuP_small_vtavs = sum(TabNeuP_small_vtavs{:,:});
SumNeuP_small_sync = sum(TabNeuP_small_sync{:,:});

SumNeuP_large_vsvta = sum(TabNeuP_large_vsvta{:,:});
SumNeuP_large_vtavs = sum(TabNeuP_large_vtavs{:,:});
SumNeuP_large_sync = sum(TabNeuP_large_sync{:,:});
Sum_Neu = sum(TabCount_Neu{:,:});
% other way to do that
% for i=1: size(count_Neu,2)
% Sum_NinP(i) = sum(TLabels_NinP{find(count_Neu(:,i)>0),i});
% Sum_Neu(i) = sum(TLabels_Neu{find(count_Neu(:,i)>0),i});
% 
% end


%%% Considering ALL ANIMAL NO MEAN
%%% small Bins
PercNeuP_small_vsvta_Tog = SumNeuP_small_vsvta./Sum_Neu; % Here I'm considering everything omiting nan results Take the sum and no distiction for animal
PercNeuP_small_vtavs_Tog = SumNeuP_small_vtavs./Sum_Neu; % Here I'm considering everything omiting nan results Take the sum and no distiction for animal
PercNeuP_small_sync_Tog = SumNeuP_small_sync./Sum_Neu; % Here I'm considering everything omiting nan results Take the sum and no distiction for animal

%%%%%% large bins
PercNeuP_large_vsvta_Tog = SumNeuP_large_vsvta./Sum_Neu; % Here I'm considering everything omiting nan results Take the sum and no distiction for animal
PercNeuP_large_vtavs_Tog = SumNeuP_large_vtavs./Sum_Neu; % Here I'm considering everything omiting nan results Take the sum and no distiction for animal
PercNeuP_large_sync_Tog = SumNeuP_large_sync./Sum_Neu; % Here I'm considering everything omiting nan results Take the sum and no distiction for animal


%%%%% EACH ANIMAL DISTINCT and Then the MEAN!!!
PercNeuP_small_vsvta_Dist = CountNeuP_small_vsvta./count_Neu;  % divided for animal, then I can take the mean %%% !!!!!This one is coerent
PercNeuP_small_vtavs_Dist = CountNeuP_small_vtavs./count_Neu;  % divided for animal, then I can take the mean %%% !!!!!This one is coerent
PercNeuP_small_sync_Dist = CountNeuP_small_sync./count_Neu;  % divided for animal, then I can take the mean %%% !!!!!This one is coerent



PercNeuP_large_vsvta_Dist = CountNeuP_large_vsvta./count_Neu;  % divided for animal, then I can take the mean %%% !!!!!This one is coerent
PercNeuP_large_vtavs_Dist = CountNeuP_large_vtavs./count_Neu;  % divided for animal, then I can take the mean %%% !!!!!This one is coerent
PercNeuP_large_sync_Dist = CountNeuP_large_sync./count_Neu;  % divided for animal, then I can take the mean %%% !!!!!This one is coerent

MeanPerNeuP_small_vsvta = mean(PercNeuP_small_vsvta_Dist,1,'omitnan');
MeanPerNeuP_small_vtavs = mean(PercNeuP_small_vtavs_Dist,1,'omitnan');
MeanPerNeuP_small_sync = mean(PercNeuP_small_sync_Dist,1,'omitnan');


MeanPerNeuP_large_vsvta = mean(PercNeuP_large_vsvta_Dist,1,'omitnan');
MeanPerNeuP_large_vtavs = mean(PercNeuP_large_vtavs_Dist,1,'omitnan');
MeanPerNeuP_large_sync = mean(PercNeuP_large_sync_Dist,1,'omitnan');


StdPerNeuP_small_vsvta = std(PercNeuP_small_vsvta_Dist,0,1,'omitnan')/sqrt(size(PercNeuP_small_vsvta_Dist,1));
StdPerNeuP_small_vtavs = std(PercNeuP_small_vtavs_Dist,0,1,'omitnan')/sqrt(size(PercNeuP_small_vtavs_Dist,1));
StdPerNeuP_small_sync = std(PercNeuP_small_sync_Dist,1,'omitnan')/sqrt(size(PercNeuP_small_sync_Dist,1));


StdPerNeuP_large_vsvta = std(PercNeuP_large_vsvta_Dist,1,'omitnan')/sqrt(size(PercNeuP_large_vsvta_Dist,1));
StdPerNeuP_large_vtavs = std(PercNeuP_large_vtavs_Dist,1,'omitnan')/sqrt(size(PercNeuP_large_vtavs_Dist,1));
StdPerNeuP_large_sync = std(PercNeuP_large_sync_Dist,1,'omitnan')/sqrt(size(PercNeuP_large_sync_Dist,1));

ZeroErr=zeros(size(StdPerNeuP_small_vsvta));
%% %%%%%%%%%% Plot Neurons types
Text=strcat('Post Pruned---',nameT);
% Text=strcat(name1);
type{1}=strcat('MSN');type{2}=strcat('FSI');type{3}=strcat('CIN');type{4}= strcat('NoId-Vs');type{5}= strcat('type1');type{6}=strcat('type2');
type{7}=strcat('type3');type{8}=strcat('NoId-Vta');
 ysup=max(max([(MeanPerNeuP_small_vsvta + StdPerNeuP_small_vsvta);PercNeuP_small_vsvta_Tog;...
     (MeanPerNeuP_small_vtavs + StdPerNeuP_small_vtavs);PercNeuP_small_vtavs_Tog;... 
    (MeanPerNeuP_small_sync + StdPerNeuP_small_sync);PercNeuP_small_sync_Tog;...
    (MeanPerNeuP_large_vsvta + StdPerNeuP_large_vsvta);PercNeuP_large_vsvta_Tog;...
     (MeanPerNeuP_large_vtavs + StdPerNeuP_large_vtavs);PercNeuP_large_vtavs_Tog;... 
    (MeanPerNeuP_large_sync + StdPerNeuP_large_sync);PercNeuP_large_sync_Tog]))+0.2;

figure(11);hold on;
% boxplot(Perc_PwN_pAn)

 subplot(3,2,1)
for i=1:size(PercNeuP_small_vsvta_Tog,2)
    text(i-0.1,max(PercNeuP_small_vsvta_Tog(i),(MeanPerNeuP_small_vsvta(i)+StdPerNeuP_small_vsvta(i)))+0.1,type{i}),hold on;
end
ylim([0,ysup])
xlim([0,size(PercNeuP_small_vsvta_Tog,2)+1])
conf=[MeanPerNeuP_small_vsvta;PercNeuP_small_vsvta_Tog];
ErrConf = [StdPerNeuP_small_vsvta; ZeroErr];
b=bar(conf','group');hold on;
b(1).FaceColor='c';hold on;
b(2).FaceColor='y';hold on;
for i=1:size(PercNeuP_small_vsvta_Tog,2)
er=errorbar(i-0.15,MeanPerNeuP_small_vsvta(:,i),StdPerNeuP_small_vsvta(:,i),'k*');hold on;
end
 legend([b(1) b(2)],'Mean','AnTog','Location','northwest');
box on;hold on;
title('Neuron types. vs->vta Small Bins');
  
  subplot(3,2,3)
for i=1:size(PercNeuP_small_vtavs_Tog,2)
    text(i-0.2,max(PercNeuP_small_vtavs_Tog(i),(MeanPerNeuP_small_vtavs(i)+StdPerNeuP_small_vtavs(i)))+0.1,type{i}), hold on;
end
ylim([0,ysup])
xlim([0,size(PercNeuP_small_vtavs_Tog,2)+1])
conf=[MeanPerNeuP_small_vtavs;PercNeuP_small_vtavs_Tog];
b=bar(conf','group');hold on;
b(1).FaceColor='c';hold on;
b(2).FaceColor='y';hold on;
for i=1:size(PercNeuP_small_vsvta_Tog,2)
er=errorbar(i-0.15,MeanPerNeuP_small_vtavs(:,i),StdPerNeuP_small_vtavs(:,i),'k*');hold on;
end
legend([b(1) b(2)],'Mean','AnTog','Location','northwest');
text(-1,0,Text,'Rotation',90,'Fontsize',14);hold on;
box on;
title('Neuron types. vta->vs Small Bins');

subplot(3,2,5)
for i=1:size(PercNeuP_small_sync_Tog,2)
    text(i-0.2,max(PercNeuP_small_sync_Tog(i),(MeanPerNeuP_small_sync(i)+StdPerNeuP_small_sync(i)))+0.1,type{i}), hold on;
end
ylim([0,ysup])
xlim([0,size(PercNeuP_small_sync_Tog,2)+1])
conf=[MeanPerNeuP_small_sync;PercNeuP_small_sync_Tog];
b=bar(conf','group');hold on;
b(1).FaceColor='c';hold on;
b(2).FaceColor='y';hold on;
for i=1:size(PercNeuP_small_vsvta_Tog,2)
er=errorbar(i-0.15,MeanPerNeuP_small_sync(:,i),StdPerNeuP_small_sync(:,i),'k*');hold on;
end
legend([b(1) b(2)],'Mean','AnTog','Location','northwest');box on;
title('Neuron types. Sync Small Bins');



subplot(3,2,2)
for i=1:size(PercNeuP_large_vsvta_Tog,2)
    text(i-0.1,max(PercNeuP_large_vsvta_Tog(i),(MeanPerNeuP_large_vsvta(i)+StdPerNeuP_large_vsvta(i)))+0.1,type{i}),hold on;
end
ylim([0,ysup])
xlim([0,size(PercNeuP_large_vsvta_Tog,2)+1])
conf=[MeanPerNeuP_large_vsvta;PercNeuP_large_vsvta_Tog];
b=bar(conf','group');hold on;
b(1).FaceColor='c';hold on;
b(2).FaceColor='y';hold on;
for i=1:size(PercNeuP_large_vsvta_Tog,2)
er=errorbar(i-0.15,MeanPerNeuP_large_vsvta(:,i),StdPerNeuP_large_vsvta(:,i),'k*');hold on;
end
legend([b(1) b(2)],'Mean','AnTog','Location','northwest'); box on;
title('Neuron types. vs->vta Large Bins');
  
  subplot(3,2,4)
for i=1:size(PercNeuP_large_vtavs_Tog,2)
    text(i-0.2,max(PercNeuP_large_vtavs_Tog(i),(MeanPerNeuP_large_vtavs(i)+StdPerNeuP_large_vtavs(i)))+0.1,type{i}), hold on;
end
ylim([0,ysup])
xlim([0,size(PercNeuP_large_vtavs_Tog,2)+1])
conf=[MeanPerNeuP_large_vtavs;PercNeuP_large_vtavs_Tog];
b=bar(conf','group');hold on;
b(1).FaceColor='c';hold on;
b(2).FaceColor='y';hold on;
for i=1:size(PercNeuP_large_vsvta_Tog,2)
er=errorbar(i-0.15,MeanPerNeuP_large_vtavs(:,i),StdPerNeuP_large_vtavs(:,i),'k*');hold on;
end
legend([b(1) b(2)],'Mean','AnTog','Location','northwest'); box on;
title('Neuron types. vta->vs Large Bins');

subplot(3,2,6)
for i=1:size(PercNeuP_large_sync_Tog,2)
    text(i-0.2,max(PercNeuP_large_sync_Tog(i),(MeanPerNeuP_large_sync(i)+StdPerNeuP_large_sync(i)))+0.1,type{i}), hold on;
end
ylim([0,ysup])
xlim([0,size(PercNeuP_large_sync_Tog,2)+1])
conf=[MeanPerNeuP_large_sync;PercNeuP_large_sync_Tog];
b=bar(conf','group');hold on;
b(1).FaceColor='c';hold on;
b(2).FaceColor='y';hold on;
for i=1:size(PercNeuP_large_sync_Tog,2)
er=errorbar(i-0.15,MeanPerNeuP_large_sync(:,i),StdPerNeuP_large_sync(:,i),'k*');hold on;
end
legend([b(1) b(2)],'Mean','AnTog','Location','northwest'); box on; hold on;
title('Neuron types. Sync Large Bins');

hold on;



%%

%%%%%%%%%%%%%% Conditional Probability = Joint Probability/ Probability to have assembly with a specific neuron-type
count_Joint_small_vsvta = zeros(TotAn, 16);
count_Joint_small_vtavs = zeros(TotAn, 16);
count_Joint_small_sync = zeros(TotAn, 16);
count_Joint_large_vsvta = zeros(TotAn, 16);
count_Joint_large_vtavs = zeros(TotAn, 16);
count_Joint_large_sync = zeros(TotAn, 16);
% count_Joint_VTA_VS_small = zeros(TotAn, 16);
% count_Joint_VTA_VS_large = zeros(TotAn, 16);


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
for k=1:TotAn
[TLabelsJoint_small_vsvta,count_Joint_small_vsvta]=CountTypeNeuPairTwoSides(struct_pair_small_vsvta{k},k,count_Joint_small_vsvta); %joint probability to have vs_k and vta_k
[TLabelsJoint_small_vtavs,count_Joint_small_vtavs]=CountTypeNeuPairTwoSides(struct_pair_small_vtavs{k},k,count_Joint_small_vtavs);
[TLabelsJoint_small_sync,count_Joint_small_sync]=CountTypeNeuPairTwoSides(struct_pair_small_sync{k},k,count_Joint_small_sync);


[TLabelsJoint_large_vsvta,count_Joint_large_vsvta]=CountTypeNeuPairTwoSides(struct_pair_large_vsvta{k},k,count_Joint_large_vsvta);
[TLabelsJoint_large_vtavs,count_Joint_large_vtavs]=CountTypeNeuPairTwoSides(struct_pair_large_vtavs{k},k,count_Joint_large_vtavs);
[TLabelsJoint_large_sync,count_Joint_large_sync]=CountTypeNeuPairTwoSides(struct_pair_large_sync{k},k,count_Joint_large_sync);
end
% Conditional Probability vs|vta
given='vs|vta';
[TabCondProb_VS_VTA_small_vsvta,CondProb_VS_VTA_small_vsvta] = CondProbAsType(count_Joint_small_vsvta,NumPwN_VTA_small_vsvta,given);
[TabCondProb_VS_VTA_small_vtavs,CondProb_VS_VTA_small_vtavs] = CondProbAsType(count_Joint_small_vtavs,NumPwN_VTA_small_vtavs,given);
[TabCondProb_VS_VTA_small_sync,CondProb_VS_VTA_small_sync] = CondProbAsType(count_Joint_small_sync,NumPwN_VTA_small_sync,given);
[TabCondProb_VS_VTA_large_vsvta,CondProb_VS_VTA_large_vsvta] = CondProbAsType(count_Joint_large_vsvta,NumPwN_VTA_large_vsvta,given);
[TabCondProb_VS_VTA_large_vtavs,CondProb_VS_VTA_large_vtavs] = CondProbAsType(count_Joint_large_vtavs,NumPwN_VTA_large_vtavs,given);
[TabCondProb_VS_VTA_large_sync,CondProb_VS_VTA_large_sync] = CondProbAsType(count_Joint_large_sync,NumPwN_VTA_large_sync,given);


% Conditional probability vta|vs
% given= 'vta|vs' ;
% [TabCondProb_VTA_VS_small,CondProb_VTA_VS_small] = CondProbAsType(count_Joint_small,NumPwN_VS_small,given);
% [TabCondProb_VTA_VS_large,CondProb_VTA_VS_large] = CondProbAsType(count_Joint_large,NumPwN_VS_large,given);

% Chance level
% in case vs|vta I can compare th conditional probability {P(Vs|Vta)} with P(Vs), in
% case vta|vs we compare P(vta|vs) with P(vta): because: 
%%%%%%%%%%% Th: P(A|B)= P(A joint B)/P(B)  (Conditional=joint/2ndProb) => 
%%%%%%%%%%% => P(A joint B) = P(A|B)*P(B)
%%%%%%%%%%% When two events are indipendent we have P(A joint B) = P(A)*P(B) =>
%%%%%%%%%%% This statment is equivalent to say that in indipendent case P(A|B)=P(A)
%%%%%%%%%%% The chance level is the level that measure the probability that
%%%%%%%%%%% two elements form a couple only by chance, namely though they
%%%%%%%%%%% are indipendent, so to find our chance level is P(A)

given = ('vs|vta');
[ChanceLevPairs_VS_VTA_small_vsvta] = ChanceLevelCond(Tab_PwN_small_vsvta, given);
[ChanceLevPairs_VS_VTA_small_vtavs] = ChanceLevelCond(Tab_PwN_small_vtavs, given);
[ChanceLevPairs_VS_VTA_small_sync] = ChanceLevelCond(Tab_PwN_small_sync, given);
[ChanceLevPairs_VS_VTA_large_vsvta] = ChanceLevelCond(Tab_PwN_large_vsvta, given);
[ChanceLevPairs_VS_VTA_large_vtavs] = ChanceLevelCond(Tab_PwN_large_vtavs, given);
[ChanceLevPairs_VS_VTA_large_sync] = ChanceLevelCond(Tab_PwN_large_sync, given);

% given = ('vta|vs');
% [ChanceLevPairs_VTA_VS_small] = ChanceLevelCond(TLabels_PwN_small, given);
% [ChanceLevPairs_VTA_VS_large] = ChanceLevelCond(TLabels_PwN_large, given);



%%%%% Ratio Assemblies bidirections/ chance level
RatioPairs_VS_VTA_small_vsvta = CondProb_VS_VTA_small_vsvta./ChanceLevPairs_VS_VTA_small_vsvta;
RatioPairs_VS_VTA_small_vtavs = CondProb_VS_VTA_small_vtavs./ChanceLevPairs_VS_VTA_small_vtavs;
RatioPairs_VS_VTA_small_sync = CondProb_VS_VTA_small_vsvta./ChanceLevPairs_VS_VTA_small_sync;
RatioPairs_VS_VTA_large_vsvta = CondProb_VS_VTA_large_vsvta./ChanceLevPairs_VS_VTA_large_vsvta;
RatioPairs_VS_VTA_large_vtavs = CondProb_VS_VTA_large_vtavs./ChanceLevPairs_VS_VTA_large_vtavs;
RatioPairs_VS_VTA_large_sync = CondProb_VS_VTA_large_sync./ChanceLevPairs_VS_VTA_large_sync;

% RatioPairs_VTA_VS_small = CondProb_VTA_VS_small./ChanceLevPairs_VTA_VS_small;
% RatioPairs_VTA_VS_large = CondProb_VTA_VS_large./ChanceLevPairs_VTA_VS_large;
%%

%%%%%%%%%%%%%%%%% Plot Probability


VsVtaType{1} = strcat('M-T1'); VsVtaType{2} = strcat('M-T2'); VsVtaType{3} = strcat('M-T3'); VsVtaType{4} = strcat('M-No');
VsVtaType{5} = strcat('F-T1'); VsVtaType{6} = strcat('F-T2'); VsVtaType{7} = strcat('F-T3'); VsVtaType{8} = strcat('F-No');
VsVtaType{9} = strcat('C-T1'); VsVtaType{10} = strcat('C-T2'); VsVtaType{11} = strcat('C-T3'); VsVtaType{12} = strcat('C-No');
VsVtaType{13} = strcat('No-T1'); VsVtaType{14} = strcat('No-T2'); VsVtaType{15} = strcat('No-T3'); VsVtaType{16} = strcat('No-No');

ysup1=max([max(max(RatioPairs_VS_VTA_small_vsvta)),max(max(RatioPairs_VS_VTA_small_vtavs)),max(max(RatioPairs_VS_VTA_small_sync)),...
    max(max(RatioPairs_VS_VTA_large_vsvta)),max(max(RatioPairs_VS_VTA_large_vtavs)),max(max(RatioPairs_VS_VTA_large_sync))]);
yinf =min([min(min(RatioPairs_VS_VTA_small_vsvta)),min(min(RatioPairs_VS_VTA_small_vtavs)),min(min(RatioPairs_VS_VTA_small_sync)),...
    min(min(RatioPairs_VS_VTA_large_vsvta)),min(min(RatioPairs_VS_VTA_large_vtavs)),min(min(RatioPairs_VS_VTA_large_sync))]);
maxRatio_small_vsvta = max(RatioPairs_VS_VTA_small_vsvta);
maxRatio_small_vtavs = max(RatioPairs_VS_VTA_small_vtavs);
maxRatio_small_sync = max(RatioPairs_VS_VTA_small_sync);

maxRatio_large_vsvta = max(RatioPairs_VS_VTA_large_vsvta);
maxRatio_large_vtavs = max(RatioPairs_VS_VTA_large_vtavs);
maxRatio_large_sync = max(RatioPairs_VS_VTA_large_sync);

figure(12);hold on;

subplot(3,2,1)
for i=1:size(ChanceLevPairs_VS_VTA_small_vsvta,2)
    if ~isnan(maxRatio_small_vsvta(i))
    text(i-0.2,maxRatio_small_vsvta(i)+0.4,VsVtaType{i}, 'Fontsize',8);hold on;
    else
    text(i-0.2,ysup1,VsVtaType{i}, 'Fontsize',8);hold on;
    end
end
plot([0 size(ChanceLevPairs_VS_VTA_small_vsvta,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((RatioPairs_VS_VTA_small_vsvta)','b*');hold on;
boxplot(RatioPairs_VS_VTA_small_vsvta);hold on;
ylim([(yinf-0.5) (ysup1+0.8)]);hold on;
title('Assembly type Ratio vs->vta Small Bins');hold on; 


subplot(3,2,3)
for i=1:size(ChanceLevPairs_VS_VTA_small_vtavs,2)
    if ~isnan(maxRatio_small_vtavs(i))
    text(i-0.2,maxRatio_small_vtavs(i)+0.4,VsVtaType{i}, 'Fontsize',8);hold on;
    else
    text(i-0.2,ysup1,VsVtaType{i}, 'Fontsize',8);hold on;
    end
end
plot([0 size(ChanceLevPairs_VS_VTA_small_vtavs,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((RatioPairs_VS_VTA_small_vtavs)','b*');hold on;
boxplot(RatioPairs_VS_VTA_small_vtavs);hold on;
ylim([(yinf-0.5) (ysup1+0.8)]); hold on;
title('Assembly type Ratio vta->vs Small Bins'); hold on;
text(-1,0,Text,'Rotation',90,'Fontsize',14);hold on;

subplot(3,2,5)
for i=1:size(ChanceLevPairs_VS_VTA_small_sync,2)
    if ~isnan(maxRatio_small_sync(i))
    text(i-0.2,maxRatio_small_sync(i)+0.4,VsVtaType{i}, 'Fontsize',8);hold on;
    else
    text(i-0.2,ysup1,VsVtaType{i}, 'Fontsize',8);hold on;
    end
end
plot([0 size(ChanceLevPairs_VS_VTA_small_sync,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((RatioPairs_VS_VTA_small_sync)','b*');hold on;
boxplot(RatioPairs_VS_VTA_small_sync);hold on;
ylim([(yinf-0.5) (ysup1+0.8)]); hold on;
title('Assembly type Ratio sync Small Bins');hold on;


subplot(3,2,2)
for i=1:size(ChanceLevPairs_VS_VTA_large_vsvta,2)
    if ~isnan(maxRatio_large_vsvta(i))
    text(i-0.2,maxRatio_large_vsvta(i)+0.4,VsVtaType{i},'Fontsize',8);hold on;
    else
    text(i-0.2,ysup1,VsVtaType{i}, 'Fontsize',8);hold on;
    end
end
plot([0 size(ChanceLevPairs_VS_VTA_large_vsvta,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((RatioPairs_VS_VTA_large_vsvta)','b*');hold on;
boxplot(RatioPairs_VS_VTA_large_vsvta);hold on;
ylim([(yinf-0.5) (ysup1+0.8)]);hold on;
title('Assembly type Ratio vs->vta Large Bins');hold on;


subplot(3,2,4)
for i=1:size(ChanceLevPairs_VS_VTA_large_vtavs,2)
    if ~isnan(maxRatio_large_vtavs(i))
    text(i-0.2,maxRatio_large_vtavs(i)+0.4,VsVtaType{i},'Fontsize',8);hold on;
    else
    text(i-0.2,ysup1,VsVtaType{i}, 'Fontsize',8);hold on;
    end
end
plot([0 size(ChanceLevPairs_VS_VTA_large_vtavs,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((RatioPairs_VS_VTA_large_vtavs)','b*');hold on;
boxplot(RatioPairs_VS_VTA_large_vtavs);hold on;
ylim([(yinf-0.5) (ysup1+0.8)]); hold on;
title('Assembly type Ratio vta->vs Large Bins');hold on;



subplot(3,2,6)
for i=1:size(ChanceLevPairs_VS_VTA_large_sync,2)
 if ~isnan(maxRatio_large_sync(i))
    text(i-0.2,maxRatio_large_sync(i)+0.4,VsVtaType{i},'Fontsize',8);hold on;
 else
  text(i-0.2,ysup1,VsVtaType{i}, 'Fontsize',8);hold on;
 end
end
plot([0 size(ChanceLevPairs_VS_VTA_large_sync,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((RatioPairs_VS_VTA_large_sync)','b*');hold on;
boxplot(RatioPairs_VS_VTA_large_sync);hold on;
ylim([(yinf-0.5) (ysup1+0.8)]); hold on;
title('Assembly type Ratio sync Large Bins');hold on;

box on;


SaveName = strcat('AsNeutypeAndProbWithDirSLAll',nameT);

save(SaveName,'-v7.3')


% % % count=Perc_NinP_small_An_Tog; %neurons in pairs animal together
% % % MSN = count(:,1); FSI = count(:,2); CIN = count(:,3); NoID_Vs = count(:,4);  type1 = count(:,5); type2 =count(:,6);
% % % type3 = count(:,7); NoID_Vta = count(:,8); 
% % % TLabelsPercNeuAnTog_small = table(MSN, FSI, CIN, NoID_Vs, type1, type2, type3, NoID_Vta);
% % % clear count MSN FSI CNI NoID_Vs type1 type2 type3 NoID_Vta
% % % 
% % % count=Perc_NinP_large_AnTog; %neurons in pairs animal together
% % % MSN = count(:,1); FSI = count(:,2); CIN = count(:,3); NoID_Vs = count(:,4);  type1 = count(:,5); type2 =count(:,6);
% % % type3 = count(:,7); NoID_Vta = count(:,8); 
% % % TLabelsPercNeuAnTog_large = table(MSN, FSI, CIN, NoID_Vs, type1, type2, type3, NoID_Vta);
% % % clear count MSN FSI CNI NoID_Vs type1 type2 type3 NoID_Vta
% % % 
% % % 
% % % count=Perc_NinP_small_AnDist; %neurons in Pairs/animal distinct
% % % MSN = count(:,1); FSI = count(:,2); CIN = count(:,3); NoID_Vs = count(:,4);  type1 = count(:,5); type2 =count(:,6);
% % % type3 = count(:,7); NoID_Vta = count(:,8); 
% % % TLabelsPercNeuAnDist_small = table(MSN, FSI, CIN, NoID_Vs, type1, type2, type3, NoID_Vta);
% % % clear count MSN FSI CNI NoID_Vs type1 type2 type3 NoID_Vta
% % % 
% % % count=Perc_NinP_large_AnDist; %neurons in Pairs/animal distinct
% % % MSN = count(:,1); FSI = count(:,2); CIN = count(:,3); NoID_Vs = count(:,4);  type1 = count(:,5); type2 =count(:,6);
% % % type3 = count(:,7); NoID_Vta = count(:,8); 
% % % TLabelsPercNeuAnDist_large = table(MSN, FSI, CIN, NoID_Vs, type1, type2, type3, NoID_Vta);
% % % clear count MSN FSI CNI NoID_Vs type1 type2 type3 NoID_Vta
