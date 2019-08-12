%% Script to count whether ensembles do prefer specific cathegories of units
addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM');
% % % % % % % % The matrix new_data is the matrix that contains the info, in particular new_data.spikeT_BegEnd has the 
% % % % % % % % spiketrain cutted properly

% % name1=strcat('Lastrev1_OrRevW.mat');
% % load(name1)
% % alg= strcat('Stand_');
% % % alg=strcat('Pru_');
% % name=strcat('An_Disp_',alg,name1);
% % % name=strcat('An_Disp_',name);
% % load(name);
% % name_bins=strcat('BinLag_',alg,name1);
% % load(name_bins)
% % % MatBin=-MBinAnP;
% % nameSmallBig=strcat('SmallBig_An_Disp',alg,name1)
% % load(nameSmallBig)


% fignum=27;
% fignum1=28;
% 
% clear c
% load('classlist.mat');
% 
% 
% Region='vsvta';
% % % % % Region='vsvs';
% % % % % Region='vtavta';
% % % % 
% % % % %Dir stays for "Directionality", can assume these values
% % % % %'vs->vta Dir','vta->vs Dir','vs->vta Inv', 'vta->vs Inv'0
% % % % %  Dir='vs->vta Dir';
% % % % % Dir='vs->vta Inv';
% % % % % Dir='vs->vta AllBins';
% % % % 
% % % % %  Dir='vta->vs Dir';
% % % % %  Dir='vta->vs Inv';
% % % % %  Dir='vta->vs AllBins';
% % % % %    Dir=0;
% % % % % Dir='sync';
% % % % % TitleP=strcat(names, 'Dir',Dir);

% % % 
% % % nnr=nn(a);
% % % [new_data] = parId(new_data,a,nnr);
% % % 
% % % 
% % % for k=1:length(new_data.par)
% % %  [new_data] = labels(new_data,c,k);
% % %  end
% % %  clear tt i events licks laser info spikes map params spM spikes spike_region j k new_spM
% % % 
% % % %%%%%%%%%%%%%%%% Find number of units of the Striatum and VTA nneuS
% % % 
% % % new_spM=new_data;
% % % 
% % % % TotAn=size(new_spM.spike_regionNoId,2);
% % % TotAn=length(new_data.par);
% % % for k=1:TotAn
% % %     nneuVS{k}=0;
% % %     nneuVTA{k}=0;
% % % end
% % % 
% % % for k=1:TotAn
% % %     for j=1:size(new_spM.spike_regionNoId{k},2)
% % %         [nneuVS{k}] = countNeu(new_spM.spike_regionNoId{k}(1,j),1,nneuVS{k});
% % %         [nneuVTA{k}] = countNeu(new_spM.spike_regionNoId{k}(1,j),2,nneuVTA{k});
% % %     end
% % % 
% % % end
% % % %%%%%
% % % 
% % % 
% % % %%%%% Find Pairs in one region and shared between two differnt regions
% % % % pairs diveded by regions with all info
% % % [pairs_r_small,pairs_vsvs_small,pairs_vsvta_small,pairs_vtavta_small] = PairsInfo(As_across_smallBins,As_order_smallBins,nneuVS);
% % % [pairs_r_large,pairs_vsvs_large,pairs_vsvta_large,pairs_vtavta_large] = PairsInfo(As_across_largeBins,As_order_largeBins,nneuVS);
% % % 
% % % clear a b A B i ii j k
% % % 
% % % [struct_pair_small] = PairsOrgByPair(pairs_vsvta_small);
% % % [struct_pair_large] = PairsOrgByPair(pairs_vsvta_large);
% % % 
% % % 
% % % [struct_pair_small] = PairNeuLabelsID(struct_pair_small,new_data);
% % % [struct_pair_large] = PairNeuLabelsID(struct_pair_large,new_data);
% % % copy_struct_pair_small=struct_pair_small;
% % % copy_struct_pair_large=struct_pair_large;
% % % clear struct_pair
%%%%%% Create the struct coerently with the directionality
%%%%%%%%%% Here the directionality doesn't cointan the info about the bins
%%%%%%%%%% because the pairs are already divided for bins

% % %  Dir1='vs->vta';
% % %  Dir2='vta->vs';
% % %  Dir3='sync';
% % %  Dir='no';
% Dir='vs->vta Inv';
% Dir='vs->vta AllBins';

%  
%  Dir='vta->vs Inv';
%  Dir='vta->vs AllBins';
%    Dir=0;
% Dir='sync';

% % % path_po='/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM';
% % % cd(path_po)
% % % 
% % % for k=1:TotAn
% % %     if ~isempty(copy_struct_pair_small)
% % % [struct_pair_small_vsvta{k}] = DirAndBinPostPru(Dir1, copy_struct_pair_small,k);
% % %  [struct_pair_small_vtavs{k}] = DirAndBinPostPru(Dir2, copy_struct_pair_small,k);
% % %  [struct_pair_small_sync{k}] = DirAndBinPostPru(Dir3, copy_struct_pair_small,k);
% % %     end
% % %     if ~isempty(copy_struct_pair_large)
% % %     [struct_pair_large_vsvta{k}] = DirAndBinPostPru(Dir1, copy_struct_pair_large,k);
% % %     [struct_pair_large_vtavs{k}] = DirAndBinPostPru(Dir2, copy_struct_pair_large,k);
% % %     [struct_pair_large_sync{k}] = DirAndBinPostPru(Dir3, copy_struct_pair_large,k);
% % %     end
% % % end
% % % 
% % % % TotAn=length(new_data.par);
% count=zeros(TotAn,7);
% for k=1:TotAn
% [count] = CountTypeNeuPair(struct_pair{k},k,count);

path_p= '/zifnas/Carla/CWEM_Project_GenProg';
cd(path_p)
MergeDataGenPostPru_1;
MergeDataGenPostPru_2;

path_po='/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM';
cd(path_po)
copy_struct_pair_small=struct_pair_small;
copy_struct_pair_large=struct_pair_large;
% clear struct_pair
Dir1='vs->vta';
 Dir2='vta->vs';
 Dir3='sync';
 Dir='no';
 TotAnPerS = TotAn;
clear TotAn;
TotAn = SumTotAn;
for k=1:TotAn
    if ~isempty(copy_struct_pair_small)
% [struct_pair_small_vsvta{k}] = DirAndBinPostPru(Dir1, copy_struct_pair_small,k);
%  [struct_pair_small_vtavs{k}] = DirAndBinPostPru(Dir2, copy_struct_pair_small,k);
%  [struct_pair_small_sync{k}] = DirAndBinPostPru(Dir3, copy_struct_pair_small,k);
 [struct_pair_small_no{k}] = DirAndBinPostPru(Dir, copy_struct_pair_small,k);
    end
    if ~isempty(copy_struct_pair_large)
%     [struct_pair_large_vsvta{k}] = DirAndBinPostPru(Dir1, copy_struct_pair_large,k);
%     [struct_pair_large_vtavs{k}] = DirAndBinPostPru(Dir2, copy_struct_pair_large,k);
% [struct_pair_large_sync{k}] = DirAndBinPostPru(Dir3, copy_struct_pair_large,k);
    [struct_pair_large_no{k}] = DirAndBinPostPru(Dir, copy_struct_pair_large,k);
    end
end

%% Do the assemblies prefer specific units?
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
[~,count_PwN_small_vsvta] = CountTypeNeuPair(struct_pair_small_vsvta{k},k,count_PwN_small_vsvta); % pairs with specific units, animal per animal
 end
 
 count_PwN_small_vtavs=zeros(TotAn,8);
 for k=1:TotAn
[~,count_PwN_small_vtavs] = CountTypeNeuPair(struct_pair_small_vtavs{k},k,count_PwN_small_vtavs); % pairs with specific units, animal per animal
 end
 
 count_PwN_small_sync=zeros(TotAn,8);
 for k=1:TotAn
[~,count_PwN_small_sync] = CountTypeNeuPair(struct_pair_small_sync{k},k,count_PwN_small_sync); % pairs with specific units, animal per animal
 end
 
 %%% Large bins different directionality
count_PwN_large_vsvta=zeros(TotAn,8);
 for k=1:TotAn
[~,count_PwN_large_vsvta] = CountTypeNeuPair(struct_pair_large_vsvta{k},k,count_PwN_large_vsvta); % pairs with specific units, animal per animal
 end
 
 count_PwN_large_vtavs=zeros(TotAn,8);
 for k=1:TotAn
[~,count_PwN_large_vtavs] = CountTypeNeuPair(struct_pair_large_vtavs{k},k,count_PwN_large_vtavs); % pairs with specific units, animal per animal
 end
 
 count_PwN_large_sync=zeros(TotAn,8);
 for k=1:TotAn
[~,count_PwN_large_sync] = CountTypeNeuPair(struct_pair_large_sync{k},k,count_PwN_large_sync); % pairs with specific units, animal per animal
 end
 

%%%%%%%%%%%%%%%%% Distict Animals
countMat_Pairs_small_vsvta=repmat(count_Pairs_small_vsvta,8,1);
countMat_Pairs_small_vtavs=repmat(count_Pairs_small_vtavs,8,1);
countMat_Pairs_small_sync=repmat(count_Pairs_small_sync,8,1);

countMat_Pairs_large_vsvta=repmat(count_Pairs_large_vsvta,8,1);
countMat_Pairs_large_vtavs=repmat(count_Pairs_large_vtavs,8,1);
countMat_Pairs_large_sync=repmat(count_Pairs_large_sync,8,1);


Perc_PwN_pAn_small_vsvta=count_PwN_small_vsvta./countMat_Pairs_small_vsvta'; % # assemblies with V_{ri}/ # assemblies Animal per animal
Perc_PwN_pAn_small_vtavs=count_PwN_small_vtavs./countMat_Pairs_small_vtavs';
Perc_PwN_pAn_small_sync=count_PwN_small_sync./countMat_Pairs_small_sync';

Perc_PwN_pAn_large_vsvta=count_PwN_large_vsvta./countMat_Pairs_large_vsvta'; % # assemblies with V_{ri}/ # assemblies Animal per animal
Perc_PwN_pAn_large_vtavs=count_PwN_large_vtavs./countMat_Pairs_large_vtavs';
Perc_PwN_pAn_large_sync=count_PwN_large_sync./countMat_Pairs_large_sync';

% TotAn=length(new_data.par);
count_Neu=zeros(TotAn,8);
for k=1:TotAn
[~,count_Neu] = CountTypeNeuAll(Data_All.par{k},k,count_Neu);
end

ChLevNeu_pAn= NaN(TotAn,8);
copy_nneuVS=nneuVS_All;
copy_nneuVTA=nneuVTA_All;
for k=1:TotAn
    if copy_nneuVS{k}==0, copy_nneuVS{k}=10^8; end
    if copy_nneuVTA{k}==0, copy_nneuVTA{k}=10^8; end
    ChLevNeu_pAn(k,1:4) = count_Neu(k,1:4)./copy_nneuVS{k};
    ChLevNeu_pAn(k,5:8) = count_Neu(k,5:8)./copy_nneuVTA{k};
end

% Perc_PwN_pAn=Perc_PwN_pAn';
% ChLevNeu_pAn=ChLevNeu_pAn';
%%
% ChLevNeu_pAn(Perc_PwN_pAn==nan)=[];
% Perc_PwN_pAn(Perc_PwN_pAn==nan)=[];
% Text=strcat('Post Pruned---',name1);
type{1}=strcat('MSN');type{2}=strcat('FSI');type{3}=strcat('CIN');type{4}= strcat('NoId-Vs');type{5}= strcat('type1');type{6}=strcat('type2');
type{7}=strcat('type3');type{8}=strcat('NoId-Vta');

ysup=max([max(max(Perc_PwN_pAn_small_vsvta./ChLevNeu_pAn)),max(max(Perc_PwN_pAn_small_vtavs./ChLevNeu_pAn)),max(max(Perc_PwN_pAn_small_sync./ChLevNeu_pAn)),...
    max(max(Perc_PwN_pAn_large_vsvta./ChLevNeu_pAn)),max(max(Perc_PwN_pAn_large_vtavs./ChLevNeu_pAn)),max(max(Perc_PwN_pAn_large_sync./ChLevNeu_pAn))]);
yinf =min([min(min(Perc_PwN_pAn_small_vsvta./ChLevNeu_pAn)),min(min(Perc_PwN_pAn_small_vtavs./ChLevNeu_pAn)),min(min(Perc_PwN_pAn_small_sync./ChLevNeu_pAn)),...
    min(min(Perc_PwN_pAn_large_vsvta./ChLevNeu_pAn)),min(min(Perc_PwN_pAn_large_vtavs./ChLevNeu_pAn)),min(min(Perc_PwN_pAn_large_sync./ChLevNeu_pAn))]);
figure(6);hold on;
% boxplot(Perc_PwN_pAn)
subplot(3,2,1)
for i=1:size(ChLevNeu_pAn,2)
    text(i-0.2,ysup,type{i});hold on;
end
plot([0 size(ChLevNeu_pAn,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((Perc_PwN_pAn_small_vsvta./ChLevNeu_pAn)','b*');hold on;
boxplot(Perc_PwN_pAn_small_vsvta./ChLevNeu_pAn);hold on;
ylim([(yinf-0.5) (ysup+0.5)]);
title('Small Bins vs->vta');

subplot(3,2,2)

for i=1:size(ChLevNeu_pAn,2)
    text(i-0.2,ysup,type{i});hold on;
end
plot([0 size(ChLevNeu_pAn,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((Perc_PwN_pAn_large_vsvta./ChLevNeu_pAn)','b*');hold on;
boxplot(Perc_PwN_pAn_large_vsvta./ChLevNeu_pAn);hold on;
ylim([(yinf-0.5) (ysup+0.5)]);
title('Large Bins vs->vta')
% ylim([(yinf-0.5) (ysup+0.5)]);

subplot(3,2,3)
for i=1:size(ChLevNeu_pAn,2)
    text(i-0.2,ysup,type{i});hold on;
end
plot([0 size(ChLevNeu_pAn,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((Perc_PwN_pAn_small_vtavs./ChLevNeu_pAn)','b*');hold on;
boxplot(Perc_PwN_pAn_small_vtavs./ChLevNeu_pAn);hold on;
ylim([(yinf-0.5) (ysup+0.5)]);
% text(-0.5,0,Text,'Rotation',90,'Fontsize',14);hold on;
title('Small Bins vta->vs')
% ylim([(yinf-0.5) (ysup+0.5)]);

subplot(3,2,4)
for i=1:size(ChLevNeu_pAn,2)
    text(i-0.2,ysup,type{i});hold on;
end
plot([0 size(ChLevNeu_pAn,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((Perc_PwN_pAn_large_vtavs./ChLevNeu_pAn)','b*');hold on;
boxplot(Perc_PwN_pAn_large_vtavs./ChLevNeu_pAn);hold on;
ylim([(yinf-0.5) (ysup+0.5)]);
title('Large Bins vta->vs')


subplot(3,2,5)
for i=1:size(ChLevNeu_pAn,2)
    text(i-0.2,ysup,type{i});hold on;
end
plot([0 size(ChLevNeu_pAn,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((Perc_PwN_pAn_small_sync./ChLevNeu_pAn)','b*');hold on;
boxplot(Perc_PwN_pAn_small_sync./ChLevNeu_pAn);hold on;
ylim([(yinf-0.5) (ysup+0.5)]);
title('Small Bins sync')


subplot(3,2,6)
for i=1:size(ChLevNeu_pAn,2)
    text(i-0.2,ysup,type{i});hold on;
end
plot([0 size(ChLevNeu_pAn,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((Perc_PwN_pAn_large_sync./ChLevNeu_pAn)','b*');hold on;
boxplot(Perc_PwN_pAn_large_sync./ChLevNeu_pAn);hold on;
ylim([(yinf-0.5) (ysup+0.5)]);
title('Large Bins sync')

% % % for i=1:size(ChLevNeu_pAn,2)
% % %     text(i-0.2,3+0.05,type{i});hold on;
% % % end
% % % plot([0 size(ChLevNeu_pAn,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
box on;

 
 
% end
%% Are there typologies of neurons most probably deputed than others to form assemblies?

% functions to count how many neurons of a cathegories are part of one
% assembly (the neurons are taken only once), the first funtion eliminates
% double elements, the second one counts
struct_pair=struct_pair_small
for k=1:TotAn
    Un{k}=struct('neuID',[],'ia',[],'ic',[],'labels',[]);
    [Un{k}.neuID,Un{k}(:,1).ia(:,1),Un{k}.ic(:,1)]=unique(struct_pair{k}.neuID);
%     [Un{k}.neuID,Un{k}.ia,Un{k}.ic]=unique(struct_pair{k}.neuID(:,1));
       Un{k}.labels(:,1) = struct_pair{k}.labels(Un{k}.ia); % I create here the structure with the index unique
end


% TotAn=length(Data_All.par);
count_NinP=zeros(TotAn,8);
for k=1:TotAn
    if~isempty(Un{k}.neuID)
        for i=1:length(Un{k}.neuID)
[~,count_NinP] = CountTypeNeuGen(Un{k}.labels(i),k,count_NinP); 
        end
    end
end

% Function that count how many neurons of specific cathegories are present

% TotAn=length(new_data.par);
count_Neu=zeros(TotAn,8);
for k=1:TotAn
[~,count_Neu] = CountTypeNeuAll(Data_All.par{k},k,count_Neu);
end


%%%%% Maketables
count=count_NinP;

MSN = count(:,1); FSI = count(:,2); CIN = count(:,3); NoID_Vs = count(:,4); type1 = count(:,5); type2 =count(:,6);
type3 = count(:,7); NoID_Vta = count(:,8); 
TLabels_NinP = table(MSN, FSI, CIN, NoID_Vs, type1, type2, type3, NoID_Vta);


clear count
count=count_Neu;
MSN = count(:,1); FSI = count(:,2); CIN = count(:,3); NoID_Vs = count(:,4);  type1 = count(:,5); type2 =count(:,6);
type3 = count(:,7); NoID_Vta = count(:,8); 
TLabels_Neu = table(MSN, FSI, CIN, NoID_Vs, type1, type2, type3, NoID_Vta);

%%%%%% Make sum and ratio for each categories of neurons
%%%%%%% Goal:  #neu_c in any assembly/ #neu_c Tot 
%%%%%%% To see whether there exist typologies of neurons more deputed to
%%%%%%% take part in any assembly

Sum_NinP = sum(TLabels_NinP{:,:});
Sum_Neu = sum(TLabels_Neu{:,:});
% for i=1: size(count_Neu,2)
% Sum_NinP(i) = sum(TLabels_NinP{find(count_Neu(:,i)>0),i});
% Sum_Neu(i) = sum(TLabels_Neu{find(count_Neu(:,i)>0),i});
% 
% end


% PercNeuPP1=SumPairsUn1./SumNeu1;
Perc_NinP_AnTog=Sum_NinP./Sum_Neu; % Here I'm considering everything omiting nan results Take the sum and no distiction for animal
Perc_NinP_AnDist=count_NinP./count_Neu;  % divided for animal, then I can do the mean
% PR=sum(PercNeu1,1,'omitnan'); % Here is equivalent to assign 0 propability for a neuron that doesn't exist
% PRM=PR/TotAn;
% for i=1: size(countNeu,2)
%   PRMM(i) = mean(PercNeu1(find(countNeu(:,i)>0),i)); % when there is find>0 is equivalent to mean(PercNeu1,'omitnan')
% end


%%%%%%%%%%%%
clear count
count=Perc_NinP_AnTog; %neurons in pairs animal together
MSN = count(:,1); FSI = count(:,2); CIN = count(:,3); NoID_Vs = count(:,4);  type1 = count(:,5); type2 =count(:,6);
type3 = count(:,7); NoID_Vta = count(:,8); 
TLabelsPercNeuAnTog = table(MSN, FSI, CIN, NoID_Vs, type1, type2, type3, NoID_Vta);

clear count
count=Perc_NinP_AnDist; %neurons in Pairs/animal distinct
MSN = count(:,1); FSI = count(:,2); CIN = count(:,3); NoID_Vs = count(:,4);  type1 = count(:,5); type2 =count(:,6);
type3 = count(:,7); NoID_Vta = count(:,8); 
TLabelsPercNeuAnDist = table(MSN, FSI, CIN, NoID_Vs, type1, type2, type3, NoID_Vta);
%% %%%%%%
type{1}=strcat('MSN');type{2}=strcat('FSI');type{3}=strcat('CIN');type{4}= strcat('NoId-Vs');type{5}= strcat('type1');type{6}=strcat('type2');
type{7}=strcat('type3');type{8}=strcat('NoId-Vta');
figure(8);hold on;
% boxplot(Perc_PwN_pAn)
% boxplot(Perc_PwN_pAn./ChLevNeu_pAn);hold on;
for i=1:size(Perc_NinP_AnTog,2)
%     plot([i-0.5,i+0.5],[ChLevNeu(i), ChLevNeu(i)],'k.-','linewidth',1.5);hold on
    text(i-0.2,Perc_NinP_AnTog(i)+0.1,type{i}),
end
ylim([0,max(max([mean(Perc_NinP_AnDist,1,'omitnan');Perc_NinP_AnTog]))+0.2])
% title(TitleP);
conf=[mean(Perc_NinP_AnDist,1,'omitnan');Perc_NinP_AnTog]
  bar(conf','group');hold on;
% bar(mean(Perc_NinP_AnDist,1,'omitnan'),'r');
box on;

