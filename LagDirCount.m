addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM');
% The matrix new_data is the matrix that contains the info, in particular new_data.spikeT_BegEnd has the 
% spiketrain cutted properly in databasePrepro.m


% load('StableOr_rev4.mat')
% load('AsStableOr_rev4.mat')
% load('DispAnStableOrRev4.mat')
name1=strcat('Lastrev1_OrRevW.mat');
load(name1)
names=strcat('An_Disp_Pru_',name1);
load(names)

fignum=25;
fignum1=26;
% load('Lastrev1_OrRevW.mat')
% load('AsStandOrW_lastrev1.mat')
% load('AsPruOrW_lastrev1.mat')
% % load('DispLastrev1OrW.mat')
% load('DispLastrev1OrWPru.mat')


%%%%% Originalreversal Analysis
% load('Lastrev1_OrRevW.mat')
% load('AsPruOrRevW_lastrev1_rl2.mat')
% load('DispLastrev1_Rl2_OrRevWPru.mat');

% load('Rev4_OrRevW.mat')
% load('AsPruOrRevW_rev4_rl2.mat')
% load('DispRev4_Rl2_OrRevWPru.mat');



clear c
load('classlist.mat');


Region='vsvta';
% Region='vsvs';
% Region='vtavta';

%Dir stays for "Directionality", can assume these values
%'vs->vta Dir','vta->vs Dir','vs->vta Inv', 'vta->vs Inv'0
%  Dir='vs->vta Dir';
% Dir='vs->vta Inv';
% Dir='vs->vta AllBins';

%  Dir='vta->vs Dir';
%  Dir='vta->vs Inv';
%  Dir='vta->vs AllBins';
%    Dir=0;
% Dir='sync';
TitleP=strcat(names, 'Dir',Dir);


nnr=nn(a);
[new_data] = parId(new_data,a,nnr);


for k=1:length(new_data.par)
 [new_data] = labels(new_data,c,k);
 end
 clear tt i events licks laser info spikes map params spM spikes spike_region j k

 
 

%%%%%%%%%%%%%%%% Find number of units of the Striatum and VTA nneuS

new_spM=new_data;
As_across_bins=As_acr_bins_pru;
As_order=As_order_pru;
TotAn=size(new_spM.spike_regionNoId,2);
for k=1:TotAn
    nneuVS{k}=0;
    nneuVTA{k}=0;
end

for k=1:TotAn
    for j=1:size(new_spM.spike_regionNoId{k},2)
        [nneuVS{k}] = countNeu(new_spM.spike_regionNoId{k}(1,j),1,nneuVS{k});
        [nneuVTA{k}] = countNeu(new_spM.spike_regionNoId{k}(1,j),2,nneuVTA{k});
    end

end


%%%%% Find Pairs in one region and shared between two differnt regions
% pairs diveded by regions with all info
[pairs_r,pairs_vsvs,pairs_vsvta,pairs_vtavta] = PairsInfo(As_across_bins,As_order,nneuVS);

clear a b A B i ii j k
%%%% # of possible pairs in the single regions and intra-region
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

[struct_pair] = PairsOrgByPair(pairs_vsvta);

[struct_pair] = PairNeuLabelsID(struct_pair,new_data);
copy_struct_pair=struct_pair;
clear struct_pair
%%%%%% Create the struct coerently with the directionality
for k=1:length(copy_struct_pair)
    if ~isempty(copy_struct_pair)
[struct_pair{k}] = DirAndBin(Dir, copy_struct_pair,k);
    end
end

% TotAn=length(new_data.par);
% count=zeros(TotAn,7);
% for k=1:TotAn
% [count] = CountTypeNeuPair(struct_pair{k},k,count);
% end
%% Are there typologies of neurons most probably deputed than others to form assemblies?

% functions to count how many neurons of a cathegories are part of one
% assembly (the neurons are taken only once), the first funtion eliminates
% double elements, the second one counts

for k=1:TotAn
    Un{k}=struct('neuID',[],'ia',[],'ic',[],'labels',[]);
    [Un{k}.neuID,Un{k}(:,1).ia(:,1),Un{k}.ic(:,1)]=unique(struct_pair{k}.neuID);
%     [Un{k}.neuID,Un{k}.ia,Un{k}.ic]=unique(struct_pair{k}.neuID(:,1));
       Un{k}.labels(:,1) = struct_pair{k}.labels(Un{k}.ia); % I create here the structure with the index unique
end


TotAn=length(new_data.par);
count_NinP=zeros(TotAn,8);
for k=1:TotAn
    if~isempty(Un{k}.neuID)
        for i=1:length(Un{k}.neuID)
[count_NinP] = CountTypeNeuGen(Un{k}.labels(i),k,count_NinP); 
        end
    end
end

% Function that count how many neurons of specific cathegories are present

TotAn=length(new_data.par);
count_Neu=zeros(TotAn,8);
for k=1:TotAn
[count_Neu] = CountTypeNeuAll(new_data.par{k},k,count_Neu);
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
figure(fignum);hold on;
% boxplot(Perc_PwN_pAn)
% boxplot(Perc_PwN_pAn./ChLevNeu_pAn);hold on;
for i=1:size(Perc_NinP_AnTog,2)
%     plot([i-0.5,i+0.5],[ChLevNeu(i), ChLevNeu(i)],'k.-','linewidth',1.5);hold on
    text(i-0.2,Perc_NinP_AnTog(i)+0.1,type{i}),
end
ylim([0,max(max([mean(Perc_NinP_AnDist,1,'omitnan');Perc_NinP_AnTog]))+0.2])
title(TitleP);
conf=[mean(Perc_NinP_AnDist,1,'omitnan');Perc_NinP_AnTog]
  bar(conf','group');hold on;
% bar(mean(Perc_NinP_AnDist,1,'omitnan'),'r');
box on;

%% Do the assemblies prefer specific units?
% Total Number of assemblies

for k=1:TotAn
   % if ~isempty(struct_pair{k}.pair)
        count_Pairs(k)=length(struct_pair{k}.pair);
   % end
end

TotAs=sum(count_Pairs);
%%%%%%%%%%%%%%
% Number of assembly with a specific class of neurons %%%%%%% !!!!! 
TotAn=length(new_data.par);
count_PwN=zeros(TotAn,8);
 for k=1:TotAn
[count_PwN] = CountTypeNeuPair(struct_pair{k},k,count_PwN); % pairs with specific units, animal per animal
 end
%%%%%%%%%%%%%%%%% Distict Animals
count_PairsM=repmat(count_Pairs,8,1);
Perc_PwN_pAn=count_PwN./count_PairsM'; % # assemblies with V_{ri}/ # assemblies Animal per animal
ChLevNeu_pAn= NaN(TotAn,8);
copy_nneuVS=nneuVS;
copy_nneuVTA=nneuVTA;
for k=1:TotAn
    if copy_nneuVS{k}==0, copy_nneuVS{k}=10^8; end
    if copy_nneuVTA{k}==0, copy_nneuVTA{k}=10^8; end
    ChLevNeu_pAn(k,1:4) = count_Neu(k,1:4)./copy_nneuVS{k};
    ChLevNeu_pAn(k,5:8) = count_Neu(k,5:8)./copy_nneuVTA{k};
end

% Perc_PwN_pAn=Perc_PwN_pAn';
% ChLevNeu_pAn=ChLevNeu_pAn';

% ChLevNeu_pAn(Perc_PwN_pAn==nan)=[];
% Perc_PwN_pAn(Perc_PwN_pAn==nan)=[];
type{1}=strcat('MSN');type{2}=strcat('FSI');type{3}=strcat('CIN');type{4}= strcat('NoId-Vs');type{5}= strcat('type1');type{6}=strcat('type2');
type{7}=strcat('type3');type{8}=strcat('NoId-Vta');
figure(fignum1);hold on;
% boxplot(Perc_PwN_pAn)
boxplot(Perc_PwN_pAn./ChLevNeu_pAn);hold on;
for i=1:size(ChLevNeu_pAn,2)
%     plot([i-0.5,i+0.5],[ChLevNeu(i), ChLevNeu(i)],'k.-','linewidth',1.5);hold on
    text(i-0.2,3+0.05,type{i});hold on;
end
plot([0 size(ChLevNeu_pAn,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
title(TitleP);
% ylim([0,max(Perc_PwN)+0.08]);
box on;

 
 
 
%% %%%%%%% All Animals Together
As_PwN = sum(count_PwN, 1);
Perc_PwN=As_PwN./TotAs;

%%%% Chance level pairs with neurons
ChLevNeu=NaN(1,8);
NeuMatVS=cellfun(@sum,nneuVS);
NeuMatVTA=cellfun(@sum,nneuVTA);
NeuMat=cellfun(@sum,Nneu);
SumAllNeuVS=sum(NeuMatVS);
SumAllNeuVTA=sum(NeuMatVTA);
SumAllNeu=sum(NeuMat);
ChLevNeu(1,1:4)=Sum_Neu(1,1:4)/SumAllNeuVS;
ChLevNeu(1,5:8)=Sum_Neu(1,5:8)/SumAllNeuVTA;
% ChLevNeu(1,7)=Sum_Neu(1,7)/SumAllNeu;



type{1}=strcat('MSN');type{2}=strcat('FSI');type{3}=strcat('CIN');type{4}= strcat('NoId-Vs');type{5}= strcat('type1');type{6}=strcat('type2');
type{7}=strcat('type3');type{8}=strcat('NoId-Vta');
figure(16);hold on;
bar(Perc_PwN);hold on;

for i=1:length(ChLevNeu)
    plot([i-0.5,i+0.5],[ChLevNeu(i), ChLevNeu(i)],'k.-','linewidth',1.5);hold on
    text(i-0.2,Perc_PwN(i)+0.05,type{i}),
end
title(TitleP);
ylim([0,max(Perc_PwN)+0.08]);
box on;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % filename = 'Mtx_LabelsAll_Lasterev1.xlsx';
% % % Write the table
% % writetable(TableLabelsAll,filename,'Sheet',1);
% % filename = 'Mtx_LabelsPairs_Lasterev1.xlsx';
% % % Write the table
% % writetable(TableLabelsPairs,filename,'Sheet',1);
