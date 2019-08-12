%% Script to count whether ensembles do prefer specific cathegories of units
addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM');
% % % % % % % % The matrix new_data is the matrix that contains the info, in particular new_data.spikeT_BegEnd has the 
% % % % % % % % spiketrain cutted properly

name1=strcat('Rev2_OrRevW.mat');
load(name1)
alg= strcat('Stand_');
% alg=strcat('Pru_');
name=strcat('An_Disp_',alg,name1);
% name=strcat('An_Disp_',name);
load(name);
name_bins=strcat('BinLag_',alg,name1);
load(name_bins)
% MatBin=-MBinAnP;
nameSmallBig=strcat('SmallBig_An_Disp',alg,name1)
load(nameSmallBig)
%%
fignum=29;
fignum1=30;

clear c
load('classlist.mat');


Region='vsvta';
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


nnr=nn(a);
[new_data] = parId(new_data,a,nnr);


for k=1:length(new_data.par)
 [new_data] = labels(new_data,c,k);
 end
 clear tt i events licks laser info spikes map params spM spikes spike_region j k new_spM

%%%%%%%%%%%%%%%% Find number of units of the Striatum and VTA nneuS

new_spM=new_data;

% TotAn=size(new_spM.spike_regionNoId,2);
TotAn=length(new_data.par);
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
%%%%%


%%%%% Find Pairs in one region and shared between two differnt regions
% pairs diveded by regions with all info
[pairs_r_small,pairs_vsvs_small,pairs_vsvta_small,pairs_vtavta_small] = PairsInfo(As_across_smallBins,As_order_smallBins,nneuVS);
[pairs_r_large,pairs_vsvs_large,pairs_vsvta_large,pairs_vtavta_large] = PairsInfo(As_across_largeBins,As_order_largeBins,nneuVS);

clear a b A B i ii j k

[struct_pair_small] = PairsOrgByPair(pairs_vsvta_small);
[struct_pair_large] = PairsOrgByPair(pairs_vsvta_large);


[struct_pair_small] = PairNeuLabelsID(struct_pair_small,new_data);
[struct_pair_large] = PairNeuLabelsID(struct_pair_large,new_data);
copy_struct_pair_small=struct_pair_small;
copy_struct_pair_large=struct_pair_large;
clear struct_pair
%%%%%% Create the struct coerently with the directionality
%%%%%%%%%% Here the directionality doesn't cointan the info about the bins
%%%%%%%%%% because the pairs are already divided for bins

 Dir1 = 'vs->vta';
 Dir2 = 'vta->vs';
 Dir3 = 'sync';
 Dir4 = 'no';
% Dir='vs->vta Inv';
% Dir='vs->vta AllBins';

%  
%  Dir='vta->vs Inv';
%  Dir='vta->vs AllBins';
%    Dir=0;
% Dir='sync';



for k=1:TotAn
    if ~isempty(copy_struct_pair_small)
[struct_pair_small{k}] = DirAndBinPostPru(Dir4, copy_struct_pair_small,k);
 
    end
    if ~isempty(copy_struct_pair_large)
    [struct_pair_large{k}] = DirAndBinPostPru(Dir4, copy_struct_pair_large,k);
    
    end
end

% TotAn=length(new_data.par);
% count=zeros(TotAn,7);
% for k=1:TotAn
% [count] = CountTypeNeuPair(struct_pair{k},k,count);



%% Do the assemblies prefer specific units?
% Total Number of assemblies

for k=1:TotAn
   % if ~isempty(struct_pair{k}.pair)
        count_Pairs_small(k)=length(struct_pair_small{k}.pair);
        count_Pairs_large(k)=length(struct_pair_large{k}.pair);
       
   % end
end

TotAs_small=sum(count_Pairs_small);
TotAs_large=sum(count_Pairs_large);

%%%%%%%%%%%%%%
%%%%%%%%%% Number of assembly with a specific class of neurons %%%%%%% !!!!! 

%%% Small bins different directionality
count_PwN_small=zeros(TotAn,8);
 for k=1:TotAn
[TLabels_PwN_small,count_PwN_small] = CountTypeNeuPair(struct_pair_small{k},k,count_PwN_small); % pairs with specific units, animal per animal
 end
 
 %%% Large bins different directionality
count_PwN_large=zeros(TotAn,8);
 for k=1:TotAn
[TLabels_PwN_large,count_PwN_large] = CountTypeNeuPair(struct_pair_large{k},k,count_PwN_large); % pairs with specific units, animal per animal
 end


%%%%%%%%%%%%%%%%% For each Animal (Then i will plot in Box plot)
countMat_Pairs_small=repmat(count_Pairs_small,8,1);
countMat_Pairs_large=repmat(count_Pairs_large,8,1);

Perc_PwN_pAn_small=count_PwN_small./countMat_Pairs_small'; % # assemblies with V_{ri}/ # assemblies Animal per animal
Perc_PwN_pAn_large=count_PwN_large./countMat_Pairs_large'; % # assemblies with V_{ri}/ # assemblies Animal per animal

% TotAn=length(new_data.par);
count_Neu=zeros(TotAn,8);
for k=1:TotAn
[TLabelsNeu,count_Neu] = CountTypeNeuAll(new_data.par{k},k,count_Neu);
end

ChLevNeu_pAn= NaN(TotAn,8);
copy_nneuVS=nneuVS;
copy_nneuVTA=nneuVTA;
for k=1:TotAn
    if copy_nneuVS{k}==0, copy_nneuVS{k}=10^8; end
    if copy_nneuVTA{k}==0, copy_nneuVTA{k}=10^8; end
    ChLevNeu_pAn(k,1:4) = count_Neu(k,1:4)./copy_nneuVS{k};
    ChLevNeu_pAn(k,5:8) = count_Neu(k,5:8)./copy_nneuVTA{k};
end

%% Conditional Probability = Joint Probability/ Probability to have assembly with a specific neuron-type
count_Joint_small = zeros(TotAn, 16);
count_Joint_large = zeros(TotAn, 16);
% count_Joint_VTA_VS_small = zeros(TotAn, 16);
% count_Joint_VTA_VS_large = zeros(TotAn, 16);


% Probability to have assembly with a specific neuron-type (already calculated but here repmat)
Region = 'vs';
[TabNumPwN_VS_small,NumPwN_VS_small] = PairsWSpecNeu_TabNum(TLabels_PwN_small,Region);  % Number of pairs with specific neurons repeated in such a way to do the ratio with the joint
[TabNumPwN_VS_large,NumPwN_VS_large] = PairsWSpecNeu_TabNum(TLabels_PwN_large,Region);
Region = 'vta';
[TabNumPwN_VTA_small,NumPwN_VTA_small] = PairsWSpecNeu_TabNum(TLabels_PwN_small,Region);
[TabNumPwN_VTA_large,NumPwN_VTA_large] = PairsWSpecNeu_TabNum(TLabels_PwN_large,Region);


% Joint Probability
for k=1:TotAn
[TLabelsJoint_small,count_Joint_small]=CountTypeNeuPairTwoSides(struct_pair_small{k},k,count_Joint_small); %joint probability to have vs_k and vta_k
[TLabelsJoint_large,count_Joint_large]=CountTypeNeuPairTwoSides(struct_pair_large{k},k,count_Joint_large);
end
% Conditional Probability vs|vta
given='vs|vta';
[TabCondProb_VS_VTA_small,CondProb_VS_VTA_small] = CondProbAsType(count_Joint_small,NumPwN_VTA_small,given);
[TabCondProb_VS_VTA_large,CondProb_VS_VTA_large] = CondProbAsType(count_Joint_large,NumPwN_VTA_large,given);


% Conditional probability vta|vs
given= 'vta|vs' ;
[TabCondProb_VTA_VS_small,CondProb_VTA_VS_small] = CondProbAsType(count_Joint_small,NumPwN_VS_small,given);
[TabCondProb_VTA_VS_large,CondProb_VTA_VS_large] = CondProbAsType(count_Joint_large,NumPwN_VS_large,given);

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
[ChanceLevPairs_VS_VTA_small] = ChanceLevelCond(TLabels_PwN_small, given);
[ChanceLevPairs_VS_VTA_large] = ChanceLevelCond(TLabels_PwN_large, given);

given = ('vta|vs');
[ChanceLevPairs_VTA_VS_small] = ChanceLevelCond(TLabels_PwN_small, given);
[ChanceLevPairs_VTA_VS_large] = ChanceLevelCond(TLabels_PwN_large, given);



%%%%% Ratio Assemblies bidirections/ chance level
RatioPairs_VS_VTA_small = CondProb_VS_VTA_small./ChanceLevPairs_VS_VTA_small;
RatioPairs_VS_VTA_large = CondProb_VS_VTA_large./ChanceLevPairs_VS_VTA_large;
RatioPairs_VTA_VS_small = CondProb_VTA_VS_small./ChanceLevPairs_VTA_VS_small;
RatioPairs_VTA_VS_large = CondProb_VTA_VS_large./ChanceLevPairs_VTA_VS_large;


%% Are there typologies of neurons most probably deputed than others to form assemblies?

% functions to count how many neurons of a cathegories are part of one
% assembly (the neurons are taken only once), the first funtion eliminates
% double elements, the second one counts

for k=1:TotAn
    Un_small{k}=struct('neuID',[],'ia',[],'ic',[],'labels',[]);
    [Un_small{k}.neuID,Un_small{k}(:,1).ia(:,1),Un_small{k}.ic(:,1)]=unique(struct_pair_small{k}.neuID);
%     [Un{k}.neuID,Un{k}.ia,Un{k}.ic]=unique(struct_pair{k}.neuID(:,1));
    Un_small{k}.labels(:,1) = struct_pair_small{k}.labels(Un_small{k}.ia); % I create here the structure with the index unique
       
       
       
    Un_large{k}=struct('neuID',[],'ia',[],'ic',[],'labels',[]);
    [Un_large{k}.neuID,Un_large{k}(:,1).ia(:,1),Un_large{k}.ic(:,1)]=unique(struct_pair_large{k}.neuID);
    Un_large{k}.labels(:,1) = struct_pair_large{k}.labels(Un_large{k}.ia); % I create here the structure with the index unique
end


count_NinP_small=zeros(TotAn,8);
for k=1:TotAn
    if~isempty(Un_small{k}.neuID)
        for i=1:length(Un_small{k}.neuID)
[TabCount_NinP_small,count_NinP_small] = CountTypeNeuGen(Un_small{k}.labels(i),k,count_NinP_small); 
        end
    end
end


count_NinP_large=zeros(TotAn,8);
for k=1:TotAn
    if~isempty(Un_large{k}.neuID)
        for i=1:length(Un_large{k}.neuID)
[TabCount_NinP_large,count_NinP_large] = CountTypeNeuGen(Un_large{k}.labels(i),k,count_NinP_large); 
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

Sum_NinP_small = sum(TabCount_NinP_small{:,:});
Sum_NinP_large = sum(TabCount_NinP_large{:,:});
Sum_Neu = sum(TabCount_Neu{:,:});
% other way to do that
% for i=1: size(count_Neu,2)
% Sum_NinP(i) = sum(TLabels_NinP{find(count_Neu(:,i)>0),i});
% Sum_Neu(i) = sum(TLabels_Neu{find(count_Neu(:,i)>0),i});
% 
% end



Perc_NinP_small_AnTog=Sum_NinP_small./Sum_Neu; % Here I'm considering everything omiting nan results Take the sum and no distiction for animal
Perc_NinP_small_AnDist=count_NinP_small./count_Neu;  % divided for animal, then I can take the mean %%% !!!!!This one is coerent

Perc_NinP_large_AnTog=Sum_NinP_large./Sum_Neu; % Here I'm considering everything omiting nan results Take the sum and no distiction for animal
Perc_NinP_large_AnDist=count_NinP_large./count_Neu;  % divided for animal, then I can take the mean %%% !!!!!This one is coerent
% PR=sum(PercNeu1,1,'omitnan'); % Here is equivalent to assign 0 propability for a neuron that doesn't exist
% PRM=PR/TotAn;
% for i=1: size(countNeu,2)
%   PRMM(i) = mean(PercNeu1(find(countNeu(:,i)>0),i)); % when there is find>0 is equivalent to mean(PercNeu1,'omitnan')
% end


%%%%%%%%%%%%

count=Perc_NinP_small_AnTog; %neurons in pairs animal together
MSN = count(:,1); FSI = count(:,2); CIN = count(:,3); NoID_Vs = count(:,4);  type1 = count(:,5); type2 =count(:,6);
type3 = count(:,7); NoID_Vta = count(:,8); 
TLabelsPercNeuAnTog_small = table(MSN, FSI, CIN, NoID_Vs, type1, type2, type3, NoID_Vta);
clear count MSN FSI CNI NoID_Vs type1 type2 type3 NoID_Vta

count=Perc_NinP_large_AnTog; %neurons in pairs animal together
MSN = count(:,1); FSI = count(:,2); CIN = count(:,3); NoID_Vs = count(:,4);  type1 = count(:,5); type2 =count(:,6);
type3 = count(:,7); NoID_Vta = count(:,8); 
TLabelsPercNeuAnTog_large = table(MSN, FSI, CIN, NoID_Vs, type1, type2, type3, NoID_Vta);
clear count MSN FSI CNI NoID_Vs type1 type2 type3 NoID_Vta


count=Perc_NinP_small_AnDist; %neurons in Pairs/animal distinct
MSN = count(:,1); FSI = count(:,2); CIN = count(:,3); NoID_Vs = count(:,4);  type1 = count(:,5); type2 =count(:,6);
type3 = count(:,7); NoID_Vta = count(:,8); 
TLabelsPercNeuAnDist_small = table(MSN, FSI, CIN, NoID_Vs, type1, type2, type3, NoID_Vta);
clear count MSN FSI CNI NoID_Vs type1 type2 type3 NoID_Vta

count=Perc_NinP_large_AnDist; %neurons in Pairs/animal distinct
MSN = count(:,1); FSI = count(:,2); CIN = count(:,3); NoID_Vs = count(:,4);  type1 = count(:,5); type2 =count(:,6);
type3 = count(:,7); NoID_Vta = count(:,8); 
TLabelsPercNeuAnDist_large = table(MSN, FSI, CIN, NoID_Vs, type1, type2, type3, NoID_Vta);
clear count MSN FSI CNI NoID_Vs type1 type2 type3 NoID_Vta

%% %%%%%% Plot All together


Text=strcat('Post Pruned---',name1);
% Text=strcat(name1);
type{1}=strcat('MSN');type{2}=strcat('FSI');type{3}=strcat('CIN');type{4}= strcat('NoId-Vs');type{5}= strcat('type1');type{6}=strcat('type2');
type{7}=strcat('type3');type{8}=strcat('NoId-Vta');

VsVtaType{1} = strcat('M-T1'); VsVtaType{2} = strcat('M-T2'); VsVtaType{3} = strcat('M-T3'); VsVtaType{4} = strcat('M-No');
VsVtaType{5} = strcat('F-T1'); VsVtaType{6} = strcat('F-T2'); VsVtaType{7} = strcat('F-T3'); VsVtaType{8} = strcat('F-No');
VsVtaType{9} = strcat('C-T1'); VsVtaType{10} = strcat('C-T2'); VsVtaType{11} = strcat('C-T3'); VsVtaType{12} = strcat('C-No');
VsVtaType{13} = strcat('No-T1'); VsVtaType{14} = strcat('No-T2'); VsVtaType{15} = strcat('No-T3'); VsVtaType{16} = strcat('No-No');



ysup=max([max(max(Perc_PwN_pAn_small./ChLevNeu_pAn)),max(max(Perc_PwN_pAn_large./ChLevNeu_pAn))]);
yinf =min([min(min(Perc_PwN_pAn_small./ChLevNeu_pAn)),min(min(Perc_PwN_pAn_large./ChLevNeu_pAn))]);
figure(fignum1);hold on;
% boxplot(Perc_PwN_pAn)

subplot(3,2,1)
for i=1:size(Perc_NinP_small_AnTog,2)
    text(i-0.2,Perc_NinP_small_AnTog(i)+0.1,type{i}),
end
ylim([0,max(max([mean(Perc_NinP_small_AnDist,1,'omitnan');Perc_NinP_small_AnTog]))+0.2])
xlim([0,size(Perc_NinP_small_AnTog,2)+1])
% title(TitleP);
conf=[mean(Perc_NinP_small_AnDist,1,'omitnan');Perc_NinP_small_AnTog];
  b=bar(conf','group');hold on;
%   b(1).color='b';
%   b(2).color='y';
  legend([b(1) b(2)],'Mean','AnTog');
  title('Percentage of Neurons in Assembly - Small Bins');
  
  subplot(3,2,2)
for i=1:size(Perc_NinP_large_AnTog,2)
    text(i-0.2,Perc_NinP_large_AnTog(i)+0.1,type{i}),
end
ylim([0,max(max([mean(Perc_NinP_large_AnDist,1,'omitnan');Perc_NinP_large_AnTog]))+0.2])
xlim([0,size(Perc_NinP_small_AnTog,2)+1])
% title(TitleP);
conf=[mean(Perc_NinP_large_AnDist,1,'omitnan');Perc_NinP_large_AnTog];
  b=bar(conf','group');hold on;
%   b(1).color='b';
%   b(2).color='y';
  legend([b(1) b(2)],'Mean','AnTog');
  title('Percentage of Neurons in Assembly - Large Bins');
  
% bar(mean(Perc_NinP_AnDist,1,'omitnan'),'r');
subplot(3,2,3)
for i=1:size(ChLevNeu_pAn,2)
    text(i-0.2,ysup,type{i});hold on;
end
plot([0 size(ChLevNeu_pAn,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((Perc_PwN_pAn_small./ChLevNeu_pAn)','b*');hold on;
boxplot(Perc_PwN_pAn_small./ChLevNeu_pAn);hold on;
ylim([(yinf-0.5) (ysup+0.5)]);
title('Assembly type (one side) Ratio - Small Bins');
text(-0.5,0,Text,'Rotation',90,'Fontsize',14);hold on;


subplot(3,2,4)
for i=1:size(ChLevNeu_pAn,2)
    text(i-0.2,ysup,type{i});hold on;
end
plot([0 size(ChLevNeu_pAn,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((Perc_PwN_pAn_large./ChLevNeu_pAn)','b*');hold on;
boxplot(Perc_PwN_pAn_large./ChLevNeu_pAn);hold on;
ylim([(yinf-0.5) (ysup+0.5)]);
title('Assembly type (one side) Ratio - Large Bins');
% ylim([(yinf-0.5) (ysup+0.5)]);

subplot(3,2,5)
for i=1:size(ChanceLevPairs_VS_VTA_small,2)
    text(i-0.2,3+0.1,VsVtaType{i}, 'Fontsize',8);hold on;
end
plot([0 size(ChanceLevPairs_VS_VTA_small,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((RatioPairs_VS_VTA_small)','b*');hold on;
boxplot(RatioPairs_VS_VTA_small);hold on;
% plot((RatioPairs_VTA_VS_small)','m*');hold on;
% boxplot(RatioPairs_VTA_VS_small);hold on;
ylim([(yinf-0.5) (ysup+0.5)]);
% legend([box(1) box(2)],'P(Vs|Vta)/P(Vs)','P(Vta|Vs)/P(Vta)');
title('Assembly type Ratio - Small Bins -')


subplot(3,2,6)
for i=1:size(ChanceLevPairs_VS_VTA_large,2)
    text(i-0.2,3+0.1,VsVtaType{i},'Fontsize',8);hold on;
end
plot([0 size(ChanceLevPairs_VS_VTA_large,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
plot((RatioPairs_VS_VTA_large)','b*');hold on;
boxplot(RatioPairs_VS_VTA_large);hold on;
% plot((RatioPairs_VTA_VS_large)','m*');hold on;
% boxplot(RatioPairs_VTA_VS_large,'Color','m');hold on;
ylim([(yinf-0.5) (ysup+0.5)]);
% legend([box(1) box(2)],'P(Vs|Vta)/P(Vs)','P(Vta|Vs)/P(Vta)');
title('Assembly type Ratio - Large Bins - ')


SaveName = strcat('PlotAsNeutypeNoDirSL',name1);

save(SaveName)
% % % for i=1:size(ChLevNeu_pAn,2)
% % %     text(i-0.2,3+0.05,type{i});hold on;
% % % end
% % % plot([0 size(ChLevNeu_pAn,2)+1 ],[1,1],'k.-','Linewidth',1.5);hold on;
% box on;
% 
%  type{1}=strcat('MSN');type{2}=strcat('FSI');type{3}=strcat('CIN');type{4}= strcat('NoId-Vs');type{5}= strcat('type1');type{6}=strcat('type2');
% type{7}=strcat('type3');type{8}=strcat('NoId-Vta');
% figure(fignum);hold on;
% % boxplot(Perc_PwN_pAn)
% % boxplot(Perc_PwN_pAn./ChLevNeu_pAn);hold on;
% for i=1:size(Perc_NinP_AnTog,2)
% %     plot([i-0.5,i+0.5],[ChLevNeu(i), ChLevNeu(i)],'k.-','linewidth',1.5);hold on
%     text(i-0.2,Perc_NinP_AnTog(i)+0.1,type{i}),
% end
% ylim([0,max(max([mean(Perc_NinP_AnDist,1,'omitnan');Perc_NinP_AnTog]))+0.2])
% title(TitleP);
% conf=[mean(Perc_NinP_AnDist,1,'omitnan');Perc_NinP_AnTog]
%   b=bar(conf','group');hold on;
%   b(1).color='b';
%   b(2).color='y';
%   legend([b(1) b(2)],'Mean','AnTog');
% % bar(mean(Perc_NinP_AnDist,1,'omitnan'),'r');
% box on;
 
% end
%%%%%%%%

% % % % % for i=1:NumDiv
% % % % % for time=0:3
% % % % % TypeI(:,i) = (NumPVtaWN_large(:,1+(time*NumDiv)));
% % % % % TypeII(:,i) = (NumPVtaWN_large(:,2+(time*NumDiv)));
% % % % % TypeIII(:,i) = (NumPVtaWN_large(:,3+(time*NumDiv)));
% % % % % NoIdVTA(:,i) = (NumPVtaWN_large(:,4+(time*NumDiv)));
% % % % % end
% % % % % end

% Joint Probability (for the join probability I would not need the specific given, but I specify in this function
% a voice to put the column of the joint matrix in right order in such a way that I can have the right Conditional probability. The elements of the Joint are the same for the two direction but the columns have different order)

% given = ('vs|vta')
% for k=1:TotAn
% [TLabelsJoint_VS_VTA_small,count_Joint_VS_VTA_small]=CountTypeNeuPairTwoSides(struct_pair_small{k},k,count_Joint_VS_VTA_small,given); %joint probability to have vs_k and vta_k
% [TLabelsJoint_VS_VTA_large,count_Joint_VS_VTA_large]=CountTypeNeuPairTwoSides(struct_pair_large{k},k,count_Joint_VS_VTA_large,given);
% end
%Joint Probability 
% given = ('vta|vs')
% for k=1:TotAn
% [TLabelsJoint_VTA_VS_small,count_Joint_VTA_VS_small] = CountTypeNeuPairTwoSides(struct_pair_small{k},k,count_Joint_VTA_VS_small,given); %joint probability to have vs_k and vta_k
% [TLabelsJoint_VTA_VS_large,count_Joint_VTA_VS_large] = CountTypeNeuPairTwoSides(struct_pair_large{k},k,count_Joint_VTA_VS_large,given);
% end



%%%%% Maketables
% % % % count=count_NinP_small;
% % % % 
% % % % MSN = count(:,1); FSI = count(:,2); CIN = count(:,3); NoID_Vs = count(:,4); type1 = count(:,5); type2 =count(:,6);
% % % % type3 = count(:,7); NoID_Vta = count(:,8); 
% % % % TLabels_NinP_small = table(MSN, FSI, CIN, NoID_Vs, type1, type2, type3, NoID_Vta);
% % % % clear count MSN FSI CNI NoID_Vs type1 type2 type3 NoID_Vta
% % % % 
% % % % 
% % % % count=count_NinP_large;
% % % % MSN = count(:,1); FSI = count(:,2); CIN = count(:,3); NoID_Vs = count(:,4); type1 = count(:,5); type2 =count(:,6);
% % % % type3 = count(:,7); NoID_Vta = count(:,8); 
% % % % TLabels_NinP_large = table(MSN, FSI, CIN, NoID_Vs, type1, type2, type3, NoID_Vta);
% % % % clear count MSN FSI CNI NoID_Vs type1 type2 type3 NoID_Vta
% % % % 
% % % % count=count_Neu;
% % % % MSN = count(:,1); FSI = count(:,2); CIN = count(:,3); NoID_Vs = count(:,4);  type1 = count(:,5); type2 =count(:,6);
% % % % type3 = count(:,7); NoID_Vta = count(:,8); 
% % % % TLabels_Neu = table(MSN, FSI, CIN, NoID_Vs, type1, type2, type3, NoID_Vta);
% % % % clear count MSN FSI CNI NoID_Vs type1 type2 type3 NoID_Vta