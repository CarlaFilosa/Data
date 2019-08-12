
% nameMerge{1} = strcat('Rev2_OrRevW.mat');
% nameMerge{2} = strcat('Lastrev1_OrRevW.mat');
% nameMerge{3} = strcat('Rev3_OrRevW.mat');
% nameMerge{4} = strcat('Revconc_OrRevW.mat');

nameMerge{1} = strcat('Lastrev1_OrRevW.mat');
nameMerge{2} = strcat('Revconc_OrRevW.mat');
nameMerge{3} = strcat('Succession_2PhW.mat');
nameT = strcat('Lastrev1; Revconc; Succession_2Ph');


alg= strcat('Stand_');

fignum =[1,2,3,4,5,6,7,8];

for i =1: length(nameMerge)
%     nameSmallBigMerge{i} = strcat('SmallBig_An_Disp',alg,nameMerge{i})
nameSmallBigMerge{i} = strcat('PostPru_BinRed_An_Disp',alg,nameMerge{i})
    load(nameSmallBigMerge{i})
clear c
load('classlist.mat');
%  nnn(i,:)=nn(a);
%  aa(i,:)=a;
[new_data] = parId(new_data,a,nn);


 for k=1:length(new_data.par)
 [new_data] = labels(new_data,c,k);
 end
 
 very_data{i}=new_data;

As_across_smallBinsT{i} = As_across_smallBins;
As_order_smallBinsT{i} = As_order_smallBins;
As_across_largeBinsT{i} = As_across_largeBins;
As_order_largeBinsT{i} = As_order_largeBins; 
clear new_data a nn As_across_smallBins As_order_smallBins As_across_largeBins As_order_largeBins
end
%%
 clearvars -except very_data As_across_smallBinsT As_order_smallBinsT As_across_largeBinsT As_order_largeBinsT c nameMerge nameT fignum
%%
for i =1:length(very_data)
    TotAn(i)=length(very_data{i}.par);
for k=1:TotAn(i)
    nneuVS{i}{k}=0;
    nneuVTA{i}{k}=0;
end
end

SumTotAn=sum(TotAn);

for i =1:length(very_data)
for k=1:TotAn(i)
    for j=1:size(very_data{i}.spike_regionNoId{k},2)
        [nneuVS{i}{k}] = countNeu(very_data{i}.spike_regionNoId{k}(1,j),1,nneuVS{i}{k});
        [nneuVTA{i}{k}] = countNeu(very_data{i}.spike_regionNoId{k}(1,j),2,nneuVTA{i}{k});
    end

end
end



As_across_smallBins_All=As_across_smallBinsT{1};
As_across_largeBins_All=As_across_largeBinsT{1};
As_order_smallBins_All=As_order_smallBinsT{1};
As_order_largeBins_All=As_order_largeBinsT{1};
nneuVS_All = nneuVS{1};
nneuVTA_All = nneuVTA{1};
Data_All.par=very_data{1}.par;
for i=1:length(very_data)-1
    [As_across_smallBins_All]=[As_across_smallBins_All,As_across_smallBinsT{i+1}];
    [As_across_largeBins_All]=[As_across_largeBins_All,As_across_largeBinsT{i+1}];
    [As_order_smallBins_All]=[As_order_smallBins_All,As_order_smallBinsT{i+1}];
    [As_order_largeBins_All]=[As_order_largeBins_All,As_order_largeBinsT{i+1}];
    [nneuVS_All]=[nneuVS_All,nneuVS{i+1}];
    [nneuVTA_All]=[nneuVTA_All,nneuVTA{i+1}];
    [Data_All.par]=[Data_All.par,very_data{i+1}.par];
end


%%
%%%%% Find Pairs in one region and shared between two differnt regions
% pairs diveded by regions with all info
[pairs_r_small,pairs_vsvs_small,pairs_vsvta_small,pairs_vtavta_small] = PairsInfo(As_across_smallBins_All,As_order_smallBins_All,nneuVS_All,SumTotAn);
[pairs_r_large,pairs_vsvs_large,pairs_vsvta_large,pairs_vtavta_large] = PairsInfo(As_across_largeBins_All,As_order_largeBins_All,nneuVS_All,SumTotAn);

clear a b A B i ii j k

[struct_pair_small] = PairsOrgByPair(pairs_vsvta_small);
[struct_pair_large] = PairsOrgByPair(pairs_vsvta_large);


[struct_pair_small] = PairNeuLabelsID(struct_pair_small,Data_All);
[struct_pair_large] = PairNeuLabelsID(struct_pair_large,Data_All);
copy_struct_pair_small=struct_pair_small;
copy_struct_pair_large=struct_pair_large;
clear struct_pair



%%%%%% Create the struct coerently with the directionality
%%%%%%%%%% Here the directionality doesn't cointan the info about the bins
%%%%%%%%%% because the pairs are already divided for bins

% % % % %Dir stays for "Directionality", can assume these values

 Dir1 = 'vs->vta';
 Dir2 = 'vta->vs';
 Dir3 = 'sync';
 Dir4 = 'no';


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



%% %%%%%% Do the assemblies prefer specific units?
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
count_Neu=zeros(SumTotAn,8);
for k=1:SumTotAn
[TabCount_Neu,count_Neu] = CountTypeNeuAll(Data_All.par{k},k,count_Neu);
end


ChLevNeu_pAn= NaN(SumTotAn,8);
copy_nneuVS_All=nneuVS_All;
copy_nneuVTA_All=nneuVTA_All;
for k=1:SumTotAn
    if copy_nneuVS_All{k}==0, copy_nneuVS_All{k}=nan; end
    if copy_nneuVTA_All{k}==0, copy_nneuVTA_All{k}=nan; end
    ChLevNeu_pAn(k,1:4) = count_Neu(k,1:4)./copy_nneuVS_All{k};
    ChLevNeu_pAn(k,5:8) = count_Neu(k,5:8)./copy_nneuVTA_All{k};
end

% Ratio with chance level 
RatioAs_small_vsvta=Perc_PwN_pAn_small_vsvta./ChLevNeu_pAn;
RatioAs_small_vtavs=Perc_PwN_pAn_small_vtavs./ChLevNeu_pAn;
RatioAs_small_sync=Perc_PwN_pAn_small_sync./ChLevNeu_pAn;


RatioAs_large_vsvta=Perc_PwN_pAn_large_vsvta./ChLevNeu_pAn;
RatioAs_large_vtavs=Perc_PwN_pAn_large_vtavs./ChLevNeu_pAn;
RatioAs_large_sync=Perc_PwN_pAn_large_sync./ChLevNeu_pAn;


 MeanRatioAs_small_vsvta = mean(RatioAs_small_vsvta,'omitnan');
 MeanRatioAs_small_vtavs = mean(RatioAs_small_vtavs,'omitnan');
 MeanRatioAs_small_sync = mean(RatioAs_small_sync,'omitnan');
 
 MeanRatioAs_large_vsvta = mean(RatioAs_large_vsvta,'omitnan');
 MeanRatioAs_large_vtavs = mean(RatioAs_large_vtavs,'omitnan');
 MeanRatioAs_large_sync = mean(RatioAs_large_sync,'omitnan');
 
 MeanRatioAs_small_vsvta(find(MeanRatioAs_small_vsvta==0))=nan;
 MeanRatioAs_small_vtavs(find(MeanRatioAs_small_vtavs==0))=nan;
 MeanRatioAs_small_sync(find(MeanRatioAs_small_sync==0))=nan;
 
 
 MeanRatioAs_large_vsvta(find(MeanRatioAs_large_vsvta==0))=nan;
 MeanRatioAs_large_vtavs(find(MeanRatioAs_large_vtavs==0))=nan;
 MeanRatioAs_large_sync(find(MeanRatioAs_large_sync==0))=nan;
 
 
 
StdRatioAs_small_vsvta = std(RatioAs_small_vsvta,0,1,'omitnan')/sqrt(SumTotAn);
StdRatioAs_small_vtavs = std(RatioAs_small_vtavs,0,1,'omitnan')/sqrt(SumTotAn);
StdRatioAs_small_sync = std(RatioAs_small_sync,0,1,'omitnan')/sqrt(SumTotAn);

StdRatioAs_large_vsvta = std(RatioAs_large_vsvta,0,1,'omitnan')/sqrt(SumTotAn);
StdRatioAs_large_vtavs = std(RatioAs_large_vtavs,0,1,'omitnan')/sqrt(SumTotAn);
StdRatioAs_large_sync = std(RatioAs_large_sync,0,1,'omitnan')/sqrt(SumTotAn);


%% Plot As Types

Text=strcat('Post Pruned---',nameT);
type{1}=strcat('MSN');type{2}=strcat('FSI');type{3}=strcat('CIN');type{4}= strcat('NoId-Vs');type{5}= strcat('type1');type{6}=strcat('type2');
type{7}=strcat('type3');type{8}=strcat('NoId-Vta');

ysupM = [max(MeanRatioAs_small_vsvta);max(MeanRatioAs_small_vtavs); max(MeanRatioAs_small_sync);...
       max(MeanRatioAs_large_vsvta);max(MeanRatioAs_large_vtavs); max(MeanRatioAs_large_sync)];
   
yinfM = [min(MeanRatioAs_small_vsvta);min(MeanRatioAs_small_vtavs); min(MeanRatioAs_small_sync);...
       min(MeanRatioAs_large_vsvta);min(MeanRatioAs_large_vtavs); min(MeanRatioAs_large_sync)];

figure(fignum(1));
subplot(3,2,1)
b=bar(MeanRatioAs_small_vsvta,'BaseValue',1);hold on;
b.FaceColor='c';hold on;
err=errorbar(MeanRatioAs_small_vsvta,StdRatioAs_small_vsvta,'k*');hold on;
for i=1:length(MeanRatioAs_small_vsvta)
    text(i-0.2,ysupM(1)+0.5,type{i});hold on;
end
ylim([(yinfM(1)-0.7) (ysupM(1)+0.9)]); hold on;
title('Small Bins vs->vta'); box on; hold on;

subplot(3,2,3)
b=bar(MeanRatioAs_small_vtavs,'BaseValue',1);hold on;
b.FaceColor='c';hold on;
err=errorbar(MeanRatioAs_small_vtavs,StdRatioAs_small_vtavs,'k*');hold on;
for i=1:length(MeanRatioAs_small_vsvta)
    text(i-0.2,ysupM(2)+0.5,type{i});hold on;
end
ylim([(yinfM(2)-0.7) (ysupM(2)+0.9)]); hold on;
text(-0.9,0,Text,'Rotation',90,'Fontsize',14);hold on;
title('Small Bins vta->vs');box on; hold on;

subplot(3,2,5)
b=bar(MeanRatioAs_small_sync,'BaseValue',1);hold on;
b.FaceColor='c';hold on;
err=errorbar(MeanRatioAs_small_sync,StdRatioAs_small_sync,'k*');hold on;
for i=1:length(MeanRatioAs_small_vsvta)
    text(i-0.2,ysupM(3)+0.5,type{i});hold on;
end
ylim([(yinfM(3)-0.7) (ysupM(3)+0.9)]); hold on;
title('Small Bins sync');box on; hold on;

subplot(3,2,2)
b=bar(MeanRatioAs_large_vsvta,'BaseValue',1);hold on;
b.FaceColor='c';hold on;
err=errorbar(MeanRatioAs_large_vsvta,StdRatioAs_large_vsvta,'k*');hold on;
for i=1:length(MeanRatioAs_small_vsvta)
    text(i-0.2,ysupM(4)+0.5,type{i});hold on;
end
ylim([(yinfM(4)-0.7) (ysupM(4)+0.9)]); hold on;
title('Large Bins vs->vta');box on; hold on;

subplot(3,2,4)
b=bar(MeanRatioAs_large_vtavs,'BaseValue',1);hold on;
b.FaceColor='c';hold on;
err=errorbar(MeanRatioAs_large_vtavs,StdRatioAs_large_vtavs,'k*');hold on;
for i=1:length(MeanRatioAs_small_vsvta)
    text(i-0.2,ysupM(5)+0.5,type{i});hold on;
end
ylim([(yinfM(5)-0.7) (ysupM(5)+0.9)]); hold on;
title('Large Bins vta->vs');box on; hold on;

subplot(3,2,6)
b=bar(MeanRatioAs_large_sync,'BaseValue',1);hold on;
b.FaceColor='c';hold on;
err=errorbar(MeanRatioAs_large_sync,StdRatioAs_large_sync,'k*');hold on;
for i=1:length(MeanRatioAs_small_vsvta)
    text(i-0.2,ysupM(6)+0.5,type{i});hold on;
end
ylim([(yinfM(6)-0.7) (ysupM(6)+0.9)]); hold on;
title('Large Bins sync');box on; hold on;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot AS with Directionality
%% Plot As Types
Text=strcat('Post Pruned---',nameT);
type{1}=strcat('MSN');type{2}=strcat('FSI');type{3}=strcat('CIN');type{4}= strcat('NoId-Vs');type{5}= strcat('type1');type{6}=strcat('type2');
type{7}=strcat('type3');type{8}=strcat('NoId-Vta');

ysup=max([max(max(Perc_PwN_pAn_small_vsvta./ChLevNeu_pAn)),max(max(Perc_PwN_pAn_small_vtavs./ChLevNeu_pAn)),max(max(Perc_PwN_pAn_small_sync./ChLevNeu_pAn)),...
    max(max(Perc_PwN_pAn_large_vsvta./ChLevNeu_pAn)),max(max(Perc_PwN_pAn_large_vtavs./ChLevNeu_pAn)),max(max(Perc_PwN_pAn_large_sync./ChLevNeu_pAn))]);
yinf =min([min(min(Perc_PwN_pAn_small_vsvta./ChLevNeu_pAn)),min(min(Perc_PwN_pAn_small_vtavs./ChLevNeu_pAn)),min(min(Perc_PwN_pAn_small_sync./ChLevNeu_pAn)),...
    min(min(Perc_PwN_pAn_large_vsvta./ChLevNeu_pAn)),min(min(Perc_PwN_pAn_large_vtavs./ChLevNeu_pAn)),min(min(Perc_PwN_pAn_large_sync./ChLevNeu_pAn))]);
figure(fignum(2));hold on;
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


CountNeuP_small_vtavs=zeros(SumTotAn,8);
for k=1:SumTotAn
    if~isempty(Un_small_vtavs{k}.neuID)
        for i=1:length(Un_small_vtavs{k}.neuID)
[TabNeuP_small_vtavs,CountNeuP_small_vtavs] = CountTypeNeuGen(Un_small_vtavs{k}.labels(i),k,CountNeuP_small_vtavs); 
        end
    end
end



CountNeuP_small_sync=zeros(SumTotAn,8);
for k=1:SumTotAn
    if~isempty(Un_small_sync{k}.neuID)
        for i=1:length(Un_small_sync{k}.neuID)
[TabNeuP_small_sync,CountNeuP_small_sync] = CountTypeNeuGen(Un_small_sync{k}.labels(i),k,CountNeuP_small_sync); 
        end
    end
end


%%%%%% Large Bins


CountNeuP_large_vsvta=zeros(SumTotAn,8);
for k=1:SumTotAn
    if~isempty(Un_large_vsvta{k}.neuID)
        for i=1:length(Un_large_vsvta{k}.neuID)
[TabNeuP_large_vsvta,CountNeuP_large_vsvta] = CountTypeNeuGen(Un_large_vsvta{k}.labels(i),k,CountNeuP_large_vsvta); 
        end
    end
end


CountNeuP_large_vtavs=zeros(SumTotAn,8);
for k=1:SumTotAn
    if~isempty(Un_large_vtavs{k}.neuID)
        for i=1:length(Un_large_vtavs{k}.neuID)
[TabNeuP_large_vtavs,CountNeuP_large_vtavs] = CountTypeNeuGen(Un_large_vtavs{k}.labels(i),k,CountNeuP_large_vtavs); 
        end
    end
end



CountNeuP_large_sync=zeros(SumTotAn,8);
for k=1:SumTotAn
    if~isempty(Un_large_sync{k}.neuID)
        for i=1:length(Un_large_sync{k}.neuID)
[TabNeuP_large_sync,CountNeuP_large_sync] = CountTypeNeuGen(Un_large_sync{k}.labels(i),k,CountNeuP_large_sync); 
        end
    end
end


% Function that count how many neurons of specific cathegories are present

% count_Neu=zeros(SumTotAn,8);
% for k=1:SumTotAn
% [TabCount_Neu,count_Neu] = CountTypeNeuAll(Data_All.par{k},k,count_Neu);
% end

CountType_All = zeros(SumTotAn,8);
for k =1: SumTotAn
    for i= 1:length(Data_All.par{k})
[TLabNeu,CountType_All] = CountTypeNeuGen(Data_All.par{k}(i).labels,k,CountType_All);
    end
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
%%%%%%%%%%%% Plot Neurons types

StdPerNeuP_small_vsvta = std(PercNeuP_small_vsvta_Dist,0,1,'omitnan')/sqrt(size(PercNeuP_small_vsvta_Dist,1));
StdPerNeuP_small_vtavs = std(PercNeuP_small_vtavs_Dist,0,1,'omitnan')/sqrt(size(PercNeuP_small_vtavs_Dist,1));
StdPerNeuP_small_sync = std(PercNeuP_small_sync_Dist,1,'omitnan')/sqrt(size(PercNeuP_small_sync_Dist,1));


StdPerNeuP_large_vsvta = std(PercNeuP_large_vsvta_Dist,1,'omitnan')/sqrt(size(PercNeuP_large_vsvta_Dist,1));
StdPerNeuP_large_vtavs = std(PercNeuP_large_vtavs_Dist,1,'omitnan')/sqrt(size(PercNeuP_large_vtavs_Dist,1));
StdPerNeuP_large_sync = std(PercNeuP_large_sync_Dist,1,'omitnan')/sqrt(size(PercNeuP_large_sync_Dist,1));

ZeroErr=zeros(size(StdPerNeuP_small_vsvta));
%% Plot Units 
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

figure(fignum(3));hold on;
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
 legend([b(1) b(2)],'Mean','AnTog','Location','northwest'); box on; hold on;
title('Neuron types. vs->vta Small Bins'); hold on;
  
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
legend([b(1) b(2)],'Mean','AnTog','Location','northwest'); hold on;
text(-1,0,Text,'Rotation',90,'Fontsize',14);box on; hold on;
title('Neuron types. vta->vs Small Bins'); hold on;

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
legend([b(1) b(2)],'Mean','AnTog','Location','northwest');box on; hold on; 
title('Neuron types. Sync Small Bins'); hold on;



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
legend([b(1) b(2)],'Mean','AnTog','Location','northwest'); box on;hold on;
title('Neuron types. vs->vta Large Bins'); hold on;
  
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
legend([b(1) b(2)],'Mean','AnTog','Location','northwest'); box on; hold on;
title('Neuron types. vta->vs Large Bins'); hold on;

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
legend([b(1) b(2)],'Mean','AnTog','Location','northwest'); box on;hold on;
title('Neuron types. Sync Large Bins'); hold on;


box on;

%%
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



%%
MeanRatioPairs_VS_VTA_small_vsvta = mean(RatioPairs_VS_VTA_small_vsvta,'omitnan');
MeanRatioPairs_VS_VTA_small_vtavs = mean(RatioPairs_VS_VTA_small_vtavs,'omitnan');
MeanRatioPairs_VS_VTA_small_sync = mean(RatioPairs_VS_VTA_small_sync,'omitnan');

MeanRatioPairs_VS_VTA_large_vsvta = mean(RatioPairs_VS_VTA_large_vsvta,'omitnan');
MeanRatioPairs_VS_VTA_large_vtavs = mean(RatioPairs_VS_VTA_large_vtavs,'omitnan');
MeanRatioPairs_VS_VTA_large_sync = mean(RatioPairs_VS_VTA_large_sync,'omitnan');


MeanRatioPairs_VS_VTA_small_vsvta(find(MeanRatioPairs_VS_VTA_small_vsvta==0)) = nan;
MeanRatioPairs_VS_VTA_small_vtavs(find(MeanRatioPairs_VS_VTA_small_vtavs==0))= nan;
MeanRatioPairs_VS_VTA_small_sync(find(MeanRatioPairs_VS_VTA_small_sync==0)) = nan;

MeanRatioPairs_VS_VTA_large_vsvta(find(MeanRatioPairs_VS_VTA_large_vsvta==0)) = nan;
MeanRatioPairs_VS_VTA_large_vtavs(find(MeanRatioPairs_VS_VTA_large_vtavs==0))= nan;
MeanRatioPairs_VS_VTA_large_sync(find(MeanRatioPairs_VS_VTA_large_sync==0)) = nan;


StdRatioPairs_VS_VTA_small_vsvta = std(RatioPairs_VS_VTA_small_vsvta,0,1,'omitnan')/sqrt(SumTotAn);
StdRatioPairs_VS_VTA_small_vtavs = std(RatioPairs_VS_VTA_small_vtavs,0,1,'omitnan')/sqrt(SumTotAn);
StdRatioPairs_VS_VTA_small_sync = std(RatioPairs_VS_VTA_small_sync,0,1,'omitnan')/sqrt(SumTotAn);

StdRatioPairs_VS_VTA_large_vsvta = std(RatioPairs_VS_VTA_large_vsvta,0,1,'omitnan')/sqrt(SumTotAn);
StdRatioPairs_VS_VTA_large_vtavs = std(RatioPairs_VS_VTA_large_vtavs,0,1,'omitnan')/sqrt(SumTotAn);
StdRatioPairs_VS_VTA_large_sync = std(RatioPairs_VS_VTA_large_sync,0,1,'omitnan')/sqrt(SumTotAn);
%% Plot As types two sides
Text=strcat('Post Pruned---',nameT);

VsVtaType{1} = strcat('M-T1'); VsVtaType{2} = strcat('M-T2'); VsVtaType{3} = strcat('M-T3'); VsVtaType{4} = strcat('M-No');
VsVtaType{5} = strcat('F-T1'); VsVtaType{6} = strcat('F-T2'); VsVtaType{7} = strcat('F-T3'); VsVtaType{8} = strcat('F-No');
VsVtaType{9} = strcat('C-T1'); VsVtaType{10} = strcat('C-T2'); VsVtaType{11} = strcat('C-T3'); VsVtaType{12} = strcat('C-No');
VsVtaType{13} = strcat('No-T1'); VsVtaType{14} = strcat('No-T2'); VsVtaType{15} = strcat('No-T3'); VsVtaType{16} = strcat('No-No');

for i= 1: length(MeanRatioPairs_VS_VTA_small_vsvta)
ysupLoc(i,:) = [(MeanRatioPairs_VS_VTA_small_vsvta(i) + StdRatioPairs_VS_VTA_small_vsvta(i));...
           (MeanRatioPairs_VS_VTA_small_vtavs(i) + StdRatioPairs_VS_VTA_small_vtavs(i));...
           (MeanRatioPairs_VS_VTA_small_sync(i) + StdRatioPairs_VS_VTA_small_sync(i));...
           (MeanRatioPairs_VS_VTA_large_vsvta(i) + StdRatioPairs_VS_VTA_large_vsvta(i));...
           (MeanRatioPairs_VS_VTA_large_vtavs(i) + StdRatioPairs_VS_VTA_large_vtavs(i));...
           (MeanRatioPairs_VS_VTA_large_sync(i) + StdRatioPairs_VS_VTA_large_sync(i))];
       
yinfLoc(i,:) = [(MeanRatioPairs_VS_VTA_small_vsvta(i) - StdRatioPairs_VS_VTA_small_vsvta(i));...
           (MeanRatioPairs_VS_VTA_small_vtavs(i) - StdRatioPairs_VS_VTA_small_vtavs(i));...
           (MeanRatioPairs_VS_VTA_small_sync(i) - StdRatioPairs_VS_VTA_small_sync(i));...
           (MeanRatioPairs_VS_VTA_large_vsvta(i) - StdRatioPairs_VS_VTA_large_vsvta(i));...
           (MeanRatioPairs_VS_VTA_large_vtavs(i) - StdRatioPairs_VS_VTA_large_vtavs(i));...
           (MeanRatioPairs_VS_VTA_large_sync(i) - StdRatioPairs_VS_VTA_large_sync(i))];
       
end

ysupMA=[max(MeanRatioPairs_VS_VTA_small_vsvta + StdRatioPairs_VS_VTA_small_vsvta),...
       max(MeanRatioPairs_VS_VTA_small_vtavs + StdRatioPairs_VS_VTA_small_vtavs),...
       max(MeanRatioPairs_VS_VTA_small_sync + StdRatioPairs_VS_VTA_small_sync),...
       max(MeanRatioPairs_VS_VTA_large_vsvta + StdRatioPairs_VS_VTA_large_vsvta),...
       max(MeanRatioPairs_VS_VTA_large_vtavs + StdRatioPairs_VS_VTA_large_vtavs),...
       max(MeanRatioPairs_VS_VTA_large_sync + StdRatioPairs_VS_VTA_large_sync)];
   
yinfMA=[min(MeanRatioPairs_VS_VTA_small_vsvta - StdRatioPairs_VS_VTA_small_vsvta),...
       min(MeanRatioPairs_VS_VTA_small_vtavs - StdRatioPairs_VS_VTA_small_vtavs),...
       min(MeanRatioPairs_VS_VTA_small_sync - StdRatioPairs_VS_VTA_small_sync),...
       min(MeanRatioPairs_VS_VTA_large_vsvta - StdRatioPairs_VS_VTA_large_vsvta),...
       min(MeanRatioPairs_VS_VTA_large_vtavs - StdRatioPairs_VS_VTA_large_vtavs),...
       min(MeanRatioPairs_VS_VTA_large_sync - StdRatioPairs_VS_VTA_large_sync)];


   figure(fignum(4)); hold on;
   subplot(3,2,1)
   b = bar(MeanRatioPairs_VS_VTA_small_vsvta,'BaseValue',1); hold on;
   b.FaceColor='c'; hold on;
   err = errorbar(MeanRatioPairs_VS_VTA_small_vsvta,StdRatioPairs_VS_VTA_small_vsvta,'k*'); hold on;
   for i=1:length(MeanRatioPairs_VS_VTA_small_vsvta)
    if ~isnan(MeanRatioPairs_VS_VTA_small_vsvta(i)) 
        if(MeanRatioPairs_VS_VTA_small_vsvta(i) >=1)
        text(i-0.2,ysupLoc(i,1)+0.2,VsVtaType{i}, 'Fontsize',8);hold on;
        elseif (MeanRatioPairs_VS_VTA_small_vsvta(i) < 1)
        text(i-0.2,yinfLoc(i,1)-0.2,VsVtaType{i}, 'Fontsize',8);hold on;
        end
    else
        text(i-0.2,ysupMA(1),VsVtaType{i}, 'Fontsize',8);hold on;
    end
   end
   ylim([yinfMA(1)-0.5; ysupMA(1)+0.5]);hold on;
   title('As Type Ratio Two Sides, vs->vta Small Bins'); hold on;
   box on;
   
   subplot(3,2,3)
   b = bar(MeanRatioPairs_VS_VTA_small_vtavs,'BaseValue',1); hold on;
   b.FaceColor='c'; hold on;
   err = errorbar(MeanRatioPairs_VS_VTA_small_vtavs,StdRatioPairs_VS_VTA_small_vtavs,'k*'); hold on;
   for i=1:length(MeanRatioPairs_VS_VTA_small_vtavs)
    if ~isnan(MeanRatioPairs_VS_VTA_small_vtavs(i)) 
        if(MeanRatioPairs_VS_VTA_small_vtavs(i) >=1)
        text(i-0.2,ysupLoc(i,2)+0.2,VsVtaType{i}, 'Fontsize',8);hold on;
        elseif (MeanRatioPairs_VS_VTA_small_vtavs(i) < 1)
        text(i-0.2,yinfLoc(i,2)-0.2,VsVtaType{i}, 'Fontsize',8);hold on;
        end
    else
        text(i-0.2,ysupMA(2),VsVtaType{i}, 'Fontsize',8);hold on;
    end
   end
   text(-1.5,0,Text,'Rotation',90,'Fontsize',14);box on; hold on;
    title('As Type Ratio Two Sides, vta->vs Small Bins'); hold on;
   ylim([yinfMA(2)-0.5; ysupMA(2)+0.5]);hold on;
   box on;  
   
   
   subplot(3,2,5)
   b = bar(MeanRatioPairs_VS_VTA_small_sync,'BaseValue',1); hold on;
   b.FaceColor='c'; hold on;
   err = errorbar(MeanRatioPairs_VS_VTA_small_sync,StdRatioPairs_VS_VTA_small_sync,'k*'); hold on;
   for i=1:length(MeanRatioPairs_VS_VTA_small_vtavs)
    if ~isnan(MeanRatioPairs_VS_VTA_small_sync(i)) 
        if(MeanRatioPairs_VS_VTA_small_sync(i) >=1)
        text(i-0.2,ysupLoc(i,3)+0.2,VsVtaType{i}, 'Fontsize',8);hold on;
        elseif (MeanRatioPairs_VS_VTA_small_sync(i) < 1)
        text(i-0.2,yinfLoc(i,3)-0.2,VsVtaType{i}, 'Fontsize',8);hold on;
        end
    else
        text(i-0.2,ysupMA(3),VsVtaType{i}, 'Fontsize',8);hold on;
    end
   end
   ylim([yinfMA(3)-0.5; ysupMA(3)+0.5]); hold on;
    title('As Type Ratio Two Sides, sync Small Bins'); hold on;
   box on;
   
   
   subplot(3,2,2)
   b = bar(MeanRatioPairs_VS_VTA_large_vsvta,'BaseValue',1); hold on;
   b.FaceColor='c'; hold on;
   err = errorbar(MeanRatioPairs_VS_VTA_large_vsvta,StdRatioPairs_VS_VTA_large_vsvta,'k*'); hold on;
    for i=1:length(MeanRatioPairs_VS_VTA_small_vtavs)
    if ~isnan(MeanRatioPairs_VS_VTA_large_vsvta(i)) 
        if(MeanRatioPairs_VS_VTA_large_vsvta(i) >=1)
        text(i-0.2,ysupLoc(i,4)+0.2,VsVtaType{i}, 'Fontsize',8);hold on;
    elseif (MeanRatioPairs_VS_VTA_large_vsvta(i) < 1)
        text(i-0.2,yinfLoc(i,4)-0.2,VsVtaType{i}, 'Fontsize',8);hold on;
        end
    else
        text(i-0.2,ysupMA(4),VsVtaType{i}, 'Fontsize',8);hold on;
    end
    end
    ylim([yinfMA(4)-0.5; ysupMA(4)+0.5]); hold on;
     title('As Type Ratio Two Sides, vs->vta Large Bins'); hold on;
   box on;
    
    
   subplot(3,2,4)
   b = bar(MeanRatioPairs_VS_VTA_large_vtavs,'BaseValue',1); hold on;
   b.FaceColor='c'; hold on;
   err = errorbar(MeanRatioPairs_VS_VTA_large_vtavs,StdRatioPairs_VS_VTA_large_vtavs,'k*'); hold on;
    for i=1:length(MeanRatioPairs_VS_VTA_small_vtavs)
    if ~isnan(MeanRatioPairs_VS_VTA_large_vtavs(i)) 
        if (MeanRatioPairs_VS_VTA_large_vtavs(i) >=1)
        text(i-0.2,ysupLoc(i,5)+0.2,VsVtaType{i}, 'Fontsize',8);hold on;
        elseif (MeanRatioPairs_VS_VTA_large_vtavs(i) <1)
        text(i-0.2,yinfLoc(i,5)-0.2,VsVtaType{i}, 'Fontsize',8);hold on;
        end
    else
        text(i-0.2,ysupMA(5),VsVtaType{i}, 'Fontsize',8);hold on;
    end
    end
    ylim([yinfMA(5)-0.5; ysupMA(5)+0.5]); hold on;
     title('As Type Ratio Two Sides, vta->vs Large Bins'); hold on;
   box on;
    
    
   subplot(3,2,6)
   b = bar(MeanRatioPairs_VS_VTA_large_sync,'BaseValue',1); hold on;
   b.FaceColor='c'; hold on;
   err = errorbar(MeanRatioPairs_VS_VTA_large_sync,StdRatioPairs_VS_VTA_large_sync,'k*'); hold on;
   for i=1:length(MeanRatioPairs_VS_VTA_small_vtavs)
    if ~isnan(MeanRatioPairs_VS_VTA_large_sync(i))
        if (MeanRatioPairs_VS_VTA_large_sync(i) >=1)
        text(i-0.2,ysupLoc(i,6)+0.2,VsVtaType{i}, 'Fontsize',8);hold on;
        elseif (MeanRatioPairs_VS_VTA_large_sync(i) <1)
        text(i-0.2,yinfLoc(i,6)-0.2,VsVtaType{i}, 'Fontsize',8);hold on;
        end
    else
        text(i-0.2,ysupMA(6),VsVtaType{i}, 'Fontsize',8);hold on;
    end
   end
   ylim([yinfMA(6)-0.5; ysupMA(6)+0.5]);hold on;
    title('As Type Ratio Two Sides, sync Large Bins'); hold on;
   box on;
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

figure(fignum(6));hold on;

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


SaveName = strcat('PlotAsNeutypeAndProbWithDirSL',nameT);

save(SaveName)



%%

% % % % % % name1=strcat('Rev2_OrRevW.mat');
% % % % % % 
% % % % % % % alg=strcat('Pru_');
% % % % % % nameSmallBig = strcat('SmallBig_An_Disp',alg,name1)
% % % % % % load(nameSmallBig)
% % % % % % 
% % % % % % clear c
% % % % % % load('classlist.mat');
% % % % % % 
% % % % % % nn1=nn(a);
% % % % % % a1=a;
% % % % % % [new_data] = parId(new_data,a,nn1);
% % % % % % 
% % % % % % 
% % % % % % for k=1:length(new_data.par)
% % % % % %  [new_data] = labels(new_data,c,k);
% % % % % % end
% % % % % %  
% % % % % % new_data1=new_data;
% % % % % % As_across_smallBins1 = As_across_smallBins;
% % % % % % As_order_smallBins1 = As_order_smallBins;
% % % % % % As_across_largeBins1 = As_across_largeBins;
% % % % % % As_order_largeBins1 = As_order_largeBins;
% % % % % % 
% % % % % % clearvars -except a1 new_data1 As_across_smallBins1 As_order_smallBins1 As_across_largeBins1 As_order_largeBins1
% % % % % % 
% % % % % % name1=strcat('Lastrev1_OrRevW.mat');
% % % % % % alg= strcat('Stand_');
% % % % % % % alg=strcat('Pru_');
% % % % % % nameSmallBig=strcat('SmallBig_An_Disp',alg,name1)
% % % % % % load(nameSmallBig)
% % % % % % clear c
% % % % % % load('classlist.mat');
% % % % % % 
% % % % % % nn2=nn(a);
% % % % % % a2=a;
% % % % % % [new_data] = parId(new_data,a,nn2);
% % % % % % 
% % % % % % for k=1:length(new_data.par)
% % % % % %  [new_data] = labels(new_data,c,k);
% % % % % % end
% % % % % %  
% % % % % % new_data2=new_data;
% % % % % % As_across_smallBins2 = As_across_smallBins;
% % % % % % As_order_smallBins2 = As_order_smallBins;
% % % % % % As_across_largeBins2 = As_across_largeBins;
% % % % % % As_order_largeBins2 = As_order_largeBins;
% % % % % % 
% % % % % % clearvars -except a1 new_data1 As_across_smallBins1 As_order_smallBins1 As_across_largeBins1 As_order_largeBins1 a2 new_data2 As_across_smallBins2 As_order_smallBins2 As_across_largeBins2 As_order_largeBins2
% % % % % % 
% % % % % % 
% % % % % % name12{1}=strcat('Rev2_OrRevW.mat');name12{2}=strcat('Lastrev1_OrRevW.mat');
% % % % % % new_new=[new_data1,new_data2]
