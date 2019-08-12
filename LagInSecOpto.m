%  MergeDataGenPostPru_1
%  MergeDataGenPostPru_2 %postpruned 
% MergeDataGenPru  %pruned
% load('AfterFriedPairs_VeryAllRevPara_Rev_VsVta_wOp.mat')
%% %%%%%%%%%%%%%%%% Lag in seconds for optogenitic
path_d = '/zifnas/Carla/CWEM_Project_AfterFried/Data/AftAnFrPairs';
path_f = '/zifnas/Carla/CWEM_Project_BInAndLag/PLOTS';
path_p = '/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM';
cd(path_d)
load('AllRevPairs_All_SqRescale_2Ph.mat')
% load('AfterAnovaPairs_Late_Or_vsvta_wOp')
BinSizes=[0.01 0.015 0.03 0.05 0.08  0.12  0.25  0.35  0.5 0.6];
copy_pair_small1=pairs_small;
 clear pairs_small;
pairs_small=newnew_pairs_small;
% pairs_small=PairsBefFr{1};
% 
pairs_small=copy_pair_small1;

%%%%% VS classification %%%%%%%%%
%%%%% 1=MSN ; 2.1=FSI_low; 2.2=FSI_High; 3=CIN; 10 = not classified VS

%%%%%% VTA classification %%%%%%%%%
%%%%%% 4= type I; 5= typeII; 6= typeIII; 7.1= optoVTA_ex; 7.2=optoVS_ex; %
%%%%%% 7.3=optoVTA_inh; 7.4=optoVTA_inh; 20=Not classifiedVTA
cd(path_p)
cvs= [1, 2.1, 2.2, 3]; % VS neurons typologies
cvta = [7.1,4,5,6,0,7.3,20];%VTA neurons typologies
addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM')
for i =1:length(cvs) 
    for j=1:length(cvta)
    
%%%%% MSN with all identified VTA types
[~,~,~,LagInSec_small{i}{j},LagInSecNoZero_small{i}{j},LagInSecZero_small{i}{j}] = LagCountTypeMatPairs(pairs_small,BinSizes,cvs(i),cvta(j));
% [~,~,~,LagInSec_large_MT{i},LagInSecNoZero_large_MT{i},LagInSecZero_large_MT{i}] = LagCountTypeMatPairs(pairs_large,BinSizes,c1,cvta);

    end
end
clear i j
%%  Lag distribution for all VS cathegories and for VTA cathegories included lasert.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd(path_f)
% Title_Title = {'MSN + Type I';'MSN + Type II'; 'MSN + Type III';'MSN + Lasertag VTA'; 'MSN + Not.Lasertag VTA';...
%                'FSI Low + Type I';'FSI Low + Type II';'FSI Low + Type III';'FSI Low + Lasertag VTA'; 'FSI Low + Not.Lasertag VTA';...
%                'FSI High + Type I';'FSI High + Type II'; 'FSI High + Type III';'FSI High + Lasertag VTA';'FSI High + Not.Lasertag VTA';...
%                'CIN + Type I';'CIN + Type II'; 'CIN + Type III';'CIN + Lasertag VTA'; 'CIN + Not.Lasertag VTA'};
           Title_Title = {'MSN + Lasertag VTA';'MSN + Type I';'MSN + Type II'; 'MSN + Type III';...
               'FSI Low + Lasertag VTA';'FSI Low + Type I';'FSI Low + Type II';'FSI Low + Type III';...
               'FSI High + Lasertag VTA';'FSI High + Type I';'FSI High + Type II'; 'FSI High + Type III';...
               'CIN + Lasertag VTA';'CIN + Type I';'CIN + Type II'; 'CIN + Type III'};
           
nameT = {'All Rev Paradigmas'};
tick=-2:1:2;
x_hlag{1} = -1.8:0.09:1.8;
x_hlag_zero{1} = -0.045:0.09:0.045;
% x_hlag{2} = -3.5:0.2:3.5;
% x_hlag_zero{2} = -0.1:0.2:0.1;
len_cvta=length(cvta)-3;
figure(96);
for i = 0: length(cvs)-1
    for j = 1:len_cvta  %% not considering not classified and laser inh
        subplot(length(cvs),len_cvta,j+(i*len_cvta))
        h_NoZero = histogram(LagInSecNoZero_small{i+1}{j},x_hlag{1},'FaceAlpha',0.5); hold on;
        h_Zero = histogram(LagInSecZero_small{i+1}{j},x_hlag_zero{1},'FaceAlpha',0.5); hold on;
        h_NoZero.FaceColor ='b';
        h_Zero.FaceColor = 'g';
        title(Title_Title{j+(i*len_cvta)});hold on;
        ylabel('# Pairs');hold on;
        if (j+(i*len_cvta))>12
        xticks(tick);hold on;
        xticklabels(tick);hold on;
        else xticks([]);hold on; xticklabels([]);hold on; end
        ylim([0 25])
        ax=gca;
        ax.XColor='k';
        ax.YColor='k';
% xlim([-1.2 1.2])
if (i==0) & (j==2) 
    text(0.8,10,'\Delta \in [0.01, 0.25] sec','Fontsize',16);
end
if (i==3) & (j==4)
text(3, 1, nameT,'Rotation',270,'Fontsize',14);
end
        set(gca,'Fontsize',14)
    end
end
%% % Only four types MSN+type I (SPN+DAN)/MSN+TypeII (SPN+GABA)/FSILow+TypeI (FS+DAN)/FSIlow+TypeII(FS+GABA)

cvs 
cvta

vsn=[1,2.1]; vtan=[4,5]

vsi(1)=find(cvs==vsn(1)); vsi(2)=find(cvs==vsn(2)); vtai(1)=find(cvta==vtan(1)); vtai(2)=find(cvta==vtan(2));

TitleText={'SPN+DAN';'SPN+GABA';'FS Low+DAN';'FS Low+GABA'}
tick=-2:1:2;
x_hlag{1} = -1.8:0.09:1.8;
x_hlag_zero{1} = -0.045:0.09:0.045;
len_vtai=length(vtai);
figure(97);
for i = 0: length(vsi)-1
    for j = 1:len_vtai  %% not considering not classified and laser inh
        subplot(length(vsi),len_vtai,j+(i*len_vtai))
        h_NoZero = histogram(LagInSecNoZero_small{vsi(i+1)}{vtai(j)},x_hlag{1},'FaceAlpha',0.5); hold on;
        h_Zero = histogram(LagInSecZero_small{vsi(i+1)}{vtai(j)},x_hlag_zero{1},'FaceAlpha',0.5); hold on;
        h_NoZero.FaceColor ='b';
        h_Zero.FaceColor = 'g';
        title(TitleText{j+(i*len_vtai)});hold on;
        ylabel('# Pairs');hold on;
        if (j+(i*len_vtai))>2
        xticks(tick);hold on;
        xticklabels(tick);hold on;
        else xticks([]);hold on; xticklabels([]);hold on; end
        ylim([0 25])
        ax=gca;
        ax.XColor='k';
        ax.YColor='k';
% xlim([-1.2 1.2])
if (i==0) & (j==2) 
    text(0.8,10,'\Delta \in [0.01, 0.25] sec','Fontsize',16);
end
        set(gca,'Fontsize',14)
    end
end



%% %%%%%%  Only for four types:
Title_Text={'MSN + Type I';'MSN + lasertag VTA';'FSI High + Type I';'FSI High + lasertag VTA';'FSI Low + TypeI';'FSI Low + lasertag VTA'}
% nameT={'Lastrev1; Rev2; Rev3'}
nameT={'Very All Rev Paradigmas'}


tick=-1.5:0.5:1.5
x_hlag{1} = -1.8:0.1:1.8;
x_hlag_zero{1} = -0.05:0.1:0.05;
x_hlag{2} = -3.5:0.2:3.5;
x_hlag_zero{2} = -0.1:0.2:0.1;

figure(79); hold on;
subplot(3,2,1)
h_NoZero = histogram(LagInSecNoZero_small{1}{2},x_hlag{1},'FaceAlpha',0.5); hold on;
h_Zero = histogram(LagInSecZero_small{1}{2},x_hlag_zero{1},'FaceAlpha',0.5); hold on;
h_NoZero.FaceColor ='b';
h_Zero.FaceColor = 'g';
title(Title_Text{1});hold on;
% xlabel('Lag (sec)');hold on;
ylabel('# Pairs');hold on;
xticks(tick);hold on;
xticklabels(tick);hold on;
ylim([0 20])
xlim([-1.2 1.2])
text(0.8,10,'\Delta \in [0.01, 0.25] sec', 'Fontsize',12)
set(gca,'Fontsize',12)

subplot(3,2,2)
h_NoZero = histogram(LagInSecNoZero_small{1}{1},x_hlag{1},'FaceAlpha',0.5); hold on;
h_Zero = histogram(LagInSecZero_small{1}{1},x_hlag_zero{1},'FaceAlpha',0.5); hold on;
h_NoZero.FaceColor ='b';
h_Zero.FaceColor = 'g';
title(Title_Text{2});hold on;
% xlabel('Lag (sec)');hold on;
ylabel('# Pairs');hold on;
xticks(tick);hold on;
xticklabels(tick);hold on;
ylim([0 20])
xlim([-1.2 1.2])
% text(0.8,10,'\Delta \in [0.01, 0.25] sec')
set(gca,'Fontsize',12)

subplot(3,2,3)
h_NoZero = histogram(LagInSecNoZero_small{3}{2},x_hlag{1},'FaceAlpha',0.5); hold on;
h_Zero = histogram(LagInSecZero_small{3}{2},x_hlag_zero{1},'FaceAlpha',0.5); hold on;
h_NoZero.FaceColor ='b';
h_Zero.FaceColor = 'g';
title(Title_Text{3});hold on;
% xlabel('Lag (sec)'); hold on;
ylabel('# Pairs');hold on;
xticks(tick);hold on;
xticklabels(tick);hold on;
ylim([0 20])
xlim([-1.2 1.2])
% text(0.8,10,'\Delta \in [0.01, 0.25] sec')
set(gca,'Fontsize',12)
% text(-3, 8, nameT,'Rotation',90,'Fontsize',12)

subplot(3,2,4)
h_NoZero = histogram(LagInSecNoZero_small{3}{1},x_hlag{1},'FaceAlpha',0.5); hold on;
h_Zero = histogram(LagInSecZero_small{3}{1},x_hlag_zero{1},'FaceAlpha',0.5); hold on;
h_NoZero.FaceColor ='b';
h_Zero.FaceColor = 'g';
title(Title_Text{4});hold on;
% xlabel('Lag (sec)');hold on;
ylabel('# Pairs');hold on;
xticks(tick);hold on;
xticklabels(tick);hold on;
ylim([0 20])
xlim([-1.2 1.2])
% text(0.8,10,'\Delta \in [0.01, 0.25] sec')
set(gca,'Fontsize',12)

subplot(3,2,5)
h_NoZero = histogram(LagInSecNoZero_small{2}{2},x_hlag{1},'FaceAlpha',0.5); hold on;
h_Zero = histogram(LagInSecZero_small{2}{2},x_hlag_zero{1},'FaceAlpha',0.5); hold on;
h_NoZero.FaceColor ='b';
h_Zero.FaceColor = 'g';
title(Title_Text{5});hold on;
xlabel('Lag (sec)'); hold on;
ylabel('# Pairs');hold on;
xticks(tick);hold on;
xticklabels(tick);hold on;
ylim([0 20])
xlim([-1.2 1.2])
% text(0.8,10,'\Delta \in [0.01, 0.25] sec')
set(gca,'Fontsize',12)
% text(-3, 8, nameT,'Rotation',90,'Fontsize',12)

subplot(3,2,6)
h_NoZero = histogram(LagInSecNoZero_small{2}{1},x_hlag{1},'FaceAlpha',0.5); hold on;
h_Zero = histogram(LagInSecZero_small{2}{1},x_hlag_zero{1},'FaceAlpha',0.5); hold on;
h_NoZero.FaceColor ='b';
h_Zero.FaceColor = 'g';
title(Title_Text{6});hold on;
xlabel('Lag (sec)');hold on;
ylabel('# Pairs');hold on;
xticks(tick);hold on;
xticklabels(tick);hold on;
ylim([0 20])
xlim([-1.2 1.2])
% text(0.8,10,'\Delta \in [0.01, 0.25] sec')
text(1.5, 8, nameT,'Rotation',270,'Fontsize',12)
set(gca,'Fontsize',12)

