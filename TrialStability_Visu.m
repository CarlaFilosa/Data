%%% Visualization activity and trial



% for k= 1:SumTotAn
% [FindInd_ArRew{k}] = FindIndicesForAssem(Data_All,binref.B,start,stop,k,'reward',0.7,0.7);  %around reward
% end
% for k= 1:SumTotAn
% [FindInd_AftRew{k}] = FindIndicesForAssem(Data_All,binref.B,start,stop,k,'reward',0,0.5);  %around reward
% end
% 
% for k= 1:SumTotAn
% [FindInd_Base{k}] = FindIndicesForAssem(Data_All,binref.B,start,stop,k,'HitRew',1.4,-0.9); % supposed rewarded, not rewarded
% end
% %%% all 
% [asTime_all] = AsOnInd(AsTime,FindInd_all);

% %%% Odour 
% [asTime_AftStim] = AsOnInd(AsTime,FindInd_AftStim);
%%% ar Reward
% [asTime_ArRew] = AsOnInd(AsTime,FindInd_ArRew);
% %%%%%% Hit Not rewarded in the interval odour plus
% % [asTime_AftRew] = AsOnInd(AsTime,FindInd_AftRew);
% % %%%% around long reward interval
% % % [asTime_US_RewT] = AsOnInd(AsTime,FindIndRewT);
% % %%%% Around the first Lick
% % [asTime_Base] = AsOnInd(AsTime,FindInd_Base);
% 
% for k=1:length(pairs)
%     if ~isempty(pairs{k})
% [as_sel_ArRew{k},MeanTr_ArRew{k},MNormTr_ArRew{k},StETr_ArRew{k},MaxA_ArRew{k},AsForMean_ArRew{k}] = As_TaskCodePairsMatPP(Dir,pairs{k},Data_All,asTime_ArRew,hit,start(k),k);
%     else as_sel_ArRew{k}=[]; MeanTr_ArRew{k}=[]; MNormTr_ArRew{k}=[]; StETr_ArRew{k}=[]; MaxA_ArRew{k}=[]; AsForMean_ArRew{k}=[];
%     end
% end
% hhh=0;
% figure(8);
% for j=1:size(AsForMean_ArRew{2}{1},2)
%     hhh=hhh+max(AsForMean_ArRew{2}{1}(:,j));
% plot(1:size(AsForMean_ArRew{2}{1},1),AsForMean_ArRew{2}{1}(:,j),'b.-'); hold on;
% end
% plot(1:size(AsForMean_ArRew{2}{1},1), MeanTr_ArRew{2}{1},'k.-', 'Linewidth',2)
%%   Visualization activity in the periods we chose 
clear all

% load('AfterFriedPairs_Late_Or_VsVta.mat')
% SumTotAn=35;
for k=1:SumTotAn
% EpochMatrix{k}=[];
Base{k} = [];
AftStim{k} = [];
PreRew{k} =[];
AftRew{k} =[];
MBA{k}=[];
MAS{k}=[];
MPR{k}=[];
MAR{k}=[];
MTOT{k}=[];
end
for k=1:SumTotAn
    if ~isempty(AsForMean_Base{k})
    for i=1: length(AsForMean_Base{k})
MBA{k}{i}= max(AsForMean_Base{k}{i},[],2);
MAS{k}{i}= max(AsForMean_AftStim{k}{i},[],2);
MPR{k}{i}= max(AsForMean_PreRew{k}{i},[],2);
MAR{k}{i}= max(AsForMean_AftRew{k}{i},[],2);
MTOT{k}{i}= max([MBA{k}{i}(1:50),MAS{k}{i}(1:50),MPR{k}{i}(1:50),MAR{k}{i}(1:50)],[],2);
if MBA{k}{i}~=0
% Base{k}{i} =AsForMean_Base{k}{i}./MBA{k}{i};
Base{k}{i} =AsForMean_Base{k}{i}(1:50,:)./MTOT{k}{i}(1:50);
else Base{k}{i} =AsForMean_Base{k}{i}(1:50,:)
end
Base{k}{i}=Base{k}{i}';
Base{k}{i}(end+1,:)=nan;
if MAS{k}{i}~=0
% AftStim{k}{i} =AsForMean_AftStim{k}{i}./MAS{k}{i};
AftStim{k}{i} =AsForMean_AftStim{k}{i}(1:50,:)./MTOT{k}{i}(1:50);
else AftStim{k}{i} =AsForMean_AftStim{k}{i}(1:50,:)
end
AftStim{k}{i}=AftStim{k}{i}'
AftStim{k}{i}(end+1,:)=nan;
if MPR{k}{i}~=0
% PreRew{k}{i} =AsForMean_PreRew{k}{i}./MPR{k}{i};
PreRew{k}{i} =AsForMean_PreRew{k}{i}(1:50,:)./MTOT{k}{i}(1:50);
else PreRew{k}{i} =AsForMean_PreRew{k}{i}(1:50,:)
end
PreRew{k}{i}=PreRew{k}{i}'
PreRew{k}{i}(end+1,:)=nan;
if MAR{k}{i}~=0
% AftRew{k}{i} =AsForMean_AftRew{k}{i}./MAR{k}{i};
AftRew{k}{i} =AsForMean_AftRew{k}{i}(1:50,:)./MTOT{k}{i}(1:50);
else AftRew{k}{i} =AsForMean_AftRew{k}{i}(1:50,:)
end
AftRew{k}{i}=AftRew{k}{i}';
AftRew{k}{i}(end+1,:)=nan;
    end
    end
end
clear i k

%%

pmulti=pmulFM;

for k=1:length(pmulti)
    
    PvTableMT{k}=[];
    PvTableFT{k}=[];
end

for k=1:length(pmulti)
    for i=1:length(pmulti{k})
    PvTableMT{k}{i}=[];
    PvTableFT{k}{i}=[];
    end
end
colpval=6;neu1=1;neu2=4;
for k=1:length(pmulti)
    if ~isempty(pmulti{k})
        for i=1:length(pmulti{k})
            if ~isempty(pmulti{k}{i})
                if pmulti{k}{i}(1,9)==neu1 && pmulti{k}{i}(1,10)==neu2
[PvTableMT{k}{i}] = pval_TableGroups(pmulti{k}{i},colpval);
                end
            end
        end
    end
end

colpval=6;neu1=2;neu2=4;
for k=1:length(pmulti)
    if ~isempty(pmulti{k})
        for i=1:length(pmulti{k})
            if ~isempty(pmulti{k}{i}) 
                if pmulti{k}{i}(1,9)==neu1 && pmulti{k}{i}(1,10)==neu2
[PvTableFT{k}{i}] = pval_TableGroups(pmulti{k}{i},colpval);
                end
            end
        end
    end
end

PvTable=PvTableMT;
%%
fpath= '/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM/Signif_Assemblies';
ticks=[1 10 20 30 40 50]
Xlabels1 = [-1.4, -1.3, -1.2, -1.1,-1.0,-0.9];
Xlabels2 = [0, 0.1, 0.2, 0.3, 0.4, 0.5];
Xlabels3 = [-1.4, -1.3, -1.2, -1.1,-1.0,-0.9];
Xlabels4 = [0, 0.1, 0.2, 0.3, 0.4, 0.5];
temp=0;
%  for k=1:SumTotAn
%      if~isempty(Base{k})
%          for i=1:length(Base{k})
%%%
%%% Visualization assembly for typologies and task moments

fileT=strcat('MT_WI_Or_BasePreRew');
gg=20;
ggg=0;
for k=1:length(PvTable)
    if~isempty(PvTable{k})
        for i=1:length(PvTable{k})
            if ~isempty(PvTable{k}{i})
                 if ~isempty(PvTable{k}{i}{2})
                     Neu1=mat2str(PvTable{k}{i}{2}(1,7));
                     Neu2=mat2str(PvTable{k}{i}{2}(1,8));
                     Neurons=strcat('vs =',Neu1,'vta=',Neu2);
%                   if isempty(PvTable{k}{i}{1}) && isempty(PvTable{k}{i}{2}) && isempty(PvTable{k}{i}{3}) && isempty(PvTable{k}{i}{4}) && isempty(PvTable{k}{i}{5}) && isempty(PvTable{k}{i}{6}) && isempty(PvTable{k}{i}{7})&& isempty(PvTable{k}{i}{8})&& isempty(PvTable{k}{i}{9})&& isempty(PvTable{k}{i}{10})&& isempty(PvTable{k}{i}{11})&& isempty(PvTable{k}{i}{12})&& isempty(PvTable{k}{i}{13})
                gg=gg+1;
                ggg=ggg+1;
                disp(k) 
                disp(i)
                figure(gg)
                
subplot(1,4,1);
p1=pcolor(Base{k}{i});hold on;
set(p1,'EdgeColor','none')
title('Baseline')
xlabel('Time(sec)','FontSize',12)
ylabel('Trial','FontSize',12)
xticks(ticks);
xticklabels(Xlabels1)
set(gca,'FontSize',12)
% text(-4,size(Base{k}{i},1)/2,Neurons,'Fontsize',12,'Rotation',90);
subplot(1,4,2);
p2=pcolor(AftStim{k}{i}); hold on;
set(p2,'EdgeColor','none')
title('PostStimulus')
xlabel('Time(sec)','FontSize',12)
xticks(ticks);
xticklabels(Xlabels2)
set(gca,'FontSize',12)
subplot(1,4,3);
p3=pcolor(PreRew{k}{i}); hold on;
set(p3,'EdgeColor','none')
title('PreRew')
xlabel('Time(sec)','FontSize',12)
xticks(ticks);
xticklabels(Xlabels3)
set(gca,'FontSize',12)
subplot(1,4,4);
p4=pcolor(AftRew{k}{i}); hold on;
set(p4,'EdgeColor','none')
title('PostRew')
xlabel('Time(sec)','FontSize',12)
xticks(ticks);
xticklabels(Xlabels4)
set(gca,'FontSize',12)
filenum=mat2str(ggg);
filename=strcat(fileT,filenum,'.fig');
filename1=strcat(fileT,filenum);
savefig(fullfile(fpath,filename))
saveas(gcf,fullfile(fpath,filename1),'jpg')

                end
            end
        end
    end
end
%% Matrix of count
% pmulti=pmulFM;
for k=1:length(pmulti)
    
    Base_Mtx{k}=[];
    AftStim_Mtx{k}=[];
    PreRew_Mtx{k}=[];
    AftRew_Mtx{k}=[];
    
end
for k=1:length(pmulti)
    for i=1:length(pmulti{k})
    Base_Mtx{k}{i}=[];
    AftStim_Mtx{k}{i}=[];
    PreRew_Mtx{k}{i}=[];
    AftRew_Mtx{k}{i}=[];
    end
end
neu1=1;neu2=4;
for k=1:length(pmulti)
    if ~isempty(pmulti{k})
        for i=1:length(pmulti{k})
            if ~isempty(pmulti{k}{i})
                if pmulti{k}{i}(1,9)==neu1 && pmulti{k}{i}(1,10)==neu2
               Base_Mtx{k}{i}=Base{k}{i};
               AftStim_Mtx{k}{i}=AftStim{k}{i};
               PreRew_Mtx{k}{i}=PreRew{k}{i};
               AftRew_Mtx{k}{i}=AftRew{k}{i};
                end
            end
        end
    end
end




Copy_PvTable_All=PvTable{1};
Copy_Base_All=Base_Mtx{1};
Copy_AftStim_All=AftStim_Mtx{1};
Copy_PreRew_All=PreRew_Mtx{1};
Copy_AftRew_All=AftRew_Mtx{1};
for k=1:length(PvTable)-1
    Copy_PvTable_All=[Copy_PvTable_All,PvTable{k+1}];
    Copy_Base_All =[Copy_Base_All,Base_Mtx{k+1}];
    Copy_AftStim_All=[Copy_AftStim_All,AftStim_Mtx{k+1}];
    Copy_PreRew_All=[Copy_PreRew_All,PreRew_Mtx{k+1}];
    Copy_AftRew_All =[Copy_AftRew_All,AftRew_Mtx{k+1}];
    
end
% Copy_PvTable_All=PvTable_All;
% Copy_Base_All=Base_All;
% Copy_AftStim_All=AftStim_All;
% Copy_PreRew_All=PreRew_All;
% Copy_AftRew_All=AftRew_All;


% clear PvTable_All Base_All AftStim_All
hh=0;
for i=1:length(Copy_PvTable_All)
    if ~isempty(Copy_PvTable_All{i})
        hh=hh+1;
        PvTable_All{hh}=Copy_PvTable_All{i};
        Base_All{hh}=Copy_Base_All{i};
        AftStim_All{hh}=Copy_AftStim_All{i};
        PreRew_All{hh}=Copy_PreRew_All{i};
        AftRew_All{hh}=Copy_AftRew_All{i};
    end
end
%%
save SignAsForMatr_MT_Rev_Late.mat PvTable pmulFM Base_Mtx AftStim_Mtx PreRew_Mtx AftRew_Mtx 


DiffMatx=nan(length(PvTable_All),42);
for i=1:length(PvTable_All)
    if~isempty(PvTable_All{i}{1})
        DiffMatx(i,1:6)=DiffMtx(Base_All{i},AftStim_All{i});
    elseif~isempty(PvTable_All{i}{2})
        DiffMatx(i,7:12)=DiffMtx(Base_All{i},PreRew_All{i});
    elseif~isempty(PvTable_All{i}{3})
        DiffMatx(i,13:18)=DiffMtx(Base_All{i},AftRew_All{i});
    elseif~isempty(PvTable_All{i}{4})
        DiffMatx(i,19:21)=DiffMtx(Base_All{i},AftStim_All{i});
        DiffMatx(i,22:24)=DiffMtx(Base_All{i},PreRew_All{i});
    elseif~isempty(PvTable_All{i}{5})
        DiffMatx(i,25:27)=DiffMtx(Base_All{i},AftStim_All{i});
        DiffMatx(i,28:30)=DiffMtx(Base_All{i},AftRew_All{i});
    elseif~isempty(PvTable_All{i}{6})
        DiffMatx(i,31:33)=DiffMtx(Base_All{i},PreRew_All{i});
        DiffMatx(i,34:36)=DiffMtx(Base_All{i},AftRew_All{i});
        elseif~isempty(PvTable_All{i}{7})
        DiffMatx(i,37:38)=DiffMtx(Base_All{i},AftStim_All{i});
        DiffMatx(i,39:40)=DiffMtx(Base_All{i},PreRew_All{i});
        DiffMatx(i,41:42)=DiffMtx(Base_All{i},AftRew_All{i});
    end
end

%%%
% % yellow	[1,1,0]
% % magenta	[1,0,1]
% % cyan	[0,1,1]
% % red	    [1,0,0]
% % green	[0,1,0]
% % blue	[0,0,1]
% % white	[1,1,1]
% % black	[0,0,0]

DiffMatx(:,end+1)=nan;
DiffPlot= -DiffMatx;
% DiffMatx(end+1,:)=nan;
map=[1, 0, 0;0,0,1];
TitleString = strcat('FSI Type I --- Significative Assemblies Reversal--- Baseline vs: ');
Text_T = {'PostStim';'PreRew';'PostRew';'PostStim/PreRew';'PostStim/PostRew';'PreRew/PostRew';'PostStim/PreRew/PostRew'}
figure(8);hold on
p=pcolor(DiffPlot);
colormap(map)
set(p,'EdgeColor','none');hold on;
for i=1:size(DiffMatx,1)
plot([0 size(DiffMatx,2)],[i,i],'k-');hold on;
end
jj=[0,6,12,18,24,30,36,42];
for i=1:length(jj)
plot([jj(i)+1 jj(i)+1],[0,size(DiffMatx,2)],'k-');hold on;
end
xtick=4:6:40;
xticklab=1:7;
for i=1:length(xtick)
text(xtick(i)-3,size(DiffMatx,1)+1,Text_T(i),'Fontsize',14);hold on;
end
xlim([1,size(DiffMatx,2)]);
ylim([1,size(DiffMatx,1)]);
xticks(xtick);
xticklabels(xticklab)
title(TitleString)
set(gca,'Fontsize',12)

% for k=1:SumTotAn
%     if ~isempty(AsForMean_Base{k})
%     for i=1: length(AsForMean_Base{k})
% % MBA{k}{i}= max(AsForMean_Base{k}{i},[],2);
% % MAS{k}{i}= max(AsForMean_AftStim{k}{i},[],2);
% % MPR{k}{i}= max(AsForMean_PreRew{k}{i},[],2);
% % MAR{k}{i}= max(AsForMean_AftRew{k}{i},[],2);
% % minBAS{k}{i}= min([size(AsForMean_AftStim{k}{i},1),size(Base{k}{i},1)]);
% % minBPR{k}{i}= min([size(AsForMean_PreRew{k}{i},1),size(Base{k}{i},1)]);
% % minBAR{k}{i}= min([size(AsForMean_AftRew{k}{i},1),size(Base{k}{i},1)]);
% %  Base{k}{i} =AsForMean_Base{k}{i};
% % AftStim{k}{i} =AsForMean_AftStim{k}{i}-mean(Base{k}{i},'omitnan');
% % PreRew{k}{i} =AsForMean_PreRew{k}{i}-mean(Base{k}{i},'omitnan');
% % AftRew{k}{i} =AsForMean_AftRew{k}{i}-mean(Base{k}{i},'omitnan');
%     AftStim{k}{i} =AftStim{k}{i}-mean(Base{k}{i},'omitnan');
% PreRew{k}{i} =PreRew{k}{i}-mean(Base{k}{i},'omitnan');
% AftRew{k}{i} =AftRew{k}{i}-mean(Base{k}{i},'omitnan');
%     
% %     Base{k}{i}=Base{k}{i}';
% %     Base{k}{i}(end+1,:)=nan;
% %     AftStim{k}{i}=AftStim{k}{i}'
% %     AftStim{k}{i}(end+1,:)=nan;
% %     PreRew{k}{i}=PreRew{k}{i}'
% %     PreRew{k}{i}(end+1,:)=nan;
% %     AftRew{k}{i}=AftRew{k}{i}';
% %     AftRew{k}{i}(end+1,:)=nan;
%     end
%     end
% 
% end
% 

% for k=1:sumTotAn
%     if ~isempty(Base{k})
%         for i =1:length(Base{k})
%          EpochMatrix{k}{i} = [Base{k}{i}]   
%             
%         end
%         
%     end
% end
%          end
%      end
%      temp=temp+length(Base{k});
%  end
