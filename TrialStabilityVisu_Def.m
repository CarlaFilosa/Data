
resumeAftAnFr;
pVal_OrgTyp;
%% PLOT Mean ACtivity ASSEMBLIES IN HEAT PLOT 

ticks=[1 11 21 31 41]; ticks1 = 52+ticks; ticks2 = 52+ticks1; ticks3 = 52+ticks2;
Ticks =[ticks,ticks1,ticks2,ticks3];
Xlabels1 = [-1,-0.9, -0.8, -0.7,-0.6];  %%%%% the other baseline go from -1.0 to -0.5
Xlabels2 = [0, 0.1, 0.2, 0.3, 0.4];
Xlabels3 = [-0.5, -0.4, -0.3, -0.2,-0.1];
Xlabels4 = [0, 0.1, 0.2, 0.3, 0.4];
XLabels = [Xlabels1,Xlabels2,Xlabels3,Xlabels4];
%%%%%%% Chose the pairs that you want to see
ForTypeBase = []; ForTypeAftStim = []; ForTypePreRew = []; ForTypeAftRew = []; ForTypeId = [];
%%%%%%%%%% Phase Choice, Typologie choice
Ph=2;
Typology = 'CT';

if Ph==1
    TextPh = {'Original'};
elseif Ph==2
    TextPh = {'Reversal'};
else disp('no valid input for the phase');
end

MeanFriedTypoSwitch;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Text details Needed for the plots

 for m=1:length(cvta)
     if ~isempty(ForTypeBase{m})
for i=1:length(ForTypeBase{m})
    
%     AA=ForTypeId{m}(ForTypeBase{m}(i),1)
%   TextID_Base{m}{i} = strcat('vs = ', mat2str(ForTypeId{m}(ForTypeBase{m}(i),1)),'vta =', mat2str(ForTypeId{m}(ForTypeBase{m}(i),2)));
  TextID_{1}{m}{i} = strcat('vs = ', mat2str(ForTypeId{m}(ForTypeAftStim{m}(i),1)),'vta =', mat2str(ForTypeId{m}(ForTypeAftStim{m}(i),2)));
  TextID_{2}{m}{i} = strcat('vs = ', mat2str(ForTypeId{m}(ForTypePreRew{m}(i),1)),'vta =', mat2str(ForTypeId{m}(ForTypePreRew{m}(i),2)));
  TextID_{3}{m}{i} = strcat('vs = ', mat2str(ForTypeId{m}(ForTypeAftRew{m}(i),1)),'vta =', mat2str(ForTypeId{m}(ForTypeAftRew{m}(i),2)));
end
     end
 end
clear i k m

 TextOrder = {'AftStim Sorted'; 'PreRew Sorted'; 'AftRew Sorted'};
 
 %%%%%%%%%%%%%%%%%%%
 %%%%%% FILE NAMES AND FILE PATH
 fpath= '/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM/withOpto/HeatPlotNewNorm'
 
 %%%%%%%%%%%%%%%%%%%%%% 
 for m=1:length(cvta)
    for l=1:length(BaseMean)
        TotMean{l}{m}=[];
    end
 end

for m=1:length(cvta)
    for l=1:length(BaseMean)
    BaseMean{l}{m}(:,end+1:end+2) = NaN;
    AftStimMean{l}{m}(:,end+1:end+2) = NaN;
    PreRewMean{l}{m}(:,end+1:end+2) = NaN;
    AftRewMean{l}{m}(:,end+1) = NaN;
    BaseMean{l}{m}(end+1,:) = NaN;
    AftStimMean{l}{m}(end+1,:) = NaN;
    PreRewMean{l}{m}(end+1,:) = NaN;
    AftRewMean{l}{m}(end+1,:) = NaN;
    TotMean{l}{m} = [BaseMean{l}{m},AftStimMean{l}{m},PreRewMean{l}{m},AftRewMean{l}{m}]
        
    end
end
clear l m

Text_Type_VTA = {'Type I','Type II', 'Type III', 'Laser', 'Non Laser', 'Laser-Inh', 'Not Class'};
numr=length(BaseMean);
g=0;
for m = 1:length(cvta)
    g=g+1;
    for l= 1:numr
        YMax = size(TotMean{l}{m},1);
        Xtext1 = round(max(ticks)/3);
        Xtext2 = Xtext1+max(ticks);
        Xtext3 = Xtext1+max(ticks1);
        Xtext4 = Xtext1+max(ticks2);
        
        figure(m);
        subplot(numr,1,l)
        pp =pcolor(TotMean{l}{m});
        set(pp,'EdgeColor','none');
        xticks(Ticks)
        xticklabels(XLabels)
        yticks(1:YMax-1);
        yticklabels(TextID_{l}{m});
%         yt = get(gca, 'YTick');
%         set(gca,'FontSize',10);
        text(Xtext1,YMax+2,'Baseline','Fontsize',12);
        text(Xtext2,YMax+2,'After Stimulus','Fontsize',12);
        text(Xtext3,YMax+2,'Pre Reward','Fontsize',12);
        text(Xtext4,YMax+2,'After Reward','Fontsize',12);
        if l==3
        text(Xtext1,-5,'Time (sec)','Fontsize',12);
        text(Xtext2,-5,'Time (sec)','Fontsize',12);
        text(Xtext3,-5,'Time (sec)','Fontsize',12);
        text(Xtext4,-5,'Time (sec)','Fontsize',12);
        end
        if l==2
            text(max(ticks3)+25,YMax/2,TextPh,'Rotation',270,'Fontsize',12)
        end
        text(max(ticks3)+15, YMax, TextOrder{l},'Rotation',270,'Fontsize',12);
        text(max(ticks3)+20, YMax, [Text_Type_VS,Text_Type_VTA{m}],'Rotation',270,'Fontsize',12);
        
    end
        filenum=mat2str(g);
        filename=strcat(fileT,filenum,'.fig');
        filename1=strcat(fileT,filenum);
        savefig(fullfile(fpath,filename))
end

clear l m
% pp=pcolor(TotMean{3}{1});
% set(pp, 'EdgeColor','none');

%% PLOT TRIAL STABILITY (Assembly per assembly)%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Assembly per assembly plot. It is possible to select for typologies
ticks=[1 10 20 30 40 50]
Xlabels1 = [-1,-0.9,-0.8,-0.7,-0.6,-0.5];  %%%%% the other baseline go from -1.0 to -0.5
Xlabels2 = [0, 0.1, 0.2, 0.3, 0.4, 0.5];
Xlabels3 = [-0.5, -0.4, -0.3, -0.2,-0.1, 0];
Xlabels4 = [0, 0.1, 0.2, 0.3, 0.4, 0.5];
PvTable=PvTableFlowT{1};
Sign=4; Couple =strcat('FSI + Type I --- Aft sign aft Rew');
Base = MeanFried.FlowT.Base{1}.ForMeanNormA{1};
AftStim = MeanFried.FlowT.AftStim{1}.ForMeanNormA{1};
PreRew = MeanFried.FlowT.PreRew{1}.ForMeanNormA{1};
AftRew = MeanFried.FlowT.AftRew{1}.ForMeanNormA{1};
StableID = MeanFried.FlowT.AftRew{1}.StableID{1};
g=3;
for k=1:length(Base)
    if ~isempty(PvTable{k})
        for i=1:length(Base{k})
            if~isempty(PvTable{k}{i})
                if~isempty(PvTable{k}{i}{Sign})
                Neu1 = mat2str(StableID{k}{i}(1));
                Neu2 = mat2str(StableID{k}{i}(2));
                Neurons = strcat('vs = ', Neu1,'----', 'vta = ', Neu2);
                g=g+1;
                figure(g);
                
                subplot(1,4,1)
                pp=pcolor(Base{k}{i}')
                set(pp, 'EdgeColor', 'none');
                grid on
                xticks(ticks)
                xticklabels(Xlabels1)
                xlabel('Time (sec)');
                ylabel('Trials');
                title('Baseline')
                
                subplot(1,4,2)
                pp=pcolor(AftStim{k}{i}')
                set(pp, 'EdgeColor', 'none');
                grid on
                xticks(ticks)
                xticklabels(Xlabels2)
                xlabel('Time (sec)');
                title('After Stimulus')
                
                subplot(1,4,3)
                pp=pcolor(PreRew{k}{i}')
                set(pp, 'EdgeColor', 'none');
                grid on
                xticks(ticks)
                xticklabels(Xlabels3)
                xlabel('Time (sec)');
                title('Pre Reward')
                
                subplot(1,4,4)
                pp=pcolor(AftRew{k}{i}')
                set(pp, 'EdgeColor', 'none');
                grid on
                xticks(ticks);
                xticklabels(Xlabels4);
                xlabel('Time (sec)');
                title('After Reward');
                text(60,size(AftRew{k}{i}',1),{Couple;Neurons},'Rotation',270, 'Fontsize',12);
               
%                 filenum=mat2str(g);
%                 filename=strcat(fileT,filenum,'.fig');
%                 filename1=strcat(fileT,filenum);
%                 savefig(fullfile(fpath,filename))
                end
            end
        end
    end
    
end



%% OLD PLOT>>>> NO!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%
numr=3;numc=4;
TextOrder = {'AftStim Sorted'; 'PreRew Sorted'; 'AftRew Sorted'};
for m=1:length(cvta)
    for l=1:length(BaseMean)
    BaseMean{l}{m} = BaseMean{l}{m};
    AftStimMean{l}{m} = AftStimMean{l}{m};
    PreRewMean{l}{m} = PreRewMean{l}{m};
    AftRewMean{l}{m} = AftRewMean{l}{m};
    
    BaseMean{l}{m}(:,end+1) = NaN;
    AftStimMean{l}{m}(:,end+1) = NaN;
    PreRewMean{l}{m}(:,end+1) = NaN;
    AftRewMean{l}{m}(:,end+1) = NaN;
    figure(m);
    %%%%%%%% Baseline
    subplot(numr,numc,((l-1)*numc)+1);
    pp = pcolor(BaseMean{l}{m});
    set(pp,'EdgeColor','none');
    if l==3
    xticks(ticks)
    xticklabels(Xlabels1)
    else xticks([])
    end
    yticks(1:size(BaseMean{l}{m},1))
    yticklabels(TextID_{l}{m})
    if l==1
    title('Baseline')
    end
    
    %%%%% AfterStimulus
    subplot(numr,numc,((l-1)*numc)+2);
    pp = pcolor(AftStimMean{l}{m});
    set(pp,'EdgeColor','none');
    if l==1
    title('AfterStimulus')
    end
    yticks([])
    if l==3
    xticks(ticks)
    xticklabels(Xlabels2)
    else xticks([])
    end
    
    %%%%% PreReward
    subplot(numr,numc,((l-1)*numc)+3);
    pp = pcolor(PreRewMean{l}{m});
    set(pp,'EdgeColor','none');
    if l==1
    title('PreReward')
    end
    yticks([])
    if l==3
    xticks(ticks)
    xticklabels(Xlabels3)
    else xticks([])
    end
    
    %%%%% AfterReward
    subplot(numr,numc,((l-1)*numc)+4);
    pp = pcolor(AftRewMean{l}{m});
    set(pp,'EdgeColor','none');
    if l==1
    title('AfterReward')
    end
    yticks([])
    if l==3
    xticks(ticks)
    xticklabels(Xlabels4)
    else xticks([])
    end
    text(54, size(AftRewMean{l}{m},1), TextOrder{l},'Rotation',270,'FontSize',11);
    
    end
end



% pcolor(MeanFried.MT.AftStim{1}.MeanMat{1}')


