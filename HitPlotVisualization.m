%% Sorting and text for label y



for i=1:size(pairs_Mat,1)
    l=Sort_O1(i);
    TextIndO{i}=strcat('vs = ',mat2str(pairs_Mat(l,7)),' vta = ', mat2str(pairs_Mat(l,8)))
    j=Sort_R1(i);
    TextIndR{i}=strcat('vs = ',mat2str(pairs_Mat(j,7)),' vta = ', mat2str(pairs_Mat(j,8)))
    z=Sort_L1(i);
    TextIndL{i}=strcat('vs = ',mat2str(pairs_Mat(z,7)),' vta = ', mat2str(pairs_Mat(z,8)))
%     mat2str((pairs(j,7)))
end
clear j i z l
%% Interval
gvect0=-50:50;
gvect=-51:49;
gvect1= -51:124;
tick=[find(gvect0==-25),find(gvect0==0),find(gvect0==25),...
      size(Act_Lick,2)+size(NoVect,2)+find(gvect==-25),...
      size(Act_Lick,2)+size(NoVect,2)+find(gvect==0),...
      size(Act_Lick,2)+size(NoVect,2)+find(gvect==25)...
      2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect1==-25),...
      2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect1==0),...
      2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect1==25),...
      2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect1==50),...
      2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect1==75),...
      2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect1==100)];%,...
%       2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect1==250),...
%       2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect1==300),...
%       2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect1==350),...
%       2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect1==400)];
% ticklabel =[-0.5,0,0.5,-0.5,0,0.5,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4];
ticklabel =[-0.25,0,0.25,-0.25,0,0.25,-0.25,0,0.25,0.5,0.75,1];
tickY = 1: size(Act1,1);
ticklabelYO=TextIndO;
ticklabelYR=TextIndR;
ticklabelYL=TextIndL;


%% Three plot
%  addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla')

switch Dir
    case 'small'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
 case 'smallSync'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
    case 'large'
BinDirText=strcat(Dir,'--BinSize Range: [0.35, 0.6]');
  case 'largeSync'
BinDirText=strcat(Dir,'--BinSize Range: [0.35, 0.6]');
    case 'vs->vta'
       switch binwidthType
           case 'small'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
           case 'large'
BinDirText=strcat(Dir,'--BinSize Range: [0.35, 0.6]');
           case 'all'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.6]');
       end
    case 'vta->vs'
        switch binwidthType
           case 'small'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
           case 'large'
BinDirText=strcat(Dir,'--BinSize Range: [0.35, 0.6]');
            case 'all'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.6]');
        end
    otherwise
        BinDirText=strcat('Region--',Region,'--All Bin');
end


switch TrialsT
    case 'Hit'
        TrialsT = strcat('Hit');
    case 'Rej'
        TrialsT = strcat('Rej');
    case 'False Alarm'
        TrialsT =strcat('False Alarm');
    case 'Miss'
        TrialsT =strcat('Miss');
end


if (ismember(start,startOr) & ismember(stop,stopOr))
    PhaseStr = strcat('Original');
elseif (ismember(start,startRev) & ismember(stop,stopRev))
    PhaseStr = strcat('Reversal');
elseif (ismember(start,startOr) & ismember(stop,stopRev))
    PhaseStr = strcat('Or and Rev');   
end

if pairs_Mat(:,9)==1 & pairs_Mat(:,10)==4
   Couple = strcat('MSN +Type I');
elseif pairs_Mat(:,9)==1 & pairs_Mat(:,10)==5
   Couple = strcat('MSN +Type II');
elseif pairs_Mat(:,9)==1 & pairs_Mat(:,10)==6
   Couple = strcat('MSN +Type III');
elseif pairs_Mat(:,9)==1 & pairs_Mat(:,10)==20
   Couple = strcat('MSN +NoId');
elseif pairs_Mat(:,9)==2 & pairs_Mat(:,10)==4
   Couple = strcat('FSI +Type I');
elseif pairs_Mat(:,9)==2 & pairs_Mat(:,10)==5
   Couple = strcat('FSI +Type II');
elseif pairs_Mat(:,9)==2 & pairs_Mat(:,10)==6
   Couple = strcat('FSI +Type III');
elseif pairs_Mat(:,9)==2 & pairs_Mat(:,10)==20
   Couple = strcat('FSI +NoId');
else Couple = strcat('All');
end
Act1(end+1,:)=nan;
Act2(end+1,:)=nan;
Act3(end+1,:)=nan;
Act4(end+1,:)=nan;
figure
 set(gcf,'numbertitle','off','name',[ 'paradigm: ',nameT ,'bin:',BinDirText,'couple:', Couple,'phase:',PhaseStr]);hold on;
subplot(3,1,1)
pp=pcolor(Act1);hold on;
xlim([1,size(Act1,2)])
ylim([1,size(Act1,1)])
xticks(tick)
xticklabels(ticklabel)
yticks(tickY)
yticklabels(TextIndO)

% text(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-10),-18,'Time(sec)')
% text(find(gvect==-10),-18,'Time(sec)')
% text(2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==-10),-18,'Time(sec)')
text(find(gvect==-10)-3,size(Act_Lick,1)+1,'Odour Onset','Fontsize',12)
% text(size(Act_Lick,2)+size(NoVect,2)+find(gvect==-10)-3,size(Act_Lick,1)+1,'First Lick LWO','Fontsize',12)
text(size(Act_Lick,2)+size(NoVect,2)+find(gvect==-10),size(Act_Lick,1)+1,'First abs Lick','Fontsize',12)
% text(size(Act_Lick,2)+size(NoVect,2)+find(gvect==-10),size(Act_Lick,1)+1,'Second abs Lick','Fontsize',12)
text(2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect==-10)-3,size(Act_Lick,1)+1,'Reward','Fontsize',12)
plot([size(Act_Lick,2)+size(NoVect,2)+find(gvect1==0),size(Act_Lick,2)+size(NoVect,2)+find(gvect1==0)],[0,size(Act2,1)],'-.r','Linewidth',1);
plot([2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect1==0),2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect1==0)],[0,size(Act2,1)],'-.r','Linewidth',1);
plot([find(gvect==0),find(gvect==0)],[0,size(Act2,1)],'-.r','Linewidth',1);
% text(find(gvect==-70)-10,size(ActH_Lick,1)/2,'Sorted on Odour','Rotation',90)
% vline(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==0),'-.r');
% vline(2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==0),'-.r');
% vline(find(gvect==0),'-.r');
set(pp, 'EdgeColor', 'none');
grid on;
ylabel('Pairs Actvity')
% title({TrialsT,'- Odour Sorting -'})
text(2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect==-10)+5,size(Act_Lick,1)+1,[TrialsT,'- Odour Sorting -'],'Fontsize',11);

subplot(3,1,2)
pp=pcolor(Act2);hold on;
xlim([1,size(Act2,2)])
ylim([1,size(Act2,1)])
xticks(tick)
xticklabels(ticklabel)
yticks(tickY)
yticklabels(TextIndL)
% text(find(gvect==-10),-18,'Time(sec)')
% text(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==-10),-18,'Time(sec)')
% text(2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==-10),-18,'Time(sec)')
text(find(gvect==-10),size(Act_Lick,1)+1,'Odour Onset','Fontsize',12)
% text(size(Act_Lick,2)+size(NoVect,2)+find(gvect==-10),size(Act_Lick,1)+1,'First Lick LWO','Fontsize',12)
% text(size(Act_Lick,2)+size(NoVect,2)+find(gvect==-10),size(Act_Lick,1)+1,'Second abs Lick','Fontsize',12)
text(size(Act_Lick,2)+size(NoVect,2)+find(gvect==-10),size(Act_Lick,1)+1,'First abs Lick','Fontsize',12)
text(2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect==-10),size(Act_Lick,1)+1,'Reward','Fontsize',12)
plot([size(Act_Lick,2)+size(NoVect,2)+find(gvect1==0),size(Act_Lick,2)+size(NoVect,2)+find(gvect1==0)],[0,size(Act2,1)],'-.r','Linewidth',1);
plot([2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect1==0),2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect1==0)],[0,size(Act2,1)],'-.r','Linewidth',1);
plot([find(gvect==0),find(gvect==0)],[0,size(Act2,1)],'-.r','Linewidth',1);
% vline(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==0),'-.r');
% vline(2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==0),'-.r');
% vline(find(gvect==0),'-.r');
set(pp, 'EdgeColor', 'none');
grid on;
ylabel('Pairs Actvity')
text(2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect==-10)+5,size(Act_Lick,1)+1,[TrialsT,'- First Lick Sorting -'],'Fontsize',11);
% title({TrialsT,'- First Lick Sorting -'})
text(size(Act1,2)+20,size(Act_Lick,1),{nameT, BinDirText,PhaseStr,Couple},'rotation',270,'Fontsize',12);

subplot(3,1,3)
pp=pcolor(Act3);hold on;
xlim([1,size(Act3,2)])
ylim([1,size(Act3,1)])
xticks(tick)
xticklabels(ticklabel)
yticks(tickY)
yticklabels(TextIndR)
text(find(gvect==-10),-1,'Time(sec)','Fontsize',12)
text(size(Act_Lick,2)+size(NoVect,2)+find(gvect==-10),-1,'Time(sec)','Fontsize',12)
text(2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect==-10),-1,'Time(sec)','Fontsize',12)
text(find(gvect==-10),size(Act_Lick,1)+1,'Odour Onset','Fontsize',12)
% text(size(Act_Lick,2)+size(NoVect,2)+find(gvect==-10),size(Act_Lick,1)+1,'First Lick LWO','Fontsize',12)
% text(size(Act_Lick,2)+size(NoVect,2)+find(gvect==-10),size(Act_Lick,1)+1,'Second abs Lick','Fontsize',12)
text(size(Act_Lick,2)+size(NoVect,2)+find(gvect==-10),size(Act_Lick,1)+1,'First abs Lick','Fontsize',12)
text(2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect==-10),size(Act_Lick,1)+1,'Reward','Fontsize',12)
plot([size(Act_Lick,2)+size(NoVect,2)+find(gvect1==0),size(Act_Lick,2)+size(NoVect,2)+find(gvect1==0)],[0,size(Act2,1)],'-.r','Linewidth',1);
plot([2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect1==0),2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect1==0)],[0,size(Act2,1)],'-.r','Linewidth',1);
plot([find(gvect==0),find(gvect==0)],[0,size(Act2,1)],'-.r','Linewidth',1);
% text(find(gvect==-70)-10,size(ActH_Lick,1)/2,'Sorted on Odour','Rotation',90)
% vline(size(ActH_Lick,2)+size(NoVect,2)+find(gvect==0),'-.r');
% vline(2*size(ActH_Lick,2)+2*size(NoVect,2)+find(gvect==0),'-.r');
% vline(find(gvect==0),'-.r');
set(pp, 'EdgeColor', 'none');
grid on;
ylabel('Pairs Actvity')
text(2*size(Act_Lick,2)+2*size(NoVect,2)+find(gvect==-10)+5,size(Act_Lick,1)+1,[TrialsT,'- Reward Sorting -'],'Fontsize',11);

%title({TrialsT,'- Reward Sorting -'})