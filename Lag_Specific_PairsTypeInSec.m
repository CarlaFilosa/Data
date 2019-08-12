%  MergeDataGenPostPru
%  MergeDataGenPostPru_2 %postpruned 
% MergeDataGenPru  %pruned
%% %%%%%%%%%%%%%%%%%%% lag in seconds for general typologies
copy_pair_small=pairs_small;
clear pairs_small;
pairs_small=new_pairs_small;
[~,~,~,LagInSec_small,LagInSecNoZero_small,LagInSecZero_small] = LagCountTypeMatPairs(pairs_small,BinSizes);
[~,~,~,LagInSec_large,LagInSecNoZero_large,LagInSecZero_large] = LagCountTypeMatPairs(pairs_large,BinSizes);

x_hlag{1} = -1.8:0.12:1.8;
x_hlag_zero{1} = -0.06:0.12:0.06;
x_hlag{2} = -3:0.25:3;
x_hlag_zero{2} = -0.125:0.25:0.125;

% x_hlag{1} = -1:0.1:1;
% x_hlag_zero{1} = -0.05:0.1:0.05;
% x_hlag{2} = -3.6:0.6:3.6;
% x_hlag_zero{2} = -0.3:0.6:0.3;
% 
%%%%%%% Lag Indices (NonZeros and Zeros) 
LagIndNoZero_small= find(LagInSec_small(:,1)~=0);  % all small
LagIndZero_small= find(LagInSec_small(:,1)==0);

LagIndNoZero_large= find(LagInSec_large(:,1)~=0);  % all large
LagIndZero_large= find(LagInSec_large(:,1)==0);


% LagInSec_All= [LagInSec_small,LagInSec_large];
% LagIndNoZero_All = [LagIndNoZero_small,LagIndNoZero_large];
% LagIndZero_All = [LagIndZero_small LagIndZero_large];



figure(13); hold on;
subplot(1,2,1)
h_NoZero = histogram(LagInSecNoZero_small,x_hlag{1},'FaceAlpha',0.5); hold on;
h_Zero = histogram(LagInSecZero_small,x_hlag_zero{1},'FaceAlpha',0.5); hold on;
h_NoZero.FaceColor ='b';
h_Zero.FaceColor = 'g';
title(strcat('\Delta \in [0.01; 0.25]'));hold on;
xlabel('Lag (sec)');hold on;
ylabel('# Pairs');hold on;
% xticks(tick);hold on;
% xticklabels(tick);hold on;
%ylim([0 13])
%text(0.8,10,'\Delta \in [0.01, 0.25] sec')

subplot(1,2,2)
h_NoZero = histogram(LagInSecNoZero_large,x_hlag{2},'FaceAlpha',0.5); hold on;
h_Zero = histogram(LagInSecZero_large,x_hlag_zero{2},'FaceAlpha',0.5); hold on;
h_NoZero.FaceColor ='b';
h_Zero.FaceColor = 'g';
title(strcat('\Delta \in [0.35; 0.6]'));hold on;
xlabel('Lag (sec)');hold on;
ylabel('# Pairs');hold on;
% xticks(tick);hold on;
% xticklabels(tick);hold on;
%ylim([0 13])
%text(0.8,10,'\Delta \in [0.01, 0.25] sec')

clear LagIndNoZero_small LagIndZero_small LagIndNoZero_large LagIndZero_large
%% %%%%%%%%%%%%%%%% Lag in seconds for putative pairs typologies

copy_pair_small=pairs_small;
clear pairs_small;
pairs_small=newnew_pairs_small;
%%
c1=1;% MSN 
c2=2; %FSI 
c3=3; %CIN
for i =1:3 
    cvta=i+3;  % Type I, Type, Type III
    
%%%%% MSN with all identified VTA types
[~,~,~,LagInSec_small_MT{i},LagInSecNoZero_small_MT{i},LagInSecZero_small_MT{i}] = LagCountTypeMatPairs(pairs_small,BinSizes,c1,cvta);
% [~,~,~,LagInSec_large_MT{i},LagInSecNoZero_large_MT{i},LagInSecZero_large_MT{i}] = LagCountTypeMatPairs(pairs_large,BinSizes,c1,cvta);

%%%%%%%% FSI with all identified vta types
[~,~,~,LagInSec_small_FT{i},LagInSecNoZero_small_FT{i},LagInSecZero_small_FT{i}] = LagCountTypeMatPairs(pairs_small,BinSizes,c2,cvta);
% [~,~,~,LagInSec_large_FT{i},LagInSecNoZero_large_FT{i},LagInSecZero_large_FT{i}] = LagCountTypeMatPairs(pairs_large,BinSizes,c2,cvta);

%%%%%%%% CIN with all identified vta types
[~,~,~,LagInSec_small_CT{i},LagInSecNoZero_small_CT{i},LagInSecZero_small_CT{i}] = LagCountTypeMatPairs(pairs_small,BinSizes,c3,cvta);
% [~,~,~,LagInSec_large_CT{i},LagInSecNoZero_large_CT{i},LagInSecZero_large_CT{i}] = LagCountTypeMatPairs(pairs_large,BinSizes,c3,cvta);
end

cvta =20;


%%%%% MSN with no class vta
[~,~,~,LagInSec_small_MT{4},LagInSecNoZero_small_MT{4},LagInSecZero_small_MT{4}] = LagCountTypeMatPairs(pairs_small,BinSizes,c1,cvta);
% [~,~,~,LagInSec_large_MT{4},LagInSecNoZero_large_MT{4},LagInSecZero_large_MT{4}] = LagCountTypeMatPairs(pairs_large,BinSizes,c1,cvta);

%%%%%%%% FSI with no class vta
[~,~,~,LagInSec_small_FT{4},LagInSecNoZero_small_FT{4},LagInSecZero_small_FT{4}] = LagCountTypeMatPairs(pairs_small,BinSizes,c2,cvta);
% [~,~,~,LagInSec_large_FT{4},LagInSecNoZero_large_FT{4},LagInSecZero_large_FT{4}] = LagCountTypeMatPairs(pairs_large,BinSizes,c2,cvta);

%%%%%%%% CIN with no class vta
[~,~,~,LagInSec_small_CT{4},LagInSecNoZero_small_CT{4},LagInSecZero_small_CT{4}] = LagCountTypeMatPairs(pairs_small,BinSizes,c3,cvta);
% [~,~,~,LagInSec_large_CT{4},LagInSecNoZero_large_CT{4},LagInSecZero_large_CT{4}] = LagCountTypeMatPairs(pairs_large,BinSizes,c3,cvta);

%%
% struct_pair = [struct_pair_small,struct_pair_large];
c1=1;% MSN 
c2=2; %FSI 
c3=3; %CIN
for i =1:3  %%%% Type of neuorns on the column
  
cvta=i+3;  % Type I, Type, Type III
[LagInSec_small_MT{i}] = LagCountLabel(struct_pair_small,c1,cvta);  %MSN small with typeI,TypeII,TypeIII
[LagInSec_large_MT{i}] = LagCountLabel(struct_pair_large,c1,cvta);  % MSN large with t1, t2, t3

 [LagInSec_small_FT{i}] = LagCountLabel(struct_pair_small,c2,cvta);  %FSi small with t1, t2, t3
 [LagInSec_large_FT{i}] = LagCountLabel(struct_pair_large,c2,cvta);  %FSI large with t1,t2,t3

 [LagInSec_small_CT{i}] = LagCountLabel(struct_pair_small,c3,cvta);  %CIN small with t1,t2,t3
 [LagInSec_large_CT{i}] = LagCountLabel(struct_pair_large,c3,cvta);  % CIN large with t1 t2 t3
end
cvta=20;
[LagInSec_small_MT{4}] = LagCountLabel(struct_pair_small,c1,cvta);   %MSN small with No id VTA:NoT
[LagInSec_large_MT{4}] = LagCountLabel(struct_pair_large,c1,cvta);   %MSN large with NoT

[LagInSec_small_FT{4}] = LagCountLabel(struct_pair_small,c2,cvta);  % FSI small with No Id VTA
[LagInSec_large_FT{4}] = LagCountLabel(struct_pair_large,c2,cvta);   % FSI large with NoT

[LagInSec_small_CT{4}] = LagCountLabel(struct_pair_small,c3,cvta);  % CIN small with No Id VTA
[LagInSec_large_CT{4}] = LagCountLabel(struct_pair_large,c3,cvta);   % CIN large with noT
% % % 
% % % % hLagNoNorm=cell(1,2);
% % % % x_hlag=cell(1,2);
% % % %  
% % % 
x_hlag{1} = -2.4:0.2:2.4;
x_hlag_zero{1} = -0.1:0.2:0.1;
x_hlag{2} = -3.5:0.2:3.5;
x_hlag_zero{2} = -0.1:0.2:0.1;

% x_hlag{3} = -2:0.05:2;
% x_hlag_zero{3} = -0.025:0.05:0.025;
% x_hlag{4} = -4:0.2:4;
% x_hlag_zero{4} = -0.1:0.2:0.1;
% maxY(i)=max(hhhNoNorm{i}(1,:));

% hhlagNoNorm{1} = histcounts(LagSecTwoScales{1}(:,1),x_hlag{1});
NTypes= length(LagInSec_small_MT);
for i =1: NTypes
    
   %%%%%%% MSN Indices (NonZeros and Zeros) 
if ~isempty(LagInSec_small_MT{i})
LagIndNoZero_small{i}= find(LagInSec_small_MT{i}(:,1)~=0);  % MSN with all small
LagIndZero_small{i}= find(LagInSec_small_MT{i}(:,1)==0);
end
% if ~isempty(LagInSec_large_MT{i})
% LagIndNoZero_large{i}= find(LagInSec_large_MT{i}(:,1)~=0);  % MSN with all large
% LagIndZero_large{i}= find(LagInSec_large_MT{i}(:,1)==0);
% end

   %%%%%% FSI Indices (NonZeros and Zeros)
if ~isempty(LagInSec_small_FT{i})
LagIndNoZero_small{i+length(LagInSec_small_MT)}= find(LagInSec_small_FT{i}(:,1)~=0);% FSI with all small
LagIndZero_small{i+length(LagInSec_small_MT)}= find(LagInSec_small_FT{i}(:,1)==0);
end
% if ~isempty(LagInSec_large_FT{i})
% LagIndNoZero_large{i+NTypes}= find(LagInSec_large_FT{i}(:,1)~=0);  % FSI with all large
% LagIndZero_large{i+NTypes}= find(LagInSec_large_FT{i}(:,1)==0);
% end


%%%%%% CIN Indices (NonZeros and Zeros)

if ~isempty(LagInSec_small_CT{i})
LagIndNoZero_small{i+(2*NTypes)}= find(LagInSec_small_CT{i}(:,1)~=0);  % CIN Indices with all small
LagIndZero_small{i+(2*NTypes)}= find(LagInSec_small_CT{i}(:,1)==0);
end
% if ~isempty(LagInSec_large_CT{i})
% LagIndNoZero_large{i+(2*NTypes)}= find(LagInSec_large_CT{i}(:,1)~=0);  % CIN Indices with all large
% LagIndZero_large{i+(2*NTypes)}= find(LagInSec_large_CT{i}(:,1)==0);
% end

end

%%

LagInSec_All= [LagInSec_small_MT,LagInSec_small_FT,LagInSec_small_CT]%,LagInSec_large_MT,LagInSec_large_FT,LagInSec_large_CT];
LagIndNoZero_All = [LagIndNoZero_small]%,LagIndNoZero_large];
LagIndZero_All = [LagIndZero_small]% LagIndZero_large];

tick=-1.5:0.5:1.5
x_hlag{1} = -1.8:0.1:1.8;
x_hlag_zero{1} = -0.05:0.1:0.05;
x_hlag{2} = -3.5:0.2:3.5;
x_hlag_zero{2} = -0.1:0.2:0.1;


%% %%%%% Lags in Histograms in second for putative pairs types
fignum=[1;2;3;4;5;6;7;8;9;10;11;12]

Title_String = {'MSN-TypeI','MSN-TypeII','MSN-TypeIII','MSN-NoT',...
                'FSI-TypeI','FSI-TypeII','FSI-TypeIII','FSI-NoT',...
                'CIN-TypeI','CIN-TypeII','CIN-TypeIII','CIN-NoT',...
                'MSN-TypeI','MSN-TypeII','MSN-TypeIII','MSN-NoT',...
                'FSI-TypeI','FSI-TypeII','FSI-TypeIII','FSI-NoT',...
                'CIN-TypeI','CIN-TypeII','CIN-TypeIII','CIN-NoT'}

Title_Text = strcat(Title_String);
figure(fignum(2));hold on;
Count=0;
ccc=0;
Ncol=6;
for j =1:Ncol
    
    for i = 0:NTypes-1
if ~isempty(LagInSec_All{ccc+(i+1)})
    subplot(NTypes,Ncol,(j)+(Ncol*i))
    if (j <=3)
    h_NoZero = histogram(LagInSec_All{ccc+(i+1)}(LagIndNoZero_All{ccc+(i+1)},1),x_hlag{1},'FaceAlpha',0.5); hold on;
    h_Zero = histogram(LagInSec_All{ccc+(i+1)}(LagIndZero_All{ccc+(i+1)},1),x_hlag_zero{1},'FaceAlpha',0.5); hold on;
    h_NoZero.FaceColor ='b';
    h_Zero.FaceColor = 'g';
    elseif (j>3)
         h_NoZero = histogram(LagInSec_All{ccc+(i+1)}(LagIndNoZero_All{ccc+(i+1)},1),x_hlag{2},'FaceAlpha',0.5); hold on;
    h_Zero = histogram(LagInSec_All{ccc+(i+1)}(LagIndZero_All{ccc+(i+1)},1),x_hlag_zero{2},'FaceAlpha',0.5); hold on;
    h_NoZero.FaceColor ='b';
    h_Zero.FaceColor = 'g';
    end
    xlabel('Lag (sec)','Fontsize',10)
if (j==1) | (j==4)
    ylabel('# Pairs','Fontsize',10)
end
    if (i~=0) & (j~=2) & (j~=5)
    title(Title_Text{ccc+(i+1)})
    elseif (i==0) & (j~=2) & (j~=5)
    title(Title_Text{ccc+(i+1)})
    elseif (i~=0) & (j==2) 
    title(Title_Text{ccc+(i+1)})
    elseif (i~=0) & (j==5) 
    title(Title_Text{ccc+(i+1)})
    elseif (i==0) & (j==2)
    title({'Small Bins', Title_Text{ccc+(i+1)}})
    elseif (i==0) & (j==5)
        title({'Large Bins', Title_Text{ccc+(i+1)}})
    end 
end

Count=Count+1;
    end
   ccc=Count
end
%% %%%%%%  Only for four types:
tick=-1.5:0.5:1.5
x_hlag{1} = -1.8:0.1:1.8;
x_hlag_zero{1} = -0.05:0.1:0.05;
x_hlag{2} = -3.5:0.2:3.5;
x_hlag_zero{2} = -0.1:0.2:0.1;
figure(fignum(2)); hold on;
subplot(2,2,1)
h_NoZero = histogram(LagInSecNoZero_small_MT{1},x_hlag{1},'FaceAlpha',0.5); hold on;
h_Zero = histogram(LagInSecZero_small_MT{1},x_hlag_zero{1},'FaceAlpha',0.5); hold on;
h_NoZero.FaceColor ='b';
h_Zero.FaceColor = 'g';
title(Title_Text{1});hold on;
xlabel('Lag (sec)');hold on;
ylabel('# Pairs');hold on;
xticks(tick);hold on;
xticklabels(tick);hold on;
ylim([0 13])
text(0.8,10,'\Delta \in [0.01, 0.25] sec')

subplot(2,2,2)
h_NoZero = histogram(LagInSecNoZero_small_FT{1},x_hlag{1},'FaceAlpha',0.5); hold on;
h_Zero = histogram(LagInSecZero_small_FT{1},x_hlag_zero{1},'FaceAlpha',0.5); hold on;
h_NoZero.FaceColor ='b';
h_Zero.FaceColor = 'g';
title(Title_Text{5});hold on;
xlabel('Lag (sec)');hold on;
ylabel('# Pairs');hold on;
xticks(tick);hold on;
xticklabels(tick);hold on;
ylim([0 13])
text(0.8,10,'\Delta \in [0.01, 0.25] sec')

subplot(2,2,3)
h_NoZero = histogram(LagInSecNoZero_small_MT{2},x_hlag{1},'FaceAlpha',0.5); hold on;
h_Zero = histogram(LagInSecZero_small_MT{2},x_hlag_zero{1},'FaceAlpha',0.5); hold on;
h_NoZero.FaceColor ='b';
h_Zero.FaceColor = 'g';
title(Title_Text{2});hold on;
xlabel('Lag (sec)'); hold on;
ylabel('# Pairs');hold on;
xticks(tick);hold on;
xticklabels(tick);hold on;
ylim([0 13])
text(0.8,10,'\Delta \in [0.01, 0.25] sec')
% text(-3, 8, nameT,'Rotation',90,'Fontsize',12)

subplot(2,2,4)
h_NoZero = histogram(LagInSecNoZero_small_FT{2},x_hlag{1},'FaceAlpha',0.5); hold on;
h_Zero = histogram(LagInSecZero_small_FT{2},x_hlag_zero{1},'FaceAlpha',0.5); hold on;
h_NoZero.FaceColor ='b';
h_Zero.FaceColor = 'g';
title(Title_Text{6});hold on;
xlabel('Lag (sec)');hold on;
ylabel('# Pairs');hold on;
xticks(tick);hold on;
xticklabels(tick);hold on;
ylim([0 13])
text(0.8,10,'\Delta \in [0.01, 0.25] sec')
%% %%%% %%%%%

for l=1: length(LagInSec_All)
Lag_All_Pos{l} = LagInSec_All{l}(find(LagInSec_All{l}>0));
Lag_All_Neg{l} = LagInSec_All{l}(find(LagInSec_All{l}<0));
end

for l= 1:length(Lag_All_Pos)
    Mean_Lag_All_Pos(l)=nan;
    Std_Lag_All_Pos(l)=nan;
end

for l= 1:length(Lag_All_Pos)
    if ~isempty(Lag_All_Pos{l})
    Mean_Lag_All_Pos(l) = mean(Lag_All_Pos{l});
    Std_Lag_All_Pos(l) = std(Lag_All_Pos{l})/sqrt(length(Lag_All_Pos{l}));
    end
end

for l= 1:length(Lag_All_Neg)
    Mean_Lag_All_Neg(l)=nan;
    Std_Lag_All_Neg(l)=nan;
end

for l= 1:length(Lag_All_Neg)
    if ~isempty(Lag_All_Neg{l})
    Mean_Lag_All_Neg(l) = mean(Lag_All_Neg{l});
    Std_Lag_All_Neg(l) = std(Lag_All_Neg{l})/sqrt(length(Lag_All_Neg{l}));
    end
end








%% %%%%%%%%  Histogram means lag
ysupM_small = max([max(Mean_Lag_All_Pos(1:12)+Std_Lag_All_Pos(1:12)),...
              max(Mean_Lag_All_Neg(1:12)+Std_Lag_All_Neg(1:12))]);

yinfM_small = min([min(Mean_Lag_All_Pos(1:12)-Std_Lag_All_Pos(1:12)),...
              min(Mean_Lag_All_Neg(1:12)-Std_Lag_All_Neg(1:12))]);

ysupM_large = max([max(Mean_Lag_All_Pos(13:24)+Std_Lag_All_Pos(13:24)),...
              max(Mean_Lag_All_Neg(13:24)+Std_Lag_All_Neg(13:24))]);
   
yinfM_large = min([min(Mean_Lag_All_Pos(13:24)-Std_Lag_All_Pos(13:24)),...
              min(Mean_Lag_All_Neg(13:24)-Std_Lag_All_Pos(13:24))]);
          
          
figure(fignum(3));
hold on;
subplot(2,1,1)
b_plot(1) = barh(1:4,Mean_Lag_All_Pos(1:4));hold on;
b_plot(2) = barh(5:8, Mean_Lag_All_Pos(5:8));hold on;
b_plot(3) = barh(9:12,Mean_Lag_All_Pos(9:12));hold on;
b_plot(4) = barh(1:4,Mean_Lag_All_Neg(1:4));hold on;
b_plot(5) = barh(5:8, Mean_Lag_All_Neg(5:8));hold on;
b_plot(6) = barh(9:12,Mean_Lag_All_Neg(9:12));hold on;
err(1) = errorbar(Mean_Lag_All_Pos(1:4),1:4,Std_Lag_All_Pos(1:4),'horizontal','k*');hold on;
err(2) = errorbar(Mean_Lag_All_Pos(5:8),5:8,Std_Lag_All_Pos(5:8),'horizontal','k*');hold on;
err(3) = errorbar(Mean_Lag_All_Pos(9:12),9:12,Std_Lag_All_Pos(9:12),'horizontal','k*');hold on;
err(4) = errorbar(Mean_Lag_All_Neg(1:4),1:4,Std_Lag_All_Neg(1:4),'horizontal','k*');hold on;
err(5) = errorbar(Mean_Lag_All_Neg(5:8),5:8,Std_Lag_All_Neg(5:8),'horizontal','k*');hold on;
err(6) = errorbar(Mean_Lag_All_Neg(9:12),9:12,Std_Lag_All_Neg(9:12),'horizontal','k*');hold on;
b_plot(1).FaceColor = 'b';
b_plot(2).FaceColor = 'm';
b_plot(3).FaceColor = 'g';
b_plot(4).FaceColor = 'c';
b_plot(5).FaceColor = 'r';
b_plot(6).FaceColor = 'y';
for i=1:length(Mean_Lag_All_Pos)/2
    if ~isnan(Mean_Lag_All_Pos(i))
    text(Mean_Lag_All_Pos(i)+Std_Lag_All_Pos(i)+0.05,i,strcat(Title_String{i}));hold on;
    else text(0.1,i,strcat(Title_String{i}));hold on;
    end
end
xlim([yinfM_small-0.3,ysupM_small+0.6]); hold on;
title('Small Bins');hold on;
xlabel('Mean Lag Delay (sec)');hold on;
ylabel('Pairs Type'); hold on;
xlim([-2.2 2.2]); hold on;
yticks([]);hold on;
box on;
subplot(2,1,2)
b_plot(1) = barh(1:4,Mean_Lag_All_Pos(13:16));hold on;
b_plot(2) = barh(5:8, Mean_Lag_All_Pos(17:20));hold on;
b_plot(3) = barh(9:12,Mean_Lag_All_Pos(21:24));hold on;
b_plot(4) = barh(1:4,Mean_Lag_All_Neg(13:16));hold on;
b_plot(5) = barh(5:8, Mean_Lag_All_Neg(17:20));hold on;
b_plot(6) = barh(9:12,Mean_Lag_All_Neg(21:24));hold on;
err(1) = errorbar(Mean_Lag_All_Pos(13:16),1:4,Std_Lag_All_Pos(13:16),'horizontal','k*');hold on;
err(2) = errorbar(Mean_Lag_All_Pos(17:20),5:8,Std_Lag_All_Pos(17:20),'horizontal','k*');hold on;
err(3) = errorbar(Mean_Lag_All_Pos(21:24),9:12,Std_Lag_All_Pos(21:24),'horizontal','k*');hold on;
err(4) = errorbar(Mean_Lag_All_Neg(13:16),1:4,Std_Lag_All_Neg(13:16),'horizontal','k*');hold on;
err(5) = errorbar(Mean_Lag_All_Neg(17:20),5:8,Std_Lag_All_Neg(17:20),'horizontal','k*');hold on;
err(6) = errorbar(Mean_Lag_All_Neg(21:24),9:12,Std_Lag_All_Neg(21:24),'horizontal','k*');hold on;
b_plot(1).FaceColor = 'b';
b_plot(2).FaceColor = 'm';
b_plot(3).FaceColor = 'g';
b_plot(4).FaceColor = 'c';
b_plot(5).FaceColor = 'r';
b_plot(6).FaceColor = 'y';
for i=1:length(Mean_Lag_All_Pos)/2
    i
    if ~isnan(Mean_Lag_All_Pos(i+length(Mean_Lag_All_Pos)/2))
    text(Mean_Lag_All_Pos(i+length(Mean_Lag_All_Pos)/2)+Std_Lag_All_Pos(i+length(Mean_Lag_All_Pos)/2)+0.05,i,strcat(Title_String{i}));hold on;
    else text(0.1,i,strcat(Title_String{i}));hold on;
    end
end
xlim([yinfM_large-0.3,ysupM_large+0.6]); hold on;
title('Large Bins');hold on;
xlabel('Mean Lag Delay (sec)');hold on;
ylabel('Pairs Type');hold on;
xlim([-3.6 3.6])
yticks([]);hold on;
box on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Histograms per neurons typologies Only identified typologies (No Id Neurons)
NeuType = 6;
for i =1: NeuType    %%%% Type of neurons on the column
 
       cn=i;  
      [LagInSec_small_NeuType{i}] = LagCountLabel(struct_pair_small,cn);  % typologies of neuros small vs
      [LagInSec_large_NeuType{i}] = LagCountLabel(struct_pair_large,cn);  % typologies of neuros large vs
end


for l = 1:NeuType
    if~isempty(LagInSec_small_NeuType)
    minLagSmall(l) = min(LagInSec_small_NeuType{l});
    maxLagSmall(l) = max(LagInSec_small_NeuType{l});
    
    minLagLarge(l) = min(LagInSec_large_NeuType{l});
    maxLagLarge(l) = max(LagInSec_large_NeuType{l});
    end
end


MINLagSmall = min(minLagSmall);
MAXLagSmall = max(maxLagSmall);

MINLagLarge = min(minLagLarge);
MAXLagLarge = max(maxLagLarge);

h_xlag_NoZero_neutype{1} =  MINLagSmall:0.05:MAXLagSmall;
h_xlag_NoZero_neutype{2} =  MINLagLarge:0.2:MAXLagLarge;

h_xlag_Zero_neutype{1} =  -0.025:0.05:0.025;
h_xlag_Zero_neutype{2} =  -0.1:0.2:0.1;

for l=1:NeuType
LagNoZero_small_NeuType{l} = LagInSec_small_NeuType{l}(find(LagInSec_small_NeuType{l})~=0)
LagZero_small_NeuType{l} = LagInSec_small_NeuType{l}(find(LagInSec_small_NeuType{l})==0)

LagNoZero_large_NeuType{l} = LagInSec_large_NeuType{l}(find(LagInSec_large_NeuType{l})~=0)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear Ncol Nrow
Ncol = 4;
Nrow = 3;
%%
%%% Keeping only the categories of units interest MSN FSI CIN

% % figure(fignum(3));hold on;
% % ccc=0;
% % for j=1:Ncol
% %     for i=1:Nrow
% % subplot(Nrow,Ncol,j+(Ncol*i))
% % if j<=2
% %     h_NoZeroNeuType = histogram(LagInSec_small_NeuType{i+ccc}(find(LagInSec_small_NeuType{i+ccc})~=0),h_xlag_NoZero_neutype{1});
% %     
% % elseif j>2
% % end
% %     end
% %     ccc=Ncol;
% % end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % %% vertical graph
% % % % % % figure(fignum(7));
% % % % % % hold on;
% % % % % % subplot(2,1,1)
% % % % % % b_plot(1) = bar(1:4,Mean_Lag_All_Pos(1:4));hold on;
% % % % % % b_plot(2) = bar(5:8, Mean_Lag_All_Pos(5:8));hold on;
% % % % % % b_plot(3) = bar(9:12,Mean_Lag_All_Pos(9:12));hold on;
% % % % % % b_plot(4) = bar(1:4,Mean_Lag_All_Neg(1:4));hold on;
% % % % % % b_plot(5) = bar(5:8, Mean_Lag_All_Neg(5:8));hold on;
% % % % % % b_plot(6) = bar(9:12,Mean_Lag_All_Neg(9:12));hold on;
% % % % % % err(1) = errorbar(Mean_Lag_All_Pos(1:4),1:4,Std_Lag_All_Pos(1:4),'k*');hold on;
% % % % % % err(2) = errorbar(Mean_Lag_All_Pos(5:8),5:8,Std_Lag_All_Pos(5:8),'k*');hold on;
% % % % % % err(3) = errorbar(Mean_Lag_All_Pos(9:12),9:12,Std_Lag_All_Pos(9:12),'k*');hold on;
% % % % % % err(4) = errorbar(Mean_Lag_All_Neg(1:4),1:4,Std_Lag_All_Neg(1:4),'k*');hold on;
% % % % % % err(5) = errorbar(Mean_Lag_All_Neg(5:8),5:8,Std_Lag_All_Neg(5:8),'k*');hold on;
% % % % % % err(6) = errorbar(Mean_Lag_All_Neg(9:12),9:12,Std_Lag_All_Neg(9:12),'k*');hold on;
% % % % % % b_plot(1).FaceColor = 'b';
% % % % % % b_plot(2).FaceColor = 'm';
% % % % % % b_plot(3).FaceColor = 'g';
% % % % % % b_plot(4).FaceColor = 'c';
% % % % % % b_plot(5).FaceColor = 'r';
% % % % % % b_plot(6).FaceColor = 'y';
% % % % % % for i=1:length(Mean_Lag_All_Pos)/2
% % % % % %     if ~isnan(Mean_Lag_All_Pos(i))
% % % % % %     text(i-0.25,Mean_Lag_All_Pos(i)+Std_Lag_All_Pos(i)+0.05,strcat(Title_String{i}));hold on;
% % % % % %     else text(0.1,i,strcat(Title_String{i}));hold on;
% % % % % %     end
% % % % % % end
% % % % % % 
% % % % % % title('Small Bins');hold on;
% % % % % % ylabel('Mean Lag Delay (sec)');hold on;
% % % % % % xlabel('Pairs Type'); hold on;
% % % % % % ylim([-2.75 2.75]); hold on;
% % % % % % box on;
% % % % % % subplot(2,1,2)
% % % % % % b_plot(1) = bar(1:4,Mean_Lag_All_Pos(13:16));hold on;
% % % % % % b_plot(2) = bar(5:8, Mean_Lag_All_Pos(17:20));hold on;
% % % % % % b_plot(3) = bar(9:12,Mean_Lag_All_Pos(21:24));hold on;
% % % % % % b_plot(4) = bar(1:4,Mean_Lag_All_Neg(13:16));hold on;
% % % % % % b_plot(5) = bar(5:8, Mean_Lag_All_Neg(17:20));hold on;
% % % % % % b_plot(6) = bar(9:12,Mean_Lag_All_Neg(21:24));hold on;
% % % % % % err(1) = errorbar(Mean_Lag_All_Pos(13:16),1:4,Std_Lag_All_Pos(13:16),'k*');hold on;
% % % % % % err(2) = errorbar(Mean_Lag_All_Pos(17:20),5:8,Std_Lag_All_Pos(17:20),'k*');hold on;
% % % % % % err(3) = errorbar(Mean_Lag_All_Pos(21:24),9:12,Std_Lag_All_Pos(21:24),'k*');hold on;
% % % % % % err(4) = errorbar(Mean_Lag_All_Neg(13:16),1:4,Std_Lag_All_Neg(13:16),'k*');hold on;
% % % % % % err(5) = errorbar(Mean_Lag_All_Neg(17:20),5:8,Std_Lag_All_Neg(17:20),'k*');hold on;
% % % % % % err(6) = errorbar(Mean_Lag_All_Neg(21:24),9:12,Std_Lag_All_Neg(21:24),'k*');hold on;
% % % % % % b_plot(1).FaceColor = 'b';
% % % % % % b_plot(2).FaceColor = 'm';
% % % % % % b_plot(3).FaceColor = 'g';
% % % % % % b_plot(4).FaceColor = 'c';
% % % % % % b_plot(5).FaceColor = 'r';
% % % % % % b_plot(6).FaceColor = 'y';
% % % % % % for i=1:length(Mean_Lag_All_Pos)/2
% % % % % %     i
% % % % % %     if ~isnan(Mean_Lag_All_Pos(i+length(Mean_Lag_All_Pos)/2))
% % % % % %     text(i-0.25,Mean_Lag_All_Pos(i+length(Mean_Lag_All_Pos)/2)+Std_Lag_All_Pos(i+length(Mean_Lag_All_Pos)/2)+0.05,strcat(Title_String{i}));hold on;
% % % % % %     else text(i-0.25,0.2,strcat(Title_String{i}));hold on;
% % % % % %     end
% % % % % % end
% % % % % % 
% % % % % % title('Large Bins');hold on;
% % % % % % ylabel('Mean Lag Delay (sec)');hold on;
% % % % % % xlabel('Pairs Type');hold on;
% % % % % % ylim([-3.6 3.6])
% % % % % % box on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Only two Types

% % % % % figure(37);hold on;
% % % % % 
% % % % %  subplot(2,2,1);hold on;
% % % % %  h_smallNozero = histogram(LagInSec_small_MT1(LagNoZeroInd{1},1),x_hlag{1},'FaceAlpha',0.5); hold on;
% % % % %  h_smallZero = histogram(LagInSec_small_MT1(LagZeroInd{1},1),x_hlag_zero{1},'FaceAlpha',0.5); hold on;
% % % % %  h_smallNozero.FaceColor='b'; hold on;
% % % % %  h_smallZero.FaceColor = 'g'; hold on;
% % % % %  
% % % % %  xlabel('Lag (sec)','Fontsize',12);hold on; 
% % % % %  ylabel('# pairs','Fontsize',12); hold on;
% % % % % %  ylim([0 12]);hold on;
% % % % %  title('Small Bins M-T1')
% % % % %  box on;
% % % % %  
% % % % %  subplot(2,2,3);hold on;
% % % % %  h_smallNozero2 = histogram(LagInSec_small_FT2(LagNoZeroInd{3},1),x_hlag{3},'FaceAlpha',0.5); hold on;
% % % % %  h_smallZero2 = histogram(LagInSec_small_FT2(LagZeroInd{3},1),x_hlag_zero{3},'FaceAlpha',0.5); hold on;
% % % % %  h_smallNozero2.FaceColor='m'; 
% % % % %  h_smallZero2.FaceColor = 'y';
% % % % %  xlabel('Lag (sec)','Fontsize',12);hold on; 
% % % % %  ylabel('# pairs','Fontsize',12); hold on;
% % % % % %  ylim([0 12]);hold on;
% % % % %  title('Small Bins F-T2')
% % % % %  box on;
% % % % %  
% % % % % 
% % % % %  subplot(2,2,2);hold on;
% % % % %  
% % % % %  h_largeNozero = histogram(LagInSec_large_MT1(LagNoZeroInd{2},1),x_hlag{2},'FaceAlpha',0.5); hold on;
% % % % %  h_largeZero = histogram(LagInSec_large_MT1(LagZeroInd{2},1),x_hlag_zero{2},'FaceAlpha',0.5); hold on;
% % % % %  h_largeNozero.FaceColor='b';
% % % % %  h_largeZero.FaceColor = 'g';
% % % % %  xlabel('Lag (sec)','Fontsize',12);hold on; 
% % % % %  ylabel('# pairs','Fontsize',12); hold on;
% % % % % %  ylim([0 12]);hold on;
% % % % %  title('Large Bins M-T1')
% % % % %  box on;
% % % % %  
% % % % %  subplot(2,2,4);hold on;
% % % % %  h_largeNozero2 = histogram(LagInSec_large_FT2(LagNoZeroInd{4},1),x_hlag{4},'FaceAlpha',0.5); hold on;
% % % % %  h_largeZero2 = histogram(LagInSec_large_FT2(LagZeroInd{4},1),x_hlag_zero{4},'FaceAlpha',0.5); hold on;
% % % % %  h_largeNozero2.FaceColor='m';
% % % % %  h_largeZero2.FaceColor = 'y';
% % % % %  xlabel('Lag (sec)','Fontsize',12);hold on; 
% % % % %  ylabel('# pairs','Fontsize',12); hold on;
% % % % % %  ylim([0 12]);hold on;
% % % % %  title('Large Bins F-T2')
% % % % %  box on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Histograms per neurons typologies with no id

% % % % % % % % for i =1:8  %%%% Type of neuorns on the column
% % % % % % % %     if (i<=3) % Vs neurons (they are assigned to the indices 1,2,3)
% % % % % % % %        cn=i;  
% % % % % % % %       [LagInSec_small_temp{i}] = LagCountLabel(struct_pair_small,cn);  % typologies of neuros small vs
% % % % % % % %       [LagInSec_large_temp{i}] = LagCountLabel(struct_pair_large,cn);  % typologies of neuros large vs
% % % % % % % %       
% % % % % % % %        elseif (i==4)
% % % % % % % %         cn=i+6; %%% No id Vs index is cn=10 and I want to put it in the 4th column
% % % % % % % %         [LagInSec_small_temp{i}] = LagCountLabel(struct_pair_small,cn);  % no id vs neurons small 
% % % % % % % %         [LagInSec_large_temp{i}] = LagCountLabel(struct_pair_large,cn);  % no id vs neusrons large
% % % % % % % %     
% % % % % % % %        elseif (i==8)
% % % % % % % %         cn=i+12; %%% No id Vta index is cn=20 and I want to put it in the 8th column
% % % % % % % %         [LagInSec_small_temp{i}] = LagCountLabel(struct_pair_small,cn);  % no id vta neurons small
% % % % % % % %         [LagInSec_large_temp{i}] = LagCountLabel(struct_pair_large,cn);  % no id vta neusrons large
% % % % % % % %         
% % % % % % % %        elseif (i>4) && (i<8) % Vs neurons (they are assigned to the indices 4, 5,6 and i want to assign them to the column 5,6,7)
% % % % % % % %         cn = i-1;
% % % % % % % %         [LagInSec_small_temp{i}] = LagCountLabel(struct_pair_small,cn);  % typologies of neurons small vta
% % % % % % % %         [LagInSec_large_temp{i}] = LagCountLabel(struct_pair_large,cn);  % typologies of neuorns large vta
% % % % % % % %     end
% % % % % % % % end


%%%%%%% if I want only the indentified typologies
% % % % % % % % for l = 1: 3
% % % % % % % % LagInSec_small_NeuType{l} = LagInSec_small_temp{l};    
% % % % % % % % LagInSec_small_NeuType{l+3} = LagInSec_small_temp{l+4};  
% % % % % % % % LagInSec_large_NeuType{l} = LagInSec_large_temp{l};    
% % % % % % % % LagInSec_large_NeuType{l+3} = LagInSec_large_temp{l+4}; 
% % % % % % % % end
% % % % % % % % 
% % % % % % % % NeuType = length(LagInSec_small_NeuType);
% % % % % % % % for l= 1:NeuType
% % % % % % % % LagIndNoZero_small_NeuType{l} = find(LagInSec_small_NeuType{l}~=0)
% % % % % % % % LagIndZero_small_NeuType{l} = find(LagInSec_small_NeuType{l}==0)
% % % % % % % % LagIndNoZero_large_NeuType{l} = find(LagInSec_large_NeuType{l}~=0)
% % % % % % % % LagIndZero_large_NeuType{l} = find(LagInSec_large_NeuType{l}==0)
% % % % % % % % end