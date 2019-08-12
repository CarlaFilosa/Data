%%%% Programs that does the lag distribution for each specific couple of
%%%% units

nameMerge{1} = strcat('Rev2_OrRevW.mat');
nameMerge{2} = strcat('Lastrev1_OrRevW.mat');
nameMerge{3} = strcat('Rev3_OrRevW.mat');
nameMerge{4} = strcat('Revconc_OrRevW.mat');
nameT = strcat('Rev2; Lastrev1; Rev3; Revconc');


alg= strcat('Stand_');
algPru= strcat('Pru_');

fignum =[1,2,3,4,5,6,7,8];
IdType='ID';
for i =1: length(nameMerge)

% nameNotPruMerge{i} = strcat('An_Disp_',alg,nameMerge{i})     % Not pruned Algorithm
% load(nameNotPruMerge{i})

namePruMerge{i} = strcat('An_Disp_',algPru,nameMerge{i})     % prepruned Algorithm
load(namePruMerge{i})
clear c
load('classlist.mat');
%  nnn(i,:)=nn(a);
%  aa(i,:)=a;
[new_data] = parId(new_data,a,nn);


 for k=1:length(new_data.par)
 [new_data] = labels(new_data,c,k,IdType);
 end
 
 very_data{i}=new_data;

% As_across_binsT{i} = As_across_bins;  %%Not pruned or postpruned
% As_order_binsT{i} = As_order;         

As_across_binsT{i} = As_acr_bins_pru;     %% Prepruned algorithm
As_order_binsT{i} = As_order_pru;

clear new_data a nn As_across_bins As_order As_acr_bins_pru As_order_pru
end
%%
clearvars -except very_data As_across_binsT As_order_binsT c nameMerge nameT fignum BinSizes MaxLags
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


As_across_bins_All=As_across_binsT{1};
As_order_bins_All=As_order_binsT{1};

nneuVS_All = nneuVS{1};
nneuVTA_All = nneuVTA{1};
Data_All.par=very_data{1}.par;

for i=1:length(very_data)-1
    [As_across_bins_All]=[As_across_bins_All,As_across_binsT{i+1}];
    [As_order_bins_All]=[As_order_bins_All,As_order_binsT{i+1}];
    [nneuVS_All]=[nneuVS_All,nneuVS{i+1}];
    [nneuVTA_All]=[nneuVTA_All,nneuVTA{i+1}];
    [Data_All.par]=[Data_All.par,very_data{i+1}.par];
    
end


[pairs_r,pairs_vsvs,pairs_vsvta,pairs_vtavta] = PairsInfo(As_across_bins_All,As_order_bins_All,nneuVS_All,SumTotAn);

clear a b A B i ii j k
[pairs_vsvta] = PairMatNeuLabelsID(pairs_vsvta,Data_All);

[struct_pairAll] = PairsOrgByPair(pairs_vsvta);

[struct_pairAll] = PairNeuLabelsID(struct_pairAll,Data_All); % Labels of neurons in pairs
clear struct_pair
struct_pair =struct_pairAll;
%%
pairs =pairs_vsvta{1};
for i = 1: length(pairs_vsvta)-1
    if ~isempty(pairs_vsvta{i+1})
pairs=[pairs;pairs_vsvta{i+1}];
    end
end

%%
c1=1;% MSN 
c2=2; %FSI 
c3=3; %CIN
 for i =1:3  %%%% Type of neuorns on the column
  
cvta=i+3;  % 
[LagPerBin_MT{i},LagPerBinNoZero_MT{i},LagPerBinZero_MT{i}] = LagCountTypeMatPairs(pairs,BinSizes,c1,cvta);
[LagPerBin_FT{i},LagPerBinNoZero_FT{i},LagPerBinZero_FT{i}] = LagCountTypeMatPairs(pairs,BinSizes,c2,cvta);
[LagPerBin_CT{i},LagPerBinNoZero_CT{i},LagPerBinZero_CT{i}] = LagCountTypeMatPairs(pairs,BinSizes,c3,cvta);
end

i=4; cvta=20;
[LagPerBin_MT{i},LagPerBinNoZero_MT{i},LagPerBinZero_MT{i}] = LagCountTypeMatPairs(pairs,BinSizes,c1,cvta);
[LagPerBin_FT{i},LagPerBinNoZero_FT{i},LagPerBinZero_FT{i}] = LagCountTypeMatPairs(pairs,BinSizes,c2,cvta);
[LagPerBin_CT{i},LagPerBinNoZero_CT{i},LagPerBinZero_CT{i}] = LagCountTypeMatPairs(pairs,BinSizes,c3,cvta);


%% %%%%%%%%%
% BinSizes_Small = [0.01, 0.015, 0.03, 0.05, 0.08, 0.12, 0.25];
for i = 1: length(BinSizes)
H_x_lag{i} = -MaxLags(i)-0.5:MaxLags(i)+0.5;
end
%% Find no zero indeces
% Ind can assume the values NoZero or Zero depending if I want to find the
% zero or non zero indices
% Ind='NoZero';
% [LagIndNoZero_MT] = FindIndLagHisto(LagPerBin_MT,Ind);
% [LagIndNoZero_FT] = FindIndLagHisto(LagPerBin_FT,Ind);
% [LagIndNoZero_CT] = FindIndLagHisto(LagPerBin_CT,Ind);
% for j= 1:length(LagPerBin_MT)
%     for i =1:length(BinSizes)
% LagPerBin_MT_NoZero{j}{i} = LagPerBin_MT{j}{i}(LagIndNoZero_MT{j}{i});
% LagPerBin_FT_NoZero{j}{i} = LagPerBin_FT{j}{i}(LagIndNoZero_FT{j}{i});
% LagPerBin_CT_NoZero{j}{i} = LagPerBin_CT{j}{i}(LagIndNoZero_CT{j}{i});
%     end
% end
clear BinSizes 
BinSizes =[0.01, 0.015, 0.03, 0.05, 0.08, 0.12,0.25,0.35, 0.5, 0.6]
TextId = {'T1'; 'T2'; 'T3'; 'NoT'}
for i =1:length(TextId)
TextIdMT{i} = strcat('MT-',TextId{i})
end

MAXY=12;
figure(8);hold on;
for i = 1: length(BinSizes)
  
    subplot(2,5,i)
    if ~isempty(LagPerBinNoZero_MT{1}{i})
%      h_LagPerBinMT(1)=histogram(LagPerBin_small_MT{1}{i}(LagIndNoZero_small_MT{1}{i}),H_x_lag{i},'FaceAlpha',0.7);hold on;
    h_LagPerBinMT(1)=histogram(LagPerBinNoZero_MT{1}{i}(:,1),H_x_lag{i},'FaceAlpha',0.7);hold on;
    h_LagPerBinMT(1).FaceColor = 'g'; hold on;
    if i<=10
        text(H_x_lag{i}(2), MAXY-0.4, 'M-T1','Color','g')
    else
        text(H_x_lag{i}(1), MAXY-0.4, 'M-T1','Color','g')
    end
    end
     ylim([0 MAXY])
    ylabel('# Pairs')
    xlabel('Lag')
    box on;
    
    if ~isempty(LagPerBinNoZero_MT{2}{i})
     h_LagPerBinMT(2)=histogram(LagPerBinNoZero_MT{2}{i}(:,1),H_x_lag{i},'FaceAlpha',0.7);hold on;
     h_LagPerBinMT(2).FaceColor = shade([0 1 0],0.5); hold on;
     if i<=10
         text(H_x_lag{i}(2), MAXY-0.8, 'M-T2','Color',shade([0 1 0],0.5))
     else
         text(H_x_lag{i}(1), MAXY-0.8, 'M-T2','Color',shade([0 1 0],0.5))
     end
    end
    ylim([0 MAXY])
    ylabel('# Pairs')
    xlabel('Lag')
    box on;
    
    if ~isempty(LagPerBinNoZero_MT{3}{i})
    h_LagPerBinMT(3)=histogram(LagPerBinNoZero_MT{3}{i}(:,1),H_x_lag{i},'FaceAlpha',0.7);hold on;
    h_LagPerBinMT(3).FaceColor = 'c'; hold on;
    if i<=10
        text(H_x_lag{i}(2), MAXY-1.2, 'M-T3','Color','c')
    else
        text(H_x_lag{i}(1), MAXY-1.2, 'M-T3','Color','c')
    end
    end
     ylim([0 MAXY])
    ylabel('# Pairs')
    xlabel('Lag')
    box on;
    
    if ~isempty(LagPerBinNoZero_MT{4}{i})
    h_LagPerBinMT(4)=histogram(LagPerBinNoZero_MT{4}{i}(:,1),H_x_lag{i},'FaceAlpha',0.7);hold on;
    h_LagPerBinMT(4).FaceColor = shade([0 1 1], 0.5); hold on;
    if i<=10
        text(H_x_lag{i}(2), MAXY-1.6, 'M-NoT','Color',shade([0 1 1], 0.5))
    else
        text(H_x_lag{i}(1), MAXY-1.6, 'M-NoT','Color',shade([0 1 1], 0.5))
    end
     ylim([0 MAXY])
    end

    if ~isempty(LagPerBinNoZero_FT{1}{i})
     h_LagPerBinFT(1)=histogram(LagPerBinNoZero_FT{1}{i}(:,1),H_x_lag{i},'FaceAlpha',0.7);hold on;
    h_LagPerBinFT(1).FaceColor = 'r'; hold on;
    if i<6
       text(H_x_lag{i}(18), MAXY-0.4, 'F-T1','Color','r')
    elseif i==6 || i==7
        text(H_x_lag{i}(8), MAXY-0.4, 'F-T1','Color','r')
    elseif i==8 || i==9 || i==10
        text(H_x_lag{i}(6), MAXY-0.4, 'F-T1','Color','r')
    elseif i==11 || i==12
        text(H_x_lag{i}(3), MAXY-0.4, 'F-T1','Color','r')
    end
    end
    ylim([0 MAXY])
    ylabel('# Pairs')
    xlabel('Lag')
    box on;
    
    if ~isempty(LagPerBinNoZero_FT{2}{i})
     h_LagPerBinFT(2)=histogram(LagPerBinNoZero_FT{2}{i}(:,1),H_x_lag{i},'FaceAlpha',0.7);hold on;
     h_LagPerBinFT(2).FaceColor = shade([1 0 0],0.5); hold on;
%     l(2)=legend(TextIdMT{2},'Location','NorthWest');hold on;
    if i<6
    text(H_x_lag{i}(18), MAXY-0.8, 'F-T2','Color',shade([1 0 0],0.5))
    elseif i==6 || i==7
     text(H_x_lag{i}(8), MAXY-0.8, 'F-T2','Color',shade([1 0 0],0.5))
     elseif i==8 || i==9 || i==10
     text(H_x_lag{i}(6), MAXY-0.8, 'F-T2','Color',shade([1 0 0],0.5))
     elseif i==11 || i==12
     text(H_x_lag{i}(3), MAXY-0.8, 'F-T2','Color',shade([1 0 0],0.5))
    end
    end
    ylim([0 MAXY])
    ylabel('# Pairs')
    xlabel('Lag')
    box on;
    
    if ~isempty(LagPerBinNoZero_FT{3}{i})
    h_LagPerBinFT(3)=histogram(LagPerBinNoZero_FT{3}{i}(:,1),H_x_lag{i},'FaceAlpha',0.7);hold on;
    h_LagPerBinFT(3).FaceColor = 'm'; hold on;
    if i<6
        text(H_x_lag{i}(18), MAXY-1.2, 'F-T3','Color','m')
    elseif i==6 || i==7
        text(H_x_lag{i}(8), MAXY-1.2, 'F-T3','Color','m')
    elseif i== 8 || i==9 || i==10
        text(H_x_lag{i}(6), MAXY-1.2, 'F-T3','Color','m')
    elseif i== 11 || i==12 
        text(H_x_lag{i}(3), MAXY-1.2, 'F-T3','Color','m')
    end
    end
     ylim([0 MAXY])
    ylabel('# Pairs')
    xlabel('Lag')
    
    if ~isempty(LagPerBinNoZero_FT{4}{i})
    h_LagPerBinFT(4)=histogram(LagPerBinNoZero_FT{4}{i}(:,1),H_x_lag{i},'FaceAlpha',0.7);hold on;
    h_LagPerBinFT(4).FaceColor = shade([1 0 1], 0.5); hold on;
    if i<6 
        text(H_x_lag{i}(18), MAXY-1.6, 'F-NoT','Color',shade([1 0 1],0.5))
    elseif i==6 || i==7
        text(H_x_lag{i}(8), MAXY-1.6, 'F-NoT','Color',shade([1 0 1],0.5))
    elseif i== 8 || i==9 || i==10
        text(H_x_lag{i}(6), MAXY-1.6, 'F-NoT','Color',shade([1 0 1],0.5))
    elseif i== 11 || i==12 
        text(H_x_lag{i}(3), MAXY-1.6, 'F-NoT','Color',shade([1 0 1],0.5))
    end
    ylim([0 MAXY])
    end
    ylabel('# Pairs')
    xlabel('Lag')
    
    if ~isempty(LagPerBinNoZero_CT{1}{i})
     h_LagPerBinCT(1)=histogram(LagPerBinNoZero_CT{1}{i}(:,1),H_x_lag{i},'FaceAlpha',0.7);hold on;
    h_LagPerBinCT(1).FaceColor = 'b'; hold on;
    if i <6
        text(H_x_lag{i}(34), MAXY-0.4, 'C-T1','Color','b')
    elseif i==6 || i==7
        text(H_x_lag{i}(14), MAXY-0.4, 'C-T1','Color','b')
    elseif i== 8 || i ==9 || i==10
        text(H_x_lag{i}(10), MAXY-0.4, 'C-T1','Color','b')
    elseif i== 11 || i ==12 
        text(H_x_lag{i}(4), MAXY-0.4, 'C-T1','Color','b')
    end
    end
     ylim([0 MAXY])
    ylabel('# Pairs')
    xlabel('Lag')
    box on;
    
    if ~isempty(LagPerBinNoZero_CT{2}{i})
     h_LagPerBinCT(2)=histogram(LagPerBinNoZero_CT{2}{i}(:,1),H_x_lag{i},'FaceAlpha',0.7);hold on;
     h_LagPerBinCT(2).FaceColor = shade([0 0 1],0.5); hold on;
     if i<6
         text(H_x_lag{i}(34), MAXY-0.8, 'C-T2','Color',shade([0 0 1],0.5))
     elseif i==6 || i==7
         text(H_x_lag{i}(14), MAXY-0.8, 'C-T2','Color',shade([0 0 1],0.5))
     elseif i== 8 || i ==9 || i==10
        text(H_x_lag{i}(10), MAXY-0.8, 'C-T2','Color',shade([0 0 1],0.5))
     elseif i== 11 || i ==12 
       text(H_x_lag{i}(4), MAXY-0.8, 'C-T2','Color',shade([0 0 1],0.5))
     end
%     l(2)=legend(TextIdMT{2},'Location','NorthWest');hold on;
    ylim([0 MAXY])
    end
    ylabel('# Pairs')
    xlabel('Lag')
    
    box on;
    if ~isempty(LagPerBinNoZero_CT{3}{i})
    h_LagPerBinCT(3)=histogram(LagPerBinNoZero_CT{3}{i}(:,1),H_x_lag{i},'FaceAlpha',0.7);hold on;
    h_LagPerBinCT(3).FaceColor = 'k'; hold on;
    if i<6
         text(H_x_lag{i}(34), MAXY-1.2, 'C-T3','Color','k')
    elseif i==6 || i==7
        text(H_x_lag{i}(14), MAXY-1.2, 'C-T3','Color','k')
    elseif i== 8 || i ==9 || i==10
        text(H_x_lag{i}(10), MAXY-1.2, 'C-T3','Color','k')
    elseif i== 11 || i ==12 
       text(H_x_lag{i}(4), MAXY-1.2, 'C-T3','Color','k')
    end
   ylim([0 MAXY])
    ylabel('# Pairs')
    xlabel('Lag');
    box on;
    
    if ~isempty(LagPerBinNoZero_CT{4}{i})
    h_LagPerBinCT(4)=histogram(LagPerBinNoZero_CT{4}{i}(:,1),H_x_lag{i},'FaceAlpha',0.7);hold on;
    h_LagPerBinCT(4).FaceColor = shade([0 0 0], 0.5); hold on;
    if i<6
         text(H_x_lag{i}(34), MAXY-1.6, 'C-NoT','Color',shade([0 0 0], 0.5))
    elseif i==6 || i==7
        text(H_x_lag{i}(14), MAXY-1.6, 'C-NoT','Color',shade([0 0 0], 0.5))
    elseif i== 8 || i ==9 || i==10
        text(H_x_lag{i}(10), MAXY-1.6, 'C-NoT','Color',shade([0 0 0], 0.5))
    elseif i== 11 || i ==12 
       text(H_x_lag{i}(4), MAXY-1.6, 'C-NoT','Color',shade([0 0 0], 0.5))     
    end
    end
    end
    ylim([0 MAXY])
    ylabel('# Pairs')
    xlabel('Lag')
    
    title({'BinSize:',mat2str(BinSizes(i))});
    box on;
    
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bin Distribution (For the moment all together)
BinSizes =[0.01, 0.015, 0.03, 0.05, 0.08, 0.12,0.25,0.35, 0.5, 0.6, 0.8, 1.6]

for k=1:SumTotAn
    NP_vsvs{k}=[];
    NP_vtavta{k}=[];
    NP_vsvta{k}=[];
end

for k=1:SumTotAn
NP_vsvs{k}=(nneuVS_All{k}*(nneuVS_All{k}+1))/2;
NP_vtavta{k}=(nneuVTA_All{k}*(nneuVTA_All{k}+1))/2;
NP_vsvta{k}=(nneuVS_All{k}*nneuVTA_All{k});
end




BinW=cell(size(struct_pair)); 
WW=cell(size(struct_pair));

for k=1:size(struct_pair,2)
    
    for i=1:size(struct_pair{k}.bin,2)
        for j=1: size(struct_pair{k}.bin{i},1)
    BinW{k}(i,j)=struct_pair{k}.bin{i}(j);
        end
    end
    if size(struct_pair{k}.bin)~=0
    WW{k}=BinW{k}(find(BinW{k}));
    if size(WW{k},2)>1
        WW{k}=WW{k}'
    end
    
    end
end
%%
% % % % % % % % % % % %%
NP_XXX=NP_vsvta;
 An=zeros(size(WW,2),size(BinSizes,2));
for k=1:size(WW,2)
   
    for i=1:size(BinSizes,2)
        count=find(WW{k}(:)==BinSizes(i));
        if size(count,1)~=0
        An(k,i)=size(count,1);
        end
    end
%    AnT(k,:)=An(k,:)/sum(An(k,:));
   AnP(k,:)=An(k,:)/NP_XXX{k}(:);
end
% % % % % % % % % % % 
% MAnT = mean(AnT,'omitnan');
% STAnT = std(AnT,0,1,'omitnan')/sqrt(size(AnT,1));
MAnP = mean(AnP,'omitnan');
STAnP = std(AnP,0,1,'omitnan')/sqrt(size(AnP,1));
% % % % % % % % % % % %%
BinSizesT= mat2str(BinSizes);
% StrTS = strcat(d_sonew.clust_params(1).session_tag);
%% %%%%%%%%%%%%%%%% PLOT BIN DISTRIBUTION

figure(40);hold on;
bar(MAnP,'r');
errorbar(MAnP,STAnP,'k*');
% title({'Mean Bin distribution-VS/VS';'Norm Poss Pairs';StrTS;BinSizesT})
% title({'Mean Bin distribution-VTA/VTA';'Norm Poss Pairs';StrTS;BinSizesT})
%  title({'Mean Bin distribution';'Norm Poss Pairs';StrTS})
 xticks(1:12)
 xticklabels(BinSizes)
 xlabel('Bin Width','Fontsize',14)
 ylabel('# Pairs averaged on animals/ # Possible Pairs')
 
 box on

 
 % % % % % % figure(55);hold on;
% % % % % % bar(MAnT);
% % % % % % errorbar(MAnT,STAnT/2,'c*');
% % % % % % % title({'Mean Bin distribution-VS/VS';'Norm Detected Pairs';StrTS;BinSizesT})
% % % % % % title({'Mean Bin distribution-VTA/VTA';'Norm Detected Pairs';StrTS;BinSizesT})
% % % % % % box on
% % % % % % % title({'Mean Bin distribution-VS/VTA';'Norm Detected Pairs';StrTS;BinSizesT})