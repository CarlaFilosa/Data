%%%% Programs that does the lag distribution for each specific couple of
%%%% units

nameMerge{1} = strcat('Rev2_OrRevW.mat');
nameMerge{2} = strcat('Lastrev1_OrRevW.mat');
nameMerge{3} = strcat('Rev3_OrRevW.mat');
nameMerge{4} = strcat('Revconc_OrRevW.mat');
nameT = strcat('Rev2; Lastrev1; Rev3; Revconc');


alg= strcat('Stand_');

fignum =[1,2,3,4,5,6,7,8];

for i =1: length(nameMerge)

nameNotPruMerge{i} = strcat('An_Disp_',alg,nameMerge{i})     % Not pruned Algorithm
load(nameNotPruMerge{i})
clear c
load('classlist.mat');
%  nnn(i,:)=nn(a);
%  aa(i,:)=a;
[new_data] = parId(new_data,a,nn);


 for k=1:length(new_data.par)
 [new_data] = labels(new_data,c,k);
 end
 
 very_data{i}=new_data;

As_across_binsT{i} = As_across_bins;
As_order_binsT{i} = As_order;

clear new_data a nn As_across_bins As_order 
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

[struct_pairAll] = PairsOrgByPair(pairs_vsvta);

[struct_pairAll] = PairNeuLabelsID(struct_pairAll,Data_All); % Labels of neurons in pairs
clear struct_pair
struct_pair =struct_pairAll;
%%


% struct_pair = [struct_pair_small,struct_pair_large];
c1=1;% MSN 
c2=2; %FSI 
c3=3; %CIN
for i =1:3  %%%% Type of neuorns on the column
  
cvta=i+3;  % Type I, Type, Type III
[LagPerBin_MT{i}] = LagCountPerBinType(struct_pair,BinSizes,c1,cvta);  %MSN small with typeI,TypeII,TypeIII
% [LagPerBin_small_FT{i}] = LagCountPerBinType(struct_pair,BinSizes,c2,cvta);  %FSi small with t1, t2, t3
%  [LagPerBin_small_CT{i}] = LagCountPerBinType(struct_pair,BinSizes,c3,cvta);  %CIN small with t1,t2,t3
 
end
cvta=20;
[LagPerBin_MT{4}] = LagCountPerBinType(struct_pair,BinSizes,c1,cvta);   %MSN small with No id VTA:NoT

[LagPerBin_FT{4}] = LagCountPerBinType(struct_pair,BinSizes,c2,cvta);  % FSI small with No Id VTA

[LagPerBin_CT{4}] = LagCountPerBinType(struct_pair,BinSizes,c3,cvta);  % CIN small with No Id VTA


%% %%%%%%%%%
BinSizes_Small = [0.01, 0.015, 0.03, 0.05, 0.08, 0.12, 0.25];
for i = 1: length(BinSizes_Small)
H_x_lag{i} = -MaxLags(i)-0.5:MaxLags(i)+0.5;
end
%% Find no zero indeces
% Ind can assume the values NoZero or Zero depending if I want to find the
% zero or non zero indices
Ind='NoZero';
[LagIndNoZero_MT] = FindIndLagHisto(LagPerBin_MT,Ind);
[LagIndNoZero_FT] = FindIndLagHisto(LagPerBin_FT,Ind);
[LagIndNoZero_CT] = FindIndLagHisto(LagPerBin_CT,Ind);
for j= 1:length(LagPerBin_MT)
    for i =1:length(BinSizes)
LagPerBin_MT_NoZero{j}{i} = LagPerBin_MT{j}{i}(LagIndNoZero_MT{j}{i});
LagPerBin_FT_NoZero{j}{i} = LagPerBin_FT{j}{i}(LagIndNoZero_FT{j}{i});
LagPerBin_CT_NoZero{j}{i} = LagPerBin_CT{j}{i}(LagIndNoZero_CT{j}{i});
    end
end

TextId = {'T1'; 'T2'; 'T3'; 'NoT'}
for i =1:length(TextId)
TextIdMT{i} = strcat('MT-',TextId{i})
end
figure(6);hold on;
for i = 1: length(BinSizes)
  
    subplot(3,4,i)
    if ~isempty(LagPerBin_MT_NoZero{1}{i})
%      h_LagPerBinMT(1)=histogram(LagPerBin_small_MT{1}{i}(LagIndNoZero_small_MT{1}{i}),H_x_lag{i},'FaceAlpha',0.7);hold on;
    h_LagPerBinMT(1)=histogram(LagPerBin_MT_NoZero{1}{i},H_x_lag{i},'FaceAlpha',0.7);hold on;
    h_LagPerBinMT(1).FaceColor = 'g'; hold on;
    text(H_x_lag{i}(2), 7.8, 'M-T1','Color','g')
    end
    ylim([0 8])
    ylabel('# Pairs')
    xlabel('Lag')
    box on;
    
    if ~isempty(LagPerBin_MT_NoZero{2}{i})
     h_LagPerBinMT(2)=histogram(LagPerBin_MT_NoZero{2}{i},H_x_lag{i},'FaceAlpha',0.7);hold on;
     h_LagPerBinMT(2).FaceColor = shade([0 1 0],0.5); hold on;
     text(H_x_lag{i}(2), 7.4, 'M-T2','Color',shade([0 1 0],0.5))
    end
    ylim([0 8])
    ylabel('# Pairs')
    xlabel('Lag')
    box on;
    
    if ~isempty(LagPerBin_MT_NoZero{3}{i})
    h_LagPerBinMT(3)=histogram(LagPerBin_MT_NoZero{3}{i},H_x_lag{i},'FaceAlpha',0.7);hold on;
    h_LagPerBinMT(3).FaceColor = 'c'; hold on;
    text(H_x_lag{i}(2), 7, 'M-T3','Color','c')
    end
    ylim([0 8])
    ylabel('# Pairs')
    xlabel('Lag')
    box on;
    
    if ~isempty(LagPerBin_MT_NoZero{4}{i})
    h_LagPerBinMT(4)=histogram(LagPerBin_MT_NoZero{4}{i},H_x_lag{i},'FaceAlpha',0.7);hold on;
    h_LagPerBinMT(4).FaceColor = shade([0 1 1], 0.5); hold on;
    text(H_x_lag{i}(2), 6.6, 'M-NoT','Color',shade([0 1 1], 0.5))
    ylim([0 8])
    end

    if ~isempty(LagPerBin_FT_NoZero{1}{i})
     h_LagPerBinFT(1)=histogram(LagPerBin_FT_NoZero{1}{i},H_x_lag{i},'FaceAlpha',0.7);hold on;
    h_LagPerBinFT(1).FaceColor = 'r'; hold on;
    if i<6
       text(H_x_lag{i}(12), 7.8, 'F-T1','Color','r')
    else
        text(H_x_lag{i}(8), 7.8, 'F-T1','Color','r')
    end
    end
    ylim([0 8])
    ylabel('# Pairs')
    xlabel('Lag')
    box on;
    
    if ~isempty(LagPerBin_FT_NoZero{2}{i})
     h_LagPerBinFT(2)=histogram(LagPerBin_FT_NoZero{2}{i},H_x_lag{i},'FaceAlpha',0.7);hold on;
     h_LagPerBinFT(2).FaceColor = shade([1 0 0],0.5); hold on;
%     l(2)=legend(TextIdMT{2},'Location','NorthWest');hold on;
    if i<6
    text(H_x_lag{i}(12), 7.4, 'F-T2','Color',shade([1 0 0],0.5))
    else
     text(H_x_lag{i}(8), 7.4, 'F-T2','Color',shade([1 0 0],0.5))
    end
    end
    ylim([0 8])
    ylabel('# Pairs')
    xlabel('Lag')
    box on;
    
    if ~isempty(LagPerBin_FT_NoZero{3}{i})
    h_LagPerBinFT(3)=histogram(LagPerBin_FT_NoZero{3}{i},H_x_lag{i},'FaceAlpha',0.7);hold on;
    h_LagPerBinFT(3).FaceColor = 'm'; hold on;
    if i<6
        text(H_x_lag{i}(12), 7, 'F-T3','Color','m')
    else
        text(H_x_lag{i}(8), 7, 'F-T3','Color','m')
    end
    end
    ylim([0 8])
    ylabel('# Pairs')
    xlabel('Lag')
    
    if ~isempty(LagPerBin_FT_NoZero{4}{i})
    h_LagPerBinFT(4)=histogram(LagPerBin_FT_NoZero{4}{i},H_x_lag{i},'FaceAlpha',0.7);hold on;
    h_LagPerBinFT(4).FaceColor = shade([1 0 1], 0.5); hold on;
    if i<6 
        text(H_x_lag{i}(12), 6.6, 'F-NoT','Color',shade([1 0 1],0.5))
    else
        text(H_x_lag{i}(8), 6.6, 'F-NoT','Color',shade([1 0 1],0.5))
    end
    ylim([0 8])
    end
    ylabel('# Pairs')
    xlabel('Lag')
    
    if ~isempty(LagPerBin_CT_NoZero{1}{i})
     h_LagPerBinCT(1)=histogram(LagPerBin_CT_NoZero{1}{i},H_x_lag{i},'FaceAlpha',0.7);hold on;
    h_LagPerBinCT(1).FaceColor = 'b'; hold on;
    if i <6
    text(H_x_lag{i}(22), 7.8, 'C-T1','Color','b')
    else
        text(H_x_lag{i}(14), 7.8, 'C-T1','Color','b')
    end
    end
    ylim([0 8])
    ylabel('# Pairs')
    xlabel('Lag')
    box on;
    
    if ~isempty(LagPerBin_CT_NoZero{2}{i})
     h_LagPerBinCT(2)=histogram(LagPerBin_CT_NoZero{2}{i},H_x_lag{i},'FaceAlpha',0.7);hold on;
     h_LagPerBinCT(2).FaceColor = shade([0 0 1],0.5); hold on;
     if i<6
         text(H_x_lag{i}(22), 7.4, 'C-T2','Color',shade([0 0 1],0.5))
     elseif i==6 || i==7
         text(H_x_lag{i}(14), 7.4, 'C-T2','Color',shade([0 0 1],0.5))
     end
%     l(2)=legend(TextIdMT{2},'Location','NorthWest');hold on;
    ylim([0 8])
    end
    ylabel('# Pairs')
    xlabel('Lag')
    
    box on;
    if ~isempty(LagPerBin_CT_NoZero{3}{i})
    h_LagPerBinCT(3)=histogram(LagPerBin_CT_NoZero{3}{i},H_x_lag{i},'FaceAlpha',0.7);hold on;
    h_LagPerBinCT(3).FaceColor = 'k'; hold on;
    if i<6
         text(H_x_lag{i}(22), 7, 'C-T3','Color','k')
    else
        text(H_x_lag{i}(14), 7, 'C-T3','Color','k')
    end
    ylim([0 8])
    ylabel('# Pairs')
    xlabel('Lag');
    box on;
    
    if ~isempty(LagPerBin_CT_NoZero{4}{i})
    h_LagPerBinCT(4)=histogram(LagPerBin_CT_NoZero{4}{i},H_x_lag{i},'FaceAlpha',0.7);hold on;
    h_LagPerBinCT(4).FaceColor = shade([0 0 0], 0.5); hold on;
    if i<6
         text(H_x_lag{i}(22), 6.6, 'C-NoT','Color',shade([0 0 0], 0.5))
    else
        text(H_x_lag{i}(14), 6.6, 'C-NoT','Color',shade([0 0 0], 0.5))
    end
    end
    end
    ylim([0 8])
    ylabel('# Pairs')
    xlabel('Lag')
    
    title({'BinSize:',mat2str(BinSizes(i))});
    box on;
    
    
end
    
% c1=2;% MSN 
% c2=20; 
% 
% BinSizes_Small = [0.01, 0.015, 0.03, 0.05, 0.08, 0.12, 0.25]
% 
% BinSizes_Large = [0.35,0.5,0.6]

% hlagNoNorm=cell(size(BinSizes));
% x_hlag=cell(size(BinSizes));
%  for i=1:length(toArrayLagType)
% if ~isempty(toArrayLagType{i})
% %    if ~isempty(LagPerBin{i})
% hlagNoNorm{i}=histcounts(toArrayLagType{i}(1,:));
% x_hlag{i}= min(toArrayLagType{i}(1,:)):max(toArrayLagType{i}(1,:))
% % x_hlag{i}=min(LagPerBin{i}(:,1)):max(LagPerBin{i}(:,1))
% maxYY(i)=max(hlagNoNorm{i}(1,:));
%    end
%   end
% 
% %%%%
% x_lag = -20:20;
% figure(35);hold on;
% histogram(toArrayLagType{5}(1,find(toArrayLagType{5}~=0)),x_lag)
% 
% figure(35);hold on;
% for i =1: length(toArrayLagType)
% if ~isempty(hlagNoNorm{i})
% subplot(2,4,i)
% Hlag=bar(x_hlag{i},hlagNoNorm{i})
% xlim([-20 20])
% end
% end
