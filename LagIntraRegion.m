%%% Lag Distribution within the regions 
%%%%%% 

MergeDataGenPostPru_1;
%%
MergeDataGenPostPru_Intra2;
%%
for k=1:SumTotAn
pairs_vsvs_small_high{k}=[];
end
for k=1:SumTotAn
    jj=0;
    if ~isempty(pairs_vsvs_small{k}) 
        for i=1:size(pairs_vsvs_small{k},1)
            if (pairs_vsvs_small{k}(i,5)<=0.05)
                jj=jj+1;
                pairs_vsvs_small_high{k}(jj,:)= pairs_vsvs_small{k}(i,:);
        end
        
    end
end
end

pairs_small =pairs_vsvs_small{1};
for i = 1: length(pairs_vsvs_small)-1
    if ~isempty(pairs_vsvs_small{i+1})
pairs_small=[pairs_small;pairs_vsvs_small{i+1}];
    end
end
clear i

pairs_large =pairs_vsvs_large{1};
for i = 1: length(pairs_vsvs_large)-1
    if ~isempty(pairs_vsvs_large{i+1})
pairs_large=[pairs_large;pairs_vsvs_large{i+1}];
    end
end
pairs_small_high =pairs_vsvs_small_high{1};
for i = 1: length(pairs_vsvs_small_high)-1
    if ~isempty(pairs_vsvs_small_high{i+1})
pairs_small_high=[pairs_small_high;pairs_vsvs_small_high{i+1}];
    end
end
clear i

% pairs_small =pairs_vtavta_small{1};
% for i = 1: length(pairs_vtavta_small)-1
%     if ~isempty(pairs_vtavta_small{i+1})
% pairs_small=[pairs_small;pairs_vtavta_small{i+1}];
%     end
% end
% clear i
% 
% pairs_large =pairs_vtavta_large{1};
% for i = 1: length(pairs_vtavta_large)-1
%     if ~isempty(pairs_vtavta_large{i+1})
% pairs_large=[pairs_large;pairs_vtavta_large{i+1}];
%     end
% end


clear i
%%%%%%%% Pairs MSN FSI low firing small
[pairs_MFlow_small] = pairs_NeuSelection(pairs_small,1,2.1);
copy_pairs_intra_small=pairs_MFlow_small;

 [aaa, bbb]=sort(pairs_MFlow_small(:,9:10),2);
  pairs_MFlow_small(:,9:10)=aaa;
  pairs_MFlow_small(bbb(:,1)==2,3)=-pairs_MFlow_small(bbb(:,1)==2,3);
  pairs_MFlow_small(bbb(:,1)==2,1)=pairs_MFlow_small(bbb(:,1)==2,2);
  pairs_MFlow_small(bbb(:,1)==2,2)=pairs_MFlow_small(bbb(:,1)==2,1);
  
  %%%%%%%% Pairs MSN FSI low firing large
% [pairs_MFlow_large] = pairs_NeuSelection(pairs_large,1,2.1);
% copy_pairs_intra_large=pairs_MFlow_large;
% 
%  [aaa, bbb]=sort(pairs_MFlow_large(:,9:10),2);
%   pairs_MFlow_large(:,9:10)=aaa;
%   pairs_MFlow_large(bbb(:,1)==2,3)=-pairs_MFlow_large(bbb(:,1)==2,3);
%   pairs_MFlow_large(bbb(:,1)==2,1)=pairs_MFlow_large(bbb(:,1)==2,2);
%   pairs_MFlow_large(bbb(:,1)==2,2)=pairs_MFlow_large(bbb(:,1)==2,1);
%    
%   
  %%%%%%%%%% Pairs MSn FSI high firing small
[pairs_MFhigh_small] = pairs_NeuSelection(pairs_small_high,1,2.2);
copy_pairs_intra_small=pairs_MFhigh_small;

 [aaa, bbb]=sort(pairs_MFhigh_small(:,9:10),2);
  pairs_MFhigh_small(:,9:10)=aaa;
  pairs_MFhigh_small(bbb(:,1)==2,3)=-pairs_MFhigh_small(bbb(:,1)==2,3);
  pairs_MFhigh_small(bbb(:,1)==2,1)=pairs_MFhigh_small(bbb(:,1)==2,2);
  pairs_MFhigh_small(bbb(:,1)==2,2)=pairs_MFhigh_small(bbb(:,1)==2,1);
  
  
  %%%%%%%%%% Pairs MSn FSI high firing large
% [pairs_MFhigh_large] = pairs_NeuSelection(pairs_large,1,2.2);
% copy_pairs_intra_large=pairs_MFhigh_large;
% 
%  [aaa, bbb]=sort(pairs_MFhigh_large(:,9:10),2);
%   pairs_MFhigh_large(:,9:10)=aaa;
%   pairs_MFhigh_large(bbb(:,1)==2,3)=-pairs_MFhigh_large(bbb(:,1)==2,3);
%   pairs_MFhigh_large(bbb(:,1)==2,1)=pairs_MFhigh_large(bbb(:,1)==2,2);
%   pairs_MFhigh_large(bbb(:,1)==2,2)=pairs_MFhigh_large(bbb(:,1)==2,1);
  
  [~,~,~,LagInSec_low_small,LagInSecNoZero_low_small,LagInSecZero_low_small] = LagCountTypeMatPairs(pairs_MFlow_small,BinSizes);
%   [~,~,~,LagInSec_low_large,LagInSecNoZero_low_large,LagInSecZero_low_large] = LagCountTypeMatPairs(pairs_MFlow_large,BinSizes);
  [~,~,~,LagInSec_high_small,LagInSecNoZero_high_small,LagInSecZero_high_small] = LagCountTypeMatPairs(pairs_MFhigh_small,BinSizes);
%   [~,~,~,LagInSec_high_large,LagInSecNoZero_high_large,LagInSecZero_high_large] = LagCountTypeMatPairs(pairs_MFhigh_large,BinSizes);
  

  
x_hlag{1} = -1:0.02:1;
x_hlag_zero{1} = -0.01:0.02:0.01
x_hlag{2} = -3.6:0.1:3.6;
x_hlag_zero{2} = -0.05:0.1:0.05;

%%%%%%% Lag Indices (NonZeros and Zeros) 
LagIndNoZero_low_small= find(LagInSec_low_small(:,1)~=0);  % all small
LagIndZero_low_small= find(LagInSec_low_small(:,1)==0);

% LagIndNoZero_low_large= find(LagInSec_low_large(:,1)~=0);  % all large
% LagIndZero__low_large= find(LagInSec_low_large(:,1)==0);

%%%%%%% Lag Indices (NonZeros and Zeros) 
LagIndNoZero_high_small= find(LagInSec_high_small(:,1)~=0);  % all small
LagIndZero_high_small= find(LagInSec_high_small(:,1)==0);

% LagIndNoZero_high_large= find(LagInSec_high_large(:,1)~=0);  % all large
% LagIndZero_high_large= find(LagInSec_high_large(:,1)==0);
%%

figure(9); hold on;
 subplot(1,2,1)
h_NoZero = histogram(LagInSecNoZero_low_small,x_hlag{1},'FaceAlpha',0.5); hold on;
h_Zero = histogram(LagInSecZero_low_small,x_hlag_zero{1},'FaceAlpha',0.5); hold on;
h_NoZero.FaceColor ='b';
h_Zero.FaceColor = 'g';
title(strcat('\Delta \in [0.01; 0.25] - FSIs low firing rate -'));hold on;
xlabel('Lag (sec)');hold on;
ylabel('# Pairs');hold on;
set(gca,'FontSize',12)
text(0.5, 60, 'MSN -> FSI','FontSize',12);
text(-1, 60, 'MSN <- FSI','FontSize',12);
% title('Lag Distribution')
 
subplot(1,2,2)
h_NoZero = histogram(LagInSecNoZero_high_small,x_hlag{1},'FaceAlpha',0.5); hold on;
h_Zero = histogram(LagInSecZero_high_small,x_hlag_zero{1},'FaceAlpha',0.5); hold on;
h_NoZero.FaceColor ='b';
h_Zero.FaceColor = 'g';
title(strcat('\Delta \in [0.01; 0.05] - FSIs high firing rate -'));hold on;
set(gca,'FontSize',12)
text(0.5, 7, 'MSN -> FSI','FontSize',12);
text(-1, 7, 'MSN <- FSI','FontSize',12);
ylim([0 8])
% text(0.1, 7.7, 'FSIs high firing rate','FontSize',12);
% title(strcat('\Delta \in [0.35; 0.6]'));hold on;
xlabel('Lag (sec)');hold on;
ylabel('# Pairs');hold on;
text(1.2, 7, nameT,'FontSize',12,'Rotation',270);
% tick=-1:0.2:1
% % xticks(tick);hold on;
% xticklabels(tick);hold on;

%ylim([0 13])
%text(0.8,10,'\Delta \in [0.01, 0.25] sec')
  
   %%
%%%% only low

figure(12); hold on;
 subplot(1,2,1)
h_NoZero = histogram(LagInSecNoZero_low_small,x_hlag{1},'FaceAlpha',0.5); hold on;
h_Zero = histogram(LagInSecZero_low_small,x_hlag_zero{1},'FaceAlpha',0.5); hold on;
h_NoZero.FaceColor ='b';
h_Zero.FaceColor = 'g';
title(strcat('\Delta \in [0.01; 0.25]'));hold on;
xlabel('Lag (sec)');hold on;
ylabel('# Pairs');hold on;
subplot(1,2,2)
h_NoZero = histogram(LagInSecNoZero_low_large,x_hlag{2},'FaceAlpha',0.5); hold on;
h_Zero = histogram(LagInSecZero_low_large,x_hlag_zero{2},'FaceAlpha',0.5); hold on;
h_NoZero.FaceColor ='b';
h_Zero.FaceColor = 'g';
title(strcat('\Delta \in [0.35; 0.6]'));hold on;
xlabel('Lag (sec)');hold on;
ylabel('# Pairs');hold on;
% xticks(tick);hold on;
% xticklabels(tick);hold on;
%ylim([0 13])
%text(0.8,10,'\Delta \in [0.01, 0.25] sec')
  

%%%% only high

figure(13); hold on;
 subplot(1,2,1)
h_NoZero = histogram(LagInSecNoZero_high_small,x_hlag{1},'FaceAlpha',0.5); hold on;
h_Zero = histogram(LagInSecZero_high_small,x_hlag_zero{1},'FaceAlpha',0.5); hold on;
h_NoZero.FaceColor ='b';
h_Zero.FaceColor = 'g';
title(strcat('\Delta \in [0.01; 0.25]'));hold on;
xlabel('Lag (sec)');hold on;
ylabel('# Pairs');hold on;
subplot(1,2,2)
h_NoZero = histogram(LagInSecNoZero_high_large,x_hlag{2},'FaceAlpha',0.5); hold on;
h_Zero = histogram(LagInSecZero_high_large,x_hlag_zero{2},'FaceAlpha',0.5); hold on;
h_NoZero.FaceColor ='b';
h_Zero.FaceColor = 'g';
title(strcat('\Delta \in [0.35; 0.6]'));hold on;
xlabel('Lag (sec)');hold on;
ylabel('# Pairs');hold on;
  