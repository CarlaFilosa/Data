%% Assembly Activity alligned to the reward

addpath('/home/carla.filosa/Tactbox_Eleonora/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder');
% The matrix new_spM is built in the program MaxDataPrepro.m
 load('Everything and then some.mat')
%%
% load('NewDataLastrev1.mat')
% load('AnalysisLastrev1.mat')
%  load('DispAnLastrev1Clust.mat')

% load('NewDataFirstrev1.mat')
% load('AnalysisFirstrev1.mat')
% load('DispAnFirstrev1Clust.mat')

% load('NewDataEndrev.mat')
% load('AnalysisEndrev.mat')
%  load('DispAnEndrevClust.mat')
%  
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Only the original Phase
load('NewDataLastrev1.mat')
load('AnalysisLastrev1_Orig.mat')
 load('DispAnLastrev1_Orig_Clust.mat')
%  
%  load('NewDataFirstrev1.mat')
% load('AnalysisFirstrev1_Orig.mat')
% load('DispAnFirstrev1_Orig_Clust.mat')

% load('NewDataEndrev.mat')
% load('AnalysisEndrev_Orig.mat')
% load('DispAnEndrev_Orig_Clust.mat')
% Add licks file and info and parameters



spM=d_sonew.spiketrains;
params=d_sonew.clust_params;
licks=d_sonew.licks;
info = d_sonew.info;
laser= d_sonew.laser;

for i=1:size(spM,2)
nn(i) = size(spM{i},1);
tt(i)= size(spM{i},2);
end
a=find(nn);

nn(a);
for i=1:length(a)
    new_spM.licks{i}=licks{a(i)};
    new_spM.laser{i}=laser{a(i)};
    new_spM.info(i)=info(a(i));
    new_spM.par.red(i)=params(a(i)); 
    
end
new_spM.par.all=d_sonew.clust_params;
clear nn tt i a

%% Make each trial centered on the reward, traslating times of timestamps of stimuli and licks

% Select only the rewarded events:

Copy_new_spM=new_spM;

new_spM.events=[];

for k=1:size(Copy_new_spM.events,2)
    l=0;
    for j=1:size(Copy_new_spM.events{k},1)
        if Copy_new_spM.events{k}(j).rew_code==1 & size(Copy_new_spM.events{k}(j).reward_time==1)
           l=l+1; 
           new_spM.events{k}(l)=Copy_new_spM.events{k}(j);
        end
    end
end
%% Allign to the reward

minInt=0.01;
for k= 1:size(new_spM.events,2)
    for l = 1: length(new_spM.events{k})
        new_spM.events{k}(l).fv_onRes=new_spM.events{k}(l).fv_on-new_spM.events{k}(l).reward_time;
        new_spM.events{k}(l).fv_offRes=new_spM.events{k}(l).fv_off-new_spM.events{k}(l).reward_time;
        new_spM.events{k}(l).Rew_cmin=new_spM.events{k}(l).reward_time-1-minInt/2;
        new_spM.events{k}(l).Rew_cplus=new_spM.events{k}(l).reward_time+5+minInt/2;
        if ~isempty(new_spM.events{k}(l).licks) & ~isempty(new_spM.events{k}(l).reward_time)
            new_spM.events{k}(l).licksRes = new_spM.events{k}(l).licks(1:end)-new_spM.events{k}(l).reward_time;
        end
    end
end
% in this way I will know that the 100th bin is equal to zero
clear k l
%%
for k= 1:size(new_spM.events,2)
    for l = 1: length(new_spM.events{k})
       rew{k}(l).res =(new_spM.events{k}(l).fv_offRes-new_spM.events{k}(l).fv_onRes)/2;
       rew{k}(l).cent= new_spM.events{k}(l).Rew_cplus - new_spM.events{k}(l).Rew_cmin;
    end
end
%% %%%% Assembly activity around the reward interval


% k = index for animal
% i = index for pairs for each animal
% l = index for trials

for k= 1:size(new_spM.events,2)
for l = 1: length(new_spM.events{k})-1
   for i =1: size(assembly_activity{k},1)
%        Dist_fvP(k,l)=new_spM.events{k}(l+1).fv_on-new_spM.events{k}(l).fv_on;
%        Dist_fv(k,l)=new_spM.events{k}(l).fv_off-new_spM.events{k}(l).fv_on;
       assA{k}{i}{l}=assembly_activity{k}{i}(find(assembly_activity{k}{i}(:,1)>=new_spM.events{k}(l).Rew_cmin & assembly_activity{k}{i}(:,1)<=new_spM.events{k}(l).Rew_cplus),:);
%        assA{k}{i}{l}(:,3)=assA{k}{i}{l}(:,1)-new_spM.events{k}(l).reward_time;
   end
end
end



% %% %%%% Assembly activity in the times of stimulus activation
% %%%%%%%%%%%%%%%%%%%%%%%%%% Two version:
% 
% %%%%%%%%%First one: without rescaling time 
% 
% % k = index for animal
% % i = index for pairs for each animal
% % l = index for trials
% 
% for k= 1:size(new_spM.events,2)
% for l = 1: length(new_spM.events{k})-1
%    for i =1: size(assembly_activity{k},1)
%        Dist_fvP(k,l)=new_spM.events{k}(l+1).fv_on-new_spM.events{k}(l).fv_on;
%        Dist_fv(k,l)=new_spM.events{k}(l).fv_off-new_spM.events{k}(l).fv_on;
%        assA{k}{i}{l}=assembly_activity{k}{i}(find(assembly_activity{k}{i}(:,1)>=new_spM.events{k}(l).fv_on & assembly_activity{k}{i}(:,1)<=new_spM.events{k}(l).fv_on+Dist_fvP(k,l)),:);
%        assA{k}{i}{l}(:,3)=assA{k}{i}{l}(:,1)-new_spM.events{k}(l).reward_time;
%    end
% end
% end
% 
% 
% clear k i l
% 


%% %%%%%%%%%%%%%%%%%%%%%%%%% Second version rescaling the time and the activity of the assembly on a time binned
%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Rescale on the same scale the assemblies
%   addpath('/home/carla.filosa/Tactbox_Eleonora/Cell_assembly_detection/Programs_and_data-end');

for k=1:size(assembly_activity,2)
binref.B{k}=new_spM.events{k}(1).fv_on:minInt:max(max(new_spM.spikeT_BegEnd{k}));
 for i=1:size(assembly_activity{k},1)
    
%       [AsTime{k}{i}(:,2)]=time_rescaleE(assembly_activity{k}{i}(:,2),assembly_activity{k}{i}(:,1), binref.B{k});
%      [AsTime{k}{i}(:,1)]=time_rescaleE(assembly_activity{k}{i}(:,1),assembly_activity{k}{i}(:,1), binref.B{k});
      [AsTime{k}{i}(:,1)]=TimeRescaleMin(assembly_activity{k}{i}(:,1),assembly_activity{k}{i}(:,1), binref.B{k});
      [AsTime{k}{i}(:,2)]=TimeRescaleMin(assembly_activity{k}{i}(:,2),assembly_activity{k}{i}(:,1), binref.B{k});
%       [AsTime{k}{i}(:,3)]=TimeRescaleMin(assembly_activity{k}{i}(:,3),assembly_activity{k}{i}(:,1), binref.B{k});
 end
end


% Bin references for each trial:
%%
for k=1:size(assembly_activity,2)
    for l=1:length(new_spM.events{k})-1
%         Dist_fvP(k,l)=new_spM.events{k}(l+1).fv_on-new_spM.events{k}(l).fv_on;
%         Dist_fv(k,l)=new_spM.events{k}(l).fv_off-new_spM.events{k}(l).fv_on;        
     FindInd{k}{l}= find(binref.B{k}>=new_spM.events{k}(l).Rew_cmin & binref.B{k}<=new_spM.events{k}(l).Rew_cplus);
    end
end

for k= 1: size(assembly_activity,2)
    for l=1:length(new_spM.events{k})-1
        for i =1: size(AsTime{k},2)
        asRes{k}{i}{l}(:,1)=AsTime{k}{i}(FindInd{k}{l},1);
        asRes{k}{i}{l}(:,2)=AsTime{k}{i}(FindInd{k}{l},2);
%         asRes{k}{i}{l}(:,3)=asRes{k}{i}{l}(:,1)-new_spM.events{k}(l).reward_time;
        end
    end
end



clear i j k l n
%% I want to plot the pairs VS/VTA VS/VS and VTA/VTA separately
%%%%% I build my usual struct_pair
% Find number of units of the Striatum and VTA nneuS

for k=1:size(new_spM.spike_regionNoId,2)
    nneuVS{k}=0;
end

for k=1:size(new_spM.spike_regionNoId,2)
    for j=1:size(new_spM.spike_regionNoId{k},2)
    if new_spM.spike_regionNoId{k}(1,j)==1
    nneuVS{k}=nneuVS{k}+1;
    end
    end
    nneuVTA{k}=size(new_spM.spike_regionNoId{k},2)-nneuVS{k};
end

%%%%%%%%% Find Pairs in one region and shared between two different regions
 for k= 1:size(As_across_bins,2)
pairs_lag{k}=[];
 pairs_lag_vsvs{k}=[];
 pairs_lag_vsvta{k} =[];
pairs_lag_vtavta{k} =[];
 end

for k= 1:size(As_across_bins,2)
pairs_lag{k}=nan(length(As_across_bins{k}),6);
for i=1:length(As_across_bins{k})
     ii=As_order{k}(i);
     pairs_lag{k}(i,6)=As_order{k}(i);    
%     ii=i;
%      pairs_lag(i,6)=i;    
    pairs_lag{k}(i,1:2)=As_across_bins{k}{ii}.elements;   
    pairs_lag{k}(i,3)=As_across_bins{k}{ii}.lag(end);
    pairs_lag{k}(i,4)=As_across_bins{k}{ii}.pr;
    pairs_lag{k}(i,5)=As_across_bins{k}{ii}.bin;    
    
    
end

[a{k}, b{k}]=sort(pairs_lag{k}(:,1:2),2);

 pairs_lag{k}(:,1:2)=a{k};
 pairs_lag{k}(b{k}(:,1)==2,3)=-pairs_lag{k}(b{k}(:,1)==2,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% divide pairs
% 
% 
A{k}=pairs_lag{k}(:,1)<=nneuVS{k};
B{k}=pairs_lag{k}(:,2) > nneuVS{k};
pairs_lag_vsvta{k}=pairs_lag{k}(logical(A{k}.*B{k}),:);
      

A{k}=pairs_lag{k}(:,1)>nneuVS{k};
B{k}=pairs_lag{k}(:,2)>nneuVS{k};
pairs_lag_vtavta{k}=pairs_lag{k}(logical(A{k}.*B{k}),:);


A{k}=pairs_lag{k}(:,1)<=nneuVS{k};
B{k}=pairs_lag{k}(:,2)<=nneuVS{k};
pairs_lag_vsvs{k}=pairs_lag{k}(logical(A{k}.*B{k}),:);       
end
%%
clear a b A B i ii j k
%% # of possible pairs in the single regions and intra-region
for k=1:size(nneu,2)
    NP_vsvs{k}=[];
    NP_vtavta{k}=[];
    NP_vsvta{k}=[];
end

for k=1:size(nneu,2)
NP_vsvs{k}=(nneuVS{k}*(nneuVS{k}+1))/2;
NP_vtavta{k}=(nneuVTA{k}*(nneuVTA{k}+1))/2;
NP_vsvta{k}=(nneuVS{k}*nneuVTA{k});
end

%% %%%%%%%%%%%%%%%%%%%%
for k=1:size(pairs_lag,2) % k= index for number of animals
pairs_lag_XXX{k}=[];
NP_XXX{k}=[];
end

for k=1:size(pairs_lag,2)
% pairs_lag_XXX{k}=pairs_lag_vsvs{k};
% NP_XXX{k}=NP_vsvs{k};
% pairs_lag_XXX{k}=pairs_lag_vtavta{k};
% NP_XXX{k}=NP_vtavta{k};
pairs_lag_XXX{k}=pairs_lag_vsvta{k};
NP_XXX{k}=NP_vsvta{k};
% struct for pair to overlap the bin in the figure

struct_pair{k} = struct ('pair',[], 'lag',[], 'pv', [],'bin',[],'ind',[],'new_ind',[],'imean',[]);
j=1; % j is uncremented in the second elseif
   
for i = 1: size(pairs_lag_XXX{k},1)
     if i == 1, 
         struct_pair{k}.pair{j}(1,1:2) = pairs_lag_XXX{k}(1,1:2);
         struct_pair{k}.lag{j} = pairs_lag_XXX{k}(1,3);
         struct_pair{k}.pv{j} = pairs_lag_XXX{k}(1,4);
         struct_pair{k}.bin{j} = pairs_lag_XXX{k}(1,5);
         struct_pair{k}.ind{j} = pairs_lag_XXX{k}(1,6);
         struct_pair{k}.new_ind{j}(i,:) =j;
         struct_pair{k}.imean{j}(i,:) =i;
    elseif i>1 & ismember(pairs_lag_XXX{k}(i,1:2),pairs_lag_XXX{k}(i-1,1:2)) == [1 ,1],
        struct_pair{k}.pair{j}(i,1:2) = pairs_lag_XXX{k}(i,1:2);
        struct_pair{k}.lag{j}(i,:) = pairs_lag_XXX{k}(i,3);
        struct_pair{k}.pv{j}(i,:) = pairs_lag_XXX{k}(i,4);
        struct_pair{k}.bin{j}(i,:) = pairs_lag_XXX{k}(i,5);
        struct_pair{k}.ind{j}(i,:) = pairs_lag_XXX{k}(i,6);
        struct_pair{k}.new_ind{j}(i,:) =j;
        struct_pair{k}.imean{j}(i,:) =i;
    elseif i>1 & (ismember(pairs_lag_XXX{k}(i,1:2),pairs_lag_XXX{k}(i-1,1:2)) == [0 ,1] | ismember(pairs_lag_XXX{k}(i,1:2),pairs_lag_XXX{k}(i-1,1:2)) == [1,0] | ismember(pairs_lag_XXX{k}(i,1:2),pairs_lag_XXX{k}(i-1,1:2)) == [0,0])
        j=j+1;
       struct_pair{k}.pair{j}(i,1:2) = pairs_lag_XXX{k}(i,1:2);
       struct_pair{k}.lag{j}(i,:) = pairs_lag_XXX{k}(i,3);
       struct_pair{k}.pv{j}(i,:) = pairs_lag_XXX{k}(i,4);
       struct_pair{k}.bin{j}(i,:) = pairs_lag_XXX{k}(i,5);
       struct_pair{k}.ind{j}(i,:) = pairs_lag_XXX{k}(i,6);
       struct_pair{k}.new_ind{j}(i,:) =j;
       struct_pair{k}.imean{j}(i,:) =i;
     end
        
   %j=j+1;
end

end

for k=1:size(struct_pair,2)
    for j=1:size(struct_pair{k}.pair,2)
n= 0; 
for i = size(struct_pair{k}.pair{j},1) - n:-1:1
    if struct_pair{k}.pair{j}(i,1)==0, 
        struct_pair{k}.pair{j}(i,:)=[]; 
        struct_pair{k}.lag{j}(i,:)=[]; 
        struct_pair{k}.pv{j}(i,:)=[];
        struct_pair{k}.bin{j}(i,:)=[];
        struct_pair{k}.ind{j}(i,:)=[];
        struct_pair{k}.new_ind{j}(i,:)=[]; 
        struct_pair{k}.imean{j}(i,:) =[];
        n=n+1;
    end
end
    end
end

clear i j k n l




%% %%%%%%%%%% mean across trials

% k=4;
for k=1:size(struct_pair,2)
count=0;
for i= 1: size(struct_pair{k}.pair,2)
cluster_by_lag{i}=find(struct_pair{k}.lag{i}<0);
a=size(cluster_by_lag{i},1);
if a~=0
    count=count+1;
for iii= 1:a
fff{i}=struct_pair{k}.ind{i};
aaa = fff{i}(iii);
in=0;
oo=0;
   for j=2:size(asRes{k}{aaa},2)
       
in=in+2;
oo=oo+1;
maxA=max(asRes{k}{fff{i}(iii)}{j}(:,2));
 if maxA==0 
     maxA=1; end

%  maxT(oo)=max(asRes{k}{aaa}{j}(:,3));
 
 for ggg=1:size(asRes{k}{aaa}{j},1)
   SzA{k}{count}(ggg,1) = size(asRes{k}{aaa}{j}(:,2),1);
   AsForMeanA{k}{count}(ggg,oo)=asRes{k}{aaa}{j}(ggg,2)/maxA;
 end
%end
    end
%  maxTtot=max(maxT);
MeanTrA{k}{count}=mean(AsForMeanA{k}{count}(1:min(SzA{k}{count}),:),2);
StETrA{k}{count}=std(AsForMeanA{k}{count}(1:min(SzA{k}{count}),:),1,2)/sqrt(size(AsForMeanA{k}{count},2));

end
end
end
end

%%
clear i ii j jj k iii in l ggg count fff aaa oo

%%
l=0;

for k=1:size(MeanTrA,2)
    if size(MeanTrA{k},2)~=0
        l=l+1;
        for i=1:size(MeanTrA{k},2)
        SZMA(l)=size(MeanTrA{k}{i},1);
        end
    end
end

% l=0;
% for k=1:size(MeanTrB,2)
%     if size(MeanTrB{k},2)~=0
%         l=l+1;
%         for i=1:size(MeanTrB{k},2)  
%         SZMB(l)=size(MeanTrB{k}{i},1);
%         end
%     end
% end

clear l
%%
l=0;
count=0;
SZCell=[];
summa=[];
MatMeanTrA=[];
for k=1:size(MeanTrA,2)
    if size(MeanTrA{k},2)~=0
        l=l+1;
        SZCell(l)=size(MeanTrA{k},2);
        summa(l)= sum(SZCell);
        for i=1:size(MeanTrA{k},2)
           if l==1
              MatMeanTrA(1:min(SZMA),i)= MeanTrA{k}{i}(1:min(SZMA));
           elseif l~=1

            MatMeanTrA(1:min(SZMA),i+summa(l-1))= MeanTrA{k}{i}(1:min(SZMA));

            end
        end
    end
end


MeanTotTrA=mean(MatMeanTrA,2);
StrTotTrA=std(MatMeanTrA,0,2)/sqrt(size(MatMeanTrA,2));
% plot(MeanTotTrA)
%%
vect=-1:0.01:5;
%
%%
minv=min(MeanTotTrA-1.5*StrTotTrA);
maxv=max(MeanTotTrA+1.5*StrTotTrA);
figure(1);hold on;
box on;
[l,p]=boundedline(vect,MeanTotTrA,StrTotTrA,'-b','alpha');hold on;
%plot([0 0],[minv,maxv],'k.-');hold on;
vline(0,'k-.','reward')
ylim([minv maxv])
title('Negative-Lags Pairs Mean Hit')
ylabel('Mean activity')
xlabel('Time interval (sec)')
% plot(vect,MatMeanTrA(:,50))

%% Mean single animal

for k=1:size(MeanTrA,2)
    for i=1:size(MeanTrA{k},2)
    ForMeanAnH{k}(:,i)=MeanTrA{k}{i};
    end

end
% Mean for animal HiT
for k=1:size(ForMeanAnH,2)
     if size(ForMeanAnH{k},2)~=0
         MeanAnH{k}=mean(ForMeanAnH{k},2);
         StErAnH{k}=std(ForMeanAnH{k},0,2)/sqrt(size(ForMeanAnH{k},2));
     end
end

% for k=1:size(MeanTrB,2)
%     for i=1:size(MeanTrB{k},2)
%     ForMeanAnCR{k}(:,i)=MeanTrB{k}{i};
%     end
% 
% end
% %Mean for animal Correct rejection
% for k=1:size(ForMeanAnCR,2)
%      if size(ForMeanAnCR{k},2)~=0
%          MeanAnCR{k}=mean(ForMeanAnCR{k},2);
%          StErAnCR{k}=std(ForMeanAnCR{k},0,2)/sqrt(size(ForMeanAnCR{k},2));
%      end
% end
% %%
l=0;
for k=1:size(MeanAnH,2)
    if size(MeanAnH{k},2)~=0
    l=l+1;
    vecH{k}=-1:0.01:5;
%     vecCR{k}=0.01:0.01:(size(MeanAnCR{k},1))*0.01;
    figure(l);hold on;
%     subplot(1,2,1);hold on;
    [mh,ph]=boundedline(vecH{k},MeanAnH{k},StErAnH{k},'-b','alpha');hold on;
%  plot([median(Dist_fv(k,:)) median(Dist_fv(k,:))], [min(MeanAnH{k})-0.05 max(MeanAnH{k})+0.05], 'k-.'); hold on;
%     subplot(1,2,2);hold on
%     [mc,pc]=boundedline(vecCR{k},MeanAnCR{k},StErAnCR{k},'-m','alpha');hold on;
%      plot([median(Dist_fv(k,:)) median(Dist_fv(k,:))], [min(MeanAnCR{k})-0.05 max(MeanAnCR{k})+0.05], 'k-.'); hold on;
    end
end

