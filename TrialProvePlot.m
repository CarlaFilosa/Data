clear all
addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla')
[ActAll,MaxInd,SortMaxInd] = ActAsAllSameLengthInt( SumTotAn, StETr_HLick)

figure(7);
boundedline(1:140,Act_HLick(find(pairs_Mat(:,8)==1261 & pairs_Mat(:,7)==1254),:),ActAll(find(pairs_Mat(:,8)==1261 & pairs_Mat(:,7)==1254),:))

As_HLick_Tr = AsForMean_HLick{1};
for i=1:TotAn-1
    As_HLick_Tr = [As_HLick_Tr, AsForMean_HLick{i+1}];
    
end

figure(11);
plot(-0.695:0.01:0.695,As_HLick_Tr{find(pairs_Mat(:,8)==1261 & pairs_Mat(:,7)==1252)}); hold on;
boundedline(-0.695:0.01:0.695,Act_HLick(find(pairs_Mat(:,8)==1261 & pairs_Mat(:,7)==1252),:),ActAll(find(pairs_Mat(:,8)==1261 & pairs_Mat(:,7)==1252),:))