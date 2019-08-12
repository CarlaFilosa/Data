%% postpruned algorithm
% addpath('/home/carla.filosa/Tactbox_Eleonora/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM');
% load('Everything and then some.mat');

addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM');
% % % % % % % % The matrix new_data is the matrix that contains the info, in particular new_data.spikeT_BegEnd has the 
% % % % % % % % spiketrain cutted properly

name1=strcat('Rev1_OrRevW.mat');
load(name1)
alg= strcat('Stand_');
% alg=strcat('Pru_');
name=strcat('An_Disp_',alg,name1);
% name=strcat('An_Disp_',name);
load(name);
clear name1
name1=strcat('Rev1_OrRevW.mat');
% name_bins=strcat('BinLag_',alg,name1);
% load(name_bins)
% clear name1
% name1=strcat('Lastrev1_OrRevW.mat');
namePostPru=strcat('PostPru_All_An_Disp_',name1);
%%

%%%%%%%%%%%%%%%%

clear k i


addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data','-end');
ref_lag=2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% POST DETECTION PRUNING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% 
 
% for k=1:TotAn
%    Xassembly{k}.bin=as_pairs{k}; 
%  
% end
% 
% %%%%%%% ASSEMBLY REORDERING  small and large
% for k=1:TotAn
%     if ~isempty(Xassembly{k}.bin)
% [As_across_bins{k},As_across_bins_index{k}]=assemblies_across_bins(Xassembly{k},BinSizes);
% 
%  act_count = 'full'; % act_count = 'partial'; % act_count = 'combined';
%  lagChoice = 'duration'; % lagChoice = 'beginning';
%  [assembly_activity{k}]=Assembly_activity_function(As_across_bins{k}, Xassembly{k}, spM_Cut{k}, BinSizes,lagChoice, act_count);
%  
%     end
% end

criteria='biggest';
for k=1:TotAn
[As_across_bins{k}, As_across_bins_index{k}]=pruning_across_bins(As_across_bins{k},As_across_bins_index{k},Nneu{k},criteria);
end

%%%%%%%%%%%%%%%%% Save files Post Pruned

%%%%%%%%%%%%%%%% Display Post Pruned
addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data','-end');
new_spM.spikeT_BegEnd= new_data.spikeT_BegEnd;

StrT1= strcat('Session',new_data.params(1).session_tag)
j=0;
 for k=1:TotAn
     if ~isempty(As_across_bins{k})
         j=j+1;
 figure(j);hold on;
 display='raw';
%  display='clustered';
[Amatrix{k},Binvector{k},Unit_order{k},As_order{k}]=assembly_assignment_matrix(As_across_bins{k}, Nneu{k}, BinSizes, display);
     end
 end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear i j k 



save(namePostPru)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
