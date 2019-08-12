%% postpruned algorithm
% addpath('/home/carla.filosa/Tactbox_Eleonora/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM');
% load('Everything and then some.mat');

addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM');
% % % % % % % % The matrix new_data is the matrix that contains the info, in particular new_data.spikeT_BegEnd has the 
% % % % % % % % spiketrain cutted properly
path_ad = '/zifnas/Carla/CWEM_Project_DATA/Cutted_Data';
cd(path_ad)
nameNew=strcat('Rev4_OrRevW.mat');
load(nameNew)
alg= strcat('Stand_');
% alg=strcat('Pru_');
path_ad1 = '/zifnas/Carla/CWEM_Project_DATA/Assemblies_Data';
cd(path_ad1)
name=strcat('An_Disp_',alg,nameNew);
load(name);
clear name1
% name_bins=strcat('BinLag_',alg,nameNew);
% load(name_bins)
% MatBin=-MBinAnP;

nameSmallBig=strcat('PostPru_BinRed_An_Disp_Prova',alg,nameNew)

% [pks,locs] = findpeaks(MatBin); %% MAnP contains the bins distributions
%%%%% 
    as_pairs_smallBins = cell(size(as_pairs));
    as_pairs_largeBins = cell(size(as_pairs));
    
    TotAn = size(new_data.par,2);
for k = 1:TotAn
    if ~isempty(as_pairs{k})
        for i=1:length(BinSizes)
            if isempty(as_pairs{k}{i}.n)
                as_pairs_smallBins{k} = cell(1);
                as_pairs_largeBins{k} = cell(1);
            elseif ~isempty(as_pairs{k}{i}.n)
           if i<=7 %locs(1)
        as_pairs_smallBins{k}{i}=as_pairs{k}{i};
           elseif i> 7 %locs(1)
               as_pairs_smallBins{k}{i}=[];
           end
           if (i > 7) && (i< 11) %locs(1) % pruning without 0.8, 1.6
        as_pairs_largeBins{k}{i}=as_pairs{k}{i};
        elseif i>= 11 %locs(1)
               as_pairs_largeBins{k}{i}=[];
           end
       
        end
        end
    end
end

%%%%%%%%%%%%%%%%

clear k i

%%
addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data','-end');
ref_lag=2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% POST DETECTION PRUNING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% 
 
for k=1:TotAn
for l = 1:length(as_pairs_smallBins{k})
 if ~isempty(as_pairs_smallBins{k}{l}) 
   Xassembly_smallBins{k}.bin=as_pairs_smallBins{k}; 
 end
end
end
for k=1:TotAn
    for l= 1:length(as_pairs_largeBins{k})
    if ~isempty(as_pairs_largeBins{k}{l}) 
   Xassembly_largeBins{k}.bin=as_pairs_largeBins{k};
    end
    end
end
%%
%%%%%%% ASSEMBLY REORDERING  small and large
 for k=1:TotAn

    if ~isempty(Xassembly_smallBins{k}) 
[As_across_smallBins{k},As_across_smallBins_index{k}]=assemblies_across_bins(Xassembly_smallBins{k},BinSizes);

%  act_count = 'full'; % act_count = 'partial'; % act_count = 'combined';
%  lagChoice = 'duration'; % lagChoice = 'beginning';
%  [assembly_activity_smallBins{k}]=Assembly_activity_function(As_across_smallBins{k}, Xassembly_smallBins{k}, spM_Cut{k}, BinSizes,lagChoice, act_count);
%  
    end

     if ~isempty(Xassembly_largeBins{k})  
         [As_across_largeBins{k},As_across_largeBins_index{k}]=assemblies_across_bins(Xassembly_largeBins{k},BinSizes);

%  act_count = 'full'; % act_count = 'partial'; % act_count = 'combined';
%  lagChoice = 'duration'; % lagChoice = 'beginning';
%  [assembly_activity_largeBins{k}]=Assembly_activity_function(As_across_largeBins{k}, Xassembly_largeBins{k}, spM_Cut{k}, BinSizes,lagChoice, act_count);

     end
 end

criteria='biggest';
for k=1:TotAn
[As_across_smallBins_pr{k}, As_across_smallBins_pr_index{k}]=pruning_across_bins(As_across_smallBins{k},As_across_smallBins_index{k},Nneu{k},criteria);
[As_across_largeBins_pr{k}, As_across_largeBins_pr_index{k}]=pruning_across_bins(As_across_largeBins{k},As_across_largeBins_index{k},Nneu{k},criteria);
end

%%%%%%%%%%%%%%%%% Save files Post Pruned

%%%%%%%%%%%%%%%% Display Post Pruned
addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data','-end');
new_spM.spikeT_BegEnd= new_data.spikeT_BegEnd;

StrT1= strcat('Session',new_data.params(1).session_tag)
j=0;
 for k=1:TotAn
     if ~isempty(As_across_smallBins{k})
         j=j+1;
 figure(j);hold on;
 display='raw';
%  display='clustered';
[Amatrix_smallBins_pr{k},Binvector_smallBins_pr{k},Unit_order_smallBins_pr{k},As_order_smallBins_pr{k}]=assembly_assignment_matrix(As_across_smallBins_pr{k}, Nneu{k}, BinSizes, display);
     end
 end
 for k=1:TotAn
     if ~isempty(As_across_largeBins{k})
         j=j+1;
 figure(j);hold on;
 display='raw';
%  display='clustered';
[Amatrix_largeBins{k},Binvector_largeBins{k},Unit_order_largeBins{k},As_order_largeBins{k}]=assembly_assignment_matrix(As_across_largeBins{k}, Nneu{k}, BinSizes, display);
     end
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear i j k 



save(nameSmallBig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
