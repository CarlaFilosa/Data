addpath('/zifnas/Carla/Tactbox/Cell_assembly_detection/Programs_and_data/Carla/Carla_MaxFolder/CWEM')
% addpath([pwd '/CWEM/']);
load('AfterFriedPairs_VeryAllRevPara_Rev_VsVta_wOp.mat','newnew_pairs_small');

% Pairs_Aft_Fried = struct('lag',[],'bin',[],'S_ID_neu1',[],'S_ID_neu2',[]);
Pairs_Aft_Fried = nan(size(newnew_pairs_small,1),4);
Pairs_Aft_Fried(:,1) = newnew_pairs_small(:,3);
Pairs_Aft_Fried(:,2) = newnew_pairs_small(:,5);
Pairs_Aft_Fried(:,3:4) = newnew_pairs_small(:,7:8);
String = {'Lag','Bin','StableID_neu1','StableID_neu2'}
AftFried_Rev_Table= array2table(Pairs_Aft_Fried, 'VariableName',String);
filename = strcat('AftFried_Rev_ForMax.mat');
filepath = [pwd,'/ForMax/']
% fullfile(filename, filepath)
save(fullfile(filepath,filename),'AftFried_Rev_Table')
