
clear all;
%%% Non Parametric version of a repeated measures Anova
MergeDataGenPostPru_1
MergeDataGenPostPru_2
%%
ChapOr=1; % these are the phases Original
ChapRev=2;

%%%% GivInEnd gives me the indices in a particular chapter and differently
%%%% for the phases 'stable', 'unstable' and 'all'
[cinOr, cendOr, ccrOr, inOr, finOr] = GiveInEnd( Data_All,ChapOr, 'all' );  % indices for the whole phases 
[cinRev,cendRev, ccrRev, inRev, finRev] = GiveInEnd( Data_All,ChapRev, 'all' );

startOr = cinOr;
stopOr = cendOr;

startRev = cinRev;
stopRev = cendRev;

%%%%%%%%%%%%%%%%%%
start = startOr+1; % because of the baseline I can't take the first trial of the original phase
stop = stopOr;
% start = startRev;
% stop = stopRev;

%%%%% For small bins
As_activity_All = As_activity_smallBins_All;
struct_pair =struct_pair_small;


Region='vsvta';
Dir='small';
hit=1;

for k=1:length(As_activity_All)
   
for i =1: 4
    
[pairs_MT{i}{k}] = [];
[pairs_FT{i}{k}] = [];
[pairs_CT{i}{k}] = [];
end

end


cvs1=1;  cvs2=2; cvs3=3;
for k=1:length(As_activity_All)
    if ~isempty(pairs_vsvta_small{k})
for i =1: 3
    cvta = i+3;
[pairs_MT{i}{k}] = pairs_NeuSelection(pairs_vsvta_small{k},cvs1,cvta);
[pairs_FT{i}{k}] = pairs_NeuSelection(pairs_vsvta_small{k},cvs2,cvta);
[pairs_CT{i}{k}] = pairs_NeuSelection(pairs_vsvta_small{k},cvs3,cvta);
end
[pairs_MT{4}{k}] = pairs_NeuSelection(pairs_vsvta_small{k},cvs1,20);
[pairs_FT{4}{k}] = pairs_NeuSelection(pairs_vsvta_small{k},cvs2,20);
[pairs_CT{4}{k}] = pairs_NeuSelection(pairs_vsvta_small{k},cvs3,20);
    end
end
pairs =pairs_vsvta_small;

FindIndicesAnovaScript;  %%% This script find indices and assembly activity in the ranges I need
%%%%% The Mean across the time AsForMean.. are the matrices in which I am interested to
%%%%% implement the Anova. The time is on the rows and the trials are on
%%%%% the columns


%%%%%Those matrices have trial on the columns and
%%%%% point in time on the rows

for k=1: SumTotAn
    Y.Base{k}=[];Y.AftStim{k}=[];Y.PreRew{k}=[]; Y.AftRew{k}=[];
    AnMat.Base{k} =  []; AnMat.AftRew{k} =  []; AnMat.PreRew{k} =  []; AnMat.AftStim{k} =  [];
end
for k=1: SumTotAn
    if ~isempty(AsForMean_Base{k})
     for i =1:length(AsForMean_Base{k})
AnMat.Base{k}{i} =  (mean(AsForMean_Base{k}{i},1))';
AnMat.AftStim{k}{i} =  (mean(AsForMean_AftStim{k}{i},1))';
AnMat.PreRew{k}{i} =  (mean(AsForMean_PreRew{k}{i},1))';
AnMat.AftRew{k}{i} =  (mean(AsForMean_AftRew{k}{i},1))';
Y.Base{k}{i} = [ones(size(AnMat.Base{k}{i}));zeros(size(AnMat.AftStim{k}{i}));zeros(size(AnMat.PreRew{k}{i}));zeros(size(AnMat.AftRew{k}{i}))];
Y.AftStim{k}{i} = [zeros(size(AnMat.Base{k}{i}));ones(size(AnMat.AftStim{k}{i}));zeros(size(AnMat.PreRew{k}{i}));zeros(size(AnMat.AftRew{k}{i}))];
Y.PreRew{k}{i} = [zeros(size(AnMat.Base{k}{i}));zeros(size(AnMat.AftStim{k}{i}));ones(size(AnMat.PreRew{k}{i}));zeros(size(AnMat.AftRew{k}{i}))];
Y.AftRew{k}{i} = [zeros(size(AnMat.Base{k}{i}));zeros(size(AnMat.AftStim{k}{i}));zeros(size(AnMat.PreRew{k}{i}));ones(size(AnMat.AftRew{k}{i}))];
     end
    else AnMat.Base{k} =  []; AnMat.AftRew{k} =  []; AnMat.PreRew{k} =  []; AnMat.AftStim{k} =  [];
    end
end


% AnMat.BaseRank = tiedrank(AnMat.Base);
% AnMat.AftStimRank= tiedrank(AnMat.AftStim);
% AnMat.PreRewRank= tiedrank(AnMat.PreRew);
% AnMat.AftRewRank= tiedrank(AnMat.AftRew);
for  k=1: SumTotAn
    AnMat.Global{k}=[];
    Y.Global{k} =[];
end
for k=1: SumTotAn
 if ~isempty(AnMat.Base{k})
     for i =1:length(AsForMean_Base{k})
AnMat.Global{k}{i} = [AnMat.Base{k}{i};AnMat.AftStim{k}{i};AnMat.PreRew{k}{i};AnMat.AftRew{k}{i}]; % data
Y.Global{k}{i} =[Y.Base{k}{i},Y.AftStim{k}{i},Y.PreRew{k}{i},Y.AftRew{k}{i}]  % Predictors
     end
 end
end
%%
k=1; i=1;

for k=1: SumTotAn
 if ~isempty(AnMat.Global{k})
     for i =1:length(AnMat.Global{k})
         if ~isempty(AnMat.Global{k}{i})
[bGLM_P{k}{i}, devGLM_P{k}{i},statsGLM_P{k}{i}] = glmfit(Y.Global{k}{i},AnMat.Global{k}{i},'poisson','link','log');
DataEst_GLM{k}{i} = glmval(bGLM_P{k}{i},Y.Global{k}{i},'log')
% [b{k}{i},~,~,~,statsREG{k}{i}] = regress(Y{k}{i},AnMat.Global{k}{i});
% [b,bint,r,rint,stats] = regress(y,X)
         end
     end
 end
end

%%

for k=1: SumTotAn
 if ~isempty(AnMat.Global{k})
     for i =1:length(AnMat.Global{k})
         if ~isempty(AnMat.Global{k}{i})
pValPoiss{k}{i}=pseudoR2(AnMat.Global{k}{i}', DataEst_GLM{k}{i}', mean(AnMat.Global{k}{i}));
         end
     end
 end
end
 