
%%%% Programs that does the lag distribution in seconds for each specific couple of
%%%% units
clear all

addpath('/zifnas/Carla/CWEM_Project_GenProg');
path_p = '/zifnas/Carla/CWEM_Project_GenProg';
path_rd = '/zifnas/Carla/CWEM_Project_DATA/max_Data';
path_dd = '/zifnas/Carla/CWEM_Project_DATA/Cutted_Data';
path_ad = '/zifnas/Carla/CWEM_Project_DATA/Assemblies_Data';
cd(path_rd);
load('methodbase.mat','d','iti')
load('rates.mat')
load('c_def.mat')
cd(path_p);
Data_Par=d.clust_params;
new_rates=p;
clear p
c_def = c;
clear c;
c_def.type1_old=c_def.type1;
c_def.type2_old=c_def.type2;
c_def.type3_old=c_def.type3;
c_def.type1 =[];
c_def.type2 = [];
c_def.type3 = [];

c_def.type1 = c_def.fix{1};
c_def.type2 = c_def.fix{2};
c_def.type3 = c_def.fix{3};
c_true=c_def;
clear c_def;

% load('classlist.mat');
% c_Stable=c;
% c.allFSIs_NoLim=c.allFSIs;
% clear c_Stable.allFSIs;
% c.allFSIs=[];
c_true.FSIs_low=[];
c_true.FSIs_high=[];
oo=0;
uu=0;
for i=1: length(c_true.allFSIs)
%     if (d.clust_params(c.allFSIs(i)).ITI_rate > 5) && (d.clust_params(c.allFSIs(i)).ITI_rate <=45)
if (new_rates.value(c_true.allFSIs(i)) > 5) && (new_rates.value(c_true.allFSIs(i)) <=45)
    oo=oo+1;
    c_true.FSIs_low(oo)=c_true.allFSIs(i);
    end
 if (new_rates.value(c_true.allFSIs(i)) > 45)   
%     if (d.clust_params(c.allFSIs(i)).ITI_rate > 45) 
        uu=uu+1;
        c_true.FSIs_high(uu)=c_true.allFSIs(i)
    end
end
[c_Stable] = ID_To_StableID(c_true,d);
 clear i oo uu
 
[c_Stable,c_true] = opto_neu(c_Stable,c_true,Data_Par);
[c_Stable,c_true] = opto_neu_New(c_Stable,c_true,Data_Par);
% type2LaserCheck;

nameMerge{1} = strcat('Firstrev1_OrRevW.mat');
nameMerge{2} = strcat('Rev1_OrRevW.mat');
nameMerge{3} = strcat('Lastrev1_OrRevW.mat');
nameMerge{4} = strcat('Revconc_OrRevW.mat');
nameMerge{5} = strcat('Endrev_OrRevW.mat');
nameMerge{6} = strcat('Rev2_OrRevW.mat');
nameMerge{7} = strcat('Rev3_OrRevW.mat');
nameMerge{8} = strcat('Rev4_OrRevW.mat');
nameMerge{9} = strcat('Admixture_OrRevW.mat');

nameT = strcat('Very all with Reversal');

%%%%% Ref lag 2
alg= strcat('Stand_');

%%%%% Reflag3 Bin reduction
% alg= strcat('Stand_10Bins_RefLag3_');

fignum =[1,2,3,4,5,6,7,8];
IdType='stableID';
for i =1: length(nameMerge)
cd(path_ad);
nameSmallBigMerge{i} = strcat('PostPru_BinRed_An_Disp_',alg,nameMerge{i})
    load(nameSmallBigMerge{i})
clear c;

cd(path_p);
 for k=1:length(new_data.par)
 [new_data] = labels(new_data,c_Stable,k,IdType);

 end
 
 very_data{i}=new_data;
As_activity_smallBinsT{i} = assembly_activity_smallBins;
As_activity_largeBinsT{i} = assembly_activity_largeBins;
As_across_smallBinsT{i} = As_across_smallBins;
As_order_smallBinsT{i} = As_order_smallBins;
As_across_largeBinsT{i} = As_across_largeBins;
As_order_largeBinsT{i} = As_order_largeBins; 
clear new_data a nn As_across_smallBins As_order_smallBins As_across_largeBins As_order_largeBins As_activity_smallBins As_activity_largeBins

end

clearvars -except Data_Par new_rates c_true very_data As_activity_smallBinsT As_activity_largeBinsT As_across_smallBinsT As_order_smallBinsT As_across_largeBinsT As_order_largeBinsT c_Stable nameMerge nameT fignum BinSizes MaxLags nameT path_ad path_rd path_dd path_p
