function [Mean_Pr] = PrDim(MTrH_U,MTrH_S,MTrR_U,MTrR_S,MTrM_U,MTrM_S,MTrF_U,MTrF_S)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
if nargin ==2
    
    LHU=[];LHS=[];


for k=1:length(MTrH_U)
    LHU(k)=length(MTrH_U{k});
    LHS(k)=length(MTrH_S{k});
    
    LMax=[sum(LHU),sum(LHS)];
end

[LMaxV,LMaxI]=max(LMax);

if LMaxI == 1 Mean_Pr = MTrH_U;
 elseif LMaxI ==2 Mean_Pr = MTrH_S;
end

elseif nargin ==4
    
    LHU=[];LHS=[];LRU=[];LRS=[];


for k=1:length(MTrH_U)
    LHU(k)=length(MTrH_U{k});
    LHS(k)=length(MTrH_S{k});
    LRU(k)=length(MTrR_U{k});
    LRS(k)=length(MTrR_S{k});
    
    LMax=[sum(LHU),sum(LHS),sum(LRU),sum(LRS)];
end

[LMaxV,LMaxI]=max(LMax);

if LMaxI == 1 Mean_Pr = MTrH_U;
 elseif LMaxI ==2 Mean_Pr = MTrH_S;
   elseif LMAxI == 3 Mean_Pr= MTrR_U;
     elseif LMaxI ==4 Mean_Pr = MTrR_S;
       
end

    
elseif nargin == 8
    
LHU=[];LHS=[];LRU=[];LRS=[];LMU=[];LMS=[];LFU=[];LFS=[];


for k=1:length(MTrH_U)
    LHU(k)=length(MTrH_U{k});
    LHS(k)=length(MTrH_S{k});
    LRU(k)=length(MTrR_U{k});
    LRS(k)=length(MTrR_S{k});
    LMU(k)=length(MTrM_U{k});
    LMS(k)=length(MTrM_S{k});
    LFU(k)=length(MTrF_U{k});
    LFS(k)=length(MTrF_S{k});
    LMax=[sum(LHU),sum(LHS),sum(LRU),sum(LRS),sum(LMU),sum(LMS),sum(LFU),sum(LFS)];
end

[LMaxV,LMaxI]=max(LMax);

if LMaxI == 1 Mean_Pr = MTrH_U;
 elseif LMaxI ==2 Mean_Pr = MTrH_S;
   elseif LMAxI == 3 Mean_Pr= MTrR_U;
     elseif LMaxI ==4 Mean_Pr = MTrR_S;
       elseif LMaxI ==5 Mean_Pr = MTrM_U;
         elseif LMaxI ==6 Mean_Pr= MTrM_S;
            elseif LMaxI ==7 Mean_Pr= MTrF_U;
                elseif LMaxI ==8 Mean_Pr= MTrF_S;
end
end
end

