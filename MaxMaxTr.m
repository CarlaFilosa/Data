function [MaxATot] = MaxMaxTr(MaxA_HS,MaxA_RS,MaxA_FS,MaxA_MS,MaxA_HU,MaxA_RU,MaxA_FU,MaxA_MU)

if nargin == 1
    MaxATot = [];
    SIZE = length(MaxA_HS);
    for i = 1:SIZE
        MaxATot(i) = max(MaxA_HS(i));
    end

elseif nargin == 2
MaxATot = [];
SIZE=max([length(MaxA_HS),length(MaxA_RS)]);
[MaxHS] = Complete(MaxA_HS,SIZE);
[MaxRS] = Complete(MaxA_RS,SIZE);

for i=1:SIZE
    MaxATot(i)=max([MaxHS(i),MaxRS(i)]);
end
    
% for k=1:size(MaxA_HS,2)
elseif nargin == 4 
MaxATot=[];
SIZE=max([length(MaxA_HS),length(MaxA_RS),length(MaxA_FS),length(MaxA_MS)]);
[MaxHS] = Complete(MaxA_HS,SIZE);
[MaxRS] = Complete(MaxA_RS,SIZE);
[MaxFS] = Complete(MaxA_FS,SIZE);
[MaxMS] = Complete(MaxA_MS,SIZE);

for i=1:SIZE
    MaxATot(i)=max([MaxHS(i),MaxRS(i),MaxFS(i),MaxMS(i)]);
end
    
elseif nargin == 8 
    
MaxATot=[];
SIZE=max([length(MaxA_HS),length(MaxA_RS),length(MaxA_FS),length(MaxA_MS),length(MaxA_HU),length(MaxA_RU),length(MaxA_FU),length(MaxA_MU)]);
[MaxHS] = Complete(MaxA_HS,SIZE);
[MaxRS] = Complete(MaxA_RS,SIZE);
[MaxFS] = Complete(MaxA_FS,SIZE);
[MaxMS] = Complete(MaxA_MS,SIZE);

[MaxHU] = Complete(MaxA_HU,SIZE);
[MaxRU] = Complete(MaxA_RU,SIZE);
[MaxFU] = Complete(MaxA_FU,SIZE);
[MaxMU] = Complete(MaxA_MU,SIZE);
% end

for i=1:SIZE
    MaxATot(i)=max([MaxHS(i),MaxRS(i),MaxFS(i),MaxMS(i),MaxHU(i),MaxRU(i),MaxFU(i),MaxMU(i)]);
end

end
end

function[MaxC] = Complete(MaxToC,SIZE)
MaxC = nan(1,SIZE);
if ~isempty(MaxToC)
    MaxC=MaxToC;
end
end
