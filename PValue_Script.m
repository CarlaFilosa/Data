for k=1:length(pmulti)
    
    PvTableMT{k}=[];
    PvTableFT{k}=[];
end

for k=1:length(pmulti)
    for i=1:length(pmulti{k})
    PvTableMT{k}{i}=[];
    PvTableFT{k}{i}=[];
    end
end
colpval=6;neu1=1;neu2=4;
for k=1:length(pmulti)
    if ~isempty(pmulti{k})
        for i=1:length(pmulti{k})
            if ~isempty(pmulti{k}{i})
                if pmulti{k}{i}(1,9)==neu1 && pmulti{k}{i}(1,10)==neu2
[PvTableMT{k}{i}] = pval_TableGroups(pmulti{k}{i},colpval);
                end
            end
        end
    end
end

colpval=6;neu1=2;neu2=4;
for k=1:length(pmulti)
    if ~isempty(pmulti{k})
        for i=1:length(pmulti{k})
            if ~isempty(pmulti{k}{i}) 
                if pmulti{k}{i}(1,9)==neu1 && pmulti{k}{i}(1,10)==neu2
[PvTableFT{k}{i}] = pval_TableGroups(pmulti{k}{i},colpval);
                end
            end
        end
    end
end

PvTable=PvTableMT;