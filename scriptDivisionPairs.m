
for k=1:length(As_activity_All)
   
for i =1: length(cvta)
    
[pairs_MT{i}{k}] = [];
[pairs_FT_low{i}{k}] = [];
[pairs_FT_high{i}{k}] = [];
[pairs_CT{i}{k}] = [];
end

end

for k=1:length(As_activity_All)
    if ~isempty(pairs_vsvta_small{k})
for i =1: length(cvta)
   
[pairs_MT{i}{k}] = pairs_NeuSelection(pairs_vsvta_small{k},cvs(1),cvta(i),col1,col2(i));
[pairs_FT_low{i}{k}] = pairs_NeuSelection(pairs_vsvta_small{k},cvs(2),cvta(i),col1,col2(i));
[pairs_FT_high{i}{k}] = pairs_NeuSelection(pairs_vsvta_small{k},cvs(3),cvta(i),col1,col2(i));
[pairs_CT{i}{k}] = pairs_NeuSelection(pairs_vsvta_small{k},cvs(4),cvta(i),col1,col2(i));
end
    end
end