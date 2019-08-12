function [as_sel1New] = Comp_RewCode(as_sel1,as_sel2)
% Comparison and completing the matrices 
% in case in which there are missing trials of reward or lick
as_sel1New=as_sel1;
for k=1:length(as_sel1New)
    if~isempty(as_sel2{k}) & ~isempty(as_sel2{k}.asAct)
for i =1: length(as_sel2{k}.asAct)
    if ~isempty(as_sel2{k}.asAct{i}) & ~isempty(as_sel1{k}.asAct{i})
    for j= 1:length(as_sel2{k}.asAct{i})
if isempty(as_sel1{k}.asAct{i}{j}) & ~isempty(as_sel2{k}.asAct{i}{j})
    as_sel1New{k}.asAct{i}{j}=nan(size(as_sel2{k}.asAct{i}{j}));
    ML(k)=0;
elseif ~isempty(as_sel1{k}.asAct{i}{j})
    ML(k) = length(as_sel1{k}.asAct{i}{j});
end
    end
    end
end
    end
end
 MML= max(ML)
%   as_sel2New=as_sel1New;
for k=1:length(as_sel1New)
if ~isempty(as_sel1New{k}) & ~isempty(as_sel1New{k}.asAct)
for i= 1:length(as_sel1New{k}.asAct)
     if ~isempty(as_sel1New{k}.asAct{i})
    for j=1:length(as_sel1New{k}.asAct{i})
        if isnan(as_sel1New{k}.asAct{i}{j})== ones(size(as_sel1New{k}.asAct{i}{j}))
           as_sel1New{k}.asAct{i}{j} = nan(MML,2);
           
    end
    end
end
end
end
 end
% for k=1:length(as_sel1New)
%     for i= 1:length(as_sel1New{k})
%         if ~isempty(as_sel1New{k}.asAct{i})
%             for j=1:length(as_sel1New{k}.asAct{i})
%                 if ~isempty(as_sel2New{k}.asAct{i}{j})
%                     as_selNew{k}.asAct{i}{j}=as_sel2New{k}.asAct{i}{j}
%                 else as_selNew{k}.asAct{i}{j}= as_sel1New{k}.asAct{i}{j}
%                 end
%             end
%         end
%     end
% end
end

