function [count] = CountTypeNeuPair(struct,k,count)
% Function to count the tyologies of neurons
% struct=struct_pair
 if~isempty(struct.pair)

    for i=1:length(struct.labels)
        if struct.labels{i}(1,1) == 1 | struct.labels{i}(1,2) == 1
            disp('ciao')
                  count(k,1)=count(k,1)+1;
        elseif struct.labels{i}(1,1) == 2 | struct.labels{i}(1,2) == 2
            disp('ciaciao')
                  count(k,2)=count(k,2)+1;
        elseif struct.labels{i}(1,1) == 3 | struct.labels{i}(1,2) == 3
                  count(k,3)=count(k,3)+1;
        elseif struct.labels{i}(1,1) == 4 | struct.labels{i}(1,2) == 4
                  count(k,4)=count(k,4)+1;
        elseif struct.labels{i}(1,1) == 5 | struct.labels{i}(1,2) == 5
                  count(k,5)=count(k,5)+1;
        elseif struct.labels{i}(1,1) == 6 | struct.labels{i}(1,2) == 6
                  count(k,6)=count(k,6)+1;
        elseif struct.labels{i}(1,1) == 0 | struct.labels{i}(1,2) == 0
                  count(k,7)=count(k,7)+1;
        end
    end
 end
end

