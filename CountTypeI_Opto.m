
Data=Data_Par;
for k=1:length(Data)
    count_typeIOpto(k)=0;
    count_typeI(k)=0;
    count_Opto(k)=0;
    for i=1:length(Data{k})
        if Data{k}(i).labels==4 
            count_typeI(k)=count_typeI(k)+1;
        end
        if Data{k}(i).labelsOp==7.1
            count_Opto(k)=count_Opto(k)+1;
        end
        if Data{k}(i).labels==4 & Data{k}(i).labelsOp==7.1
            count_typeIOpto(k)=count_typeIOpto(k)+1;
        end
    end
end
sum(count_Opto)
sum(count_typeI)
sum(count_typeIOpto)
