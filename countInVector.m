Aus=sort(Aus)
Conta = zeros(1,length(Aus))
CopyAus=Aus
for i=1: length(Aus)
    for j=1:length(Aus)
   if ismember(Aus(i),Aus(j))
       Conta(i)=Conta(i)+1;
   end
    end
end

Conta1=unique(Conta)
sum(Conta>1)

ll=0;
IN=[];
Conta(:,end+1)=nan;
for i=1:length(Conta)-1
    if ~ismember(Conta(i),Conta(i+1))
        ll=ll+1;
        IN(ll)=i
    end
    
    
end