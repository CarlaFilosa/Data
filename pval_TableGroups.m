function [PvTableClass] = pval_TableGroups(PvTable_Multi,colpval)
%%%%% The function divides the table p value (classified for animal and units) 
%%%%%% in subgroups accordigly to the neurons typologies

 for i=1:13
    PvTableClass{i}=[];
 end

% if PvTable_Multi(1,9)==neu1 && PvTable_Multi(1,10)==neu2
 for i=1:size(PvTable_Multi,1)

        %%%%%%%% Baseline vs After Stimulus
    if PvTable_Multi(i,1)==1 && PvTable_Multi(i,2)==2  && PvTable_Multi(i,colpval)<=0.05 && PvTable_Multi(i+1,colpval)>0.05 && PvTable_Multi(i+2,colpval)>0.05 
    
        PvTableClass{1}=PvTable_Multi;
        PvTableClass{1}(:,end+1)=12;
    end
        %%%%%%%% Baseline vs PreReward
    if PvTable_Multi(i,1)==1 && PvTable_Multi(i,2)==3  && PvTable_Multi(i,colpval)<=0.05 && PvTable_Multi(i+1,colpval)>0.05 && PvTable_Multi(i-1,colpval)>0.05
    
        PvTableClass{2}=PvTable_Multi;
        PvTableClass{2}(:,end+1)=13;
    end
      %%%%%%%% Baseline vs After Reward
    if PvTable_Multi(i,1)==1 && PvTable_Multi(i,2)==4  && PvTable_Multi(i,colpval)<=0.05 && PvTable_Multi(i-2,colpval)>0.05 && PvTable_Multi(i-1,colpval)>0.05 
   
        PvTableClass{3}=PvTable_Multi;
        PvTableClass{3}(:,end+1)=14;
    end
    
    %%%%%%%% Baseline vs After Stimulus and Pre Reward
    if PvTable_Multi(i,1)==1 && PvTable_Multi(i,2)==2  && PvTable_Multi(i,colpval)<=0.05 && PvTable_Multi(i+1,colpval)<=0.05 && PvTable_Multi(i+2,colpval)>0.05 
   
        PvTableClass{4}=PvTable_Multi;
        PvTableClass{4}(:,end+1)=123;
    end
    
    %%%%%%%% Baseline vs After Stimulus and After Reward
    if PvTable_Multi(i,1)==1 && PvTable_Multi(i,2)==2  && PvTable_Multi(i,colpval)<=0.05 && PvTable_Multi(i+2,colpval)<=0.05 && PvTable_Multi(i+1,colpval)>0.05 
    
        PvTableClass{5}=PvTable_Multi;
        PvTableClass{5}(:,end+1)=124;
    end
    
    %%%%%%%% Baseline vs PreReward and After Reward
    if PvTable_Multi(i,1)==1 && PvTable_Multi(i,2)==3  && PvTable_Multi(i,colpval)<=0.05 && PvTable_Multi(i+1,colpval)<=0.05 && PvTable_Multi(i-1,colpval)>0.05 
   
        PvTableClass{6}=PvTable_Multi;
        PvTableClass{6}(:,end+1)=134;
    end
    
    %%%%%%%% Baseline vs After Stimulus, Pre Reward and After Reward
    if PvTable_Multi(i,1)==1 && PvTable_Multi(i,2)==2  && PvTable_Multi(i,colpval)<=0.05 && PvTable_Multi(i+1,colpval)<=0.05 && PvTable_Multi(i+2,colpval)<=0.05  
    
        PvTableClass{7}=PvTable_Multi;
        PvTableClass{7}(:,end+1)=1234;
    end
 end
 
 for i=1:size(PvTable_Multi,1)
      %%%%%%%%%  Post Stim/ PreReward
   if isempty(PvTableClass{1}) && isempty(PvTableClass{2}) && isempty(PvTableClass{3}) && isempty(PvTableClass{4}) && isempty(PvTableClass{5}) && isempty(PvTableClass{6}) && isempty(PvTableClass{7})
 
     if PvTable_Multi(i,1)==2 && PvTable_Multi(i,2)==3 && PvTable_Multi(i,colpval)<=0.05 && PvTable_Multi(i+1,colpval)>0.05 && PvTable_Multi(i+2,colpval)>0.05
        PvTableClass{8}=PvTable_Multi;
        PvTableClass{8}(:,end+1)=23;
     end
     %%%%%%%% Post Stim/ After reward
     if PvTable_Multi(i,1)==2 && PvTable_Multi(i,2)==4 && PvTable_Multi(i,colpval)<=0.05 && PvTable_Multi(i+1,colpval)>0.05 && PvTable_Multi(i-1,colpval)>0.05
        PvTableClass{9}=PvTable_Multi;
        PvTableClass{9}(:,end+1)=24;
     end
     %%%%%%%% Pre Reward/ After reward
     if PvTable_Multi(i,1)==3 && PvTable_Multi(i,2)==4 && PvTable_Multi(i,colpval)<=0.05 && PvTable_Multi(i-2,colpval)>0.05 && PvTable_Multi(i-1,colpval)>0.05
        PvTableClass{10}=PvTable_Multi;
        PvTableClass{10}(:,end+1)=34;
     end
    %%%%%% After stimulus vs pre reward and after reward
    if PvTable_Multi(i,1)==2 && PvTable_Multi(i,2)==3 && PvTable_Multi(i,colpval)<=0.05 && PvTable_Multi(i+1,colpval)<=0.05 && PvTable_Multi(i+2,colpval)>0.05
        PvTableClass{11}=PvTable_Multi;
        PvTableClass{11}(:,end+1)=234;
    end
    %%%%%%% Pre reward vs After stimulus and after reward
    if PvTable_Multi(i,1)==2 && PvTable_Multi(i,2)==3 && PvTable_Multi(i,colpval)<=0.05 && PvTable_Multi(i+2,colpval)<=0.05 && PvTable_Multi(i+1,colpval)>0.05
        PvTableClass{12}=PvTable_Multi;
        PvTableClass{12}(:,end+1)=324;
    end
    
    %%%% After reward vs sfter stimulus and pre reward
    if PvTable_Multi(i,1)==2 && PvTable_Multi(i,2)==4 && PvTable_Multi(i,colpval)<=0.05 && PvTable_Multi(i+1,colpval)<=0.05 && PvTable_Multi(i-1,colpval)>0.05
        PvTableClass{13}=PvTable_Multi;
        PvTableClass{13}(:,end+1)=423;
    end
    
   end
 end
% end

end

