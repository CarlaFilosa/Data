function [ spM ] = LicksRewIntEq( spM,TotAn,prepost )
% Find licks Interval Time
for k=1:TotAn

    for l=1:length(spM{k})
       lu= length(spM{k}(l).licks);
       lr = length(spM{k}(l).reward_time);
      
    if lu~=0

      LickDiff=0.3;

      spM{k}(l).lickPre(1:lu)=spM{k}(l).licks(1:lu)-LickDiff/2;
      spM{k}(l).lickPost(1:lu)=spM{k}(l).licks(1:lu)+LickDiff/2;
          

          spM{k}(l).lickFirstPre = spM{k}(l).licks(1)-prepost;
          spM{k}(l).lickLastPost = spM{k}(l).licks(end)+prepost; 

       end
      %%% Reward Intervals
       if lr==1
           spM{k}(l).rewardPre = spM{k}(l).reward_time-prepost;
           spM{k}(l).rewardPost = spM{k}(l).reward_time+prepost;
       end
       
    end
    
 end


end