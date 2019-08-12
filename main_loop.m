% Main loop for TD learning

pat = pat+1;
if pat > buffer_length
  epoch = epoch + 1;
  pat = 1;
  buffer(:) = 0;
  history(:) = 0;
  V = 0;
end

 % shift new stimulus into buffer
oldbuff = buffer;
buffer(2:end,:) = buffer(1:(end-1),:);
buffer(1,:) = stimpat(pat,:);

R = reward(pat);

 % TD learning rule
oldV = V;
V = dot(W(:),buffer(:));
delta = R + Gamma*V - oldV;
history(2:end,:) = history(1:(end-1),:);
history(1,:) = [V oldV R delta];

W(:) = W(:) + eta*delta*oldbuff(:);

  % Now do the graphic display

subplot(3,1,1)
bar(buffer,'grouped')
axis([0 1+buffer_length 0 1.5])
title('TD Learning Simulation')
ylabel('Stimuli')

subplot(3,1,2)
bar(W,'grouped')
axis([0 1+buffer_length -1.2 1.2])
ylabel('Weights')

subplot(3,1,3)
bar(history)
xlabel(sprintf('Trial %d',epoch))
ylabel('History')
% h = legend('V(t)','V(t-1)','Reward','Error',1);
% p = get(h,'Position');
p(1:2) = [0.75 0.25];
% set(h,'Position',p)
axis([0 1+buffer_length -0.2 1.2])
drawnow

if pause_every_step == 1
  pause
elseif pat==buffer_length
  pause(1)
end
