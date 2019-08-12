% Simulation of TD (Temporal Difference) learning in the Schultz task.
%
% David S. Touretzky.  
% October, 1998

buffer_length = 8;
num_stimuli = 3;

buffer = zeros(buffer_length,num_stimuli);
reward = zeros(buffer_length,1);
history = zeros(buffer_length,4);  % [V, oldV, R, delta]

W = buffer;
V = 0;
Gamma = 0.9;
eta = 0.25;

% set up stimulus and reward patterns

stimpat = buffer;
stimpat(1,1) = 1;
stimpat(4,2) = 1;
stimpat(3,3) = 1;
reward(5)=1;
% reward(5:buffer_length:(buffer_length*Ntrials)/2) = 1;
% reward((buffer_length*Ntrials)/2+8:buffer_length:end) = 1;


epoch = 0;
pat = Inf;

figure(1), clf, whitebg(gcf,[0 0 0])

if ~exist('pause_every_step'), pause_every_step = 0; end

if pause_every_step
  disp('Hit the space bar to advance each timestep.')
  disp('Set pause_every_step=0 to run automatically.')
else
  disp('Set pause_every_step=1 to manually step through the simulation.')
end

if ~exist('Ntrials'), Ntrials=20; end

for i=1:(Ntrials*buffer_length)
  main_loop
end
