name1=strcat('Lastrev1_OrRevW.mat');
load(name1)

clear c
load('classlist.mat');

nnr=nn(a);
[new_data] = parId(new_data,a,nnr);


for k=1:length(new_data.par)
 [new_data] = labels(new_data,c,k);
 end
 clear tt i events licks laser info spikes map params spM spikes spike_region j k new_spM
