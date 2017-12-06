addpath ./;
addpath ../../../../matlab;
addpath ../../../../matlab/seq;
addpath ../../../../matlab/Chip;
disp('Starting compilation...')
mcc -mC SV_cluster_forBP
disp('Finished compilation.')
quit;
