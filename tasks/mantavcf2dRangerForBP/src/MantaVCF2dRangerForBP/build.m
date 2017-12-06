addpath ../../../../matlab;
addpath ../../../../matlab/seq;
addpath ../../../../matlab/mike;
disp('Starting compilation...')
mcc -mC MantaVCF2dRangerForBP
disp('Finished compilation.')
quit;
