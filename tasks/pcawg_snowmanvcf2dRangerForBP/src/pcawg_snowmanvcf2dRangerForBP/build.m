addpath ../../../../matlab;
addpath ../../../../matlab/seq;
addpath ../../../../matlab/mike;
disp('Starting compilation...')
mcc -mC pcawg_snowmanvcf2dRangerForBP
disp('Finished compilation.')
quit;
