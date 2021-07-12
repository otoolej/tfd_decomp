

addpath(genpath('\\fileserver02\staff\nathanSt\GitHub\tfd_decomp\'))

load test_signal_eeg1  % active sleep in an infant
wl = 63; 
flag1 = 1;
len1 = 128;
bw = 15; fw = 15;
M = 127; bwf = 'lossless';
decomp1 = tf_decomposition(eeg1-mean(eeg1), wl, flag1, bw, fw, len1, M, bwf);
decomp1 = decomp1(:,ceil(M/2):end); % correct for time delay of FIR filters

load bat1; bat2 = [zeros(56,1); bat ; zeros(56,1)]; 
wl = 127; 
flag1 = 0;
len1 = 200;
bw = 5; fw = 5;
M = 63; bwf = 'lossless';
decomp2 = tf_decomposition(bat2'-mean(bat2), wl, flag1, bw, fw, len1, M, bwf);
decomp2 = decomp2(:,ceil(M/2):end); % correct for time delay of FIR filters

% these two take ages as my TFDs will be N x N due to laziness
load whale1; whale1a = whale(1:4096);
wl = 81; 
flag1 = 0;
len1 = 511;
bw = 65; fw = 31;
M = 255; bwf = 'lossless';
decomp3a = tf_decomposition(whale1a'-mean(whale1a), wl, flag1, bw, fw, len1, M, bwf);
decomp3a = decomp3a(:,ceil(M/2):end); % correct for time delay of FIR filters

load whale1; whale1b = whale(end-4095:end);
wl = 81; 
flag1 = 0;
len1 = 511;
bw = 65; fw = 31;
M = 255; bwf = 'lossless';
decomp3b = tf_decomposition(whale1b'-mean(whale1b), wl, flag1, bw, fw, len1, M, bwf);
decomp3b = decomp3b(:,ceil(M/2):end); % correct for time delay of FIR filters
