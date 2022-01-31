

addpath(genpath('\\fileserver02\staff\nathanSt\GitHub\tfd_decomp\'))

load test_signal_eeg1  % active sleep in an infant
flag = 1;
len1 = 128;
decomp1 = tf_decomposition_v1(eeg1(1:2048), len1, flag); 


%decomp1 = tf_decomposition(eeg1(1:2048), wl, flag1, bw, fw, len1, M, bwf);
decomp1 = decomp1(:,ceil(M/2):end); % correct for time delay of FIR filters


load bat1; bat2 = [zeros(56,1); bat-mean(bat) ; zeros(56,1)]; 
flag = 2; len1 = 100;
% wl = 127; 
% flag1 = 0;
% len1 = 200;
% bw = 5; fw = 5;
% M = 63; bwf = 'lossless';
% len1, wx, bw, fw, M

decomp2 = tf_decomposition_v1(bat2', len1, flag); 
% flag is a parameter that gives 3 different TFD types - 
% 1 - track long duration components with slowly varying IF law
% 2 - track long duration components with rapidly varying IF law
% 3 - track transients/short duration components 
% this function also uses a non-zero threshold for signal energy - it is half the
% total TF energy divided by the discrete size of the TF plane

%decomp2 = tf_decomposition(bat2', wl, flag1, bw, fw, len1, M, bwf);
%decomp2 = decomp2(:,ceil(M/2):end); % correct for time delay of FIR filters
imfs1 = emd(bat2); %imfs1 = imfs1(il(M/2):end,:);
imfs2 = vmd(bat2); %imfs2 = imfs2(ceil(M/2):end,:);

%wsst(quadchirp,Fs,'bump')
figure; set(gcf, 'Position', [100 100 1800 600])
stp = [0.25 0.5 0.7 0.9 1.1 1.3];
A = size(decomp2);
subplot(1,3,1); hold on
for ii = 1:5
    plot(decomp2(ii,:)-stp(ii))
end
axis([1 512 -1.25 -0.05])
set(gca, 'Ytick', [-1.1 -0.9 -0.7 -0.5 -0.25], 'Yticklabel', [5 4 3 2 1], 'Position', [0.075 0.125 0.25 0.85], 'fontsize', 16)
ylabel('My Method'); xlabel('Samples')

A = size(imfs1);
subplot(1,3,2); hold on
for ii = 1:5
    plot(imfs1(:,ii)-stp(ii))
end
axis([1 512 -1.25 -0.05])
set(gca, 'Ytick', [-1.1 -0.9 -0.7 -0.5 -0.25], 'Yticklabel', [5 4 3 2 1], 'Position', [0.4 0.125 0.25 0.85], 'fontsize', 16)
ylabel('EMD'); xlabel('Samples')

A = size(imfs2);
subplot(1,3,3); hold on
for ii = 1:5
    plot(imfs2(:,ii)-stp(ii))
end
axis([1 512 -1.25 -0.05])
set(gca, 'Ytick', [-1.1 -0.9 -0.7 -0.5 -0.25], 'Yticklabel', [5 4 3 2 1], 'Position', [0.725 0.125 0.25 0.85], 'fontsize', 16)
ylabel('VMD'); xlabel('Samples')

figure; set(gcf, 'Position', [100 100 1800 600])
subplot(1,3,2); hht(imfs1)
set(gca, 'Position', [0.4 0.125 0.25 0.8], 'fontsize', 16)
title('EMD'); 
subplot(1,3,3); hht(imfs2);
set(gca, 'Position', [0.725 0.125 0.25 0.8], 'fontsize', 16)
title('VMD'); 
subplot(1,3,1); hht(decomp2'); 
set(gca, 'Position', [0.075 0.125 0.25 0.8], 'fontsize', 16)
title('My Method'); 


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
