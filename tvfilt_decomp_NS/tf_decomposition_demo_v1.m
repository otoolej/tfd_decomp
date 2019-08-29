%
%
% DEMO - active sleep, term infant
%
%

load test_signal_eeg1

N = 2048;
wl = 63; ref1 = N/2+1; 
am2 = hanning(wl)*hanning(2*wl+1)';
am_sep = zeros(N, N);
am_sep(ref1-floor(wl/2):ref1+floor(wl/2), ref1-floor(wl):ref1+floor(wl)) = am2;

amb = ambig(diff(eeg1)); % Using diff to whiten spectrum, not sure if necessary
tf_tau = ifft(ifftshift((amb.*am_sep'),2),[],2); 
tfrep = real(1./N.*fft(ifftshift(tf_tau,1), N, 1));

a = size(tfrep);
tfrep(tfrep<0) = 0;
im_bin1 = zeros(a);
for kk = 1:a(2)
    q1 = tfrep(:,kk);
    u = diff(q1);
    uu = q1;
    zc_max =[];
    count_max = 1;
    for zz = 1:(length(u)-2)
        if u(zz)>=0
           if u(zz+1)<0
            zc_max(count_max) = zz;
            count_max = count_max+1;
        end
    end
    end
    q2 = zeros(1,a(1));
    q2(zc_max+1) = 1;
    im_bin1(kk,:) = q2;
    clear q* z* c* ref
end
 
% % Set up the search region based on the size of the separable kernel.
ref = [0:32 -1:-1:-32];
[val, idx] = sort(abs(ref));
nref = ref(idx); sch = [];
% BIG SEARCH SPACE MINIMIZES COMPONENTS
for ii = 1:32
sch = [sch [ii*ones(1,length(nref)) ; nref]];
end
% % Search binary image for nonstationary components.
bw = 103;
[el1, ei1] = edge_link(im_bin1, 63, sch, 2*bw);

% Energy constraint
en = zeros(1, length(el1));
for ii = 1:length(el1)
    c1 = el1{ii}; a = zeros(1,length(c1));
    for jj = 1:length(c1)
       a(jj) = tfrep(c1(jj,2),c1(jj,1));
    end        
    en(ii) = sum(a);
    len(ii) = length(c1);
end
[rnk, idx] = sort(en, 'descend');
M = 63;
alpha = (M-1)/2; 
n = 0:M-1;
decomp1 = zeros(length(idx), N); tfp = zeros(N,N); HDF = tfp;
for jj = 1:length(idx)
    jj
    cx = el1{idx(jj)}; % take component defined by IF
    % INTERPOLATE EDGES TO ENSURE CONTINUITY
    c3 = min(cx(:,1)):max(cx(:,1));
    c2 = interp1(cx(:,1), cx(:,2), c3);
    c1 = zeros(length(c2), 2);
    c1(:,1) = c3'; c1(:,2) = c2'; 
    c1 = round(c1);
    % SETUP FILTER CHARACTERISTICS
    ife = zeros(1,N); fc1 = ife; fc2 = ife;
    ife(c1(:,1)) = c1(:,2)/(2*N);
    rng(1) = min(c1(:,1));
    rng(2) = max(c1(:,1));
    fc1(ife>0) = ife(ife>0)-0.5/16; % This bandwidth is related to the underlying filter order - need more sophisticated definition here but is doable
    fc2(ife>0) = ife(ife>0)+0.5/16;
    fc1(fc1<0) = 0;
    fc2(fc2>0.5) = 0.5;
    hd = zeros(N, M); % filter response, time domain
    hdf = zeros(N, N); % filter response, frequency domain
    % time-varying filter is simple windowed, bandpass FIR at the moment
    for ii = rng(1):rng(2)
        if mod(M,2)==1
            hd(ii,n~=alpha) = 1./(pi.*(n(n~=alpha)-alpha)).*(sin(2*pi*fc2(ii)*(n(n~=alpha)-alpha))-sin(2*pi*fc1(ii)*(n(n~=alpha)-alpha)));  % n~=alpha M = odd
            hd(ii,n==alpha) = 2*(fc2(ii)-fc1(ii));
        else
            hd(ii,:) = 1./(pi.*(n-alpha)).*(sin(2*pi*fc2(ii)*(n-alpha))-sin(2*pi*fc1(ii)*(n-alpha)));  % n~=alpha M = odd
        end
        hd(ii,:) = hd(ii,:).*hamming(M)';
        hdf(ii,:) = abs(fft(hd(ii,:), N));
        hd(ii,:) = hd(ii,:)./max(hdf(ii,:));
    end
    HDF = HDF+hdf;  % Span of the TF domain by the filter bank
    filtered = zeros(N, N); 
    % DO FILTERING (using filtfilt for no actual reason)
    for ii = 1:N
        filtered(ii,:) = filtfilt(hd(ii,:), 1, eeg1(1:N));
    end
    decomp1(jj,:) = diag(filtered);  
end

% NOTE: in general pretty good but weird missing bit in reconstructed
% signal at around 1750 samples, will look into it, best thing about it is
% that is shows temporal segmentation, not just frequency based
% segmentation
