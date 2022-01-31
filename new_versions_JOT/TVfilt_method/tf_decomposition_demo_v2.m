%
%
% DEMO - active sleep, term infant
%
%

addpath(genpath('\\fileserver02\staff\nathanSt\GitHub\tfd_decomp\'))

load test_signal_eeg1; N = 2048; wl = 63
%load bat1; bat2 = [zeros(56,1); bat ; zeros(56,1)]; N = 512; wl = 63;

ref1 = N/2+1; 
am2 = hanning(wl)*hanning(2*wl+1)';
am_sep = zeros(N, N);
am_sep(ref1-floor(wl/2):ref1+floor(wl/2), ref1-floor(wl):ref1+floor(wl)) = am2;

amb = ambig(diff(eeg1)); % Using diff to whiten spectrum, not sure if necessary
%amb = ambig(bat2');
tf_tau = ifft(ifftshift((amb.*am_sep'),2),[],2); 
tfrep = real(1./N.*fft(ifftshift(tf_tau,1), N, 1));
amb = ambig(eeg1(1:N)); % Using diff to whiten spectrum, not sure if necessary
tf_tau = ifft(ifftshift((amb.*am_sep'),2),[],2); 
tfrep1 = real(1./N.*fft(ifftshift(tf_tau,1), N, 1));
tfrep1 = tfrep;

a = size(tfrep);
tfrep(tfrep<0) = 0;
im_bin1 = zeros(a); im_bin2 = zeros(a);
for kk = 1:a(2)
    q1 = tfrep(:,kk); q3 = tfrep(kk,:);
    u = diff(q1); v = diff(q3);
    uu = q1; vv = q3;
    zc_max1 =[]; zc_max2 =[];
    count_max1 = 1; count_max2 = 1;
    for zz = 1:(length(u)-2)
        if u(zz)>=0
           if u(zz+1)<0
            zc_max1(count_max1) = zz;
            count_max1 = count_max1+1;
           end
        end
        if v(zz)>=0
           if v(zz+1)<0
            zc_max2(count_max2) = zz;
            count_max2 = count_max2+1;
           end
        end
    end
    q2 = zeros(1,a(1));
    q2(zc_max1+1) = 1;
    im_bin1(kk,:) = q2;
    q4 = zeros(1,a(1));
    q4(zc_max2+1) = 1;
    im_bin2(:,kk) = q4;

    clear q* z* c* v* ref
end
im_bin = im_bin1; %+im_bin2;
im_bin(im_bin>0)=1;

% % Search binary image for nonstationary components.
%bw = 31; fw = 31;
bw = 5; fw = 5;
[el1,ei1] = edge_link_new(im_bin, tfrep', 210, bw, fw);
%[el1, ei1] = edge_link(im_bin1, 63, sch, 2*bw);

% New version tries to get better components from the edge linking
% Sort components by signal energy - good b
en = zeros(1, length(el1));
for ii = 1:length(el1)
    c1 = el1{ii}; %a = zeros(1,length(c1));
    a = tfrep1(sub2ind(size(tfrep1), c1(:,1), c1(:,2)));
    %for jj = 1:length(c1)
    %   a(jj) = tfrep(c1(jj,2),c1(jj,1));
    %end        
    en(ii) = sum(a);
    len(ii) = length(c1);
end
[rnk, idx] = sort(en, 'descend');
M = 127;
alpha = (M-1)/2; 
n = 0:M-1;
decomp1 = zeros(length(idx), N); HDF = zeros(N,N/2); %HDF = tfp;
bwf = 'variable';
for jj = 1:length(idx)
    jj
    cx = el1{idx(jj)}; % take component defined by IF
    % INTERPOLATE EDGES TO ENSURE CONTINUITY
    %c3 = min(cx(:,1)):max(cx(:,1));
    %c2 = interp1(cx(:,1), cx(:,2), c3);
    c1 = zeros(length(cx), 2);
    c1(:,1) = cx(:,1); c1(:,2) = cx(:,2); 
    %c1 = round(c1);
    % SETUP FILTER CHARACTERISTICS
    ife = zeros(1,N); fc1 = ife; fc2 = ife;
    ife(c1(:,2)) = c1(:,1)/(2*N);
    rng(1) = min(c1(:,2));
    rng(2) = max(c1(:,2));
    switch bwf
        case 'fixed'
            fc1(ife>0) = ife(ife>0)-0.5/M; 
            fc2(ife>0) = ife(ife>0)+0.5/M;
            fc1(fc1<0) = 0;
            fc2(fc2>0.5) = 0.5;
        case 'variable'
            for qq = 1:length(c1)
                dum = ei1(:,c1(qq,2));
                ref = find(dum==1)-c1(qq,1);
                if isempty(find(ref<0,1,'last'))==0
                    fc1(c1(qq,2)) = ife(c1(qq,2))+floor(ref(find(ref<0,1,'last'))/2)/(2*N);
                else
                    fc1(c1(qq,2)) = 0;
                end
                if isempty(find(ref>0,1,'first'))==0
                    fc2(c1(qq,2)) = ife(c1(qq,2))+floor(ref(find(ref>0,1,'first'))/2)/(2*N);
                else
                    fc2(c1(qq,2)) = 0.5;
                end
            end
    end
    hd = zeros(N, M); % filter response, time domain
    hdf = zeros(N, N/2); % filter response, frequency domain
    % time-varying filter is simple windowed, bandpass FIR at the moment
    for ii = rng(1):rng(2)
        if mod(M,2)==1
            hd(ii,n~=alpha) = 1./(pi.*(n(n~=alpha)-alpha)).*(sin(2*pi*fc2(ii)*(n(n~=alpha)-alpha))-sin(2*pi*fc1(ii)*(n(n~=alpha)-alpha)));  % n~=alpha M = odd
            hd(ii,n==alpha) = 2*(fc2(ii)-fc1(ii));
        else
            hd(ii,:) = 1./(pi.*(n-alpha)).*(sin(2*pi*fc2(ii)*(n-alpha))-sin(2*pi*fc1(ii)*(n-alpha)));  % n~=alpha M = odd
        end
        hd(ii,:) = hd(ii,:).*hamming(M)';
        %hd(ii,:) = hd(ii,:)./max(hd(ii,:));
        hdum = abs(fft(hd(ii,:), N));
        hdum(hdum<max(hdum)/2)=0; hdum(hdum>0)=1;
        hdf(ii,:) = hdum(1:N/2);
    end
    
    HDF = HDF+hdf;  % Span of the TF domain by the filter bank
    filtered = zeros(N, N); 
    % DO FILTERING (using filtfilt for no actual reason)
    for ii = 1:N
        filtered(ii,:) = filter(hd(ii,:), 1, eeg1(1:N)-mean(eeg1(1:N)));
         %filtered(ii,:) = filtfilt(hd(ii,:), 1, bat2');
    end
    decomp1(jj,:) = diag(filtered);  
end

A = size(decomp1); figure; hold on;
for ii = 1:8
    plot(decomp1(ii,:)+5*ii)
end

% REDO WITH OTHER SIGNAL
%load bat1
