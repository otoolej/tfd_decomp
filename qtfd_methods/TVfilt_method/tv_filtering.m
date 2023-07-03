function decomp1 = tv_filtering(el1, ei1, signal1, tfrep1, M, bwf, N_components)
% This function applies a time-varying filter to a signal. The method is
% currently based on FIR filters (window method), IIR versions should be
% realtively trivial to implement
% 
% Inputs:
%       el1 - the IF of the time varying fitler a two column vector with
%       time and frequency - cell array of size L
%       signal1 - the signal to be filtered - size N
%       tfrep1 - the TFD of the signal to be filtered
%       M - the order of the filter - must be large enough to generate a
%       narrow enough band to extract desired components
%       bwf - filtering method either 'lossy' or lossless
%       N_components - only decompose this many components
%
% Outputs: 
%       decomp1 - a matrix containing the components extracted/filtered
%       from the signal (N x L)
%
% Nathan Stevenson
% July 2021
if(nargin < 7 || isempty(N_components)), N_components = length(el1); end


N = length(signal1);
% Sort to ensure largest TF energy components are extracted first
en = zeros(1, length(el1));
for ii = 1:length(el1)
    c1 = el1{ii}; %a = zeros(1,length(c1));
    a = tfrep1(sub2ind(size(tfrep1), c1(:,1), c1(:,2)));
    en(ii) = sum(a);
end
[~, idx] = sort(en, 'descend');

if(length(idx) < N_components)
    N_components = length(idx);
end


% Do TV filtering
alpha = (M-1)/2; 
n = 0:M-1;
decomp1 = zeros(N_components, N); %HDF = zeros(N,N/2); %HDF = tfp;
disp('Decomposing Signal'); aa_old = 11;

for jj = 1:N_components
    cx = el1{idx(jj)}; % take component defined by IF
    % DEFINE IF LAW for component jj
    c1 = zeros(length(cx), 2);
    c1(:,1) = cx(:,1); c1(:,2) = cx(:,2); 
    % SETUP FILTER CHARACTERISTICS
    ife = zeros(1,N); fc1 = ife; fc2 = ife;
    ife(c1(:,2)) = c1(:,1)/(2*N);
    rng(1) = min(c1(:,2));
    rng(2) = max(c1(:,2));
    switch bwf
        case 'lossy'
            fc1(ife>0) = ife(ife>0)-0.5/M; 
            fc2(ife>0) = ife(ife>0)+0.5/M;
            fc1(fc1<0) = 0;
            fc2(fc2>0.5) = 0.5;
        case 'lossless'
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
    %hdf = zeros(N, N/2); % filter response, frequency domain
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
     %   hdum = abs(fft(hd(ii,:), N));
     %   hdum(hdum<max(hdum)/2)=0; hdum(hdum>0)=1;
     %   hdf(ii,:) = hdum(1:N/2);
    end
    
    %HDF = HDF+hdf;  % Span of the TF domain by the filter bank
    filtered = zeros(N, N); 
    % DO FILTERING (using filtfilt for no actual reason if you wish)
    for ii = 1:N
        filtered(ii,:) = filter(hd(ii,:), 1, signal1);
    end
    decomp1(jj,:) = diag(filtered);  
    
    % aa_new = round(10*jj/length(idx)); str = [];
    % if aa_old-aa_new~=0
    %     for ii = 1:aa_new; str  = [str '*']; end
    %     aa_old = aa_new;
    %     disp(str)
    % end        
    
end
