function decomp1 = tv_filtering(if_tracks, tfd_mask, signal1, M, bwf, N_components)
% This function applies a time-varying filter to a signal. The method is
% currently based on FIR filters (window method), IIR versions should be
% realtively trivial to implement
% 
% Inputs:
%       if_tracks - the IF of the time varying fitler a two column vector with
%       time and frequency - cell array of size L
%       signal1 - the signal to be filtered - size N
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
if(nargin < 7 || isempty(N_components)), N_components = length(if_tracks); end


N_x = length(signal1);
% zero-pad signal:
signal1 = [signal1 zeros(1, M+1)];
N = length(signal1);
N_freq = size(tfd_mask, 1);

% Do TV filtering
alpha = (M-1)/2; 
n = 0:M-1;
decomp1 = zeros(N_components, N); %HDF = zeros(N,N/2); %HDF = tfp;
                                  % disp(N_components);

for jj = 1:N_components
    cx = if_tracks{jj}; % take component defined by IF
                        % DEFINE IF LAW for component jj
    c1 = zeros(length(cx), 2);
    c1(:,1) = cx(:,1); c1(:,2) = cx(:,2); 
    % SETUP FILTER CHARACTERISTICS
    ife = zeros(1,N_freq); fc1 = ife; fc2 = ife;
    ife(c1(:,2)) = c1(:,1)/(2*N_freq);
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
            dum = tfd_mask(:,c1(qq,2));
            ref = find(dum==1)-c1(qq,1);
            if isempty(find(ref<0,1,'last'))==0
                fc1(c1(qq,2)) = ife(c1(qq,2))+floor(ref(find(ref<0,1,'last'))/2)/(2*N_freq);
            else
                fc1(c1(qq,2)) = 0;
            end
            if isempty(find(ref>0,1,'first'))==0
                fc2(c1(qq,2)) = ife(c1(qq,2))+floor(ref(find(ref>0,1,'first'))/2)/(2*N_freq);
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
    end
    

    % DO FILTERING 
    for ii = rng(1):rng(2)
        filt_signal1 = filter(hd(ii,:), 1, signal1);
        decomp1(jj, ii) = filt_signal1(ii);
    end

    % for ii = rng(1):rng(2)
    %         % Element-wise multiplication using '.*' and then summing the result
    %     filtered_signal = sum(signal1(ii:ii+M-1) .* hd(ii, :));
    %     decomp1(jj, ii) = filtered_signal;  % Assign the result directly
    % end

end
