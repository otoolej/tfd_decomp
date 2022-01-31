
% FIGURE BITS AND PIECES RE: V FILTERING METHOD

load bat1; bat2 = [zeros(56,1); detrend(bat) ; zeros(56,1)]; 
% NOTE BAT SIGNAL SEEMS TO HAVE A ERALLY SLOW TREND SO USE DETREND RATHER
% THAN SUBTRACT MEAN AND TENDS TO GIVE BETTER RESULTS
len1 = 100;

N = length(bat2); M = round(sqrt(N));

wx = M; wx = wx+(1-rem(wx,2));
% ESTIMATE TFD
tfrep = estimate_tfd(bat2', wx);
% EXTRACT COMPONENTS
bw = N/wx/4; fw = N/wx/2;
bw = bw+(1-rem(bw,2)); fw = fw+(1-rem(fw,2));
[el1, ei1] = find_components_v1(tfrep, tfrep, 2*bw, 2*fw, len1);

%subplot(2,2,2); 

% DO FILTERING
M = ceil(2*M); M = M+(1-rem(M,2));        
%N = length(signal1);
% Sort to ensure largest TF energy components are extracted first
idx = [3 2 1]
%M = 127;
% Do TV filtering
alpha = (M-1)/2; 
n = 0:M-1;
decomp1 = zeros(length(idx), N); %HDF = zeros(N,N/2); %HDF = tfp;
%disp('Decomposing Signal'); 
aa_old = 11;
%bwf = 'variable';
%h1 = figure; h2 = f
cmap = [45 168 216 ; 75 61 143 ; 55 169 135]/256; figure; hold on;
fc = cell(3,2);
for jj = 1:length(idx)
%jj=1
    cx = el1{idx(jj)}; % take component defined by IF
    % DEFINE IF LAW for component jj
    c1 = zeros(length(cx), 2);
    c1(:,1) = cx(:,1); c1(:,2) = cx(:,2); 
    % SETUP FILTER CHARACTERISTICS
    ife = zeros(1,N); fc1 = ife; fc2 = ife;
    ife(c1(:,2)) = c1(:,1)/(2*N);
    rng(1) = min(c1(:,2));
    rng(2) = max(c1(:,2));
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
    
    if jj == 2; hd_save = hd; end
    
    fc{jj, 1} = fc1;
    fc{jj, 2} = fc2;
%    fill([fc1(1) fc2 fc2(N) fc1(N:-1:1)], [1 1:N N N:-1:1], cmap(jj,:), 'EdgeColor', cmap(jj,:), 'EdgeAlpha', 0)
    
    %HDF = HDF+hdf;  % Span of the TF domain by the filter bank
    filtered = zeros(N, N); 
    % DO FILTERING (using filtfilt for no actual reason if you wish)
    for ii = 1:N
        filtered(ii,:) = filter(hd(ii,:), 1, bat2);
    end
    decomp1(jj,:) = diag(filtered);  
    
    if jj == 2; filtered_save = filtered; end
    
%     aa_new = round(10*jj/length(idx)); str = [];
%     if aa_old-aa_new~=0
%         for ii = 1:aa_new; str  = [str '*']; end
%         aa_old = aa_new;
%         disp(str)
%     end        
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PLOT STUFF - colour coded
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CONTOUR PLOT OF TFREP COLOUR CODED TO COMPONENT NUMBER
figure; hold on;
%subplot(2,2,1)
n = 1:N;
slope = [-1.67 -1.47 -1.36 0];
icept = [250 550 750 500];
for ii = 1:length(idx)
    tdum = zeros(N);
    f1 = slope(ii)*n+icept(ii);
    f2 = slope(ii+1)*n+icept(ii+1);
    %plot(f1+ii, 'color', cmap(ii,:)); %plot(f2+ii, 'color', cmap(ii,:))
    for jj = 1:N
       r1 = ceil(f1(jj)); r2 = floor(f2(jj));
       if r1<1; r1 = 1; end
       if r2>N; r2 = N; end
       %r3(jj) = r1; r4(jj) = r2;
       tdum(jj,r1:r2) = tfrep(jj, r1:r2);  
    end
    [~,h] = contour(tdum',[max(max(tfrep))./(2.^[5:-0.5:1])]);
    h.LineColor = cmap(ii,:);
end
xlabel('normalized frequency'); ylabel('time');
set(gca, 'fontsize', 14, 'Xtick', [1 128 256 384 512], 'Xticklabel', [0 0.125 0.25 0.375 0.5])
axis([1 N 1 N])

% CONTOUR PLOT IF LAWS EXTRACTED (COMPLETE REPRESENTATION - i.e. extension in time on selected components)
figure; hold on;
for ii = 1:length(el1)
    plot(el1{ii}(:,1), el1{ii}(:,2), 'linewidth', 2, 'color', cmap(ii,:))
end
xlabel('normalized frequency'); ylabel('time');
set(gca, 'fontsize', 14, 'Xtick', [1 128 256 384 512], 'Xticklabel', [0 0.125 0.25 0.375 0.5])
axis([1 N 1 N])

% CONTOUR PLOT OF REGION OF INTEREST OF TV FILTERS 
figure; hold on;
for jj = 1:length(idx)
    fill(round([fc{jj,1}(1) fc{jj,2} fc{jj,2}(N) fc{jj,1}(N:-1:1)]*2*N), [1 1:N N N:-1:1], cmap(idx(jj),:), 'EdgeColor', cmap(idx(jj),:), 'EdgeAlpha', 0)
end
xlabel('normalized frequency'); ylabel('time');
set(gca, 'fontsize', 14, 'Xtick', [1 128 256 384 512], 'Xticklabel', [0 0.125 0.25 0.375 0.5])
axis([1 N 1 N])

% PLOT OF TV FILTERS - just for component 1 in terms of signal energy (the middle one in TF)
figure; 
set(gcf, 'Position', [600 400 1200 400])
subplot(3,3, [1 4 7]);
cc = colormap('gray')
colormap(1-cc);
imagesc(hd_save)
axis([1 M 1 N])
set(gca, 'Xtick', [1 12 24 36 47], 'Xticklabel', [-24 -12 0 12 23], 'fontsize', 14)
set(gca, 'Position', [0.075 0.15 0.275 0.8])
xlabel('time lag (m)'); ylabel('time (n)')
axis xy

subplot(3,3,2)
plot(hd_save(360,:), 'color', cmap(2,:));
axis([1 47 -0.375 0.6])
set(gca, 'Xtick', [1 12 24 36 47], 'Ytick', [-0.3 0 0.3], 'Xticklabel', [], 'Position', [0.425 0.72 0.25 0.24])
ylabel('h(360,m)'); 
set(gca, 'fontsize', 14)

subplot(3,3,5)
plot(hd_save(280,:), 'color', cmap(2,:));
set(gca, 'fontsize', 14)
axis([1 47 -0.375 0.6])
set(gca, 'Xtick', [1 12 24 36 47], 'Ytick', [-0.3 0 0.3], 'Xticklabel', [], 'Position', [0.425 0.435 0.25 0.24])
ylabel('h(280,m)'); 

subplot(3,3,8)
plot(hd_save(115,:), 'color', cmap(2,:));
set(gca, 'fontsize', 14)
axis([1 47 -0.375 0.6])
set(gca, 'Xtick', [1 12 24 36 47], 'Ytick', [-0.3 0 0.3], 'Xticklabel', [-24 -12 0 12 23], 'Position', [0.425 0.15 0.25 0.24])
ylabel('h(115,m)'); xlabel('time lag (m)')

HD = abs(fft(hd_save(115,:)));
subplot(3,3,9)
plot(HD(1:floor(M/2)), 'color', cmap(2,:))
set(gca, 'fontsize', 14)
axis([1 23 0 1.1])
set(gca, 'Xtick', [1 11.5 23], 'Xticklabel', [0 0.25 0.5], 'Position', [0.735 0.15 0.25 0.24])
ylabel('h(115,f)'); xlabel('normalized frequency (f)')

HD = abs(fft(hd_save(280,:)));
subplot(3,3,6)
plot(HD(1:floor(M/2)), 'color', cmap(2,:))
set(gca, 'fontsize', 14)
axis([1 23 0 1.1])
set(gca, 'Xtick', [1 11.5 23], 'Xticklabel', [], 'Position', [0.735 0.435 0.25 0.24])
ylabel('h(280,f)');

HD = abs(fft(hd_save(360,:)));
subplot(3,3,3)
plot(HD(1:floor(M/2)), 'color', cmap(2,:))
set(gca, 'fontsize', 14)
axis([1 23 0 1.1])
set(gca, 'Xtick', [1 11.5 23], 'Xticklabel', [], 'Position', [0.735 0.72 0.25 0.24])
ylabel('h(360,f)');

% APPLYING A TF FILTER h(n,m) - I need something better like a
% waterall plot, but waterfall looks crap so image it is for now
figure; set(gcf, 'Position', [600 400 800 400])
colormap(1-cc);
subplot(2,2,[1 3])
imagesc(filtered_save); hold on; plot([1 512], [1 512], 'color', [0 0 0])
axis xy
set(gca, 'fontsize', 14, 'Position', [0.08 0.15 0.4 0.75], 'Xtick', [1 128 256 384 512], 'Ytick', [1 128 256 384 512])
xlabel('time lag'); ylabel('time'); title('time-varying filtering')
subplot(2,2, 2);
plot(bat2)
set(gca, 'fontsize', 14, 'Position', [0.575 0.55 0.4 0.35], 'Xtick', [], 'Ytick', [-0.2 0 0.2])
axis([1 512 -0.25 0.2]); ylabel('input')
subplot(2,2, 4);
plot(diag(filtered_save), 'color', cmap(2,:))
set(gca, 'fontsize', 14, 'Position', [0.575 0.15 0.4 0.35], 'Xtick', [1 128 256 384 512], 'Ytick', [-0.2 0 0.2])
axis([1 512 -0.25 0.2]); ylabel('output'); xlabel('time')

% OUTPUT
decomp2 = tf_decomposition_v1(bat2', len1, 2); 
figure;
plot(bat2); hold on;
plot(decomp2(1,:)-0.4, 'color', cmap(idx(2),:))
plot(decomp2(2,:)-0.65, 'color', cmap(idx(3),:))
plot(decomp2(3,:)-0.8, 'color', cmap(idx(1),:))
axis([1 512 -0.9 0.2])
set(gca, 'fontsize', 14, 'Ytick', [-0.8 -0.7 -0.4 0], 'Yticklabel', {'Comp 3','Comp 2','Comp 1', 'Signal'})
xlabel('time')


