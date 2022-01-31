N = 512;

fs = 2000;
t = 1 / fs:1 / fs:1;
x1 = cos(2*pi*(1.5 + 150*t - 100*(t.^2) + 416*(t.^3) - 200*(t.^4)));
% x1 = x1(201:400);
% t = t(201:400);

% x1 = diff(wfbm(0.8, 201));


N = length(x1);
t = 1:N;


% short-time window:
L_bw = 20;
w = tukeywin(L_bw);

% set_figure(1); 
% plot(x1);

% buffer(L_)

y1 = zeros(1, N);
w_all = zeros(1, N);
t_bw = (1:L_bw) - 1;

i_nf_mid = floor(L_bw / 2);

N_epochs=floor( (N-L_bw) );

for k = 1:N_epochs
    nf = mod(t_bw + (k-1), N) + 1;    

    % dispVars(rem(n, 50));
    % if(rem(n, 10) == 0)
    %     p = unifrnd(0, 2 * pi, 1);
    % else
    %     p = 1;
    % end

    w = hamming(make_odd(randi([3, 11], 1)));
    w = shiftWin(padWin(shiftWin(w), L_bw)).';
    
    % set_figure(60); 
    % plot(w);
    % disp('--- paused; hit key to continue ---'); pause;


    % y1(nf) = y1(nf) + (x1(nf) .* w.');
    % y1(nf) = y1(nf) + (conv(x1(nf), w, 'same'));
    
    yc = conv(x1(nf), w .* sum(w), 'same');
    imid = nf(i_nf_mid);
    y1(imid) = yc(i_nf_mid);


    w_all(nf) = w_all(nf) + w;
end

y1 = y1 ./ sum(w_all);
set_figure(2); 
plot(t, y1);

set_figure(3);
subplot(3, 1, 1);
plot(w_all);
subplot(3, 1, 2); hold all;
plot(y1);
plot(x1);
subplot(3, 1, 3); hold all;
plot(x1 - y1);



qtfd = qtfd_sep_kern(y1, {55, 'hamm'}, {33, 'dolph', 100}, N, N);
set_figure(4); 
vtfd(qtfd); axis('tight');


% z_all=zeros(1,N); win_summed=zeros(1,N);
% for k=1:N_epochs
%     nf=mod(nw+(k-1)*L_hop,N);
%     x_epoch=est_IF(nf+1);

%     ev=median(x_epoch);

%     z_all(nf+1)=z_all(nf+1) + (ones(1,L_epoch)*ev);
%     win_summed(nf+1)=win_summed(nf+1)+ones(1,L_epoch);            
% end
% t_stat=z_all./win_summed;
