%-------------------------------------------------------------------------------
% xtfd_sep_kern: decimated cross time-frequency distribution 
%           for *separable-kernel only*
%
% NOTE: ONLY valid when Ntime>2N and Nfreq>N as the assumptions for conjugate symmetry
% that apply for real-valued TFDs are not valid here.
%
% Syntax: xtf=xtfd_sep_kern(x,y,dopp_win_params,lag_win_params,Ntime,Nfreq)
%
% INPUT: 
%      x = N-point time-domain signal
%      y = N-point time-domain signal
%
% 
%
% OUTPUT:
%      xtf = Ntime x Nfreq cross-DTFD p(nT,k/2NT)
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 27-08-2021
%
% last update: Time-stamp: <2021-08-31 18:26:52 (otoolej)>
%-------------------------------------------------------------------------------
function [xtf,R]=xtfd_sep_kern(x,y,dopp_win_params,lag_win_params,Ntime,Nfreq)
if(nargin<5 || isempty(Ntime)) Ntime=length(x); end
if(nargin<6 || isempty(Nfreq)) Nfreq=length(x)/2; end


if(length(x(:)) ~= length(y(:)))
    error('x and y must be the same length.');
end

%---------------------------------------------------------------------
% 1. get analytic signals
%---------------------------------------------------------------------
if(isreal(x)) 
    x = get_analytic_signal(x); 
end 
if(isreal(y)) 
    y = get_analytic_signal(y); 
end 
N2 = length(x);
N = N2 / 2;
Nh = ceil(Nfreq / 2);


%---------------------------------------------------------------------
% 2. compute lag and Doppler windows
%---------------------------------------------------------------------
[g2, ~, ~, Nfreq] = gen_lag_kern(lag_win_params, N, Nfreq);
[G1, ~, Ntime] = gen_Doppler_kern(dopp_win_params, N, Ntime);
g2 = padWin(g2, Nfreq);
G1 = padWin(G1, Ntime);


%---------------------------------------------------------------------
% 1. Form the time-lag function K[n, m]
%---------------------------------------------------------------------
if(Ntime>N)
    K = zeros(Ntime, Nfreq);
else
    K = zeros(N, Nfreq);
end


m = -(Nh-1):Nh;
for n = 0:N-1
    i1 = mod(n + m, N2);
    i2 = mod(n - m, N2);
    i3 = mod(m, Nfreq);

    ikeep = find(g2(i3 + 1) ~= 0);    
    K(n + 1, i3(ikeep) + 1) = x(i1(ikeep) + 1) .* conj( y(i2(ikeep) + 1) ) .* g2(i3(ikeep) + 1);
end




%---------------------------------------------------------------------
% 2. Multiply in the Doppler--lag domain
%---------------------------------------------------------------------
R = zeros(Ntime, Nfreq);


Nth = ceil(Ntime / 2);
lp = 0:Nth;
ln = 1:Nth - 1;

% transform to Doppler--lag domain and multiply by kernel:
K = fft(K);
R = K .* G1;

% and back again to time--lag domain:
R = ifft(R);



if(Ntime > N)
    R = R(1:N, :);
end


%---------------------------------------------------------------------
% 4. DFT to the time--frequency domain to get the DWVD
%---------------------------------------------------------------------
xtf = fft(R.').';
xtf = xtf ./ Nfreq;

