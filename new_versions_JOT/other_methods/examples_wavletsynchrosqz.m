%-------------------------------------------------------------------------------
% examples_wavletsynchrosqz: 
%
% Syntax: [] = examples_wavletsynchrosqz()
%
% Inputs: 
%      - 
%
% Outputs: 
%     [] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 02-09-2021
%
% last update: Time-stamp: <2021-09-02 18:20:07 (otoolej)>
%-------------------------------------------------------------------------------
function examples_wavletsynchrosqz()

%---------------------------------------------------------------------
% bat signal
%---------------------------------------------------------------------

d = load('batsignal');

bat = d.batsignal;
Fs_bat = 1 / d.DT;

y_bat_wss = synchrosqueeze_wavelet(bat, Fs_bat, 5, true, 'bump', 5);

y_bat_xtfd = extract_components_xTFD(bb, 50, [], 4, true, []);  
