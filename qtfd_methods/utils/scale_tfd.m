%-------------------------------------------------------------------------------
% scale_tfd: TFDs are scaled to uphold the total energy property, i.e.
%            Σ|x[n]|² =  1/2 Σ Σ p[n, k]
%
% Syntax: qtfd = scale_tfd(qtfd, N)
%
% Inputs: 
%     qtfd - TFD (matrix)
%     g2   - lag window
%     N    - signal length
%
% Outputs: 
%     qtfd         - scaled TFD (matrix)
%     scale_factor - scaling factor
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 01-08-2019
%
% last update: Time-stamp: <2021-08-26 18:19:05 (otoolej)>
%-------------------------------------------------------------------------------
function [qtfd, scale_factor] = scale_tfd(qtfd, g2, N)

[Ntime, Nfreq] = size(qtfd);
scale_factor = (Nfreq / Ntime) * N / sum(g2);

qtfd = qtfd .* scale_factor;

