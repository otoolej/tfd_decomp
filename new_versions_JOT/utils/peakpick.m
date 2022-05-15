%-------------------------------------------------------------------------------
% peakpick: find peaks
%
% Syntax: y = peakpick(x, method_type, k0_estimate)
%
% Inputs: 
%     x           - input signal
%     method_type - type of method, either 'harmonic_peaks' or 'diff' (default)
%     k0_estimate - estimate of fundamental component (for harmonic_peaks method only)
%
% Outputs: 
%     y - mask with peaks
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 26-08-2021
%
% last update: Time-stamp: <2021-08-26 18:37:59 (otoolej)>
%-------------------------------------------------------------------------------
function y = peakpick(x, method_type, k0_estimate)
if(nargin<2 || isempty(method_type)) method_type='diff'; end
N=length(x);
y=zeros(N,1);


switch method_type
  case 'harmonic_peaks'
    win_width=4;
    w_l=0; end_freq=0;
    
    k_lower=floor( k0_estimate-(k0_estimate/win_width) );
    k_upper=floor( k0_estimate+(k0_estimate/win_width) );
    
    
    while(end_freq==0)
        n_window=(w_l+k_lower):(w_l+k_upper);
        
        if(n_window(end)>length(x))
            end_freq=1;
        else
            [xmax,imax]=max(x(n_window));
            w_l=n_window(imax);        
            y(w_l)=1;
        end
    end
    
    
  case 'diff'
    y(2:N-1)=sign(-sign(diff(sign(diff(x)))+0.5)+1);
    
  otherwise
    error('which peak-picking method to use?');
end


