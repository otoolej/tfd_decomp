%-------------------------------------------------------------------------------
% fw_fprintf: fprintf to multiple locations
%
% Syntax: =fw_fprintf(fids,str,varargin)
%
% Inputs: 
%      - 
%
% Outputs: 
%      - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 24-03-2014
%
% last update: Time-stamp: <2014-03-24 13:01:59 (otoolej)>
%-------------------------------------------------------------------------------
function fw_fprintf(fids,str,varargin)
if(nargin<2), error('need more inputs.'); end

if(nargin>2), str=sprintf(str,varargin{:}); end

for n=1:length(fids)
    fprintf(fids(n),str);
end
