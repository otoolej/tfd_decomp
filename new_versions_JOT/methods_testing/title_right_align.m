%-------------------------------------------------------------------------------
% title_right_align: 
%
% Syntax: [] = title_right_align(title_txt, FONT_NAME, FONT_SIZE)
%
% Inputs: 
%     title_txt, FONT_NAME, FONT_SIZE - 
%
% Outputs: 
%     [] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 19-11-2021
%
% last update: Time-stamp: <2021-11-19 17:53:47 (otoolej)>
%-------------------------------------------------------------------------------
function title_right_align(title_txt, FONT_NAME, FONT_SIZE)
    
htit = title(title_txt, 'fontname', FONT_NAME, 'fontsize', FONT_SIZE, ...
             'fontweight', 'normal', 'horizontalalignment', 'right');
    
set(htit, 'horizontalAlignment', 'right');
set(htit, 'units', 'normalized');
h1 = get(htit, 'position');
set(htit, 'position', [1 h1(2) h1(3)])
