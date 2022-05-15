%-------------------------------------------------------------------------------
% set_gca_fonts: 
%
% Syntax: set_gca_fonts(FONT_NAME,FONT_SIZE,ha_all,exclude_text,ADJ_LINES)
%
% Inputs: 
%     FONT_NAME  - default 'Arial'
%     FONT_SIZE  - default 12
%     ha_all     - handles 
%     exclude_text - text not to change size
%     ADJ_LINES    - include ticks (default=1)
%
% Outputs: 
%      - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 03-04-2017
%
% last update: Time-stamp: <2017-04-03 11:28:48 (otoolej)>
%-------------------------------------------------------------------------------
function set_gca_fonts(FONT_NAME,FONT_SIZE,ha_all,exclude_text,ADJ_LINES)
if(nargin<1 || isempty(FONT_NAME)), FONT_NAME='Arial'; end
if(nargin<2 || isempty(FONT_SIZE)), FONT_SIZE=12; end
if(nargin<3 || isempty(ha_all)), ha_all=gca; end
if(nargin<4 || isempty(exclude_text)), exclude_text=[]; end
if(nargin<5 || isempty(ADJ_LINES)), ADJ_LINES=1; end


if(~isempty(exclude_text))
    if(~iscell(exclude_text)), exclude_text={exclude_text}; end
else
    exclude_text={'NS','*','**','***'};
end



L=max(size(ha_all));
for l=1:L
    ha=ha_all(l);
    
    all_text=findall(ha, 'Type', 'Text');
 
    % if want to keep some text a certain size:
    if(~isempty(all_text))
        if(~isempty(exclude_text))
            itext_remove=[];
            for p=1:length(all_text)
                if(any(find( strcmp( exclude_text, get(all_text(p),'string') ) )))
                    itext_remove=[itext_remove p];
                end
            end
            all_text(itext_remove)=[];
        end
        set(all_text,'FontName',FONT_NAME,'FontSize',FONT_SIZE);
    end
    
    set(ha,'FontName',FONT_NAME,'FontSize',FONT_SIZE);

    % axis ticks:
    if(ADJ_LINES)
        set(ha, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.015 .015] , ...
            'XMinorTick'  , 'off'      , ...
            'YMinorTick'  , 'off'      , ...
            'YGrid'       , 'off'      , ...
            'Xgrid'       , 'off'      , ...
            'LineWidth'   , 1         , ...
            'FontName'   , FONT_NAME , ...
            'FontSize'   , FONT_SIZE  ...    
            );
    end
% $$$     'XColor'      , [.3 .3 .3], ...
% $$$     'YColor'      , [.3 .3 .3], ...
end
