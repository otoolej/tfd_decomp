%-------------------------------------------------------------------------------
% plot_components: plot signal and decomposed components
%
% Syntax: [hp]=plot_components(x,components,Fs,fignum,L_components,PRINT_PREPARE)
%
% Inputs: 
%     x,components,Fs,fignum,L_components,PRINT_PREPARE - 
%
% Outputs: 
%     [hp] - 
%
% Example:
%     
%

% John M. O' Toole, University of Deusto
% Started: 24-01-2013
%-------------------------------------------------------------------------------
function [hp]=plot_components(x,components,Fs,fignum,L_components,PRINT_PREPARE)
if(nargin<3 || isempty(Fs)) Fs=1; end
if(nargin<4 || isempty(fignum)) fignum=1; end
if(nargin<5 || isempty(L_components)) L_components=[]; end
if(nargin<6 || isempty(PRINT_PREPARE)) PRINT_PREPARE=0; end


FONT_TYPE='Helvetica';
FONT_SIZE=15;

EQUAL_SPACING=0;
SPACE_FACTOR = 0.3; % / abs(max(x) - min(x));





% weed out any empty components:
irem = [];
for n = 1:length(components)
    if(all(isnan(components{n})))
        irem = [irem n];
    end
end
if(~isempty(irem))
    components(irem) = [];
end



M=length(components);
N=length(x);
max_time=N./Fs;

% if want to put a limit on the number of components:
if(~isempty(L_components) && L_components <= M)
    M=L_components;
end

t=(1:N)./Fs;


%---------------------------------------------------------------------
% calculate shift
%---------------------------------------------------------------------
mi(1)=min(x); ma(1)=max(x); me(1)=mean(x);
labs{1}='signal';
for p=1:M
    ip=p-1; if(ip<1) ip=1; end
    mi(p+1)=min(components{p})+SPACE_FACTOR*min(components{ip});
    ma(p+1)=max(components{p})+SPACE_FACTOR*max(components{ip});
    me(p+1)=mean(components{p});
    labs{p+1}=['comp. ' num2str(p)];
end
ma=ma.'; mi=mi.';


if(EQUAL_SPACING)
    shift = cumsum([0; abs(max(ma).*ones(M,1))+abs(max(mi).*ones(M,1))]);
else
    space = abs(ma(1:end-1))+abs(mi(2:end));
    % uncomment this for the overlapping signals (e.g. VMD for noise case)
    % space(space < 1) = 1.1;
    
    shift = cumsum([0; space]);    
end

% minimum shift:
% b = diff(shift);
% iz = find(diff(shift) < 1e-3);
% if(~isempty(iz))
%     for n = 1:length(iz)
%         shift(iz) = shift(iz) + shift(iz - 1);
%     end
% end
% shift



%---------------------------------------------------------------------
% plot the data
%---------------------------------------------------------------------
hFig=set_figure(fignum); clf;
hp(1)=plot(t,x+shift(1),'k'); hold on;
for p=1:M
    n_component=find(components{p}~=0);

    if(~isempty(n_component))
        hp(p+1)=plot(n_component./Fs,components{p}(n_component)+shift(p+1));
    else
        hp(p+1)=plot(t, NaN(1, N)+shift(p+1));
    end

   
end
hold off;

set(gca,'xlim',[t(1),t(end)]);
ylim([mi(1), ma(end)+shift(end)]);
set(gca,'ytick',me+shift.','yticklabel',labs);
% $$$ grid on;

set(gca,'YDir','reverse')

% $$$ return;    



if(PRINT_PREPARE)
% $$$     height=1.0/1.618; % width/golden ratio
% $$$     width=1;
% $$$     scale=800; % 3.13 inches
% $$$     xpos=150; ypos=300;
% $$$     set(gcf,'Position',[xpos ypos scale*width scale*height]);
% $$$     set(gca,'ticklength',[0,0]);
    

    xl=xlim; yl=ylim;
    ylim_offset=0.05*abs(yl(2)-yl(1));
    ylim([yl(1)-ylim_offset yl(2)+ylim_offset]);

    xl=xlim; yl=ylim;
    fcolor=get(hFig,'Color');
    set(gca,'XColor',fcolor,'box','off','TickLength',[0, 0]);    
    line([xl(2), xl(2)],[yl(1), yl(2)],'color','w');    

    set(gca,'FontName',FONT_TYPE,'FontSize',FONT_SIZE);    
    

end








