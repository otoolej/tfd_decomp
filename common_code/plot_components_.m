%-------------------------------------------------------------------------------
% plot_components: 
%
% Syntax: [hax]=plot_components(x_comp,y_gap)
%
% Inputs: 
%     x_comp,y_gap - 
%
% Outputs: 
%     [hax] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 03-08-2016
%
% last update: Time-stamp: <2019-09-06 16:36:20 (otoolej)>
%-------------------------------------------------------------------------------
function [hax,hlines]=plot_components_(x_comp,y_gap,xticks, fig_num)
if(nargin<2 || isempty(y_gap)), y_gap=0; end
if(nargin < 3 || isempty(fig_num)), fig_num = 1; end
if(nargin < 4 || isempty(xticks)), xticks=[]; end


figure(fig_num);
clf; hold all;

[L,N]=size(x_comp);
if(isempty(xticks))
    xticks=1:N;
end


y_gap=nanmean(nanstd(x_comp'));

hlines=zeros(1,L);
all_gap=0;
yl = [];
for l=1:L
    if(l>1)
        yheight=max(x_comp(l,:))-nanmean(x_comp(l,:));
        all_gap=y_gap+abs(yl)+yheight;
    else
        x_comp(1,:)=x_comp(1,:)-nanmean(x_comp(1,:));
    end
    
    hlines(l)=plot(xticks,x_comp(l,:)-all_gap);
    yl=min([get(hlines(l),'ydata') yl]);
end


xlim([xticks(1) xticks(end)]);
ylim([yl max(get(hlines(1),'ydata'))]);
