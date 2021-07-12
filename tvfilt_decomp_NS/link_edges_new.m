function [ml, imbm] = link_edges_new(imb, im, rx, ry, sch)
% Link edges 
%
% [ml, imbm] = link_edges(imb, rx, ry, sch, M1, M2)
%
%This function tracks a edge through the image via a particular
% neighbourhood search pattern
%
% INPUTS: imb - the binary image to be linked
%                  [rx, ry] - the starting point of the component search
%                  sch - the user-defined search region
%                  [M1, M2] - the size of the binary image
%
% OUTPUTS: ml - the [x,y] co-ordinates of the linked component
%                     imbm - the updated binary image with the linked
%                                   component removed.
%
% UQCCR, level 04
% Nathan Stevenson
% August 2008

% imb = zeros(64);
% imb(1:16,32)=1;
% for ii = 1:31; imb(ii+16, ii+32)=1; end
% for ii = 1:31; imb(ii+16, 32-ii)=1; end
% for ii = 1:31; imb(ii+32, 48-ii)=1; end

% [M1, M2] = size(imb);
% c1 = 1; 
% mlx = rx;
% mly = ry;
% flag = 1;
% imb(ry,rx)=0;
% M = size(imb); 
%temp = zeros(max(sch()))
%for ii = 1:length(sch);
rx = refq1(1); ry = refq2(1); 
bw = sch; qq = 1;

    cmp = elb_step(rx, ry, 15, imb, im);
    


% search and pick most energetic component
% long low energy beats short high energy
% ml = [mly-M2 ; mlx-M1];
% imbm = imb;






