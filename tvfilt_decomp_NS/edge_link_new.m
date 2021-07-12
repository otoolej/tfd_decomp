function [el,ei] = edge_link_new(imb, im, len, bw, fw)
% Edge linking for binary images.
%
%  [el,ei] = edge_link_new(imb, im, len, bw, fw)
%
% This functions links together sequences of 1's in a binary image. It is a 
% modified version of the edge linking algorithm presented in 
% Farag A, Delp E, Edge linking by sequential search. Pattern
% Recognition. 1995; 28: 611--33. The modifications include an user defined
% search region and the ability to deal with bifurcations in components -
% the component with the highest TF energy is select if bifurcations are
% present
%
% INPUT: imb - the binary image under analysis
%        im - TFD
%        len - the minimum length of an edge
%        bw - lateral search distance (i.e. frequency)
%        fw - vertical search distance (i.e. time)
%
% OUTPUT: el - a Lx1 cell array containing vectors that define the L detected 
%              edges in the image. Each cell contains the x and y co-ordinates 
%              of each linked edge. 
%         ei - is an image with all linked edge projected onto a blank image.
%
%
% Notes - lots of while loops so there is potential for instability
%
% Nathan Stevenson
% July 2021

% Initialise binary image
N = size(imb);
% Perform edge-inking procedure
[refq1, refq2] = find(imb'==1); % original values
xx1 = length(refq1);
el = cell(1); count = 1; dum1 = imb;
disp('Finding Components in TFD'); aa_old = 11;
while isempty(refq1)==0
    %count
    rx = refq1(1); ry = refq2(1);
    cmp = elb_step(rx, ry, bw, fw, dum1, im);
    AA = size(cmp);
    if AA(1)>len
        el{count} = cmp;
        count = count+1;
    end
    % reduce refq1 and refq2 do not re-estimate
    dumI = zeros(size(imb));
    dumI(sub2ind(size(imb), refq2, refq1))=1;
    dumI(sub2ind(size(imb), cmp(:,2), cmp(:,1)))=0;
    [refq1, refq2] = find(dumI'==1);
    %imagesc(dumI); drawnow; pause
    dum1 = dumI;
    aa_new = round(10*length(refq1)/xx1); str = [];
    if aa_old-aa_new~=0
        for ii = 1:aa_new; str  = [str '*']; end
        aa_old = aa_new;
        disp(str)
    end
        
end
disp([num2str(length(el)) ' Components Found'])
% Generate binary image of linked edges
ei = zeros(N);
for ii = 1:length(el)
    el1 = el{ii};
    for jj = 1:length(el1)
       ei(el1(jj,1), el1(jj,2)) =1; 
    end
    clear el1
end
