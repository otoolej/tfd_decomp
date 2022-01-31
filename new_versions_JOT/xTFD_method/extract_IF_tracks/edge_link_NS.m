function [el,ei] = edge_link_NS(imb, im, len, bw, fw)
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
% disp('Finding Components in TFD'); aa_old = 11;
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
    % aa_new = round(10*length(refq1)/xx1); str = [];
    % if aa_old-aa_new~=0
    %     for ii = 1:aa_new; str  = [str '*']; end
    %     aa_old = aa_new;
    %     disp(str)
    % end
        
end
disp([num2str(length(el)) ' Components Found'])

m1 = zeros(1,length(el)); m2 = m1;
for ii = 1:length(el)
   m1(ii) = min(el{ii}(:,2));
   m2(ii) = max(el{ii}(:,2));
end
if m1>1;
    aa = find(m1==min(m1));
    aa = aa(1);
    el{aa} = [el{aa}(1,1)*ones(min(m1),1) [1:min(m1)]' ; el{aa}];
end
if m2<N(1);
    aa = find(m2==max(m2));
    aa = aa(1);
    % dispVars(aa, m2, max(m2));
    el{aa} = [el{aa} ; el{aa}(end,1)*ones(N(1)-max(m2),1) [max(m2)+1:N(1)]'];
end


% Generate binary image of linked edges
ei = zeros(N);
for ii = 1:length(el)
    el1 = el{ii};
    for jj = 1:length(el1)
       ei(el1(jj,1), el1(jj,2)) =1; 
    end
    clear el1
end
