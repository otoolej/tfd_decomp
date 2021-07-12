function [el1, ei1] = find_components(tfrep, tfrep1, bw, fw, len1);
%
%
%
%
% dependencies edge_link_new.m

a = size(tfrep1);
tfrep1(tfrep1<0) = 0;
im_bin1 = zeros(a); im_bin2 = zeros(a);
for kk = 1:a(2)
    q1 = tfrep1(:,kk); q3 = tfrep1(kk,:);
    u = diff(q1); v = diff(q3);
    %uu = q1; vv = q3;
    zc_max1 =[]; zc_max2 =[];
    count_max1 = 1; count_max2 = 1;
    for zz = 1:(length(u)-2)
        if u(zz)>=0
           if u(zz+1)<0
            zc_max1(count_max1) = zz;
            count_max1 = count_max1+1;
           end
        end
        if v(zz)>=0
           if v(zz+1)<0
            zc_max2(count_max2) = zz;
            count_max2 = count_max2+1;
           end
        end
    end
    q2 = zeros(1,a(1));
    q2(zc_max1+1) = 1;
    im_bin1(kk,:) = q2;
    q4 = zeros(1,a(1));
    q4(zc_max2+1) = 1;
    im_bin2(:,kk) = q4;

    clear q* z* c* v* ref
end
im_bin = im_bin1; %+im_bin2;
%im_bin(im_bin>0)=1;

% % Search binary image for nonstationary components.
%bw = 31; fw = 31;
%bw = 5; fw = 5;
[el1,ei1] = edge_link_new(im_bin, tfrep', len1, bw, fw);
%[el1, ei1] = edge_link(im_bin1, 63, sch, 2*bw);