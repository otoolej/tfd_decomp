function lc1 = elb_step_v2(rx1, ry1, bw, imb, im)


% imb = zeros(64);
% imb(1:16,32)=1;
% for ii = 1:31; imb(ii+16, ii+32)=1; end
% for ii = 1:31; imb(ii+16, 32-ii)=1; end
% for ii = 1:31; imb(ii+32, 48-ii)=1; end

c3=1; [M1, M2] = size(imb); bd = 1;
bif{1} = [1 2 3]; comp = cell(1); bif_old = bif;
while isempty(bif{1})==0
    bif = cell(1); yr = M1; c1=1; c2=1; qq=1; rx = rx1; ry = ry1; %rx = 32; ry = 1; 
  %while length(qq)==1 || isempty(qq)==1
    while isempty(qq)==0 & ry<M1 & bd~=0
        r1 = rx(c1)-bw; r2 = rx(c1)+bw;
        if r1<1; r1 = 1; end
        if r2>M2; r2 = M2; end
        val = r1:r2;
        qq1 = find(imb(ry(c1)+1,val)==1);
        if length(qq1)>1
        
        
        if length(qq1)>1 && ry(c1)<yr
           bif{c2} = [val(qq1) c1]; 
           yr = ry(c1);
           c2 = c2+1; 
        end
        qq = find(imb(ry(c1)+1,val)==1);
        qq = qq(find(abs(qq-(bw+1)) == min(abs(qq-(bw+1))),1));
        %find(qq-bw)
        if length(qq)==1
            c1 = c1+1;
            rx(c1) = val(qq); 
            ry(c1) = ry(c1-1)+1;
        end
    end
    comp{c3,1} = rx; comp{c3,2} = ry; c3 = c3+1;
    if isempty(bif{1})==0
        imb(ry(bif{end}(3)+1:end), rx(bif{end}(3)+1:end)) = 0;           
    end  
    if length(bif{end})==length(bif_old)
        bd = sum(bif_old - bif{end});
    else
        bd=1;
    end
    bif_old = bif{end};
end

for ii = 1:length(comp)-1
   val = sum(im(comp{ii,2}));
end
zz = find(val==max(val),1);
lc1 = [comp{zz,1} comp(zz,2)];