function [wx, bw, fw, M] = setparameters(N, flag)

M = sqrt(N);
switch flag
    case 1  % slow IF laws
        wx = M; wx = wx+(1-rem(wx,2));
        bw = N/wx/2; fw = N/wx/4;
        bw = bw+(1-rem(bw,2)); fw = fw+(1-rem(fw,2));
        M = ceil(2*M); M = M+(1-rem(M,2));
    case 2 % fast IF laws
        wx = M; wx = wx+(1-rem(wx,2));
        bw = N/wx/4; fw = N/wx/4;
        bw = bw+(1-rem(bw,2)); fw = fw+(1-rem(fw,2));
        M = ceil(2*M); M = M+(1-rem(M,2));
    case 3 % transients
        wx = M; wx = wx+(1-rem(wx,2));
        bw = N/wx/4; fw = N/wx/2;
        bw = bw+(1-rem(bw,2)); fw = fw+(1-rem(fw,2));
        M = ceil(2*M); M = M+(1-rem(M,2));        
end