function [wx, bw, fw, M] = setparameters(N)

M = sqrt(N);
wx = M; wx = wx+(1-rem(wx,2));
bw = N/wx/2; fw = N/wx/4;
bw = bw+(1-rem(bw,2)); fw = fw+(1-rem(fw,2));
M = ceil(2*M); M = M+(1-rem(M,2));
