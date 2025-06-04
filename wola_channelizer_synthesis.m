function y = wola_channelizer_synthesis(r,fs,M,h)
%----------------------------------------------
% Synthesis channelizer using Weighted Overlap-Add (WOLA) method. Outputs
% resynthesized signal using output from analysis filterbank.
% See text by Crochiere & Rabiner, Chapter 7
%
% y = wola_channelizer_synthesis(r,fs,h)
%
% r:            channelized signal (real or complex) (k=
% fs:           sample rate (Hz)
% M:            overlap in samples
% h:            filter coeffs
%
% Author: drohm
%----------------------------------------------

[c,K]=size(r);      % number of channels K from analysis filterbank, c samples 

L = length(h);      % filter length
ind = [1:L];        % length of time series slice
blksz = L/K;   

sr = zeros(1,L);    % initialize block length
y=[];
xm=[];
m = 1;              % start of phase factor index
for m=1:c           % mth block index
    
    %--get block index m from all channels
    Xk =  r(m,:);
    
    %--multiply by phase factor
    pf = exp(sqrt(-1)*2*pi*[0:K-1]*m*M/K);
    
    
    %--do K point IFFT
    Xk = K.*ifft(Xk.*pf, K);  % mult by K since fft does 1/K
    
    %--form periodic extension
    xm=[];
    for i=1:blksz
        xm = [xm Xk];
    end
    
    %--multimply by window (filter) - i.e., convolution in freq domain
    f = xm.*h;
    
    %--add to shift reg of same length (L)
    sr = sr + f;
   
    %--shift contents to left.  replace right side with zeros for next ov-add
    sr = [sr(M+1:end) zeros(1,M)];
     
    %--overlap-add, take first M pts out
    y = [y M.*sr(1:M)];     % mult by M since interpolating
    
end

