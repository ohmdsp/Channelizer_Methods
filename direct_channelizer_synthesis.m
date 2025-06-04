function xs = direct_channelizer_synthesis(r,fs,h,ovsfact)
%----------------------------------------------
% Synthesis channelizer using Direct Down Conversion (DDC) Method.
%
% y = direct_channelizer_synthesis(r,fs,h,osvfact)
%
% r:            channelized signal matrix (KxN)
% fs:           sample rate (Hz)
% h:            filter coeffs
% ovsfact:      oversample factor of baseband channels
%
% Author: drohm
%----------------------------------------------

[K,N]=size(r);      % get size of channelized data matrix
                    % K = # of channels, N = # samples per channel
                   
fo = fs/K;          % channel frequency spacing
P = K*ovsfact;        % upsample after synthesis
Np = N*P;           % number of samples output            
x = zeros(1,Np);    % initialize vector for output

%--Upsample all baseband channels, then mix to center frequency (fo)
for i=1:K
    tmp = kron(r(i,:),[1 zeros(1,P-1)]);
    tmp_flt = P*filter(h,1,tmp);        % filter the channel data
    mix = exp(sqrt(-1)*2*pi*fo*(i-1)*1/fs*[0:Np-1]); % shift to fo
    x = x + mix.*tmp_flt;
end
xs = x;