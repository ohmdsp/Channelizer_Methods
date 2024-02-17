function r = direct_channelizer_analysis(x,fs,K,h,ovsfact)
%-----------------------------------------------
% Analysis channelizer using Direct Down Conversion (DDC) method.
% Returns K channels spaced fs/K apart
%
% r = direct_channelizer_rx(x,fs,K,h,ovsfact)
%
% x:            input signal (real or complex)
% fs:           sample rate (Hz)
% K:            number of filter bank channels
% h:            low pass filter coeffs
% ovsfact:      oversample factor
%
% Author: drohm
%-----------------------------------------------

%--TO TEST: Uncomment below, comment out function at top and run as a script
% K = 64;
% ovsfact = 1;
% BW = 100;
% fs = K*BW;
% N = 10000;
% x = sin(2*pi*310/fs*[0:N-1]);
% L = K*4;                          % filter length
% %h = nuttallwin(L)';
% h = fir1(L-1, BW/fs, kaiser(L, 4));
%--Uncomment To Here to Test

%--Setup channelizer
chan_space = fs/K;
Np = length(x);
M = K*ovsfact;

%--Write out information
disp(['Direct channelizer with K=',num2str(K),', channel spacing=',num2str(chan_space),'Hz, ovsfactor=',num2str(ovsfact)]);

Q = M/ovsfact;          % decimation factor (adjust osfact to increase)
D = floor(Np/Q);        % number of samples per channel
r = zeros(K,D);         % initialize output array (channels, samples) 
fo = fs/K;              % center freq spacing

for i=1:K
    mix = exp(-1j*2*pi*fo*(i-1)*1/fs*[0:Np-1]); % center frequencies for each channel
    x_if = mix.*x;
    x_if = filter(h,1,x_if);      % low pass filter
    r(i,:) = x_if(1:Q:Q*D);
end
