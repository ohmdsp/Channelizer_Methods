function [r, fso] = wola_channelizer_analysis(x, fs, K, h, M)
%-------------------------------------------------------------
% Analysis channelizer using the Weighted OverLap-Add (WOLA) method that 
% generates K output channels spaced fs/K apart. 
% See text by Crochiere & Rabiner, Chapter 7
%
% r = wola_channelizer_analysis(x,fs,K,h,M)
%
% Output: r output matrix [samples per channel, # channels]
%
% x:        input signal (real or complex)
% fs:       sample rate (Hz)
% K:        number of channels (FFT size)
% h:        filter coeffs (length must be div into K)
% M:        channelizer decimation/interpolation ratio
% 
% Note: oversample factor M/K = 1 (M=K ->critically sampled)
%
% Author: drohm
%-------------------------------------------------------------
%%--TO TEST: uncomment below, comment out function at top and run as a script
% M = 64;
% K = 64;             % number of frequeny channels
% BW = 100;           % channel filter bandwidth 
% fs = K*BW;          % sample rate
% N = 10000;         % number of samples for input signal
% 
% %--simulate signal with two frequencies
% x = 0.5*sin(2*pi*310/fs*[0:N-1]) + 0.5*sin(2*pi*1290/fs*[0:N-1]);
% 
% %--design channel low pass filter
% L = K*4;                  % filter length
% h = fir1(L-1, BW/fs, kaiser(L, 4)); % Matlab filter design using window method
%%--Uncomment to here

 
chan_space = fs/K;          % spacing between channels
ovsfact = M/K;
disp(['WOLA channelizer with K=',num2str(K),', channel spacing=',num2str(chan_space),'Hz, ovsfactor=',num2str(ovsfact)]);

L = length(h);              % length of FIR filter
ind = [1:L];                % create index vector for filter

%--Prepend L zeros to input. This will make outputs match with direct channelizer.
x = [zeros(1,L) x];

%--Determine overlap for input segments based on M shift
Nov = (K-M);                % number of overlap samples per block 
N = length(x);              % length of input data

%--Determine number of input data segments to compute
blks = floor( (N - L - Nov)/(K - Nov ) ) ;

Y = zeros(blks,K);          % initialize channel data matrix
 
for m=1:blks
    [Y(m,:),f] = wola(x(ind).*h,fs,K);
    Y(m,:) = Y(m,:).*exp(-sqrt(-1)*2*pi*[0:K-1]*m*M/K);
    ind = ind+M;
end

%--Wola method uses critically sampled scenereo with M=K, otherwise oversampled by K/M
%--and K/M is not constrained to be an integer
fso = chan_space * (K/M); 
r = Y;



