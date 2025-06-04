function [Y,f] = wola(x,fs,K)
%------------------------------------------------------------
% Weighted OverLap-Add (WOLA) filter bank function
% See text by Crochiere & Rabiner, Chapter 7
%
% [Y,f] = wola(x,fs,K)
%
% x:        input data
% fs:       sample rate
% K:        channels and FFT size
%
% Author: drohm
%-------------------------------------------------------------
N = length(x);          % length of input data
D = N/K;                % number of data segment blocks

if(D - floor(D) ~= 0)
    disp('block len must be int multiple of input!')
    Y=[], f=[];
    return
end
    
%--Segment input into blocks and sum
%--The sum will include x(n) + x(n+D) + x(n+2D) ... , n=0..M-1, M=N/DP
xp = reshape(x, K, D);      % reshape data into 3-D array  
xs = sum(xp');              % uses built in matlab function "sum"

%--Compute the FFT of summed blocks (len = N/D), keep fs the same
%--Note that this is a decimation in frequency domain
Y = fft(xs);
f = [0:length(Y)-1]*fs/(length(Y));
