% Script for exploring prototype filter designs use in channelizer
%
% Methods:
%  - Remez algorithm
%  - Window filter method
% 
% Author: drohm
%------------------------------------------------------------------------
%------------------------------------------------------------------------
clear all; close all; clc

%--Method used in Fred Hariss paper
%--flag=0 for flat sidelobes, flag=1 for falling sidelobes
flag = 1;
hh1=remez(169,[0 40 60 500]/500,[1 1 0 0],[1 100]);
frq=[0 40 60 99 100 149 150 199 200 249 250 299 300 349 350 399 400 449 450 500]/500; 
gn=[1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
pn= [ 1 100 140 180 220 260 300 340 380 420];
hh2=remez(169,frq,gn,pn); 
hh=hh1;
if flag==1
    hh=hh2;
end

figure
subplot(2,1,1)
plot(hh)
grid
title('Impulse Response: Prototype Filter')
xlabel('Normalized time nT/T')
ylabel('Amplitude')
subplot(2,1,2) 
plot((-0.5:1/1024:.5-1/1024)*1000,fftshift(20*log10(0.000001+abs(fft(hh,1024))))) 
grid
axis([-500 500 -90 10])
title('Frequency Response: Prototype Filter')
xlabel('Frequency (kHz)')
ylabel('Log-Magnitude (dB)')

%--Matlab's filter design using windowing method
L = 170;
BW = 100;
fs = 16*100;
h = fir1(L-1, BW/fs, kaiser(L, 4)); 

figure
subplot(2,1,1)
plot(h)
grid
title('Impulse Response: Prototype Filter')
xlabel('Normalized time nT/T')
ylabel('Amplitude')
subplot(2,1,2) 
plot((-0.5:1/1024:.5-1/1024)*1000,fftshift(20*log10(0.000001+abs(fft(h,1024))))) 
grid
axis([-500 500 -90 10])
title('Frequency Response: Prototype Filter')
xlabel('Frequency (kHz)')
ylabel('Log-Magnitude (dB)')


