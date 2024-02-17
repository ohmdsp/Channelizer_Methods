% Demo script of channelization using Direct Down Conversion (DDC)
% method. Input signals are generated as real-valued sinusoids. 
% To Do: use optimal reconstruction filter for better MSE
% 
% Author: drohm
%------------------------------------------------------------------------
%------------------------------------------------------------------------
clear all; close all

K = 64;             % number of filterbank channels
ovsfact = 1;        % oversample factor
BW = 100;           % channel filter bandwidth 
fs = K*BW;          % sample rate
%N = 10000;          % number of samples for input signal

%--Simulate signal with multiple frequencies
%w = 0.1*randn(1,N);
%x = 0.5*sin(2*pi*305/fs*[0:N-1])+0.5*sin(2*pi*1305/fs*[0:N-1]) ;  

%--Create increasing carrier freq list with small freq offsets in each band
fmax = fs/2;        % max freq of input signal
flag = 1;
freqs=[1];
step = BW *(1+0.01);
i=step;
while(flag)
    freqs = [freqs i];
    i=i+step;
    if(freqs(end) > fmax*0.95)
        flag = 0;
    end
end
%--Generate signals
N = 10*1024;
x = zeros(1,N);         % initialize output vector
for i=1:length(freqs)
    x = x + sin(2*pi*freqs(i)/fs.*[0:N-1]) ; 
end

%--Design channel low pass filter
L = K*4;                  % filter length
h = fir1(L-1, BW/fs, kaiser(L, 4)); % Matlab filter design 

%--Call the DDC filterbank function
r = direct_channelizer_analysis(x,fs,K,h,ovsfact);

%--Plot channel outputs
figure
KP = K/2;       % only half needed with real-valued signal input

for kk=1:KP
   subplot(round(sqrt(KP)), round(sqrt(KP)), kk)
   plot(real(r(kk,:)))
   hold on;
   grid  
   axis([1 fs/K -1.1 1.1])
   title(['Fc: ', num2str((kk-1)*fs/K),' Hz'])
end
sgtitle('DDC Channelizer')

%--Compute direct resynthesis of signal
xs = direct_channelizer_synthesis(r,fs,h,ovsfact);

%--Plot resynthesized signals (note: filter delay corrected)
figure
disp('Plotting original signal vs resynthsized signal')
plot(1/fs*[L:3*L],x(L:3*L),'k')
hold on
plot(1/fs*[L:3*L],real(xs(L+L/2+1:3*L+L/2+1)),'b')
hold off
xlim(1/fs*[L 3*L]);xlabel('time (sec)')
grid
legend('Original','Direct');
sgtitle('Original vs. Resynthesized (Direct) Signal')
disp(' ')


%--Compute mean-squared-error between WOLA and Direct outputs (synthesis)
mse = mean(( x(L:3*L) - real(xs(L+L/2+1:3*L+L/2+1)) ).^2);
disp(['MSE between resynthesized and original: ',num2str(mse)]);
disp(' ')

