% Demo script for channelization using polyphase filter bank
% method. Input signals are generated as sum of real-valued sinusoids. 
% 
% Author: drohm
%------------------------------------------------------------------------
%------------------------------------------------------------------------
clear all; close all

K = 64;         % # of filter bank channels
BW = 100;       % channel filter bandwidth 
fs = K*BW;        % sample rate
N = 10000;        % number of samples for input signal

%--Simulate signal with two frequencies
%x = randn(1,N);
%x = 0.5*sin(2*pi*310/fs*[0:N-1]) + 0.5*sin(2*pi*1290/fs*[0:N-1]);
%x = 0.5*sin(2*pi*210/fs*[0:N-1]);

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
%h = nuttallwin(L)';     % nuttall window of width L
h = fir1(L-1, BW/fs, kaiser(L, 4)); % Matlab filter design using window method

L = length(h);
figure,subplot(2,1,1)
stem([0:L-1],h);
xlim([0 L]);
title('The Prototype Lowpass Filter');
subplot(2,1,2)
plot(linspace(-1,1,4*L),20*log10(abs(fftshift(fft(h,4*L)))));
title('The Prototype Lowpass Filter Spectrum (dB)');
grid on;


%--Call the DDC filterbank function
ovsfact = 1;
r = polyphase_channelizer_analysis(x,fs,K,h,ovsfact);

%--Plot channel outputs
figure
Kp = K/2; 

for kk=1:Kp
   subplot(round(sqrt(Kp)), round(sqrt(Kp)), kk)
   plot(real(r(kk,:)));
   hold on;
   grid  
   axis([1 fs/K -1.1 1.1]);
   title(['Fc: ', num2str((kk-1)*fs/K),' Hz'])
end
sgtitle('Polyphase Channelizer')


%--Compute direct resynthesis of signal
xs = polyphase_channelizer_synthesis(r,fs,h);

%--Plot resynthesized signals (note: filter delay corrected)
figure
disp('Plotting original signal vs resynthsized signal')
plot(1/fs*[L:3*L],x(L:3*L),'k-')
hold on
plot(1/fs*[L:3*L],real(xs(L+L/2+1:3*L+L/2+1)),'r-')
hold off
xlim(1/fs*[L 3*L]);xlabel('time (sec)')
grid
legend('Original','Polyphase');
sgtitle('Original vs. Resynthesized (Polyphase) Signal')
disp(' ')


%--Compute mean-squared-error between outputs (synthesis)
mse = mean(( x(L:3*L) - real(xs(L+L/2+1:3*L+L/2+1)) ).^2);
disp(['MSE between resynthesized and original: ',num2str(mse)]);
disp(' ')

