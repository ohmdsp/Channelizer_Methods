% Demo script showing various methods of channelization 
%
% Methods:
%  - Direct Down Conversion (DDC)
%  - Weighted Overlap-Add (WOLA)
%  - Polyphase
% 
% Author: drohm
%------------------------------------------------------------------------
%------------------------------------------------------------------------
clear all; close all; clc

%--Generate sum of sinusoids or BPSK signals for channelizer demonstration
%car_only = 1;          % 1-Sinusoids, 0-BPSK signals
realorim = 0;           % demo uses real-valued signal as input, set flag to 0
K = 64;                 % number of channels output from channelizer
ovsfact = 1;            % oversample factor applied to output sample rate = fs/K  
                        % 1X only supported with polyphase channelizer
M = K*ovsfact; 

%--Channels are centered at n*fs/K with BW=fs/K
BW = 100;           % channel bandwidth
fs = K*BW;          % sample rate of input
chan_space = fs/K;

alphaBW = 0.2;      % decimation filter excess BW
fmax = fs/2;        % max freq of input signal

%--Create increasing carrier freq list with small freq offsets in each band
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
N = 100*1024;
x = zeros(1,N);         % initialize output vector
for i=1:length(freqs)
    x = x + sin(2*pi*freqs(i)/fs.*[0:N-1]) ; 
end


%--Design filter for channelizer
disp('Creating decimation filter')
L = K*8;                  % filter length
h = fir1(L-1, BW/fs, kaiser(L, 4)); % using Matlab's filter design

disp('Plotting filter bank frequency response')
%--Plot input signal(s) spectrum and filter bank frequency response together
Lp = 16*1024;
if(realorim == 0)
    Kp = K/2;   %only need half the channels for real signal
else
    Kp = K;
end
[Pxx,F] = pwelch(x,hamming(Lp),Lp/2,Lp,fs);
figure
plot(F, (db(Pxx)) , 'r-','LineWidth',1);
grid
[H FH] = freqz(h,1,1024,fs,'twosided');
hold on; 
for i = 1:Kp
     plot(FH-fs/2 + chan_space*(i-1), fftshift(db(H)), 'k-','LineWidth',.1);
end
hold off;
xlim([-fs/4,fs-fs/4])
sgtitle('PSD of Filter Bank and Input Signal')
ylabel('dB')
xlabel('Freq (Hz)')
legend('Input Signals','Channelizer Response');
pause


%--List frequencues were signals fall into bands, 1st chan = DC
chans = (freqs)/chan_space + 1;
disp(['Channels with signals: ',mat2str(floor(chans))]);


disp(' ')
disp('========== ANALYSIS =============================')
disp(' ')

%================================================
%======== Polyphase Channelizer =================
%================================================
tic
ovsfact = 1;
yy = polyphase_channelizer_analysis(x,fs,K,h,ovsfact);
e = toc;
disp(['Polyphase channelizer done in ',num2str(e), ' sec'])
%[K c] = size(yy);

if(realorim == 0)
    Kp = K/2;   %only need half the channels for real signal
else
    Kp = K;
end

%--The output sample rate of each band is fs/K, scale by ovsfact
fso = fs/K;

%--Plot the transmux channel outputs
figure
for kk=1:Kp
   subplot(round(sqrt(Kp)), round(sqrt(Kp)), kk)
   plot(real(yy(kk,:)))
   hold on;
   grid  
   axis([fs/K 2*fs/K -1.1 1.1])
   title(['Fc: ', num2str((kk-1)*chan_space),' Hz'])
end
sgtitle('Polyphase Channelizer - Analysis')
disp(' ')

%==================================================
%======== WOLA Analysis channelizer ===============
%==================================================
tic;
Moverlap = round(K*ovsfact);   %can do noninteger ovsfact for wola
[yw, fsow] = wola_channelizer_analysis(x,fs,K,h,Moverlap);
e2 = toc;
disp(['WOLA channelizer done in ',num2str(e), ' sec'])
[c K] = size(yw);

if(realorim == 0)
    Kp = K/2;   % only need half the channels for real signal
else
    Kp = K;
end

figure
for kk=1:Kp
   subplot(round(sqrt(Kp)), round(sqrt(Kp)), kk)
   plot(real(yw(:,kk)))
   hold on;
   grid  
   axis([fs/K 2*fs/K -1.1 1.1])
   title(['Fc: ', num2str((kk-1)*fs/K),' Hz'])
end
sgtitle('WOLA Channelizer - Analysis')
disp(' ')

%=====================================================
%======= Direct Down-Conversion Channelizer ==========
%=====================================================
tic;
yyd = direct_channelizer_analysis(x,fs,K,h,ovsfact);
e3 = toc;
disp(['Direct channelizer done in ',num2str(e2), ' sec'])

if(realorim == 0)
    Kp = K/2;   % only need half the channels for real signal
else
    Kp = K;
end
    
figure
for kk=1:Kp
   subplot(round(sqrt(Kp)), round(sqrt(Kp)), kk)
   plot(real(yyd(kk,:)),'r')
   hold on;
   grid
   axis([fs/K 2*fs/K -1.1 1.1])
   title(['Fc: ', num2str((kk-1)*fs/K),' Hz'])
end
sgtitle('DDC Channelizer - Analysis')
disp(' ')

%=====================================================
%======= Analysis Performance Comparisons ============
%=====================================================
disp(['Analysis speed of Polyphase Channelizer over Direct Channelizer  =  ',num2str(e3/e), 'X'])
disp(['Analysis speed of WOLA Channelizer over Direct Channelizer       =  ',num2str(e3/e2), 'X'])
disp(' ')

% for kk=1:K  
%    mind = min(length(yyd(1,:)), length(yy(1,:)));   
%    err1 = (yyd(kk,1:mind))-(yy(kk,1:mind));     % polyphase vs direct error
%    mind = min(length(yyd(1,:)), length(yw(1,:)));
%    err2 = (yyd(kk,1:mind)) - (yw(kk,1:mind));   % wola vs direct error 
%    disp(['Channel ',num2str(kk),' MSE (Polyphase to direct): ',num2str((abs(mean(err1)))^2),', MSE (WOLA to direct): ',num2str(abs(mean(err2))^2)])  
% end
% disp(' ')



%--Compute synthesis direction for each channelizer
disp(' ')
disp('============= SYNTHESIS ===========================')
disp('Direct channelizer resynthesis from baseband channels')
tic;
xu = direct_channelizer_synthesis(yyd,fs,h,ovsfact);
e = toc;
disp(['Direct resynth done in ',num2str(e2), ' sec']);
disp(' ')

disp('WOLA channelizer resynthesis from baseband channels')
tic;
xuw = wola_channelizer_synthesis(yw,fs,Moverlap,h);
e2 = toc;
disp(['WOLA resynth done in ',num2str(e2), ' sec']);
disp(' ')

disp('Polyphase channelizer resynthesis from baseband channels')
tic;
xup = polyphase_channelizer_synthesis(yy,fs,h);
e3 = toc;
disp(['Polyphase resynth done in ',num2str(e3), ' sec']);
disp(' ')

%--Compute mean-squared-error between original and Direct outputs (synthesis)
mse = mean(( x(2*L:3*L) - real(xu(2*L+5:3*L+5)) ).^2);
disp(['MSE between Direct synthesis and original signal: ',num2str(mse)]);

%--Compute mean-squared-error between original and WOLA outputs (synthesis)
mse = mean(( x(2*L:3*L) - real(xuw(2*L+L/2+2:3*L+L/2+2)) ).^2);
disp(['MSE between WOLA synthesis and original signal: ',num2str(mse)]);

%--Compute mean-squared-error between original and Polyphase outputs (synthesis)
mse = mean(( x(2*L:3*L) - real(xup(2*L+5:3*L+5)) ).^2);
disp(['MSE between Polyphase synthesis and original signal: ',num2str(mse)]);
disp(' ')


%--Plot resynthesized signals
disp('Plotting channelizer resynthesized signals')
figure
plot(1/fs*[L:3*L],x(L:3*L),'k')
hold on
plot(1/fs*[L:3*L],real(xu(L+5:3*L+5)),'g')   
plot(1/fs*[L:3*L],real(xuw(L+L/2+2:3*L+L/2+2)),'b')
plot(1/fs*[L:3*L],real(xup(L+5:3*L+5)),'r')  
xlim(1/fs*[2*L 3*L]);xlabel('time (sec)')
grid
legend('Original','Direct','WOLA', 'Polyphase');
sgtitle('Original vs. Resynthesized Signals')


%--Plot spectrums for resynthesized data
disp('Plotting spectrum of original signal vs resynthesized signals')
figure; 
X = fft(x(1:8192));
subplot(4,1,1)
plot([0:length(X)-1]*fs/length(X)-fs/2, fftshift(db(X)),'b')
xlim([-fs/2,fs/2]);ylim([-10 90]);
xlabel('frequency (Hz)')
title('Original Spectrum')
grid

U = fft(xu(1:8192));
subplot(4,1,2)
plot([0:length(U)-1]*fs/length(U)-fs/2, fftshift(db(U)),'r')
xlim([-fs/2,fs/2]);ylim([-10 90]);
xlabel('frequency (Hz)')
title('DDC Synthesis Spectrum')
grid

Uw = fft(xuw(1:8192));
subplot(4,1,3)
plot([0:length(Uw)-1]*fs/length(Uw)-fs/2, fftshift(db(Uw)),'r')
xlim([-fs/2,fs/2]);ylim([-10 90]);
xlabel('frequency (Hz)')
title('WOLA Synthesis Spectrum')
grid

Up = fft(xup(1:8192));
subplot(4,1,4)
plot([0:length(Up)-1]*fs/length(Up)-fs/2, fftshift(db(Up)),'r')
xlim([-fs/2,fs/2]);ylim([-10 90]);
xlabel('frequency (Hz)')
title('Polyphase Synthesis Spectrum')
grid
disp(' ')



% %=====================================================
% %======= Synthesis Performance Comparisons ===========
% %=====================================================
disp(['Synthesis speedup of Polyphase over Direct Channelizer  =  ',num2str(e/e3), 'X'])
disp(['Synthesis speedup of WOLA over Direct Channelizer       =  ',num2str(e/e2), 'X'])
disp(' ')
