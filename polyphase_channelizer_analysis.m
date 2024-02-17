function r = polyphase_channelizer_analysis(x, fs, K, h, ovsfact)
%----------------------------------------------
% Polyphase filterbank implementation that returns K channels spaced fs/K apart
% See text by Crochiere & Rabiner, Chapter 7
%
%
% y = polyphase_channelizer_analysis(x,fs,K,h,ovsfact)
%
% x:                input signal (real or complex)
% fs:               sample rate (Hz)
% K:                number of filterbank channels
% h:                filter coeffs (length must be divisable into K)
% ovsfact:          oversample factor (only 1x and 2x oversample supported)
%
% Author: drohm
%----------------------------------------------
%--TO TEST: Uncomment below, comment out function at top and run as a script
% K = 64;
% ovsfact = 1;
% BW = 100;
% fs = K*BW;
% N = 10*1024;
% x = sin(2*pi*310/fs*[0:N-1]);
% L = K*8;                          % filter length
% h = nuttallwin(L)';
% h = fir1(L-1, BW/fs, kaiser(L, 4));
%--TO TEST: Uncomment To Here

M = K;                          % decimation amount
chan_space = fs/M;

disp(['Polyphase channelizer with K=',num2str(K),', channel spacing=',num2str(chan_space),'Hz, ovsfact=',num2str(ovsfact)]);

N = length(h);
hp = reshape(h, M, (N)/M);      % reshape prototype filter into poly branches

[rb,cb] = size(hp);              % rb = channels, cb = # poly branch taps (i.e., fftsize) 
hist = zeros(rb,cb);
cnt = 1;
j = 1;

Np = length(x);
r = zeros(M,ovsfact*round(Np/M)+1);

for i=1:length(x)

  % clockwise commutator - insert new sample

  hist(cnt,1) = x(i);           % insert first sample as initialization
  cnt = cnt-1;

  %--Run filter once when commutator distributes M samples to network
  if(cnt == 0)

    %--Decimated output from all branches
    y = sum((hist.*hp)');       % compute filter on all poly branches          

    if(ovsfact == 2)

        cnt = M/2;              %start filling poly branches from M/2 up to 0    

        % 1st percolate history buffer to do the polyphase oversample
        % here shift hist buf by M/2 and repartition back to polyphase
        tmph = reshape(hist, 1, rb*cb);
        tmph = filter([zeros(1,M/2) 1], 1, tmph); 
        hist = reshape(tmph, rb, cb);

        % 2nd step is to account for phase shift out of the filter.  For 2x
        % ovsfact we can do a circular shift of N/2 before FFT, or do pi radian
        % phase correction on odd bins after the FFT

        % Here we do the time shift prior to fft (NOTE: only do on alternating
        % blocks for 2x oversample)
        if( mod(j,2) == 0)
            y = [y(M/2 + 1: end) , y(1:M/2)];
        end

    elseif(ovsfact == 1)
        cnt = M;

        % pupdate history buffer by shifting samples in poly-branches
        hist = filter([0 1], 1, hist')'; 

    end
    
    % Send branch outputs through FFT for modulation and sum
    % Use ifft so synthesis step does translation correctly
    Y= M*ifft(y, M);    % scale ifft by M since it has 1/M
    r(:,j) = Y;         % output from FFT is M channels at baseband

    j=j+1;
  end


end

