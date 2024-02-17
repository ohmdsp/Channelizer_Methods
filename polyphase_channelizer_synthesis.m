function xs = polyphase_channelizer_synthesis(r, fs, h)
%----------------------------------------------
% Polyphase filterbank implementation that resynthesizes K channels spaced fs/K apart
% See text by Crochiere & Rabiner, Chapter 7
%
% y = polyphase_channelizer_synthesis(x,fs,h)
%
% r:                input filter bank matrix (KxN)
% fs:               sample rate (Hz)
% h:                filter coeffs
%
% Author: drohm
%----------------------------------------------

[K,Nc] = size(r);               % get input matrix size (K=number of channels)
Nh = length(h);                 % length of FIR filter
hp = reshape(h, K , (Nh)/K);    % determine number of polyphase branches, 
                                % reshape into filter matrix with poly branch
                                % per column
cb = floor(Nh/K);               % cb is # of poly filter branches                 

hist = zeros(K,cb);             % initialize commutator matrix buffer 


j = 1;                          % initialize counters/indexes for loops
xs = zeros(1,round(Nc*K)+1);             % initialize output vector

for i=1:Nc

  %--Grab first sample from each channel - K samples
  tmp = flipud(r(:,i));         % load opposite direction of analysis commutator
  
  %--Compute FFT and load into Po branch
  hist(:,1) = K*fft(tmp, K); 
  %pause
  
  %--Apply filter 
  y = sum((hist.*hp)');   % apply filter to each branch 

  % update history buffer   
  hist = filter([0 1], 1, hist')';    % shift filtered data to next column
  
  
  xs(1,(j-1)*K+1:j*K) = y;
  j = j+1;

end

