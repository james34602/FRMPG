[data,fs] = audioread('Storm Unity.wav'); % Load impulse response
% N-octave smoothing, extract interest points from vast amount of impulse response
Noct = 4;
[freq,gV]=smoothSpectrumImpulseResponse(data,fs,Noct); % N-octave smoothing process
% Add 0 bin and nyquist bin handling
freq=freq./(fs/2);
freq=[0 freq']';
freq(freq>1)=1;
gV=[gV(1) gV']';
% Minimum phase frequency sampling FIR filter design
y = minphasefir2(length(data),freq,db2mag(gV));
audiowrite('Storm Unity_mps.wav',y,fs,'BitsPerSample',32);
freqz(y);