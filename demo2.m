fs = 48000.0;
freq = [0 200 800 1900 3000 5000 12000 fs/2];
freq=freq./(fs / 2);
gV = [-30 0 16 0 -8 0.1 5 0];
% Minimum phase frequency sampling FIR filter design
bhi = minphasefir2(8192,freq,db2mag(gV));
freqz(bhi);