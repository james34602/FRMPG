function [freqVector,gainPoints] = smoothSpectrumImpulseResponse(xn,fs,Noct)
% Matthes (2019). Short time FFT with octave smooth (https://www.mathworks
% .com/matlabcentral/fileexchange/19228-short-time-fft-with-octave-smooth),
%MATLAB Central File Exchange. Retrieved December 14, 2019. 
% Fung (2019) <james34602@gmail.com>
[octfreq,Px]=averfft(xn,fs,Noct,2^nextpow2(length(xn)));
[~,idx] = unique(Px);
idx=sort(idx);
freqVector=octfreq(idx);
gainPoints=Px(idx);
end
function [fc,Px_oct]=averfft(sig,fs,Noct,Nfft)
if size(sig,2)>size(sig,1)
    sig=sig';
end
if nargin < 2
    Noct=0;
end
fsdiv2=fs/2;
L=length(sig);
Px=20*log10(abs(fft(sig,Nfft)));
Px=Px(1:end/2,:);
freq=(0:fs/(Nfft-1):fsdiv2)';
%--------------------------------------------------------------------------
% octave smoothing
Noct=2*Noct;
% octave center frequencies
f1=1;
i=0;
while f1 < fsdiv2
    f1=f1*10^(3/(10*Noct));
    i=i+1;
    fc(i,:)=f1;
end
% octave edge frequencies
for i=0:length(fc)-1
    i=i+1;
    f1=10^(3/(20*Noct))*fc(i);
    fe(i,:)=f1;
end
% find nearest frequency edges
for i=1:length(fe)
    fe_p=find(freq>fe(i),1,'first');
    fe_m=find(freq<fe(i),1,'last');
    fe_0=find(freq==fe(i));
    if isempty(fe_0)==0
        fe(i)=fe_0;
    else
        p=fe_p-fe(i);
        m=fe(i)-fe_m;
        if p<m
            fe(i)=fe_p;
        else
            fe(i)=fe_m;
        end
    end
end
Px_oct = zeros(length(fe)-1,1);
for i=1:length(fe)-1
    Px_i=Px(fe(i):fe(i+1),:);
    Px_oct(i,1:size(Px,2))=mean(Px_i);
end
fc=fc(2:end);
end