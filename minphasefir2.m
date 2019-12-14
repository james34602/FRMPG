function b = minphasefir2(n, f, m)
% Copyright (C) 2019 James Fung <james34602@gmail.com>
% Copyright (C) 2000 Paul Kienzle <pkienzle@users.sf.net>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; see the file COPYING. If not, see
% <https://www.gnu.org/licenses/>.
%
% Produce an order @var{n} minimum phase FIR filter with arbitrary frequency
% response @var{m} over frequency bands @var{f}, returning the @var{n}+1 filter
% coefficients in @var{b}.  The vector @var{f} specifies the frequency band
% edges of the filter response and @var{m} specifies the magnitude response
% at each frequency.
%
% The vector @var{f} must be nondecreasing over the range [0,1], and the
% first and last elements must be 0 and 1, respectively.  A discontinuous
% jump in the frequency response can be specified by duplicating a band edge
% in @var{f} with different values in @var{m}.
t = length(f);
if t<2 || f(1)~=0 || f(t)~=1 || any(diff(f)<0)
    error ("minphasefir2: frequency must be nondecreasing starting from 0 and ending at 1");
elseif t ~= length(m)
    error ("minphasefir2: frequency and magnitude vectors must be the same length");
end
n = n + 1;
if n < 512
    grid_n = 512;
else
    grid_n = 2.^ceil(log(n)/log(2));
end
ramp_n = fix (grid_n / 25);
% Apply ramps to discontinuities
if (ramp_n > 0)
    % remember original frequency points prior to applying ramps
    basef = f(:);
    basem = m(:);
    % separate identical frequencies, but keep the midpoint
    idx = find (diff(f) == 0);
    f(idx) = f(idx) - ramp_n/grid_n/2;
    f(idx+1) = f(idx+1) + ramp_n/grid_n/2;
    f = [f(:);basef(idx)]';
    % make sure the grid points stay monotonic in [0,1]
    f(f<0) = 0;
    f(f>1) = 1;
    tmp = basef(idx);
    f = unique([f(:); tmp(:)]');
    % preserve window shape even though f may have changed
    m = interp1(basef, basem, f, 'pchip');
end
% interpolate between grid points
grid = interp1(f,m,linspace(0,1,grid_n+1), 'pchip');
H = mps([grid grid(grid_n-1:-1:2)]);
ht = real(ifft(H));
b = ht(1:(n-1));
end
function [rw] = fold(r)
[m,n] = size(r);
if m*n ~= m+n-1
    error('minphaseminphasefir2:fold.m: input must be a vector');
end
flipped = 0;
if (m > n)
    n = m;
    r = r.';
    flipped = 1;
end
if n < 3, rw = r; return;
elseif mod(n,2)==1
    nt = (n+1)/2;
    rw = [ r(1), r(2:nt) + conj(r(n:-1:nt+1)), ...
        0*ones(1,n-nt) ];
else
    nt = n/2;
    rf = [r(2:nt),0];
    rf = rf + conj(r(n:-1:nt+1));
    rw = [ r(1) , rf , 0*ones(1,n-nt-1) ];
end
if flipped
    rw = rw.';
end
end
function [clipped] = clipdb(s,cutoff)
clipped = s;
as = abs(s);
mas = max(as(:));
if mas==0, return; end
if cutoff >= 0, return; end
thresh = mas*10^(cutoff/20); % db to linear
clipped = s;
clipped(as < thresh) = thresh;
end
function [sm] = mps(s)
% [sm] = mps(s)
% create minimum-phase spectrum sm from complex spectrum s
sm = exp(fft(fold(ifft(log(clipdb(s,-140))))));
end