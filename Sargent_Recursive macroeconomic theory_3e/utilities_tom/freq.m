function s = freq(b,a,n)
% FREQ	Frequency response of a digital filter.
%	FREQ(B,A,N) returns the n-point frequency response of filter B/A:
%
%		                      -1              -(nb-1) 
%		   B(z)   b(1) + b(2)z + .... + b(nb)z
%		   ---- = ---------------------------
%		                      -1              -(na-1)
%		   A(z)    1   + a(2)z + .... + a(na)z
%
%	given the numerator and denominator coefficients in vectors B and A. 
%	If N is not a power of two it will be rounded up to the next power.
%	The frequency response is evaluated at N points equally spaced around
%	the unit circle.  The magnitude and phase can be plotted with the
%	commands:
%			s = freq(b,a,n);
%			mag = abs(s);
%			phase = imag(log(s));
%			f = pi*(0:n-1)/n;
%			semilogy(f,mag), plot(f,phase)
%
%	See also FDESIGN, FILTER, and FFT.
na = max(size(a));
nb = max(size(b));
s  = fft([b zeros(1,n-nb)]) ./ fft([a zeros(1,n-na)]);
