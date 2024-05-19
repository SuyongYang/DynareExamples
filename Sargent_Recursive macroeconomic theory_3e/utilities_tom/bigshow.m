function z = bigshow (varargin)

%  function z = bigshow (a,b)
%  function z = bigshow (a,b,tsim,timpulse)
% This function takes a ARMA process of the form:
%
% a*y_t = b*w_t 
%
%  and graphs a simulation, the impulse response,
%  and the spectrum.
%  e.g., for
%  y_t = .9 y_{t-1} + w_t 
%  the input vectors would be:
%
%  a=[1 -.9]; and b=1;
%
%  The function would then be called:
%  
%  >> bigshow (a,b)

%a=[1 -.99]; b=1;  %  near random walk
%a=[1 0 0 0 -.8]; b=1;  %
%a=[1 -.8]; b=[1 -.7];  % muth process
%a=[1 -1.3 +.7]; b=1;    %  business cycle process

% Set up the input variables.

if nargin < 2 | nargin >4
   error ('Error.  Wrong number of input arguments.')
else
   if nargin == 4
      timpulse = varargin{4};
   else
      timpulse = 26;
   end;
   if nargin >= 3;
      tsim = varargin{3};
   else
      tsim = 70;
   end
   a = varargin{1};
   b = varargin{2};
end;

% first generate the simulation (Gaussian errors)
nn=100;   % number of observations to be discarded
u=randn(tsim+nn,1);
y=filter(b,a,u);
subplot(2,2,1),
axis ([0 tsim -Inf Inf]);
plot(y(nn+1:length(y)));
title ('Simulation');

%  Now generate the impulse response

subplot(2,2,2);
nn=timpulse;
y2=dimpulse(b,a,nn);
xx=[0:nn-1];  % fix the x axis to start at zero
y2=y2(1:nn);
%axis ([0 timpulse -Inf Inf]);
plot(xx',y2)
axis ([0 timpulse -Inf Inf]);

title ('Impulse Response')

%  Now generate the spectrum

n=256;
subplot(2,2,3);
%axis ([0 pi -1 1]);
[spect2,f2]=shownew(b,a,n);
title ('Log Spectrum')
