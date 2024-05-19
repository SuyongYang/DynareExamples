%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
% Discretization of AR(1) process using the Tauchen and Hussey (1991) method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% 
% [s,p,probst,arho,asigma]=tauch_hussey(xbar,rho,sigma,n)
% x_t=rho*x_t-1+(1-rho)*xbar+u_t
% u_t normal(0,sigma^2)
% 
% Input:
% xbar: mean of the process x
% rho:  persistence
% sigma: volatility
% n:     number of nodes
%
% Output: 
% s      : n by 1 vector of discretized state space of y_t 
% p      : Transition probability
% probst : Invariant distribution
% arho   : Theoretical first order autoregression coefficient for Markov chain
% asigma : Theoretical standard deviation for Markov chain Y_t

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s,p,probst,arho,asigma]=tauch_hussey(xbar,rho,sigma,n)

[xx,wx]=gauss_herm(n);		% nodes and weights for x
s=sqrt(2)*sigma*xx+xbar;     
x=xx(:,ones(n,1));
z=x';
w=wx(:,ones(n,1))';

p=(exp(z.*z-(z-rho*x).*(z-rho*x)).*w)./sqrt(pi);
sm=sum(p')';
p=p./sm(:,ones(n,1));

% calculate the invariant distribution of Markov chain

Trans= p';
probst = (1/n)*ones(n,1); % initial distribution of states
test = 1;

  while test > 10^(-8);
      probst1 = Trans*probst;
      test=max(abs(probst1-probst));
      probst = probst1;   
   end
   s       = s';
   meanm   = s*probst;             % mean of invariant distribution
   varm    = ((s-meanm).^2)*probst;  % variance of invariant distribution
    
   midaut1 = (s-meanm)'*(s-meanm); % cross product of deviation from the
                                   % mean of y_t and y_t-1
                                   
   probmat = probst*ones(1,n);     % each column is invariant distribution   
   
   
   midaut2 = p.*probmat.*midaut1; % product of the first two terms is 
                                     % the joint distribution of (Y_t-1,Y_t)
                                                                      
   autcov1 = sum(sum(midaut2));      % first-order auto-covariance
     
   arho    = autcov1/varm;              % approximated rho
   
   asigma  = sqrt(varm)*sqrt(1-arho^2);  % approximated sigma
   
   s=s';
