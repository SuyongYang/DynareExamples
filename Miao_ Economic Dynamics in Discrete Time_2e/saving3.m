%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  THIS PROGRAM COMPUTES THE CONSUMPTION-SAVINGS PROBLEM 
%  USING THE DISCRETE STATE SPACE VALUE FUNCTION ITERATION METHOD                             
%  Add Interpolation
%      max E[\sum \beta^t u(c_t)]
% subject to 
% a_{t+1}+c_{t}=(1+r)*a_{t}+w*exp(z_t),
% a_{t} \geq -b,
% z_{t} = \rho*z_{t-1}+ \epsilon_{t}.
% 
% BY JIANJUN MIAO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all

%  set parameter values
%
gamma  =1.5;               		% risk aversion              
r      = 0.02;						% interet rate
beta   = 0.95;						% discount factor
b 	   = 0.0;					% Borrowing limit
w	   = 1;                     % wage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discretize AR(1) shock using Tauchen and Hussey (1991) method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$

nz		 = 2;
rho    = 0.8;
sigma  = 0.2;
[gridz, prob,probst,arho,asigma] = tauch_hussey(0,rho,sigma,nz);
%[gridz,prob,probst,arho,asigma]=markovappr(rho,sigma,3,nz)

income  = exp(gridz)*w;

%income=[0.1,1];
%prob=[0.5 0.5; 1-0.925 0.925]

% form asset holdings grid
na        = 100;     
nap       =2000;
amax      = 2;                     	  % maximum value of asset holdings grid  
amin      = -b;                     % borrowing constraint

a     = linspace(amin,amax,na)';	      % grid assets 
aprime = linspace(amin,amax,nap)';        % asset choice grid

%  initialize some variables
   
   v       = zeros(na,nz);
   decis   = zeros(na,nz);     
   cons    = zeros(nap,na,nz);
   util    = zeros(nap,na,nz); 
   tv      = zeros(na,nz);          % updated value function
   tdecis  = zeros(na,nz);        % updated decision rule
   
   
   %  tabulate the utility function such that for zero or negative
   %  consumption utility remains a large negative number so that
   %  such values will never be chosen as utility maximizing      

   for j=1:nz            
      cons(:,:,j) = repmat(income(j) + (1+r)*a',nap,1) - repmat(aprime,1,na);        
      util(:,:,j) = (cons(:,:,j).^(1-gamma)-1)/(1-gamma);
      utilj 		= util(:,:,j);
      Aj 			= cons(:,:,j);                 % give negative value to 
      i		 		= find( Aj <= 0);              % infeasible consumption choice
      utilj(i) 	= -10000;
      util(:,:,j) = utilj;
   end
   
   %  iterate on Bellman's equation and get the decision 
   %  rules and the value function at the optimum         
   
       test  = 1;
       iter=1;
 while test ~= 0;
          
   	for j=1:nz      
         vi		       = interp1(a,v,aprime,'spline'); %nbc by nba matrix
      [t1,t2] 		= max(util(:,:,j) + beta*repmat(vi*prob(j,:)',1,na));
      tv(:,j) 		=	t1';
      tdecis(:,j)	=	t2';
   	end
   
      test		= max(any(tdecis-decis));
      v			= tv;
      decis	   = tdecis;
      disp(sprintf('Iteration # %2d \tCriterion: %g',iter,test));
      iter=iter+1;
 end;
 
 for j=1:nz;
     ap(:,j)=aprime(decis(:,j));
end
   
   figure
                
   plot(a',ap(:,1),a',ap(:,nz),'--',a,a,'-.','linewidth',2)
   title('Policy Function')
   xlabel('asset of current period')
   ylabel('asset of next period')
   grid on
   legend('Low','High','Location','SouthEast')
   axis([amin amax amin amax])