%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear Quadratic Gaussion Model under full information with discounting
% beta
% Quadratic return function x'*Q*x+u'*R*u+2x'*S*u
% State transition equation x_t+1 = A*x_t + B*u_t
%
% Output: 
% P:  Value function x*P*x
% F:  Decision rule u = -F*x
%
% @ Written by Jianjun Miao, Boston University, Email: miaoj@bu.edu 
%          
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P,F,err,iter] = LQGRobust(A,B,C,S,Q,R,beta,theta,P0)

  maxiter=10000; 
  tol = 10^(-10);
  err = 1;
  iter = 0;
  n = size(C,2);
  while (err>tol & iter< maxiter);
    D = P0+P0*C/(theta*eye(n)-C'*P0*C)*C'*P0;
    P = Q+beta*A'*D*A-(beta*A'*D*B+S)/(R+beta*B'*D*B)*(beta*B'*D*A+S');
    err = norm(P0-P,'fro');
    P0 = P; % UpdateC
    iter=iter+1;
  end;
  F =(R+beta*B'*D*B)\(beta*B'*D*A+S');
end