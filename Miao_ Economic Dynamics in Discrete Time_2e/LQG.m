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
% When using this code, please cite 
% Miao, Jianjun and Jieran Wu, 2020, A Matlab Toolbox to Solve Dynamic Multivariate RI Problems in the LQG Framework, working paper, Boston University
% Miao, Jianjun, Jieran Wu, and Eric Young, 2020, Multivariate Rational Inattention, working paper, Boston University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P,F,err,iter] = LQG(A,B,S,Q,R,beta,P0)

  maxiter=10000; 
  tol = 10^(-10);
  err = 1;
  iter = 0;
  while (err>tol & iter< maxiter);
    P = Q+beta*A'*P0*A-(beta*A'*P0*B+S)/(R+beta*B'*P0*B)*(beta*B'*P0*A+S');
    err = norm(P0-P,'fro');
    P0 = P; % Update
    iter=iter+1;
  end;
  F =(R+beta*B'*P*B)\(beta*B'*P*A+S');
end