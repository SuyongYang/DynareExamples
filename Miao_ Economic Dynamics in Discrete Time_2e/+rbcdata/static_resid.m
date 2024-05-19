function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = rbcdata.static_resid_tt(T, y, x, params);
end
residual = zeros(9, 1);
lhs = 1/exp(y(2));
rhs = T(1)*T(4);
residual(1) = lhs - rhs;
lhs = exp(y(2))*params(2)/(1-exp(y(5)));
rhs = exp(y(7));
residual(2) = lhs - rhs;
lhs = exp(y(7));
rhs = exp(y(9))*(1-params(4))*T(5)*T(6);
residual(3) = lhs - rhs;
lhs = exp(y(2))+exp(y(4));
rhs = exp(y(1));
residual(4) = lhs - rhs;
lhs = exp(y(1));
rhs = T(3)*exp(y(9))*T(5);
residual(5) = lhs - rhs;
lhs = exp(y(4));
rhs = exp(y(3))-exp(y(3))*(1-params(3));
residual(6) = lhs - rhs;
lhs = exp(y(6));
rhs = exp(y(1))/exp(y(5));
residual(7) = lhs - rhs;
lhs = exp(y(8));
rhs = params(4)*exp(y(1))/exp(y(3));
residual(8) = lhs - rhs;
lhs = y(9);
rhs = y(9)*params(5)+x(1);
residual(9) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
