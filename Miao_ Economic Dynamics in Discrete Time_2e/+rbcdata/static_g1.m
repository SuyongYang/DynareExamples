function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
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
%   g1
%

if T_flag
    T = rbcdata.static_g1_tt(T, y, x, params);
end
g1 = zeros(9, 9);
g1(1,2)=(-exp(y(2)))/(exp(y(2))*exp(y(2)))-T(4)*params(1)*(-exp(y(2)))/(exp(y(2))*exp(y(2)));
g1(1,3)=(-(T(1)*T(3)*exp(y(9))*params(4)*exp(y(3))*getPowerDeriv(exp(y(3)),params(4)-1,1)));
g1(1,5)=(-(T(1)*T(2)*T(8)));
g1(1,9)=(-(T(1)*T(2)*T(3)));
g1(2,2)=exp(y(2))*params(2)/(1-exp(y(5)));
g1(2,5)=(-(exp(y(2))*params(2)*(-exp(y(5)))))/((1-exp(y(5)))*(1-exp(y(5))));
g1(2,7)=(-exp(y(7)));
g1(3,3)=(-(T(6)*exp(y(9))*(1-params(4))*T(7)));
g1(3,5)=(-(exp(y(9))*(1-params(4))*T(5)*exp(y(5))*getPowerDeriv(exp(y(5)),(-params(4)),1)));
g1(3,7)=exp(y(7));
g1(3,9)=(-(exp(y(9))*(1-params(4))*T(5)*T(6)));
g1(4,1)=(-exp(y(1)));
g1(4,2)=exp(y(2));
g1(4,4)=exp(y(4));
g1(5,1)=exp(y(1));
g1(5,3)=(-(T(3)*exp(y(9))*T(7)));
g1(5,5)=(-(exp(y(9))*T(5)*T(8)));
g1(5,9)=(-(T(3)*exp(y(9))*T(5)));
g1(6,3)=(-(exp(y(3))-exp(y(3))*(1-params(3))));
g1(6,4)=exp(y(4));
g1(7,1)=(-(exp(y(1))/exp(y(5))));
g1(7,5)=(-((-(exp(y(5))*exp(y(1))))/(exp(y(5))*exp(y(5)))));
g1(7,6)=exp(y(6));
g1(8,1)=(-(params(4)*exp(y(1))/exp(y(3))));
g1(8,3)=(-((-(exp(y(3))*params(4)*exp(y(1))))/(exp(y(3))*exp(y(3)))));
g1(8,8)=exp(y(8));
g1(9,9)=1-params(5);
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
