function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = rbcdata.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(9, 1);
lhs = 1/exp(y(4));
rhs = T(1)*T(4);
residual(1) = lhs - rhs;
lhs = exp(y(4))*params(2)/(1-exp(y(7)));
rhs = exp(y(9));
residual(2) = lhs - rhs;
lhs = exp(y(9));
rhs = (1-params(4))*exp(y(11))*T(5)*T(6);
residual(3) = lhs - rhs;
lhs = exp(y(4))+exp(y(6));
rhs = exp(y(3));
residual(4) = lhs - rhs;
lhs = exp(y(3));
rhs = exp(y(11))*T(5)*T(7);
residual(5) = lhs - rhs;
lhs = exp(y(6));
rhs = exp(y(5))-exp(y(1))*(1-params(3));
residual(6) = lhs - rhs;
lhs = exp(y(8));
rhs = exp(y(3))/exp(y(7));
residual(7) = lhs - rhs;
lhs = exp(y(10));
rhs = params(4)*exp(y(3))/exp(y(1));
residual(8) = lhs - rhs;
lhs = y(11);
rhs = params(5)*y(2)+x(it_, 1);
residual(9) = lhs - rhs;

end
