function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = rbcdata.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(9, 15);
g1(1,4)=(-exp(y(4)))/(exp(y(4))*exp(y(4)));
g1(1,12)=(-(T(4)*params(1)*(-exp(y(12)))/(exp(y(12))*exp(y(12)))));
g1(1,5)=(-(T(1)*T(3)*exp(y(14))*params(4)*exp(y(5))*getPowerDeriv(exp(y(5)),params(4)-1,1)));
g1(1,13)=(-(T(1)*T(2)*exp(y(13))*getPowerDeriv(exp(y(13)),1-params(4),1)));
g1(1,14)=(-(T(1)*T(2)*T(3)));
g1(2,4)=exp(y(4))*params(2)/(1-exp(y(7)));
g1(2,7)=(-(exp(y(4))*params(2)*(-exp(y(7)))))/((1-exp(y(7)))*(1-exp(y(7))));
g1(2,9)=(-exp(y(9)));
g1(3,1)=(-(T(6)*(1-params(4))*exp(y(11))*T(8)));
g1(3,7)=(-((1-params(4))*exp(y(11))*T(5)*exp(y(7))*getPowerDeriv(exp(y(7)),(-params(4)),1)));
g1(3,9)=exp(y(9));
g1(3,11)=(-((1-params(4))*exp(y(11))*T(5)*T(6)));
g1(4,3)=(-exp(y(3)));
g1(4,4)=exp(y(4));
g1(4,6)=exp(y(6));
g1(5,3)=exp(y(3));
g1(5,1)=(-(T(7)*exp(y(11))*T(8)));
g1(5,7)=(-(exp(y(11))*T(5)*exp(y(7))*getPowerDeriv(exp(y(7)),1-params(4),1)));
g1(5,11)=(-(exp(y(11))*T(5)*T(7)));
g1(6,1)=exp(y(1))*(1-params(3));
g1(6,5)=(-exp(y(5)));
g1(6,6)=exp(y(6));
g1(7,3)=(-(exp(y(3))/exp(y(7))));
g1(7,7)=(-((-(exp(y(7))*exp(y(3))))/(exp(y(7))*exp(y(7)))));
g1(7,8)=exp(y(8));
g1(8,3)=(-(params(4)*exp(y(3))/exp(y(1))));
g1(8,1)=(-((-(exp(y(1))*params(4)*exp(y(3))))/(exp(y(1))*exp(y(1)))));
g1(8,10)=exp(y(10));
g1(9,2)=(-params(5));
g1(9,11)=1;
g1(9,15)=(-1);

end
