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
    T = logRBCNK.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(12, 23);
g1(1,7)=(-params(4));
g1(1,18)=params(4);
g1(1,8)=(-1);
g1(1,19)=1;
g1(2,7)=(-params(4));
g1(2,10)=1;
g1(2,12)=(-(1/params(5)));
g1(3,8)=params(13);
g1(3,19)=(-params(13));
g1(3,20)=(-params(15));
g1(4,3)=1-params(3);
g1(4,13)=(-1);
g1(4,14)=params(3);
g1(5,6)=1;
g1(5,12)=(-(1-params(2)));
g1(5,3)=(-params(2));
g1(5,16)=(-1);
g1(6,10)=(-(1-params(2)));
g1(6,11)=(-params(2));
g1(6,15)=1;
g1(6,16)=1;
g1(7,10)=1;
g1(7,11)=(-1);
g1(7,12)=1;
g1(7,3)=(-1);
g1(8,9)=1;
g1(8,19)=(-params(1));
g1(8,15)=(-((1-params(10))*(1-params(1)*params(10))/params(10)));
g1(9,6)=params(14);
g1(9,7)=(-params(16));
g1(9,14)=(-params(17));
g1(9,17)=(-(params(14)*params(6)));
g1(10,1)=(-((1-params(7))*(-params(9))));
g1(10,6)=(-((1-params(7))*params(9)));
g1(10,2)=(-params(7));
g1(10,8)=1;
g1(10,9)=(-((1-params(7))*params(8)));
g1(10,23)=(-1);
g1(11,4)=(-params(11));
g1(11,16)=1;
g1(11,21)=(-1);
g1(12,5)=(-params(12));
g1(12,17)=1;
g1(12,22)=(-1);

end
