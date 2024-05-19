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
    T = logRBCNK.static_g1_tt(T, y, x, params);
end
g1 = zeros(12, 12);
g1(1,3)=(-1);
g1(1,4)=1;
g1(2,2)=(-params(4));
g1(2,5)=1;
g1(2,7)=(-(1/params(5)));
g1(3,3)=params(13);
g1(3,4)=(-params(13));
g1(3,6)=(-params(15));
g1(4,8)=(-(1-(1-params(3))));
g1(4,9)=params(3);
g1(5,1)=1;
g1(5,7)=(-(1-params(2)));
g1(5,8)=(-params(2));
g1(5,11)=(-1);
g1(6,5)=(-(1-params(2)));
g1(6,6)=(-params(2));
g1(6,10)=1;
g1(6,11)=1;
g1(7,5)=1;
g1(7,6)=(-1);
g1(7,7)=1;
g1(7,8)=(-1);
g1(8,4)=1-params(1);
g1(8,10)=(-((1-params(10))*(1-params(1)*params(10))/params(10)));
g1(9,1)=params(14);
g1(9,2)=(-params(16));
g1(9,9)=(-params(17));
g1(9,12)=(-(params(14)*params(6)));
g1(10,3)=1-params(7);
g1(10,4)=(-((1-params(7))*params(8)));
g1(11,11)=1-params(11);
g1(12,12)=1-params(12);

end
