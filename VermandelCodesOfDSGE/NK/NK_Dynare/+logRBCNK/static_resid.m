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
    T = logRBCNK.static_resid_tt(T, y, x, params);
end
residual = zeros(12, 1);
    residual(1) = (0) - (y(3)-y(4));
    residual(2) = (y(5)) - (1/params(5)*y(7)+params(4)*y(2));
    residual(3) = ((y(3)-y(4))*params(13)) - (params(15)*y(6));
    residual(4) = (params(3)*y(9)) - (y(8)-y(8)*(1-params(3)));
    residual(5) = (y(1)) - (y(11)+y(8)*params(2)+y(7)*(1-params(2)));
    residual(6) = (y(10)) - (y(6)*params(2)+y(5)*(1-params(2))-y(11));
    residual(7) = (y(5)+y(7)) - (y(6)+y(8));
    residual(8) = (y(4)) - (y(4)*params(1)+y(10)*(1-params(10))*(1-params(1)*params(10))/params(10));
    residual(9) = (y(1)*params(14)) - (y(2)*params(16)+y(9)*params(17)+params(14)*params(6)*y(12));
    residual(10) = (y(3)) - (y(3)*params(7)+(1-params(7))*y(4)*params(8)+x(3));
    residual(11) = (y(11)) - (y(11)*params(11)+x(1));
    residual(12) = (y(12)) - (y(12)*params(12)+x(2));

end
