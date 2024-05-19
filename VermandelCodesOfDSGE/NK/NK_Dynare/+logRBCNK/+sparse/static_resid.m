function [residual, T_order, T] = static_resid(y, x, params, T_order, T)
if nargin < 5
    T_order = -1;
    T = NaN(0, 1);
end
[T_order, T] = logRBCNK.sparse.static_resid_tt(y, x, params, T_order, T);
residual = NaN(12, 1);
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
