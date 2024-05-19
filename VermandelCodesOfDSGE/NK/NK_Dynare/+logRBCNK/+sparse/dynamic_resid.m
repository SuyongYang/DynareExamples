function [residual, T_order, T] = dynamic_resid(y, x, params, steady_state, T_order, T)
if nargin < 6
    T_order = -1;
    T = NaN(0, 1);
end
[T_order, T] = logRBCNK.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
residual = NaN(12, 1);
    residual(1) = (params(4)*(y(26)-y(14))) - (y(15)-y(28));
    residual(2) = (y(17)) - (1/params(5)*y(19)+params(4)*y(14));
    residual(3) = ((y(15)-y(28))*params(13)) - (params(15)*y(30));
    residual(4) = (params(3)*y(21)) - (y(20)-(1-params(3))*y(8));
    residual(5) = (y(13)) - (y(23)+y(8)*params(2)+y(19)*(1-params(2)));
    residual(6) = (y(22)) - (params(2)*y(18)+y(17)*(1-params(2))-y(23));
    residual(7) = (y(17)+y(19)) - (y(8)+y(18));
    residual(8) = (y(16)) - (y(28)*params(1)+y(22)*(1-params(10))*(1-params(1)*params(10))/params(10));
    residual(9) = (y(13)*params(14)) - (y(14)*params(16)+y(21)*params(17)+params(14)*params(6)*y(24));
    residual(10) = (y(15)) - (params(7)*y(3)+(1-params(7))*(y(16)*params(8)+params(9)*(y(13)-y(1)))+x(3));
    residual(11) = (y(23)) - (params(11)*y(11)+x(1));
    residual(12) = (y(24)) - (params(12)*y(12)+x(2));
end
