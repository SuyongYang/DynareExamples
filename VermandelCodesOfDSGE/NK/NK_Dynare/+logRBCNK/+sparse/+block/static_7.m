function [y, T, residual, g1] = static_7(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(5, 1);
  residual(1)=(params(3)*y(9))-(y(8)-y(8)*(1-params(3)));
  residual(2)=(y(1))-(y(11)+y(8)*params(2)+y(7)*(1-params(2)));
  residual(3)=(y(5)+y(7))-(y(6)+y(8));
  residual(4)=(y(1)*params(14))-(y(2)*params(16)+y(9)*params(17)+params(14)*params(6)*y(12));
  residual(5)=(y(5))-(1/params(5)*y(7)+params(4)*y(2));
if nargout > 3
    g1_v = NaN(12, 1);
g1_v(1)=(-(1-(1-params(3))));
g1_v(2)=(-params(2));
g1_v(3)=(-1);
g1_v(4)=1;
g1_v(5)=params(14);
g1_v(6)=(-(1-params(2)));
g1_v(7)=1;
g1_v(8)=(-(1/params(5)));
g1_v(9)=params(3);
g1_v(10)=(-params(17));
g1_v(11)=(-params(16));
g1_v(12)=(-params(4));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 5, 5);
end
end
