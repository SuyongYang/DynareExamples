function [y, T, residual, g1] = dynamic_2(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(10, 1);
  residual(1)=(y(22))-(params(2)*y(18)+y(17)*(1-params(2))-y(23));
  residual(2)=(y(16))-(y(28)*params(1)+y(22)*(1-params(10))*(1-params(1)*params(10))/params(10));
  residual(3)=(y(13)*params(14))-(y(14)*params(16)+y(21)*params(17)+params(14)*params(6)*y(24));
  residual(4)=(y(17))-(1/params(5)*y(19)+params(4)*y(14));
  residual(5)=((y(15)-y(28))*params(13))-(params(15)*y(30));
  residual(6)=(params(3)*y(21))-(y(20)-(1-params(3))*y(8));
  residual(7)=(y(13))-(y(23)+y(8)*params(2)+y(19)*(1-params(2)));
  residual(8)=(y(17)+y(19))-(y(8)+y(18));
  residual(9)=(y(15))-(params(7)*y(3)+(1-params(7))*(y(16)*params(8)+params(9)*(y(13)-y(1)))+x(3));
  residual(10)=(params(4)*(y(26)-y(14)))-(y(15)-y(28));
if nargout > 3
    g1_v = NaN(34, 1);
g1_v(1)=(-params(7));
g1_v(2)=1-params(3);
g1_v(3)=(-params(2));
g1_v(4)=(-1);
g1_v(5)=(-((1-params(7))*(-params(9))));
g1_v(6)=(-(1-params(2)));
g1_v(7)=1;
g1_v(8)=1;
g1_v(9)=1;
g1_v(10)=(-((1-params(10))*(1-params(1)*params(10))/params(10)));
g1_v(11)=(-params(17));
g1_v(12)=params(3);
g1_v(13)=(-(1/params(5)));
g1_v(14)=(-(1-params(2)));
g1_v(15)=1;
g1_v(16)=params(13);
g1_v(17)=1;
g1_v(18)=(-1);
g1_v(19)=(-1);
g1_v(20)=params(14);
g1_v(21)=1;
g1_v(22)=(-((1-params(7))*params(9)));
g1_v(23)=(-params(2));
g1_v(24)=(-1);
g1_v(25)=1;
g1_v(26)=(-((1-params(7))*params(8)));
g1_v(27)=(-params(16));
g1_v(28)=(-params(4));
g1_v(29)=(-params(4));
g1_v(30)=(-params(15));
g1_v(31)=(-params(1));
g1_v(32)=(-params(13));
g1_v(33)=1;
g1_v(34)=params(4);
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 10, 30);
end
end
