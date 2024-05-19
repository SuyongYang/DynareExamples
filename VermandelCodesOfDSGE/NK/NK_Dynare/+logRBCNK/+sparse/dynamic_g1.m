function [g1, T_order, T] = dynamic_g1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T_order, T)
if nargin < 9
    T_order = -1;
    T = NaN(0, 1);
end
[T_order, T] = logRBCNK.sparse.dynamic_g1_tt(y, x, params, steady_state, T_order, T);
g1_v = NaN(44, 1);
g1_v(1)=(-((1-params(7))*(-params(9))));
g1_v(2)=(-params(7));
g1_v(3)=1-params(3);
g1_v(4)=(-params(2));
g1_v(5)=(-1);
g1_v(6)=(-params(11));
g1_v(7)=(-params(12));
g1_v(8)=1;
g1_v(9)=params(14);
g1_v(10)=(-((1-params(7))*params(9)));
g1_v(11)=(-params(4));
g1_v(12)=(-params(4));
g1_v(13)=(-params(16));
g1_v(14)=(-1);
g1_v(15)=params(13);
g1_v(16)=1;
g1_v(17)=1;
g1_v(18)=(-((1-params(7))*params(8)));
g1_v(19)=1;
g1_v(20)=(-(1-params(2)));
g1_v(21)=1;
g1_v(22)=(-params(2));
g1_v(23)=(-1);
g1_v(24)=(-(1/params(5)));
g1_v(25)=(-(1-params(2)));
g1_v(26)=1;
g1_v(27)=(-1);
g1_v(28)=params(3);
g1_v(29)=(-params(17));
g1_v(30)=1;
g1_v(31)=(-((1-params(10))*(1-params(1)*params(10))/params(10)));
g1_v(32)=(-1);
g1_v(33)=1;
g1_v(34)=1;
g1_v(35)=(-(params(14)*params(6)));
g1_v(36)=1;
g1_v(37)=params(4);
g1_v(38)=1;
g1_v(39)=(-params(13));
g1_v(40)=(-params(1));
g1_v(41)=(-params(15));
g1_v(42)=(-1);
g1_v(43)=(-1);
g1_v(44)=(-1);
if ~isoctave && matlab_ver_less_than('9.8')
    sparse_rowval = double(sparse_rowval);
    sparse_colval = double(sparse_colval);
end
g1 = sparse(sparse_rowval, sparse_colval, g1_v, 12, 39);
end
