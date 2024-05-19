function [y, T] = dynamic_1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(23)=params(11)*y(11)+x(1);
  y(24)=params(12)*y(12)+x(2);
end
