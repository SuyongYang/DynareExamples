clear all;
load housing.data
b = housing(:,14);
% housing values in $1000s
A = housing(:,1:13);
% other factors,
A = [A ones(size(A,1),1)]; % Pad with ones y-intercept
x = regress(b,A);
plot(b,'k-o');
hold on, plot(A*x,'r-o');
[b sortind] = sort(housing(:,14)); % sorted values
plot(b,'k-o')
hold on, plot(A(sortind,:)*x,'r-o')

