x = 3;
a = [-2:.25:2]';
b = a*x + 1* randn(size(a));
plot(a,x*a,'k')
hold on, plot(a,b,'rx')
[U,S,V] = svd(a,'econ');
xtilde = V* inv(S)*U'*b;
plot(a,xtilde*a,'b--')