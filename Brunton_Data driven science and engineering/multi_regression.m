clear all;
load hald; % Load Portlant Cement dataset
A = ingredients;
b = heat;
[U,S,V] = svd(A,'econ');
x = V* inv(S)*U'*b; % Solve Ax=b using the SVD
plot(b,'k'); hold on
plot(A*x,'r-o'); % Plot data
% Plot fit
x1 = regress(b,A); % Alternative 1(regress)
x2 = pinv(A)*b; % Alternative 2(pinv)

%hold off;
%seq = 1:4;
%plot(seq,x,'gr-x');hold on
%plot(seq,x1,'k');
%plot(seq,x2,'b-o')

