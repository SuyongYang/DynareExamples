clear;
A=imread('gimono.jpeg');
X=double(rgb2gray(A)); % Convert RBG->gray, 256 bit->double.
nx = size(X,1); ny = size(X,2);

[U,S,V] = svd(X);
for r=[5 20 50 80 100 400]; % Truncation value
Xapprox = U(:,1:r)*S(1:r,1:r)*V(:,1:r)'; % Approx. image
figure, imagesc(Xapprox), axis on % axis off
title(['r=',num2str(r,'%d')]);
end

subplot(1,2,1), semilogy(diag(S), 'k')
subplot(1,2,2), plot(cumsum(diag(S))/sum(diag(S)),'k')