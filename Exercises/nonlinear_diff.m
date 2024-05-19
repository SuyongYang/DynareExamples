n=390;
u=linspace(-4,4,n);
v=linspace(-4,4,n);
[x,y] = meshgrid(u,v);
px = 0.2*ones(size(y));
py = x.^2 ./ (ones(size(y))-y.^2);
quiver(x,y,px,py,10);
xlim([-4 4]),ylim([-4 4])