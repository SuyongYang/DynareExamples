n=9
lamda=2
t=linspace(0,100,1000)
f=@(t) lamda*exp(-lamda*t)*(lamda*t)^(n-1)/factorial(n-1)
plot(t,f(t))