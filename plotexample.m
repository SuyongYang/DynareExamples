x=linspace(0.1,100,10000);
f=@(x) exp(-(log(x)-1).^2./2)*1./x;

plot(x,f(x))
hold on

f=@(x) exp(-(log(x)-2).^2./2)*1./x;
plot(x,f(x))

f=@(x) exp(-(log(x)-3).^2./2)*1./x;
plot(x,f(x))

f=@(x) exp(-(log(x)-4).^2./2)*1./x;
plot(x,f(x))

f=@(x) exp(-(log(x)-5).^2./2)*1./x;
plot(x,f(x))
