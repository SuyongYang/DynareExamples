tspan = [-3,0];
z0=[-1,0]
[t,z]=ode45(@myode,tspan,z0);
plot(z(:,1),z(:,2))
hold on
tspan = [0,3];

z0=[z(end,1),z(end,2)]
[t,z]=ode45(@myode,tspan,z0);
plot(z(:,1),z(:,2))