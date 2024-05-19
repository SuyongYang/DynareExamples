syms y(t)
[V] = odeToVectorField(diff(y, 2) == diff(y) + 2*y + 0.02*exp(-t));
M = matlabFunction(V,'vars', {'t','Y'});
sol = ode45(M,[0 20],[0 0]);

fplot(@(x)deval(sol,x,1), [0, 20])
