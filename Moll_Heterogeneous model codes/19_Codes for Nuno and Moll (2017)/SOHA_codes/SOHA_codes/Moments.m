%Moments
disp('Capital')
disp(K)
disp('Output')
disp(K^parameters.alpha)
disp('K/Y')
disp(K^(1-parameters.alpha))
disp('C')
if mode ==1
results.c(parameters.I/2:end,:) = results.c(parameters.I/2,parameters.J/2); %
end
C = sum(sum(results.c.*results.g*results.da*results.dz));
disp(C)
disp('r')
disp(100*results.r)
disp('w')
disp(results.w)
disp('Tail exponent')
if mode ==0
    disp(parameters.eta*parameters.ga/(results.r - parameters.rho))
else
    disp(parameters.eta/(parameters.eta+results.r))
end
u = results.c.^(1-parameters.ga)/(1-parameters.ga);
Welfare = sum(sum(u.*results.g*results.da*results.dz));
disp('Welfare')
disp(Welfare)