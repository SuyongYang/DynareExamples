function K = equilibrium_SS(parameters,K,lambda)
%Computes the capital in the model

maxitK = parameters.maxitK;
critK  = parameters.critK;
relax  = parameters.relax;


for n=1:maxitK
    results = HJB_SS(parameters,K,lambda);
    results = KFE_SS(parameters,results);
    da   = results.da;
    dz   = results.dz;
    g    = results.g;
    a    = results.a;
    
    S = sum(g'*a*da*dz);
    if abs(K-S)<critK
        break
    end
    %update
    K = relax*K +(1-relax)*S;           %relaxation algorithm (to ensure convergence)
    
    
%     disp('iteration = ')
%     disp(n)
    if mod(n,20)==0
        disp('capital =')
        disp(K)
    end
end

if n == maxitK
    disp('Equilibrium NOT Converged, ERROR ')
end