function lambda = lagrange_SS(parameters,lambda,K_guess)
%function error_lambda = lagrange_SS(parameters,lambda,K)
% %Computes the error in the estimation of teh SS lagrange multiplier
% K       = equilibrium_SS(parameters,K,lambda);
% results = HJB_SS(parameters,K,lambda);
% results = KFE_SS(parameters,results);
% I       = parameters.I;                % numer of points
% J       = parameters.J;                % numer of points
% alpha   = parameters.alpha;         % Production function F = K^alpha * L^(1-alpha) 
% 
% da   = results.da;
% dz   = results.dz;
% g    = results.g;
% aa   = results.aa;
% zz   = results.zz;
% V    = results.V;
% 
% Dg_a = (g(2:I,:)-g(1:I-1,:))/da;        %Forward difference
% 
% %Lagrange multipliers
% lambda1 =  sum(sum(V(1:I-1,:).*(g(1:I-1,:) + aa(1:I-1,:).*Dg_a )))*da*dz;
% lambda2 = -sum(sum(V(1:I-1,:).*zz(1:I-1,:).*Dg_a))*da*dz *K;
% error_lambda  = (lambda1 + lambda2) * alpha * (1-alpha) *K^(alpha-2) - lambda;

%Computes the error in the estimation of teh SS lagrange multiplier
I       = parameters.I;             % numer of points
alpha   = parameters.alpha;         % Production function F = K^alpha * L^(1-alpha)
maxitla = parameters.maxitla;
critla  = parameters.critla;
relaxla  = parameters.relaxla;  

%%%%%%%%%%%%%%%%%%%
K    = K_guess;

for it=1:maxitla
    K = equilibrium_SS(parameters,K,lambda);
    results   = HJB_SS(parameters,K,lambda);
    results   = KFE_SS(parameters,results);
    results.K = K;
    
    da   = results.da;
    dz   = results.dz;
    g    = results.g;
    aa   = results.aa;
    zz   = results.zz;
    V    = results.V;
    
    Dg_a = (g(2:I,:)-g(1:I-1,:))/da;        %Forward difference
    
    %Lagrange multipliers
    lambda1 =  sum(sum(V(1:I-1,:).*(g(1:I-1,:) + aa(1:I-1,:).*Dg_a )))*da*dz;
    lambda2 = -sum(sum(V(1:I-1,:).*zz(1:I-1,:).*Dg_a))*da*dz *K;
    
    lambda_next = (lambda1 + lambda2) * alpha * (1-alpha) *K^(alpha-2);
    
    
    if abs(lambda-lambda_next)<critla
        break
    end
    
    lambda   = relaxla* lambda + (1-relaxla)* lambda_next;
    
    disp('ITERATION LAMBDA = ')
    disp(it)
    disp('ERROR LAMBDA =')
    disp((abs(lambda-lambda_next)))
     disp('LAMBDA =')
    disp((lambda))
end