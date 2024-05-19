%{
=============================================================================================
Hopenhayn Model in Continuous Time - Matlab code
Filename: hopenhayn.m
Scope: main script

Author: Riccardo A. Cioffi
Date: 02/12/2017
============================================================================================= 
%}

clear all;
clc;
close all;

% Initial guess for (p,w)
p = 0.5;
w = 1;

% Numerical parameters (on productivity grid)
I    = 1000;
zmin = 0.0;
zmax = 1;
z    = linspace(zmin,zmax,I)';
dz   = z(2) - z(1);
dz2  = dz^2;
dt   = 0.001;

% Model Parameters
rho    = 0.05;
alpha  = 0.5;
eps    = 0.5;
phi    = 0.5; 
vstar  = 0.0;
c_f    = 0.05;
c_e    = 0.6;
m_bar  = 0.1;
eta    = 1000;

% Entry distribution
psi_lb                   = 0.7;
psi_lb_idx               = ceil(psi_lb*I);
psi_dist                 = zeros(I,1);
psi_dist(psi_lb_idx:end) = 1/((I-psi_lb_idx+1)*dz);

% Loop objects
relax_w  = 0.2;   %Relaxation parameter on wage loop
relax_p  = 0.001; %Relaxation parameter on price loop
max_iter = 10000;
tol      = 1e-5;

% Wage loop
iter_w     = 0;
stop_w     = 0;
while ~stop_w
    iter_w = iter_w + 1;
    
    % Price loop
    iter_p     = 0;
    stop_p     = 0;
    while ~stop_p
        iter_p = iter_p + 1;

        % Firm choices
        n  = (alpha.*z.*p./w).^(1/(1-alpha));
        f  = z.*n.^alpha;
        pi = p.*f - w.*n - c_f;

        % Scrap value of the firm
        scrap = vstar*ones(I,1);

        % diffusion parameters
        mu_bar  = -0.01;
        mu      = ones(I,1).*mu_bar;

        sig_bar = 0.01;
        sig     = z.*sig_bar;
        sig2    = sig.^2;

        %Construct matrix
        X = -min(mu,0)/dz + sig2/(2*dz2);
        Y = -max(mu,0)/dz + min(mu,0)/dz - sig2/dz2;
        Z = max(mu,0)/dz + sig2/(2*dz2);
        A = spdiags(Y,0,I,I) + spdiags(X(2:I),-1,I,I) + spdiags([0;Z(1:I-1)],1,I,I);

        %Manually adjust the boundary values at the corners.
        A(1,1) = Y(1) + X(1); %Reflecting barrier
        A(I,I) = Y(I) + Z(I); %Reflecting barrier

        %Solve option problem as LCP
        B = rho*speye(I) - A;
        q = -pi + B*scrap; 

        %TRY DIFFERENT ALGORITHMS FOR SOLVING LCP

        %Test 1: Yuval Tassa's Newton-based LCP solver, download from http://www.mathworks.com/matlabcentral/fileexchange/20952
        y0 = zeros(I,1); l = zeros(I,1); u = Inf*ones(I,1);
        y = LCP(B,q,l,u,y0,0);

        %Test 2: Andreas Almqvist pivoting (Lemke) LCP solver, download from http://www.mathworks.com/matlabcentral/fileexchange/41485
        %[w,z,retcode] = LCPSolve(B,q);

        % Check if the solution satisfies the LCP problem
        LCP_error = max(abs(y.*(B*y + q)));
        if LCP_error > 10^(-6)
            error('LCP not solved')
        end

        V_LCP = y+scrap; %calculate value function

        x_idx           = find(V_LCP-vstar,1,'first'); %exit threshold
        exit            = zeros(I,1);
        exit(1:x_idx-1) = 1;

        m = m_bar.*exp(eta.*(sum(V_LCP.*psi_dist.*dz) - c_e));

        AA = A;
        AA(:,1:x_idx-1) = 0;
        AA(1:x_idx-1,1:x_idx-1) = speye(x_idx-1);
        AAT = AA';

        gg = -AAT\(m.*psi_dist);

        % Aggregates
        N = max(sum(n.*gg.*dz),1e-4);
        Q = max(sum(f.*gg.*dz),1e-4);

        p_implied = Q^(-eps);
        
        % Convergence checks
        err_p = max(abs(p - p_implied));
        p = relax_p*p_implied + (1-relax_p)*p;
        fprintf('    p Iteration number = %4d; p = %.3f; Error = %3.0e; m = %.3g \n',iter_p,p,err_p,m);
        if err_p < tol || iter_w > max_iter
            stop_p = 1;
        end   
    end
    
    w_implied = N^phi;
    
    % Convergence checks
    err_w = max(abs(w - w_implied));
    w = relax_w*w_implied + (1-relax_w)*w;
    fprintf('w Iteration number = %4d; w = %.3f; Error = %3.0e; m = %.3g \n',iter_w,w,err_w,m);
    if err_w < tol || iter_w > max_iter
        stop_w = 1;
    end
end
       
figure()
plot(z,V_LCP,z,scrap,'--','LineWidth',2)
set(gca,'FontSize',16,'TickLabelInterpreter', 'latex')
legend({'$v(z)$','$v^*$'},'interpreter','latex','Location','NorthWest')
xlabel('$z$','interpreter','latex')

figure()
plot(z,gg*dz,'LineWidth',2)
set(gca,'FontSize',16,'TickLabelInterpreter', 'latex')
legend({'$g(z)$'},'interpreter','latex','Location','NorthWest')
xlabel('$z$','interpreter','latex')















