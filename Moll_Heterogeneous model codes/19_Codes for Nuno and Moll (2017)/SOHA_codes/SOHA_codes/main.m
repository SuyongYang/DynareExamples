%SIMULATES AN AIYAGARI ECONOMY WITH IDIOSYNCRATIC BROWNIAN MOTION AND RANDOM DEATHS

%GALO NUNO, based on codes from BENJAMIN MOLL
%Optimized for speed by SeHyoun Ahn
%The algorithm is based on a relaxation scheme to find K. The value function is
%found every round by solving the HJB equation through an upwind finite
%differences scheme. The distribution is found by also solving using finite
%differences the Fokker-Planck (Kolmogorov forward) equation

%--------------------------------------------------
%MODE
clear results parameters
%mode            = 0;              % 0 Competitive equilibrium, 1 efficient allocation

%PARAMETERS
parameters.psi   = 0.0;           % long-run growth rate
parameters.ga    = 2;             % CRRA utility with parameter gamma
parameters.alpha = 0.36;          % Production function F = K^alpha * L^(1-alpha) 
parameters.eta   = 0.02;          % death probability
parameters.delta = 0.08;          % Capital depreciation 10% - death rate
parameters.zmean = 1.038;         % mean O-U process (in levels). This parameter has to be adjusted to ensure that the mean of z (truncated gaussian) is 1.
parameters.Corr  = 1-0.6;         % persistence   O-U
parameters.sig2  = 0.2^2 * (1 - (1-parameters.Corr)^2);   % sigma^2 O-U
parameters.rho   = 0.04;          % discount rate
parameters.I     = 300;           % number of a points 
parameters.J     = 40;            % number of z points 
parameters.zmin  = 0.2;           % Range z
parameters.zmax  = 1.8;
parameters.amin  = 0;             % borrowing constraint
parameters.amax  = 100;           % range a

%simulation parameters
parameters.relax   = 0.90;         % relaxation parameter 
parameters.maxit   = 100;          % maximum number of iterations in the HJB loop
parameters.maxitK  = 1000;         % maximum number of iterations in the K loop
parameters.crit    = 10^(-6);      % criterion HJB loop
parameters.critK   = 1e-4;         %  criterion K loop
parameters.Delta   = 1000;         % delta in HJB algorithm
parameters.maxitla = 100;         % maximum number of iterations in the lambda loop
parameters.critla  = 1e-4;        %  criterion lambda loop
parameters.relaxla = 0.99;        % relaxation parameter lambda loop
K_guess            = 5;          % initial aggregate capital. It is important to guess a value close to the solution for the algorithm to converge

%Steady State
tic

%Competitive equilibrium
mode      = 0;
lambda    = 0;
K         = equilibrium_SS(parameters,K_guess,lambda);
results   = HJB_SS(parameters,K,lambda);
results   = KFE_SS(parameters,results);
results.K = K;
save competitive
g0        = results.g;
save g0 g0
plotting;
Moments
toc
%%
% %Optimal planning
mode      = 1;
lambda    = 0.0223; %Guess
lambda    = lagrange_SS(parameters,lambda,K_guess);
K         = equilibrium_SS(parameters,K_guess,lambda);
results   = HJB_SS(parameters,K,lambda);
results   = KFE_SS(parameters,results);
results.K = K;
plotting;
Moments
save planner
toc

%%
% % %Loop sensitivity
parameters.maxitla = 1;        
lambda_i = linspace(0.0,0.028,30);
lambda   = zeros(size(lambda_i)); 
for i =1:length(lambda_i)
    lambda(i) = lagrange_SS(parameters,lambda_i(i),K_guess);
end
save sensitivity
figure
plot(lambda_i,lambda,'Linewidth',2)
hold on
plot(lambda_i,lambda_i,'Linewidth',2)
xlabel('Lagrange multiplier, $\lambda$','FontSize',16,'interpreter','latex')
ylabel('$T\lambda$','FontSize',16,'interpreter','latex')
text(0.01,0.01,'Local maximum','FontSize',16,'interpreter','latex')
text(0.01,0.01,'Global maximum','FontSize',16,'interpreter','latex')

