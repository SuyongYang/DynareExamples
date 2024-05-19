
%Parametres exogenes RBC Cycles
beta 	= 0.99; %Discount Factor
delta 	= 0.025;	%Depreciation rate
alpha 	= 0.36;	%Capital share
gy 		= 0.2;    %Public spending in GDP
sigmaC 	= 2.5;
sigmaL 	= .5; %Elasticity of labor
hc 		= 0; % Consumption habit

% autoregressive roots parameters
rho_a	= 0.95;
rho_g	= 0.95;

%% Steady States
R		= 1/beta;
H		= 1/3;
Q		= 1;
Z		= (R*Q-(1-delta)*Q);
K		= H*(Z/alpha)^(1/(alpha-1));
Y		= K^alpha*H^(1-alpha);
I		= delta*K;
W		= (1-alpha)*Y/H;
C		= (1-gy)*Y;
chi		= W/((H^(1/sigmaL))*((C-hc*C)^sigmaC));

