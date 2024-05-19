% Poutineau - Vermandel (2017). A welfare analysis of macroprudential policy rules. Revue d'Economie Politique.
% Optimal policy results are slightly different than in the paper because parameters/shocks were rounded

close all;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var y_h y_f c_h c_f r_IB pi_h pi_f pic_h pic_f w_h w_f h_h h_f mc_P_h mc_P_f k_h k_f i_h i_f z_h z_f q_h q_f d_h d_f ToT omega_h omega_f n_h n_f r_K_h r_K_f
	ld_h ld_f s_h s_f k_n_h k_n_f v_h v_f ls_h ls_f r_L_h r_L_f mc_L_h mc_L_f bk_h bk_f 
	mp_h mp_f r_D_h r_D_f dd_h dd_f mu_D_h mu_D_f lb_c_h lb_c_f ac_I_h ac_I_f
	%shocks
	e_A_h e_G_h e_A_f e_G_f e_I_h e_I_f e_U_h e_U_f e_N_h e_N_f e_Q_h e_Q_f e_D_h e_D_f e_R
	eta_h eta_f ss_h ss_f d1_h d2_h delta_L_h r_L_star_h d1_f d2_f delta_L_f r_L_star_f s2_h s1_h delta_P_h pstar_h s2_f s1_f delta_P_f pstar_f f1_h f2_h delta_D_h r_D_star_h f1_f f2_f delta_D_f r_D_star_f
	pXI_h pXI_f  pX_h pX_f cg_h cg_f diffY_h diffY_f diffLs_h diffLs_f ut_h ut_f welf_h welf_f ha_h ha_f;

varexo	eta_A_h eta_G_h eta_U_h eta_I_h eta_A_f eta_G_f eta_U_f eta_I_f eta_N_h eta_N_f eta_Q_h eta_Q_f eta_D_h eta_D_f eta_R;

parameters	beta alpha delta theta_P_h theta_P_f chi_D_h chi_D_f chi_I_h chi_I_f 
			rho_A_h rho_G_h rho_A_f rho_G_f rho_U_f rho_U_h rho_I_h rho_I_f rho_N_h rho_N_f rho_Q_h rho_Q_f rho_D_h rho_D_f rho_AG_h rho_AG_f
			phi_r phi_y rho_R rho sigmaC_h sigmaC_f mu sigmaL_h sigmaL_f theta_L_h theta_L_f  epsilon_P
			alpha_I_h alpha_I_f alpha_C_h alpha_C_f R A H Z I Q K MC
			%phiP_h phiP_f phiB_h phiB_f
			Y C W xi_P_h xi_P_f hc_h hc_f
			% financial accelerator
			eta_s tau_E R_K RR_L RR R_L N_E eta_E eta_B OMEGA ka wmin w_sup 
			w_inf varkappa_h varkappa_f
			h_L_h h_L_f V_E
			mu_b epsilon_L
			% deposit 
			theta_D_h theta_D_f epsilon_D R_D
			% asymmetric
			n gy Ls_h Ls_f Ld C_h C_f MCB D_h D_f  BK_h BK_f tau_B 
			
			OA_h OA_f F1_h F1_f F2_h F2_f chi_h chi_f
			% MacroPrud
			phi_ly phi_l
			mp_cg_h mp_cg_f mp_cgy_h mp_cgy_f lambda_i
			mu_d U_h U_f WELF_h WELF_f D1_h D1_f D2_h D2_f HA_h HA_f S1_h S1_f S2_h S2_f
			S_B S 
			;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------


%% CALIBRATION
alpha   	= 0.36;							% share of capital in ouput
beta		= 1/(1+4/400);					% on the mean of the ECB rate discount factor
R_D			= 1/beta-1;						% deposit rate
R 			= R_D+1.5/400;					% taux ECB
delta 		= 0.025;						% Depreciation rate
gy 			= 0.24;    						% Public spending in GDP
chi_D_h 	= 0.07; 						% portfolio adjustment cost
chi_D_f 	= 0.07; 						% portfolio adjustment cost
mu			= 2.5; 							% subtituability final goods
alpha_I_h 	= 0.0439;						% investment goods
alpha_I_f 	= 0.0439;						% investment goods
alpha_C_h 	= 0.0922;						% final goods
alpha_C_f 	= 0.0922;						% final goods
eta_E		= .3; 							% leverage ratio N/K
eta_B		= .11;							% leverage ratio N/K bank
R_L			= R_D+4.3/400;					% ss value of Bank rates
mu_b 		= 0.12;   						% Recovery cost
epsilon_P 	= 10;							% markup

% Financial parameters
wmin		= 1-eta_E;
ka			= 1/eta_E;
eta_s 		= (1-.025)^(1/4);					% ss quarterly failure rate
% RBC ss
Q			= 1;								% asset price
A			= 1;								% productivity
P			= 1;								% price
H			= 1/3;								% share of hours worked by day
n			= .65;								% share of the home contry in EMU

% loss function
lambda_i	= 0.077; % Woodford

sigmaC_f	= 1.398962110194814;	% Consumption utility
sigmaC_h	= 1.208871515355716;	% Consumption utility
sigmaL_f	= 1.250519517923662;	% Labor utility
sigmaL_h	= 1.698886594538130;	% Labor utility
theta_P_h	= 0.709677752883044;	% Price rigidity
theta_P_f	= 0.617234185017766;	% Price rigidity
xi_P_f		= 0.341871779893139;	% Price indexation
xi_P_h		= 0.173915161239743;	% Price indexation
theta_L_f	= 0.585633061514279;	% Credit rate rigidity
theta_L_h	= 0.489723871935481;	% Credit rate rigidity
theta_D_f	= 0.699216407972913;	% Deposit rate rigidity
theta_D_h	= 0.733561373565272;	% Deposit rate rigidity
varkappa_f	= 0.071312250940474;	% External finance premium elasticity
varkappa_h	= 0.053453624488272;	% External finance premium elasticity
alpha_C_f	= 0.145873462049645;	% share of foreign goods in consumption basket
alpha_C_h	= 0.170912531682801;	% share of foreign goods in consumption basket
alpha_I_f	= 0.133144026072124;	% share of foreign goods in investment basket
alpha_I_h	= 0.107842411706345;	% share of foreign goods in investment basket
chi_D_f		= 0.107676254545201;	% cost of adjustment of deposits
chi_D_h		= 0.092859742652570;	% cost of adjustment of deposits
chi_I_f		= 5.522076893667387;	% cost of adjusting investment
chi_I_h		= 4.265681504454944;	% cost of adjusting investment
hc_f		= 0.379892538637070;	% habits formation - consumption
hc_h		= 0.172359059019540;	% habits formation - consumption
h_L_f		= 0.818331259491459;	% habits formation - loans
h_L_h		= 0.759467440973027;	% habits formation - loans
mu			= 2.439725451120586;	% substitution between home foreign consumption goods
phi_r		= 2.372897456365311;	% monetary policy inflation reaction
phi_y		= 0.161065433263715;	% monetary policy GDP growth reaction
rho			= 0.854672138986762;

rho_A_f		= 0.989939318337348;	% AR(1) home productivity
rho_A_h		= 0.948020431281802;	% AR(1) foreign productivity
rho_G_h		= 0.814950836750651;	% AR(1) home spending
rho_G_f		= 0.910778772923741;	% AR(1) foreign spending
rho_I_f		= 0.586153154740352;	% AR(1) home investment cost
rho_I_h		= 0.593646344607546;	% AR(1) foreign investment cost
rho_U_h		= 0.712453721934572;	% AR(1) home preferences
rho_U_f		= 0.644214859350958;	% AR(1) foreign preferences
rho_N_h		= 0.605987902199315;	% AR(1) home net wealth
rho_N_f		= 0.341831522760077;	% AR(1) foreign net wealth
rho_Q_h		= 0.786641123455303;	% AR(1) home efp
rho_Q_f		= 0.783607289533309;	% AR(1) foreign efp
rho_D_h		= 0.704972898799605;	% AR(1) home premium
rho_D_f		= 0.702257601151616;	% AR(1) foreign premium
rho_R		= 0.395559403909082;	% AR(1) monetary policy
rho_AG_h	= 0.815686820207272;	% home comptetitiveness effect on demand
rho_AG_f	= 0.805435186917309;	% foreign comptetitiveness effect on demand


% financial contracts ss
R_D			= 1/beta-1;
RR 			= 1+R; 						% ECB rate
mu_d		= R_D/R;					% deposit markdown
epsilon_D	= (R_D/R)/(R_D/R-1);		% substituability deposits
RR_L 		= 1 + R_L;					% Interest rate
S_B 		= (1+R_L)/(1+R);			% spread between ECB rate and FI rate
MCB			= (1+R)/(1-mu_b*(1-eta_s))-1;% Marginal cost of loans
epsilon_L	= R_L/(R_L-MCB);			% substituability loans
OMEGA		= wmin*(1+((1+R)/(1+MCB)-1)/mu_b)^(-1/ka);
R_K			= (RR_L/OMEGA)*(1-eta_E)-1;
Gw			= 1-((wmin/OMEGA)^(ka-1));
eta_d 		= 1-(wmin/OMEGA)^ka;
eta_s 		= (wmin/OMEGA)^ka;
w_sup 		= (ka/(ka-1))*OMEGA;
w_inf 		= (ka/(ka-1))*(((OMEGA^(1-ka))-(wmin^(1-ka)))/((OMEGA^(-ka))-(wmin^(-ka)))); % ou encore (ka/(ka-1))*((wmin-eta_s*OMEGA)/(1-eta_s))


% RBC ss
MC		= (epsilon_P-1)/epsilon_P;
Z		= ((1+R_K)*Q-(1-delta)*Q);
K		= H*(Z/(alpha*MC))^(1/(alpha-1));
Y		= K^alpha*H^(1-alpha);
I		= delta*K;
W		= (1-alpha)*MC*Y/H;
C		= (1-gy)*Y-I;% Consumption if n=.5
% Compute asymmetric consumption
rez		= get_consumption([Y;gy;alpha_C_h;alpha_C_f;alpha_I_h;alpha_I_f;I;n;C]);
C_h		= rez(1);
C_f		= rez(2);
% steady states financial system
N_E		= eta_E*Q*K;
Ld		= Q*K-N_E;
gamma	= w_sup^(1-varkappa_h);
V_E 	= (eta_s*OMEGA*K*(1+R_K)/(ka-1));
tau_E	= 1-N_E/(V_E);
L_RF 	= (1-eta_E-eta_B)*K;
Ls_h	= Ld;
Ls_f	= Ld;
D_h  	= .46*Ls_h;
D_f  	= .46*Ls_f;
BK_h  	= .11*Ls_h;
BK_f  	= .11*Ls_f;
tau_B	= 1-BK_h/((1-mu_b*(1-eta_s))*(1+R_L)*Ls_h-(1+R)*(Ls_h - D_h - BK_h)-(1+R_D)*D_h);

% Computing ss for 2nd order variables
DELTAC_h	= (C_h-hc_h*C_h)^-sigmaC_h;
DELTAC_f	= (C_f-hc_f*C_f)^-sigmaC_f;
chi_h		= DELTAC_h*W/(H^sigmaL_h);
chi_f		= DELTAC_f*W/(H^sigmaL_f);
S1_h		= Y/(1-beta*theta_P_h);
S1_f		= Y/(1-beta*theta_P_f);
S2_h		= Y*MC/(1-beta*theta_P_h);
S2_f		= Y*MC/(1-beta*theta_P_f);
D1_h		= Ls_h*R_L/(1-beta*theta_L_h);
D1_f		= Ls_f*R_L/(1-beta*theta_L_f);
D2_h		= Ls_h*MCB/(1-beta*theta_L_h);
D2_f		= Ls_f*MCB/(1-beta*theta_L_f);
F1_h		= D_h*R_D/(1-beta*theta_D_h);
F1_f		= D_f*R_D/(1-beta*theta_D_f);
F2_h		= D_h*R/(1-beta*theta_D_h);
F2_f		= D_f*R/(1-beta*theta_D_f);
U_h			= ((C_h-hc_h*C_h)^(1-sigmaC_h))/(1-sigmaC_h) - chi_h*(H^(1+sigmaL_h))/(1+sigmaL_h);
U_f			= ((C_f-hc_f*C_f)^(1-sigmaC_f))/(1-sigmaC_f) - chi_f*(H^(1+sigmaL_f))/(1+sigmaL_f);
WELF_h 		=  U_h/(1-beta);
WELF_f 		=  U_f/(1-beta);
S 			= (1+R_K)/(1+R_L);
OA_h 		= C_h + gy*Y - W*H - (1-MC)*Y - R_D*D_h;
OA_f 		= C_f + gy*Y - W*H - (1-MC)*Y - R_D*D_f;

% Macroprudential 
phi_l = 0;
phi_ly = 0;
mp_cg_h = 0;
mp_cg_f = 0;
mp_cgy_h = 0;
mp_cgy_f = 0;
HA_h = (H^(1+sigmaL_h))*chi_h/((1+sigmaL_h)*(1-beta));
HA_f = (H^(1+sigmaL_f))*chi_f/((1+sigmaL_f)*(1-beta));


steady_state_model;
  k_h = K;	k_f = K;
  c_h = C_h; c_f = C_f;
  y_h = Y;	y_f = Y;
  h_h = H;	h_f = H;
  w_h = W;	w_f = W;
  z_h = Z;	z_f = Z;
  i_h = I;	i_f = I;
  lb_c_h = 1; lb_c_f = 1;
  mc_P_h = MC; mc_P_f = MC;
  delta_P_f = 1; pstar_f = 1;
  delta_P_h = 1; pstar_h = 1;
  pi_h = 1;	pi_f = 1;
  pic_h = 1;pic_f = 1;
  pX_h = 1;pX_f = 1;
  q_h = 1; q_f = 1;
  ac_I_h = 0; ac_I_f = 0;
  ToT = 1; pXI_h = 1; pXI_f = 1;
  r_IB = R;
  eta_h=eta_s; eta_f=eta_s;
  r_L_h=R_L; r_L_f=R_L;
  mc_L_h = MCB; mc_L_f = MCB;
  r_L_star_h = R_L; r_L_star_f = R_L;
	mp_h = 1; mp_f = 1;
  delta_L_h = 1; delta_L_f = 1;
  r_K_h=R_K; r_K_f=R_K;
  n_h=N_E; n_f=N_E;
  v_h=V_E; v_f=V_E;
  k_n_h=1/eta_E; k_n_f=1/eta_E;
  ls_h=Ls_h; ls_f=Ls_f;
  ld_h=Ld; ld_f=Ld;
  s_h=(1+R_K)/RR_L; s_f=(1+R_K)/RR_L;
  ss_h=(1+R_K)/RR_L; ss_f=(1+R_K)/RR_L;
  omega_h=OMEGA; omega_f=OMEGA;
  bk_h = BK_h; bk_f = BK_f;
  r_D_h = R_D; r_D_f = R_D;
  r_D_star_h = R_D; r_D_star_f = R_D; 
  d_h = D_h; d_f = D_f;
  dd_h = D_h; dd_f = D_f;
  mu_D_h = mu_d; mu_D_f = mu_d;
  delta_D_h = 1; delta_D_f = 1;
  cg_h = 1; cg_f = 1;
  e_A_h = 0; e_G_h = 0; e_I_h = 0; e_U_h = 0; e_N_h = 0; e_Q_h = 0; e_D_h = 0;
  e_A_f = 0; e_G_f = 0; e_I_f = 0; e_U_f = 0; e_N_f = 0; e_Q_f = 0; e_D_f = 0;
  ut_h = U_h; ut_f = U_f;
  welf_h = WELF_h; welf_f = WELF_f;
  d1_h = D1_h; d1_f = D1_f;
  d2_h = D2_h; d2_f = D2_f;
  s2_f = S2_f; s1_f = S1_f;
  s2_h = S2_h; s1_h = S1_h;
  f2_f = F2_f; f1_f = F1_f;
  f2_h = F2_h; f1_h = F1_h;
  ha_h = HA_h; ha_f = HA_f;
  diffLs_h = 1;diffLs_f = 1;
  diffY_h = 1;diffY_f = 1;
end;

%----------------------------------------------------------------
% 3. Model (the number refers to the equation in the paper)
%----------------------------------------------------------------
model;

    %% Household
	% Euler
	((c_h(+1)-hc_h*c_h)/(c_h-hc_h*c_h(-1)))^sigmaC_h=beta*((1+r_D_h)/pic_h(+1))/(1+pX_h*chi_D_h/100*(d_h-D_h)/D_h)*(exp(e_U_h(+1))/exp(e_U_h));
	((c_f(+1)-hc_f*c_f)/(c_f-hc_f*c_f(-1)))^sigmaC_f=beta*((1+r_D_f)/pic_f(+1))/(1+pX_f*chi_D_f/100*(d_f-D_f)/D_f)*(exp(e_U_f(+1))/exp(e_U_f));	
	% hours supply
	chi_h*(h_h^sigmaL_h)*((c_h-hc_h*c_h(-1))^sigmaC_h)=w_h;
	chi_f*(h_f^sigmaL_f)*((c_f-hc_f*c_f(-1))^sigmaC_f)=w_f;
	% variations in marginal utilty of consumption
	lb_c_h = (((c_h(-1)-hc_h*c_h(-2))^sigmaC_h))/(((c_h-hc_h*c_h(-1))^sigmaC_h));
	lb_c_f = (((c_f(-1)-hc_f*c_f(-2))^sigmaC_f))/(((c_f-hc_f*c_f(-1))^sigmaC_f));
	% Deposits law of motion
	d_h*pic_h(+1) - d_h(-1) = STEADY_STATE(r_D_h)*STEADY_STATE(d_h) + w_h*h_h + (pX_h - mc_P_h)*y_h - c_h - pX_h*gy*Y*exp(e_G_h) + OA_h - pX_h*(chi_D_h/200)*((d_h-D_h)^2)/D_h;
	d_f*pic_f(+1) - d_f(-1) = STEADY_STATE(r_D_f)*STEADY_STATE(d_f) + w_f*h_f + (pX_f - mc_P_f)*y_f - c_f - pX_f*gy*Y*exp(e_G_f) + OA_f - pX_f*(chi_D_f/200)*((d_f-D_f)^2)/D_f;
	% Utility
	ut_h = exp(e_U_h)*(((c_h-hc_h*c_h(-1))^(1-sigmaC_h))/(1-sigmaC_h) - chi_h*(h_h^(1+sigmaL_h)/(1+sigmaL_h))) - lambda_i*(r_IB-R)^2;
	ut_f = exp(e_U_f)*(((c_f-hc_f*c_f(-1))^(1-sigmaC_f))/(1-sigmaC_f) - chi_f*(h_f^(1+sigmaL_f)/(1+sigmaL_f))) - lambda_i*(r_IB-R)^2;
	% Expected utility
	welf_h = ut_h + beta*welf_h(+1);
	welf_f = ut_f + beta*welf_f(+1);
	% Welfare loss auxiliary variable to convert welfare into unconditional consumption
	ha_h = exp(e_U_h)*chi_h*(h_h^(1+sigmaL_h)/(1+sigmaL_h)) - lambda_i*(r_IB-R)^2 + beta*ha_h(+1);
	ha_f = exp(e_U_f)*chi_f*(h_f^(1+sigmaL_f)/(1+sigmaL_f)) - lambda_i*(r_IB-R)^2 + beta*ha_f(+1);
	

    %% FIRMS
	% Cost minimization condition
	w_h*h_h/(1-alpha) = z_h*k_h(-1)/alpha;
	w_f*h_f/(1-alpha) = z_f*k_f(-1)/alpha;
	% Production function
	delta_P_h*y_h = exp(e_A_h)*(k_h(-1)^alpha)*(h_h^(1-alpha));
	delta_P_f*y_f = exp(e_A_f)*(k_f(-1)^alpha)*(h_f^(1-alpha));
	% Marginal cost
	mc_P_h = (1/exp(e_A_h))*(z_h/alpha)^alpha*(w_h/(1-alpha))^(1-alpha);
	mc_P_f = (1/exp(e_A_f))*(z_f/alpha)^alpha*(w_f/(1-alpha))^(1-alpha);
	% NKPC forward and backward
	s1_h = y_h*pstar_h^-epsilon_P + beta*theta_P_h*lb_c_h(+1)*((pi_h^xi_P_h/pi_h(+1))^(1-epsilon_P))*((pstar_h/pstar_h(+1))^(-epsilon_P))*s1_h(+1);
	s1_f = y_f*pstar_f^-epsilon_P + beta*theta_P_f*lb_c_f(+1)*((pi_f^xi_P_f/pi_f(+1))^(1-epsilon_P))*((pstar_f/pstar_f(+1))^(-epsilon_P))*s1_f(+1);
	s2_h = y_h*mc_P_h*pstar_h^(-epsilon_P-1) + beta*theta_P_h*lb_c_h(+1)*((pi_h^xi_P_h/pi_h(+1))^(-epsilon_P))*((pstar_h/pstar_h(+1))^(-epsilon_P-1))*s2_h(+1);
	s2_f = y_f*mc_P_f*pstar_f^(-epsilon_P-1) + beta*theta_P_f*lb_c_f(+1)*((pi_f^xi_P_f/pi_f(+1))^(-epsilon_P))*((pstar_f/pstar_f(+1))^(-epsilon_P-1))*s2_f(+1);
	s1_h = (epsilon_P/(epsilon_P-1))*s2_h;	
	s1_f = (epsilon_P/(epsilon_P-1))*s2_f;	
	% price dispersion law of motion
	delta_P_h = (1-theta_P_h)*pstar_h^-epsilon_P + theta_P_h*(pi_h/(pi_h(-1)^xi_P_h))^epsilon_P*delta_P_h(-1);
	delta_P_f = (1-theta_P_f)*pstar_f^-epsilon_P + theta_P_f*(pi_f/(pi_f(-1)^xi_P_f))^epsilon_P*delta_P_f(-1);
	% pstar - optimal price setting by firms with price signal
	1 = (1-theta_P_h)*pstar_h^(1-epsilon_P) + theta_P_h*(pi_h(-1)^xi_P_h/pi_h)^(1-epsilon_P);
	1 = (1-theta_P_f)*pstar_f^(1-epsilon_P) + theta_P_f*(pi_f(-1)^xi_P_f/pi_f)^(1-epsilon_P);
	% Resources constraints
	y_h = (1-alpha_C_h)*((pX_h)^-mu)*c_h + ((1-n)/n)*alpha_C_f*((pX_f/ToT)^-mu)*c_f 
	    + (1-alpha_I_h)*((pXI_h)^-mu)*i_h*(1+ac_I_h) + ((1-n)/n)*alpha_I_f*((pXI_f/ToT)^-mu)*i_f*(1+ac_I_f)
		+ gy*Y*exp(e_G_h) + chi_D_h/200*((d_h-D_h)^2)/D_h;
	y_f = (1-alpha_C_f)*((pX_f)^-mu)*c_f + (n/(1-n))*alpha_C_h*((pX_h*ToT)^-mu)*c_h 
	    + (1-alpha_I_f)*((pXI_f)^-mu)*i_f*(1+ac_I_f) + (n/(1-n))*alpha_I_h*((pXI_h*ToT)^-mu)*i_h*(1+ac_I_h)
		+ gy*Y*exp(e_G_f) + chi_D_f/200*((d_f-D_f)^2)/D_f;
	
	
   %%% CAPITAL PRODUCERS
	% Ex post return on capital
	(1+r_K_h)/pic_h = (z_h+(1-delta)*q_h)/q_h(-1);
	(1+r_K_f)/pic_f = (z_f+(1-delta)*q_f)/q_f(-1);
	% Optimal capital demand
	q_h = pXI_h/pX_h + q_h*(chi_I_h/2)*(1+( 3*exp(e_I_h)*i_h/i_h(-1)-4)*exp(e_I_h)*i_h/i_h(-1)) + beta*q_h(+1)*lb_c_h(+1)*chi_I_h*(1-exp(e_I_h(+1))*i_h(+1)/i_h)*exp(e_I_h(+1))*(i_h(+1)/i_h)^2;
	q_f = pXI_f/pX_f + q_f*(chi_I_f/2)*(1+( 3*exp(e_I_f)*i_f/i_f(-1)-4)*exp(e_I_f)*i_f/i_f(-1)) + beta*q_f(+1)*lb_c_f(+1)*chi_I_f*(1-exp(e_I_f(+1))*i_f(+1)/i_f)*exp(e_I_f(+1))*(i_f(+1)/i_f)^2;
	% Capital Law of motion
	(1-ac_I_h)*i_h = k_h-(1-delta)*k_h(-1);
	(1-ac_I_f)*i_f = k_f-(1-delta)*k_f(-1);
	% Investment adjustment cost
	ac_I_h = (chi_I_h/2)*((exp(e_I_h)*i_h/i_h(-1)-1)^2);
	ac_I_f = (chi_I_f/2)*((exp(e_I_f)*i_f/i_f(-1)-1)^2);

	
   %%% ENTREPRENEUR
	% Credit demand
	ld_h-h_L_h*(ld_h(-1)-STEADY_STATE(ld_h)) + n_h = q_h*k_h;
	ld_f-h_L_f*(ld_f(-1)-STEADY_STATE(ld_f)) + n_f = q_f*k_f;
	% Entrepreneur net wealth accumulation
	n_h = (1-tau_E)*v_h/exp(e_N_h);
	n_f = (1-tau_E)*v_f/exp(e_N_f);
	% one period financial gains
	v_h = ((wmin^ka)/(ka-1))*(( ((k_n_h(-1)-1)/(k_n_h(-1)*ss_h(-1)))^(1-ka))*((1+r_K_h)/pic_h)*q_h(-1)*k_h(-1));
	v_f = ((wmin^ka)/(ka-1))*(( ((k_n_f(-1)-1)/(k_n_f(-1)*ss_f(-1)))^(1-ka))*((1+r_K_f)/pic_f)*q_f(-1)*k_f(-1));
	% External finance premium
	s_h = (1+r_K_h(+1))/(1+r_L_h);
	s_f = (1+r_K_f(+1))/(1+r_L_f);
	% FOC Cost-Of-Funds
	ss_h = ((ka/(ka-1))^varkappa_h*w_sup^-1)*(1-n_h/(q_h*k_h))^varkappa_h;
	ss_f = ((ka/(ka-1))^varkappa_f*w_sup^-1)*(1-n_f/(q_f*k_f))^varkappa_f;
	% Adding the unanticipated shock
	s_h	= exp(e_Q_h)*ss_h;
	s_f	= exp(e_Q_f)*ss_f;
	% ex post threshold
	omega_h*ss_h*q_h(-1)*k_h(-1)=ld_h(-1)-h_L_h*(ld_h(-2)-STEADY_STATE(ld_h));
	omega_f*ss_f*q_f(-1)*k_f(-1)=ld_f(-1)-h_L_f*(ld_f(-2)-STEADY_STATE(ld_f));
	% capital to net wealth ratio
	k_n_h = q_h*k_h/n_h;
	k_n_f = q_f*k_f/n_f;
	% survival rate of entrepreneurs
	eta_h = (wmin/omega_h)^ka;
	eta_f = (wmin/omega_f)^ka;

	

	%%% BANKING SYSTEM
	% Net Wealth
	bk_h = (1-tau_B)*((1-mu_b*(1-eta_h))*(1+r_L_h(-1))*ls_h(-1)-(1+r_IB(-1))*(ls_h(-1) - dd_h(-1) - bk_h(-1))-(1+r_D_h(-1))*dd_h(-1));
	bk_f = (1-tau_B)*((1-mu_b*(1-eta_f))*(1+r_L_f(-1))*ls_f(-1)-(1+r_IB(-1))*(ls_f(-1) - dd_f(-1) - bk_f(-1))-(1+r_D_f(-1))*dd_f(-1));
	%% LENDING activities
	% credit market equilibrium	
	ls_h = delta_L_h*ld_h;
	ls_f = delta_L_f*ld_f;
	% marginal cost of loan
	1+mc_L_h =  (1+r_IB)*mp_h/(1-mu_b*(1-eta_h(+1)));
	1+mc_L_f =  (1+r_IB)*mp_f/(1-mu_b*(1-eta_f(+1)));
	% New Keynesian Philips Curve for credit rates
	d1_h = ld_h*((r_L_h/r_L_star_h)^epsilon_L)*r_L_star_h + beta*theta_L_h*lb_c_h(+1)*((r_L_star_h/r_L_star_h(+1))^(1-epsilon_L))*d1_h(+1);
	d1_f = ld_f*((r_L_f/r_L_star_f)^epsilon_L)*r_L_star_f + beta*theta_L_f*lb_c_f(+1)*((r_L_star_f/r_L_star_f(+1))^(1-epsilon_L))*d1_f(+1);
	d2_h = ld_h*mc_L_h*(r_L_h/r_L_star_h)^epsilon_L + beta*theta_L_h*lb_c_h(+1)*((r_L_star_h/r_L_star_h(+1))^(-epsilon_L))*d2_h(+1);
	d2_f = ld_f*mc_L_f*(r_L_f/r_L_star_f)^epsilon_L + beta*theta_L_f*lb_c_f(+1)*((r_L_star_f/r_L_star_f(+1))^(-epsilon_L))*d2_f(+1);
	d1_h = (epsilon_L/(epsilon_L-1))*d2_h;
	d1_f = (epsilon_L/(epsilon_L-1))*d2_f;
	% Law of motion of credit rates dispersion
	delta_L_h = (1-theta_L_h)*(r_L_star_h/r_L_h)^-epsilon_L + theta_L_h*((r_L_h(-1)/r_L_h)^-epsilon_L)*delta_L_h(-1);
	delta_L_f = (1-theta_L_f)*(r_L_star_f/r_L_f)^-epsilon_L + theta_L_f*((r_L_f(-1)/r_L_f)^-epsilon_L)*delta_L_f(-1);
	% Aggregate prices
	r_L_h^(1-epsilon_L) = (1-theta_L_h)*r_L_star_h^(1-epsilon_L) + theta_L_h*(r_L_h(-1)^(1-epsilon_L));
	r_L_f^(1-epsilon_L) = (1-theta_L_f)*r_L_star_f^(1-epsilon_L) + theta_L_f*(r_L_f(-1)^(1-epsilon_L));
	
	%% DEPOSITS activities
	% New Keynesian Philips Curve for deposits rates 
	f1_h = d_h*((r_D_h/r_D_star_h)^(mu_D_h/(mu_D_h-1)))*r_D_star_h + beta*theta_D_h*lb_c_h(+1)*((r_D_star_h/r_D_star_h(+1))^(1/(1-mu_D_h)))*f1_h(+1);
	f1_f = d_f*((r_D_f/r_D_star_f)^(mu_D_f/(mu_D_f-1)))*r_D_star_f + beta*theta_D_f*lb_c_f(+1)*((r_D_star_f/r_D_star_f(+1))^(1/(1-mu_D_f)))*f1_f(+1);
	f2_h = d_h*r_IB*(r_D_h/r_D_star_h)^(mu_D_h/(mu_D_h-1)) + beta*theta_D_h*lb_c_h(+1)*((r_D_star_h/r_D_star_h(+1))^(-mu_D_h/(mu_D_h-1)))*f2_h(+1);	
	f2_f = d_f*r_IB*(r_D_f/r_D_star_f)^(mu_D_f/(mu_D_f-1)) + beta*theta_D_f*lb_c_f(+1)*((r_D_star_f/r_D_star_f(+1))^(-mu_D_f/(mu_D_f-1)))*f2_f(+1);	
	f1_h = mu_D_h*f2_h;
	f1_f = mu_D_f*f2_f;
	% Law of motion of deposit rate dispersion
	delta_D_h = (1-theta_D_h)*(r_D_star_h/r_D_h)^(-mu_D_h/(mu_D_h-1)) + theta_D_h*((r_D_h(-1)/r_D_h)^(-mu_D_h/(mu_D_h-1)))*delta_D_h(-1);
	delta_D_f = (1-theta_D_f)*(r_D_star_f/r_D_f)^(-mu_D_f/(mu_D_f-1)) + theta_D_f*((r_D_f(-1)/r_D_f)^(-mu_D_f/(mu_D_f-1)))*delta_D_f(-1);
	% Aggregate markup rates
	1 = (1-theta_D_h)*(r_D_star_h/r_D_h)^(1/(1-mu_D_h)) + theta_D_h*(r_D_h(-1)/r_D_h)^(1/(1-mu_D_h));
	1 = (1-theta_D_f)*(r_D_star_f/r_D_f)^(1/(1-mu_D_h)) + theta_D_f*(r_D_f(-1)/r_D_f)^(1/(1-mu_D_h));
	% markup
	mu_D_h = epsilon_D/(epsilon_D-1) + 100*e_D_h;
	mu_D_f = epsilon_D/(epsilon_D-1) + 100*e_D_f;
	% deposit market equilibrium	
	dd_h = delta_D_h*d_h;
	dd_f = delta_D_f*d_f;


	
	% International macroeconomics
	% Terms of trades
	ToT = (pi_f/pi_h)*ToT(-1);
	% Inflation Price
	pic_h^(1-mu) = (1-alpha_C_h)*((pX_h(-1)*pi_h)^(1-mu)) + alpha_C_h*((pX_h(-1)*pi_f*ToT(-1))^(1-mu));
	pic_f^(1-mu) = (1-alpha_C_f)*((pX_f(-1)*pi_f)^(1-mu)) + alpha_C_f*((pX_f(-1)*pi_h/ToT(-1))^(1-mu));
	% Relative Price, intermediary goods
	pX_h^(mu-1) = (1-alpha_C_h) + alpha_C_h*(ToT)^(1-mu); 
	pX_f^(mu-1) = (1-alpha_C_f) + alpha_C_f*(1/ToT)^(1-mu);
	% Relative Price, investment goods
	pXI_h^(mu-1) = (1-alpha_I_h) + alpha_I_h*(ToT)^(1-mu); 
	pXI_f^(mu-1) = (1-alpha_I_f) + alpha_I_f*(1/ToT)^(1-mu);
	
	
	%%% Monetary policy
	% Pricing rule
	(1+r_IB)/(1+R) = 	(((1+r_IB(-1))/(1+R))^rho) * ((pic_h^n*pic_f^(1-n))^phi_r * ((y_h/y_h(-1))^n*(y_f/y_f(-1))^(1-n))^phi_y )^(1-rho) * exp(e_R)
						;%* (((cg_h^n)*(cg_f^(1-n)))^phi_l);
	
	% MACROPRUDENTIAL POLICY
	% uncomment one of the scenarii
	% scenario (1)
	mp_h = (((cg_h^n)*(cg_f^(1-n)))^phi_l);
	mp_f = (((cg_h^n)*(cg_f^(1-n)))^phi_l);
	% scenario (2)
	%mp_h = (((cg_h^n)*(cg_f^(1-n)))^mp_cg_h);
	%mp_f = (((cg_h^n)*(cg_f^(1-n)))^mp_cg_f);
	% scenario (3)
	%mp_h = (cg_h)^phi_l;
	%mp_f = (cg_f)^phi_l;
	% scenario (4)
	%mp_h = (cg_h^mp_cg_h);
	%mp_f = (cg_f^mp_cg_f);
	
    % Exogenous shocks
	% Home shocks
	e_A_h = rho_A_h*e_A_h(-1) + eta_A_h;
	e_G_h = rho_G_h*e_G_h(-1) + eta_G_h + rho_AG_h*eta_A_h;
	e_I_h = rho_I_h*e_I_h(-1) + eta_I_h;
	e_U_h = rho_U_h*e_U_h(-1) + eta_U_h;
	e_N_h = rho_N_h*e_N_h(-1) + eta_N_h;
	e_Q_h = rho_Q_h*e_Q_h(-1) + eta_Q_h;
	e_D_h = rho_D_h*e_D_h(-1) + eta_D_h;
	% Foreign shocks
	e_A_f = rho_A_f*e_A_f(-1) + eta_A_f;
	e_G_f = rho_G_f*e_G_f(-1) + eta_G_f + rho_AG_f*eta_A_f;
	e_I_f = rho_I_f*e_I_f(-1) + eta_I_f;
	e_U_f = rho_U_f*e_U_f(-1) + eta_U_f;
	e_N_f = rho_N_f*e_N_f(-1) + eta_N_f;
	e_Q_f = rho_Q_f*e_Q_f(-1) + eta_Q_f;
	e_D_f = rho_D_f*e_D_f(-1) + eta_D_f;
	% Common MP Shock
	e_R   = rho_R*e_R(-1) + eta_R;
	
	% growth variables
	cg_h = ls_h/ls_h(-1);
	cg_f = ls_f/ls_f(-1);
	diffY_h = y_h/y_h(-1);
	diffY_f = y_f/y_f(-1);
	diffLs_h = ls_h/ls_h(-1);
	diffLs_f = ls_f/ls_f(-1);
end;

resid(1);


%% Shock Variance-covariance matrix
%% small difference with the paper: taken at posterior mode here
shocks;
	var eta_A_h;  stderr 0.008645210580004;
	var eta_A_f;  stderr 0.012056718729987;
	corr eta_A_h, eta_A_f = 0.406000845748930;
	var eta_G_h;  stderr 0.012985597964117;
	var eta_G_f;  stderr 0.020243714020978;
	corr eta_G_h, eta_G_f = 0.057350987693517;
	var eta_U_h;  stderr 0.005887797467257;
	var eta_U_f;  stderr 0.007078713105775;
	corr eta_U_h, eta_U_f = -0.204711867818667;
	var eta_I_h;  stderr 0.019117783203897;
	var eta_I_f;  stderr 0.019423322108333;
	corr eta_I_h, eta_I_f = 0.305572351250239;
	var eta_N_h;  stderr 0.001430274684715;
	var eta_N_f;  stderr 0.010177141273495;
	corr eta_N_h, eta_N_f = 0.068370932869852;
	var eta_Q_h;  stderr 0.003660907016634;
	var eta_Q_f;  stderr 0.003945477716704;
	corr eta_Q_h, eta_Q_f = 0.002884026738987;
	var eta_D_h;  stderr 0.000422968734815;
	var eta_D_f;  stderr 0.000500165901502;
	corr eta_D_h, eta_D_f = 0.959742054036376;
	var eta_R;		stderr 0.000908680830648;
end;







if 1 == 1 % deactivate or activate the search code for optimal policies
%% OPTIMAL MONETARY POLICIES
step_size = 0.01;
% change here the optimized variable and the step size of the search grid
% Variable 1
varname1	= 'phi_l';  
interval1	= 0 : step_size : .2;
% Variable 2
varname2	= 'phi_y';  
interval2	= 0 : step_size : 0;
% Variable 3
varname3	= 'mp_cg_h';
interval3	= 0.0 : step_size : 0.0;
% Variable 4
varname4	= 'mp_cg_f';
interval4	= 0 : step_size : 0;

mytime = clock; 
disp(['I start calculating the welfare at ' datestr(datenum(mytime(1),mytime(2),mytime(3),mytime(4),mytime(5),mytime(6)))])
disp(['If I only could be running up that hill.'])

% country weight in the objective function
nn = n;		% Union average
% nn = 1;	% core only in objective
% nn = 0;	% periphery only in objective


% Initialization
val_uet=zeros(length(interval1),length(interval2),length(interval3),length(interval4));
val_uet(:,:,:,:)=NaN;
counter = 0;
countertotal = length(interval1)*length(interval2)*length(interval3)*length(interval4);
for i1 = 1:length(interval1)
	% update param 1
	set_param_value(varname1,interval1(i1));
	
	for i2 = 1:length(interval2)
		% update param 2
		set_param_value(varname2,interval2(i2));
		
		for i3 = 1:length(interval3)
			% update param 3
			set_param_value(varname3,interval3(i3));
			
			for i4 = 1:length(interval4)
				% update param 4
				set_param_value(varname4,interval4(i4));

				
				% Computing the welfare
				stoch_simul(order=2,irf=0,noprint,nodisplay,nocorr,nofunctions,nograph)  welf_h welf_f;
				val_uet(i1,i2,i3,i4) = (nn*oo_.mean(1)+(1-nn)*oo_.mean(2)) - (nn*WELF_h+(1-nn)*WELF_f) ;
				
				% updating counter
				counter = 1 + counter; 
				% getting a variation point
				if counter == 1
					startpoint = val_uet(i1,i2,i3,i4);
				end
				
				% showing in the screen the result
				mytime = clock;
				mygain = 100*((val_uet(i1,i2,i3,i4)-startpoint)/abs(startpoint));
				% Visual tips
				if mygain > 0
					mygain = ['+++' num2str(mygain)];
					jump4 = 0;
				elseif mygain < 0
					mygain = ['---' num2str(mygain)];
				else
					mygain = num2str(mygain);
				end
				
				disp([ num2str(counter) '/' num2str(countertotal) ': ' varname1 '=' num2str(interval1(i1)) ', ' varname2 '=' num2str(interval2(i2)) ', ' varname3 '=' num2str(interval3(i3)) ', ' varname4 '=' num2str(interval4(i4)) '-> welf=' num2str(val_uet(i1,i2,i3,i4)) ' (' mygain ') at ' datestr(datenum(mytime(1),mytime(2),mytime(3),mytime(4),mytime(5),mytime(6)))]);
				
				if counter > 1
					startpoint = max(val_uet(:));  %val_uet(i1,i2,i3,i4);
				end
		
			end % closing param 4 exploration
			
		end % closing param 3 exploration
			
	end % closing param 2 exploration
		
end % closing param 1 exploration
mytime = clock; 
disp(['Grid-search finished at ' datestr(datenum(mytime(1),mytime(2),mytime(3),mytime(4),mytime(5),mytime(6)))])

% Find max
[max_val, position] = max(val_uet(:)); 
[param1_id,param2_id,param3_id,param4_id] = ind2sub(size(val_uet),position);
disp([ 'local maximum : welfare of ' num2str(val_uet(param1_id,param2_id,param3_id,param4_id)) ' for ' varname1 '=' num2str(interval1(param1_id)) ' and ' varname2 '=' num2str(interval2(param2_id)) ' and ' varname3 '=' num2str(interval3(param3_id)) ' and ' varname4 '=' num2str(interval4(param4_id)) ]);
max1 = interval1(param1_id);
max2 = interval2(param2_id);
max3 = interval3(param3_id);
max4 = interval4(param4_id);

end % closing the search code



if 1 == 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% COMPUTING UTILITY LOSS
%% 1ST CASE : OPTIMAL RULE
% Optimal Monetary Policy
set_param_value('phi_r',3);	
set_param_value('phi_y',0);
% Estimated Monetary Policy
%set_param_value('phi_r',2.3729);
%set_param_value('phi_y',0.1611);
stoch_simul(order=2,irf=0,noprint,nograph) welf_h welf_f;
wwelf_h_a=oo_.mean(1);
wwelf_f_a=oo_.mean(2);
%% ANOTHER POLICY : 2ND RULE
set_param_value(varname1,max1);
set_param_value(varname2,max2);
set_param_value(varname3,max3);
set_param_value(varname4,max4);
stoch_simul(order=2,irf=0,noprint,nodisplay,nograph) welf_h welf_f ha_h ha_f;
wwelf_h_b=oo_.mean(1);
wwelf_f_b=oo_.mean(2);
hha_h_b=oo_.mean(3);
hha_f_b=oo_.mean(4);
new_welf = n*wwelf_h_b + (1-n)*wwelf_f_b;
% reinverting the vacov matrix of shocks to get good theoretical moments
for i1 = 1:size(M_.Sigma_e,1)
	for i2 = 1:size(M_.Sigma_e,1)
		M_.Sigma_e(i1,i2) = (sqrt(M_.Sigma_e(i1,i2))*100)^2;
	end
end
stoch_simul(order=1,irf=0,noprint,nodisplay,nograph,periods=1000) welf_h welf_f ha_h ha_f pi_h pi_f diffY_h diffY_f diffLs_h diffLs_f y_h y_f;
var_mat = diag(oo_.var);
% Euro Area standard deviations
varpi = (n*var_mat(5) + (1-n)*var_mat(6))^.5;
vary = (n*var_mat(7) + (1-n)*var_mat(8))^.5;
varls = (n*var_mat(9) + (1-n)*var_mat(10))^.5;
dycorrel = oo_.var(7,8)/(oo_.var(7,7)*oo_.var(8,8))^.5;
ycorrel = oo_.var(11,12)/(oo_.var(11,11)*oo_.var(12,12))^.5;

%%%% WELFARE LOSS
options = optimset('Display','off');
phi = fsolve(@(phi)n*wwelf_h_a+(1-n)*wwelf_f_a-n*((1-phi)^(1-sigmaC_h))*(wwelf_h_b+hha_h_b)-(1-n)*((1-phi)^(1-sigmaC_f))*(wwelf_f_b+hha_f_b)+n*hha_h_b+(1-n)*hha_f_b,.02,options);
phi2 = fsolve(@(phi2)wwelf_h_a-((1-phi2)^(1-sigmaC_h))*(wwelf_h_b+hha_h_b)+hha_h_b,.02,options);
phi3 = fsolve(@(phi3)wwelf_f_a-((1-phi3)^(1-sigmaC_f))*(wwelf_f_b+hha_f_b)+hha_f_b,.02,options);

disp(['		' varname1 '		' varname2 '		' varname3 '		' varname4 '		Welf		%Gain		VarPi		VarY		VarL		correl(Y,Y)'])
disp(['		' num2str(max1) '			' num2str(max2) '			' num2str(max3) '			' num2str(max4) '		' num2str(new_welf) '	' num2str(phi*100) 'o/o		' num2str(varpi) '		' num2str(vary) '		' num2str(varls) '		' num2str(dycorrel) '		' num2str(ycorrel) '		'])
disp('----------------------------------------------------')
disp('		W_u,t			W_h,t			 W_f,t')
disp(['		' num2str(phi*100) '		' num2str(phi2*100) '		' num2str(phi3*100) '		'])
disp('----------------------------------------------------')
disp('National Statistics sd')
disp(['		sd(pi_c):' num2str(var_mat(5)^.5) '			sd(Dy_c):' num2str(var_mat(7)^.5) '			sd(Dl_c):' num2str(var_mat(9)^.5)])
disp(['		sd(pi_p):' num2str(var_mat(6)^.5) '			sd(Dy_p):' num2str(var_mat(8)^.5) '			sd(Dl_p):' num2str(var_mat(10)^.5)])

end