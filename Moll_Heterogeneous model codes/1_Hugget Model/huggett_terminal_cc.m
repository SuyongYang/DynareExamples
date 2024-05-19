%Credit crunch extension by Gustavo Mellior

clear all; clc; close all;

tic;

s = 2;
rho = 0.05;
z1 = .12;
z2 = .25;
z = [z1,z2];
la1 = 1.15;
la2 = 1;
la = [la1,la2];

r0 = 0.03;
rmin = 0.001;
rmax = 0.045;

%This part is explained in huggett_initial_creditcrunch.m
I= 800;
amin = -0.15;
amax = 5;
a = linspace(amin,amax,I)';
da = (amax-amin)/(I-1);
num = 3; 

for i=num:-1:1
    first(num+1-i) = amin-i*da;
end
a = [first';a];
I = I+num;

aa = [a,a];
zz = ones(I,1)*z;

maxit= 100;
crit = 10^(-6);
Delta = 100;

dVf = zeros(I,2);
dVb = zeros(I,2);
c = zeros(I,2);

Aswitch = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)];

Ir = 40;
crit_S = 10^(-8);

%INITIAL GUESS
r = r0;
v0(:,1) = (z(1) + r.*a).^(1-s)/(1-s)/rho;
v0(:,2) = (z(2) + r.*a).^(1-s)/(1-s)/rho;

for ir=1:Ir

r_r(ir)=r;
rmin_r(ir)=rmin;
rmax_r(ir)=rmax;    

v = v0;

for n=1:maxit
    V = v;
    V_n(:,:,n)=V;
    % forward difference
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:) = (z + r.*amax).^(-s); %will never be used, but impose state constraint a<=amax just in case
    % backward difference
    dVb(num+2:I,:) = (V(num+2:I,:)-V(num+1:I-1,:))/da;
    dVb(num+1,:) = (z + r.*amin).^(-s); %state constraint boundary condition as before
    
    %Consumption, derivative of value function at steady state and new
    %constraints
    c0 = zz + r.*aa;
    for j=num:-1:1 %New state constraints
    dVb(num+1-j,:) = (z + r.*(amin-j*da)-j*da).^(-s); %moving state constraint boundary condition
    c0(num+1-j) = c0(num+1-j)-j*da;
    %In the terminal period, we force households
    %to not be in the inadmissible region. So we apply a "hard" saving rule
    %that leaves no mass to the left of the new debt limit in the terminal
    %state (notice that if they are in the inadmissable region they will have
    %to move j steps to right).
    end
    dV0 = c0.^(-s);
    
    I_concave = dVb > dVf; %indicator whether value function is concave (problems arise if this is not the case)
    
    %consumption and savings with forward difference
    cf = dVf.^(-1/s);
    ssf = zz + r.*aa - cf;
    %consumption and savings with backward difference
    cb = dVb.^(-1/s);
    ssb = zz + r.*aa - cb;
    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift    
    If = ssf > 0; %positive drift --> forward difference
    Ib = ssb < 0; %negative drift --> backward difference
    I0 = (1-If-Ib); %at steady state
    I0(I0<0)=0;
    %make sure backward difference is used at amax
    %Ib(I,:) = 1; If(I,:) = 0;
    %STATE CONSTRAINT at amin: USE BOUNDARY CONDITION UNLESS muf > 0:
    %already taken care of automatically
    
    dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0; %important to include third term
    c = dV_Upwind.^(-1/s);
    u = c.^(1-s)/(1-s);
    
    %CONSTRUCT MATRIX
    X = - min(ssb,0)/da;
    Y = - max(ssf,0)/da + min(ssb,0)/da;
    Z = max(ssf,0)/da;
    
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;

    B = (1/Delta + rho)*speye(2*I) - A;
    
    u_stacked = [u(:,1);u(:,2)];
    V_stacked = [V(:,1);V(:,2)];
    
    b = u_stacked + V_stacked/Delta;
    V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS
    
    V = [V_stacked(1:I),V_stacked(I+1:2*I)];
    
    Vchange = V - v;
    v = V;

    dist(n) = max(max(abs(Vchange)));
    if dist(n)<crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
end
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOKKER-PLANCK EQUATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%
AT = A';
b = zeros(2*I,1);

%need to fix one value, otherwise matrix is singular
i_fix = 1;
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,2*I-i_fix)];
AT(i_fix,:) = row;

%Solve linear system
gg = AT\b;
g_sum = gg'*ones(2*I,1)*da;
gg = gg./g_sum;

g = [gg(1:I),gg(I+1:2*I)];

check1 = g(:,1)'*ones(I,1)*da;
check2 = g(:,2)'*ones(I,1)*da;

g_r(:,:,ir) = g;
adot(:,:,ir) = zz + r.*aa - c;
V_r(:,:,ir) = V;

S(ir) = g(:,1)'*a*da + g(:,2)'*a*da;

%UPDATE INTEREST RATE
if S(ir)>crit_S
    disp('Excess Supply')
    rmax = r;
    r = 0.5*(r+rmin);
elseif S(ir)<-crit_S;
    disp('Excess Demand')
    rmin = r;
    r = 0.5*(r+rmax);
elseif abs(S(ir))<crit_S;
    display('Equilibrium Found, Interest rate =')
    disp(r)
    break
end

end

save huggett_terminal_creditcrunch.mat