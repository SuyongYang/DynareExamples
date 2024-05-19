function results = HJB_SS(parameters,K,lambda)
%Computes the HJB in Steady State

%PARAMETERS
psi   = parameters.psi ;          % long-run growth rate
ga    =  parameters.ga;           % CRRA utility with parameter gamma
alpha = parameters.alpha;         % Production function F = K^alpha * L^(1-alpha) 
eta   = parameters.eta;           % death probability
delta = parameters.delta;         % Capital depreciation 10% - death rate
zmean = parameters.zmean;         % mean O-U process (in levels). This parameter has to be adjusted to ensure that the mean of z (truncated gaussian) is 1.
sig2  = parameters.sig2;          % sigma^2 O-U
Corr  = parameters.Corr;          % persistence   O-U
rho   = parameters.rho;           % discount rate + death prob
I     = parameters.I;             % number of a points 
J     = parameters.J;             % number of z points 
zmin  = parameters.zmin ;         % Range z
zmax  = parameters.zmax;
amin  = parameters.amin;          % borrowing constraint
amax  = parameters.amax;          % range a

%simulation parameters
maxit = parameters.maxit;         % maximum number of iterations in the HJB loop
crit  = parameters.crit;          % criterion HJB loop
Delta = parameters.Delta ;        % delta in HJB algorithm

%ORNSTEIN-UHLENBECK IN LEVELS
the = Corr;

%--------------------------------------------------
%VARIABLES 
a  = linspace(amin,amax,I)';  %wealth vector
da = (amax-amin)/(I-1);      
z  = linspace(zmin,zmax,J);   % productivity vector
dz = (zmax-zmin)/(J-1);
dz2 = dz^2;
aa = a*ones(1,J);
zz = ones(I,1)*z;

mu = the*(zmean - z);        %DRIFT (FROM ITO'S LEMMA)
s2 = sig2.*ones(1,J);        %VARIANCE (FROM ITO'S LEMMA)

%Finite difference approximation of the partial derivatives
Vaf = zeros(I,J);             
Vab = zeros(I,J);

%CONSTRUCT MATRIX Aswitch SUMMARIZING EVOLUTION OF z
yy = - s2/dz2 - mu/dz;
chi =  s2/(2*dz2);
zeta = mu/dz + s2/(2*dz2) ;

%This will be the upperdiagonal of the matrix Aswitch
updiag=zeros(I,1); %This is necessary because of the peculiar way spdiags is defined.
for j=1:J
    updiag=[updiag;repmat(zeta(j),I,1)];
end

%This will be the center diagonal of the matrix Aswitch
centdiag=repmat(chi(1)+yy(1),I,1);
for j=2:J-1
    centdiag=[centdiag;repmat(yy(j),I,1)];
end
centdiag=[centdiag;repmat(yy(J)+zeta(J),I,1)];

%This will be the lower diagonal of the matrix Aswitch
lowdiag=repmat(chi(2),I,1);
for j=3:J
    lowdiag=[lowdiag;repmat(chi(j),I,1)];
end

%Add up the upper, center, and lower diagonal into a sparse matrix
Aswitch=spdiags(centdiag,0,I*J,I*J)+spdiags(lowdiag,-I,I*J,I*J)+spdiags(updiag,I,I*J,I*J);

%----------------------------------------------------
%INITIAL GUESS
r    = alpha     * K^(alpha-1) -delta; %interest rates
w    = (1-alpha) * K^(alpha);          %wages
v0   = (w*zz + (r  -psi +eta).*aa).^(1-ga)/(1-ga)/(rho - (1-ga)*psi + eta);
v    = v0;
dist = zeros(1,maxit);

%-----------------------------------------------------
% HAMILTON-JACOBI-BELLMAN EQUATION %
    for n=1:maxit
        V = v;
        % forward difference
        Vaf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
        Vaf(I,:) = (w*z + (r -psi +eta).*amax).^(-ga); %will never be used, but impose state constraint a<=amax just in case
        % backward difference
        Vab(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
        Vab(1,:) = (w*z + (r -psi +eta).*amin).^(-ga);  %state constraint boundary condition

        %I_concave = Vab > Vaf;              %indicator whether value function is concave (problems arise if this is not the case)

        %consumption and savings with forward difference
        cf = Vaf.^(-1/ga);
        sf = w*zz + (r -psi +eta).*aa - cf;
        %consumption and savings with backward difference
        cb = Vab.^(-1/ga);
        sb = w*zz + (r -psi +eta).*aa - cb;
        %consumption and derivative of value function at steady state
        c0 = w*zz + (r -psi +eta).*aa;
        Va0 = c0.^(-ga);

        % dV_upwind makes a choice of forward or backward differences based on
        % the sign of the drift
        If = sf > 0; %positive drift --> forward difference
        Ib = sb < 0; %negative drift --> backward difference
        I0 = (1-If-Ib); %at steady state
        %make sure backward difference is used at amax
        %     Ib(I,:) = 1; If(I,:) = 0;
        %STATE CONSTRAINT at amin: USE BOUNDARY CONDITION UNLESS sf > 0:
        %already taken care of automatically

        Va_Upwind = Vaf.*If + Vab.*Ib + Va0.*I0; %important to include third term

        c = Va_Upwind.^(-1/ga);
        u = c.^(1-ga)/(1-ga)+ lambda*(aa-K);

        %CONSTRUCT MATRIX A
        X = - min(sb,0)/da;
        Y = - max(sf,0)/da + min(sb,0)/da -eta;  %Includes eta (random death)
        Z = max(sf,0)/da;
        
        updiag=0; %This is needed because of the peculiarity of spdiags.
        for j=1:J
            updiag=[updiag;Z(1:I-1,j);0];
        end
        
        centdiag=reshape(Y,I*J,1);
        
        lowdiag=X(2:I,1);
        for j=2:J
            lowdiag=[lowdiag;0;X(2:I,j)];
        end
        
        AA=spdiags(centdiag,0,I*J,I*J)+spdiags([updiag;0],1,I*J,I*J)+spdiags([lowdiag;0],-1,I*J,I*J);
        
        A = AA + Aswitch;
        B = (1/Delta + rho - (1-ga)*psi)*speye(I*J) - A;

        u_stacked = reshape(u,I*J,1);
        V_stacked = reshape(V,I*J,1);

        b = u_stacked + V_stacked/Delta;

        V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS

        V = reshape(V_stacked,I,J);

        Vchange = V - v;
        v = V;

        dist(n) = max(max(abs(Vchange)));
        if dist(n)<crit
%             disp('Value Function Converged, Iteration = ')
%             disp(n)
            break
        end
    end
    
    if n == maxit
        disp('Value Function NOT Converged, ERROR ')
    end
    
   % clear A AA AT B
  % RESULTS
results.a  = a;
results.aa = aa;
results.z  = z;
results.zz = zz;
results.V  = V;
results.A  = A;
results.c  = c;
results.r  = r;
results.w  = w;
results.da = da;
results.dz = dz;
results.s  = w*zz + (r- psi +eta).*aa - c;