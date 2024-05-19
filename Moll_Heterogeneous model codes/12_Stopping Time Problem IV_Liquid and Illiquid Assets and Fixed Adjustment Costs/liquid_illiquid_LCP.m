% Solves the household problem in the model with a liquid and an illiquid
% asset with Poisson shocks using the finite difference upwind method with
% implicit updating.
% written by Ben Moll, Peter Maxted and Alberto Polo
%--------------------------------------------------------------------------

clear all; clc; close all;

% tic
%% PARAMETERS AND INITIAL OBJECTS:

UseNoAdjustmentAsGuess = 1;
AssumeConcavity = 0; % this is relevant in the part with adjustment
BilinW = 1; % when all neighboring points of (b',a') are non-adjustment points, use bilinear weights in dividing mass to be rebalanced 
            % across them in the KF algorithm

% Parameters:
s = 2; % risk aversion in CRRA utility
ra    = 0.0386; % return on illiquid asset
rbPos = -0.01; % return on liquid asset
rbNeg = 0.15;
rho   = 0.037; % discount rate
ka = 0.1; % cost of adjustment
          % assuming the wage w is 1

% Illiquid asset 
J   = 60; % number of grid points %original: 30
amin   = 0;
amax   = 20; % orig 6
a  = linspace(amin,amax,J); % row vector
da = a(2)-a(1);

I    = 70; % number of grid points %original: 40
bmin    = 0;
bmax    = 3;
b = linspace(bmin,bmax,I)'; % column vector
db = b(2)-b(1);

% Poisson shocks:
Nz = 2; % number of exog. states
z      = [.8,1.2]; % state values (row vector)
la = [1/15,1/15]; % Poisson intensities

% Create mesh: b in 1st dimension, a in 2nd dimension, z in 3rd dimension
bb = b*ones(1,J); % repeat b column vector J times in the 2nd dimension
aa = ones(I,1)*a; % repeat a row vector I times in the 1st dimension
zz = ones(J,1)*z; % repeat z row vector J times in the 1st dimension (will then take the 2 columns with 1st and 2nd state values separately)

bbb = zeros(I,J,Nz); aaa = zeros(I,J,Nz); zzz = zeros(I,J,Nz); % these are the meshes
bbb(:,:,1) = bb; bbb(:,:,2) = bb;
aaa(:,:,1) = aa; aaa(:,:,2) = aa;
zzz(:,:,1) = z(1);zzz(:,:,2)= z(2);

%allow for differential borrowing and lending rates (only matters when bmin<0)
Rb = rbPos.*(bbb>=0) + rbNeg.*(bbb<0); % mesh of returns on liquid asset to allow for different borrowing rate if borrowing is allowed

L = I*J*Nz; % total number of nodes

% Adjustment decision set-up:
Nx      = 600; % number of grid points for cash-in-hand when adjusting %original: 400
NaP     = 600; % number of grid points for a' when adjusting
xmin      = amin + bmin + ka; % minimum cash-in-hand
xmax      = amax + bmax; % maximum
x     = linspace(xmin,xmax,Nx);

sumGrid   = aa + bb; % mesh of cash-in-hand from (b, a) nodes
sumGrid   = sumGrid(:); % transform to column vector


% Iteration set-up:
Delta     = 1000; % step size (for time) in the implicit method (as opposed to db and db which are steps in the state dimensions)
maxit     = 200;
tol       = 1e-6;
iter      = 0;
dist      = zeros(maxit,1); % hold errors during iteration


% Initial guess
phia      = 0.8;
phib      = 0.01;
V0(:,:,1) = (1-s)^(-1)*(z(1) + phia*ra.*aa + phib.*bb).^(1-s)/rho; % assuming u(c) = c^(1 - s)/(1 - s). This is similar to value of "staying put" u(z + rb*b)/rho
                                                                  
V0(:,:,2) = (1-s)^(-1)*(z(2) + phia*ra.*aa + phib.*bb).^(1-s)/rho;
v         = V0;  

tau = 15; % from two_assets_kinked.m, this is justified as: if ra>>rb, impose tax on ra*a at high a, otherwise some households accumulate infinite illiquid wealth
          % (not needed if ra is close to - or less than - rb). It implies aDrift = ra*a*(1 - (a/(amax * 0.999))^(tau - 1))
tau0 = ra.*(amax*.999)^(1-tau); % in two_assets_kinked.m the multiplier is 1.33 instead of 0.999. That multiplier makes the tax much less aggressive, as with 
                                % 0.999 the drift of a goes down to slightly negative towards amax, while with 1.33 the drift becomes just slightly concave near amax
T = tau0.*a.^tau;
aDrift = ra.*a - T;
%plot(a,aDrift,a,zeros(J,1),a,ra.*a)

%aDrift = ra.*a;

%%

%%%%% Build matrix Bswitch summarizing evolution of (a,z) % i.e. states which are exogenous in the no adjustment case: illiquid asset and (properly) exogenous state
chi = -min(aDrift,0)/da; % this is analogous to X, Y, Z constructed below for b (but call chi instead of x since x is cash-in-hand). This (and below) has 
                         % similarities with p. 16/17 in Achdou et al. 2017 - Online Appendix
yy  = - max(aDrift,0)/da + min(aDrift,0)/da;
zeta= max(aDrift,0)/da;

%This will be the upperdiagonal of the B_switch
updiag=zeros(I,1); %This is necessary because of the peculiar way spdiags is defined (since this will be the Ith upper diagonal, the first I elements are dropped)
for j=1:J
    updiag=[updiag;repmat(zeta(j),I,1)]; % this is analogous to zeta(j) * ones(I, 1) i.e. repeating zeta(j) I times in a vector
end

%This will be the center diagonal of the B_switch
centdiag=repmat(chi(1)+yy(1),I,1);
for j=2:J-1
    centdiag=[centdiag;repmat(yy(j),I,1)];
end
centdiag=[centdiag;repmat(yy(J)+zeta(J),I,1)];

%This will be the lower diagonal of the B_switch
lowdiag=repmat(chi(2),I,1);
for j=3:J
    lowdiag=[lowdiag;repmat(chi(j),I,1)]; % this will have length I*(J-1) because the last I elements in constructing the Ith lower diagonal will be dropped anyway
end

%Add up the upper, center, and lower diagonal into a sparse matrix
Baux=spdiags(centdiag,0,I*J,I*J)+spdiags(lowdiag,-I,I*J,I*J)+spdiags(updiag,I,I*J,I*J);

Bswitch = [Baux - speye(I*J)*la(1), speye(I*J)*la(1);speye(I*J)*la(2),Baux - speye(I*J)*la(2)]; % this adds the transitions between z1 and z2 (keep in mind they have
                                                                                                % to balance along the row and sum to 0)

%% SOLVE WITHOUT ADJUSTMENT
if UseNoAdjustmentAsGuess == 1

    for n=1:maxit

        V = v;
        % forward difference for the derivative of v wrt b
        Vbf(1:I-1,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/db; % remember b changes along the 1st dimension (column vector)
        Vbf(I,:,:) = (zzz(I,:,:) + Rb(I,:,:).*bmax).^(-s); %will never be used, but impose state constraint b<=bmax just in case
                                                           % this is v'b(I,j,k)_F = u'(z(k) + rb(I) * b(I))
        % backward difference
        Vbb(2:I,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/db;
        Vbb(1,:,:) = (zzz(1,:,:) + Rb(1,:,:).*bmin).^(-s); %state constraint boundary condition: v'b(1,j,k)_B = u'(z(k) + rb(1) * b(1))

        Vbf = max(Vbf,10^(-6)); % this is to avoid dealing with too small numbers in parts where v is almost flat wrt b
        Vbb = max(Vbb,10^(-6));

        %consumption and savings with forward difference
        cf = Vbf.^(-1/s); % from FOC
        sf = zzz + Rb.*bbb - cf;
        %consumption and savings with backward difference
        cb = Vbb.^(-1/s);
        sb = zzz + Rb.*bbb - cb;
        %consumption and derivative of value function at steady state
        c0 = zzz + Rb.*bbb; % this is where s=0 i.e. db = 0

        % make a choice of forward or backward differences based on the sign of the drift (upwind scheme)
        Ib = (sb < 0); %negative drift --> backward difference. This is where we assume concavity, because we don't check sf <= 0 as well (see *** reference)
        If = (sf > 0).*(1-Ib); %positive drift --> forward difference. Note that the (1-Ib) makes sure we take only points where also sb>=0, so Ib and If are 
                               % disjoint. Concavity of v wrt b would imply that, because with concavity sf <= sb, but this enforces it in a setting
                               % where there could be some numerical instability
        I0 = (1-If-Ib); %at steady state. These have to be points where sb>=0 & sf<=0. If we didn't impose disjointness in the previous step, this could = -1 where
                        % (sb < 0 & sf > 0) and we would need to split these points (see *** reference). Instead we assume concavity i.e. these points don't arise 

        c = cf.*If + cb.*Ib + c0.*I0;
        u = c.^(1-s)/(1-s);

        %CONSTRUCT MATRIX (see p. 5 Achdou et al. 2017 - Online Appendix)
        X = -Ib.*sb/db; % lower
        Y = -If.*sf/db + Ib.*sb/db; % center
        Z = If.*sf/db; % upper

        updiag = zeros(I*J,Nz); % notice that endogenous nodes (b, a) are listed in the 1st dimension, while exogenous state z is listed in the 2nd
        lowdiag = zeros(I*J,Nz);
        centdiag = zeros(I*J,Nz);

        for i = 1:Nz
            centdiag(:,i) = reshape(Y(:,:,i),I*J,1);
        end

        % (* reference)
        lowdiag(1:I-1,:) = X(2:I,1,:); % exclude the 1st row i.e. b(1). This is only at a(1) (same for updiag below). Notice lowdiag(I,:)=0
        updiag(2:I,:) = Z(1:I-1,1,:); % exclude the last row i.e. b(I). Notice upperdiag(1,:)=0
        for j = 2:J
            lowdiag(1:j*I,:) = [lowdiag(1:(j-1)*I,:);squeeze(X(2:I,j,:));zeros(1,Nz)]; % this stacks the columns of X in the right way (adding a 0 at the end, 
                                                                                       % with size [1, Nz] to account for exogenous states)
            updiag(1:j*I,:) = [updiag(1:(j-1)*I,:);zeros(1,Nz);squeeze(Z(1:I-1,j,:))];
        end


        AA1=spdiags(centdiag(:,1),0,I*J,I*J)+spdiags(updiag(:,1),1,I*J,I*J)+spdiags(lowdiag(:,1),-1,I*J,I*J); % notice that when filling the upper diagonal
                                                                                       % the first element of the vector is dropped. When filling the lower diagonal
                                                                                       % the last element is dropped. This is consistent with how lowdiag and updiag
                                                                                       % are constructed initially (see * reference). Note: removed unneeded ; 0]
        AA2=spdiags(centdiag(:,2),0,I*J,I*J)+spdiags(updiag(:,2),1,I*J,I*J)+spdiags(lowdiag(:,2),-1,I*J,I*J);
        AA = [AA1, sparse(I*J,I*J); sparse(I*J,I*J), AA2]; % this completes the I*J matrixes to bring them to I*J*Nz

        A = AA + Bswitch; % adds the transition of "exogenous states" (in the no adjustment case illiquid asset evolves endogenously)

        B = (1/Delta + rho)*speye(I*J*Nz) - A; % this and the next 3 steps implement eq. 15 p.6 Achdou et al. 2017 - Online Appendix

        u_stacked = [reshape(u(:,:,1),I*J,1);reshape(u(:,:,2),I*J,1)]; % stack into vector I*J*Nz
        V_stacked = [reshape(V(:,:,1),I*J,1);reshape(V(:,:,2),I*J,1)];

        vec = u_stacked + V_stacked/Delta;

        %IMPLICIT UPDATING
        V_stacked_12 = B\vec; %SOLVE SYSTEM OF EQUATIONS

        V(:,:,1) = reshape(V_stacked_12(1:I*J),I,J); % bring back into mesh form
        V(:,:,2) = reshape(V_stacked_12(I*J+1:I*J*Nz),I,J);

        Vchange = V - v;
        v = V;

        dist(n) = max(max(max(abs(Vchange))));
        disp(['Initial guess, no adjustment, iter ', int2str(n), ', dist ' num2str(dist(n)) ]);
        if dist(n)<tol
            disp(['Value Function Converged, Iteration = ' int2str(n)]);
            break
        end
    end
end

%% SOLVE WITH ADJUSTMENT

for n=1:maxit

    % except for Hf, Hb, H0 and the AssumeConcavity == 0, this part (up to adjustment decision for x grid) is the same as above without adjustment
    V = v;
    % forward difference
    Vbf(1:I-1,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/db;
    Vbf(I,:,:) = (zzz(I,:,:) + Rb(I,:,:).*bmax).^(-s); %will never be used, but impose state constraint a<=amax just in case
    % backward difference
    Vbb(2:I,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/db;
    Vbb(1,:,:) = (zzz(1,:,:) + Rb(1,:,:).*bmin).^(-s); %state constraint boundary condition
    
    Vbf = max(Vbf,10^(-6));
    Vbb = max(Vbb,10^(-6));
    
    %consumption and savings with forward difference
    cf = Vbf.^(-1/s);
    sf = zzz + Rb.*bbb - cf;
    Hf = cf.^(1-s)/(1-s) +Vbf.*sf;
    
    %consumption and savings with backward difference
    cb = Vbb.^(-1/s);
    sb = zzz + Rb.*bbb - cb;
    Hb = cb.^(1-s)/(1-s) +Vbb.*sb;
    
    %consumption and derivative of value function at steady state
    c0 = zzz + Rb.*bbb;
    %hamiltonians;
    
    % make a choice of forward or backward differences based on the sign of the drift
    if AssumeConcavity == 1 % same as above
        Ib = (sb < 0); %negative drift --> backward difference
        If = (sf > 0).*(1-Ib); %positive drift --> forward difference
        I0 = (1-If-Ib); %at steady state
    else % (*** reference)
        I0 = (1-(sf>0)) .* (1-(sb<0)); % points where sf <= 0 & sb >= 0
        Iunique = (sb<0).*(1-(sf>0)) + (1-(sb<0)).*(sf>0); % points where (sf <= 0 & sb < 0) U (sf > 0 & sb >= 0). "Unique" because we know here we can easily take sb
                                                           % over (sf <= 0 & sb < 0) and sf over (sf > 0 & sb >= 0)
        Iboth = (sb<0).*(sf>0);
        Ib = Iunique.*(sb<0) + Iboth.*(Hb>Hf); % split the points where (sb < 0 & sf > 0) based on u(cb) + Vbb * sb </> u(cf) + Vbf * sf: if > then backward. Note
                                               % that the other elements in the HJB are the same between backward and forward
        If = Iunique.*(sf>0) + Iboth.*(Hb<=Hf);
    end
%       if sum(sum(sum(Iunique==0)))>0
%         disp('convex problem');         
%     end    
    
    
   
    
    c = cf.*If + cb.*Ib + c0.*I0;
    u = c.^(1-s)/(1-s);
    
    %CONSTRUCT MATRIX
    X = -Ib.*sb/db;
    Y = -If.*sf/db + Ib.*sb/db;
    Z = If.*sf/db;
    
    updiag = zeros(I*J,Nz);
    lowdiag = zeros(I*J,Nz);
    centdiag = zeros(I*J,Nz);
    
    for i = 1:Nz
        centdiag(:,i) = reshape(Y(:,:,i),I*J,1);
    end

    lowdiag(1:I-1,:) = X(2:I,1,:);
    updiag(2:I,:) = Z(1:I-1,1,:);
    for j = 2:J
        lowdiag(1:j*I,:) = [lowdiag(1:(j-1)*I,:);squeeze(X(2:I,j,:));zeros(1,Nz)];
        updiag(1:j*I,:) = [updiag(1:(j-1)*I,:);zeros(1,Nz);squeeze(Z(1:I-1,j,:))];
    end
    
    
    AA1=spdiags(centdiag(:,1),0,I*J,I*J)+spdiags(updiag(:,1),1,I*J,I*J)+spdiags(lowdiag(:,1),-1,I*J,I*J);
    AA2=spdiags(centdiag(:,2),0,I*J,I*J)+spdiags(updiag(:,2),1,I*J,I*J)+spdiags(lowdiag(:,2),-1,I*J,I*J);
    AA = [AA1, sparse(I*J,I*J); sparse(I*J,I*J), AA2];
    
    A = AA + Bswitch;
    
    
    %%%%% Adjustment decision for x grid: this constructs objects on the predetermined grid for x, before interpolating at the x values corresponding to nodes (b, a)
    vstarAux = zeros(2,Nx); % 2 stays for Nz. The Aux elements are those defined on the grid for x
    bAdjAux  = zeros(2,Nx);
    aAdjAux  = zeros(2,Nx);
    vAdj = zeros(2,NaP); % for each x in the grid, search over a grid for a' - backing out b' = x - a' - ka (see ** reference) - to find the maximum. This vadj
                         % holds the values over such grid for a'
    aPmin = amin;
    G1 = griddedInterpolant(bb,aa,V(:,:,1)); % this is the interpolant for v over the nodes (b, a) at z(1)
    G2 = griddedInterpolant(bb,aa,V(:,:,2));

    for i = 1:Nx
        aPmax = min(x(i) - ka - bmin,amax); %CONSTRAINT a'<=aMax, taking into account the minimum value for b
        aP = linspace(aPmin,aPmax,NaP);
        aP = max(aP,x(i)-ka-bmax); %CONSTRAINT b'<=bMax
        bP = x(i) - ka - aP; % (** reference)
        vAdj(1,:) = G1(bP',aP'); % collect values over the grid for a' (and implied grid for b', given x)
        vAdj(2,:) = G2(bP',aP');

        [VstarAux(1,i),idx1] = max(vAdj(1,:)); % find maximum and argmax, and store maximum
        [VstarAux(2,i),idx2] = max(vAdj(2,:));
        aAdjAux(1,i) = aP(idx1); aAdjAux(2,i) = aP(idx2); % store argmax
        bAdjAux(1,i) = bP(idx1); bAdjAux(2,i) = bP(idx2);
    end
    
    %%%%% Interpolate for all possible values of a+b:
    Vstar = [lininterp1(x',VstarAux(1,:)',sumGrid);lininterp1(x',VstarAux(2,:)',sumGrid)]; % recall that sumGrid is the vector of x=b+a over nodes (b,a)
    
  
    
    % SOLVE USING LCP  
    
    B = (1/Delta + rho)*speye(I*J*Nz) - A;

    u_stacked = [reshape(u(:,:,1),I*J,1);reshape(u(:,:,2),I*J,1)];
    V_stacked = [reshape(V(:,:,1),I*J,1);reshape(V(:,:,2),I*J,1)];
    
    vec = u_stacked + V_stacked/Delta; % from B = ... to here it's again the same as in the without adjustment case
    
    % we need to cast the problem in the LCP solver form:
    % z' * (M * z + q) = 0
    % z >= 0
    % M * z + q >= 0
    % 
    % from the discretized HJBVI we have z = v(n+1) - Vstar, M=B so that the problem would be
    % z' * (B * v(n+1) - vec) = 0
    % z >= 0
    % B * v(n+1) - vec >= 0
    % but then we can just add and subtract B*Vstar to get the desired form, implying q = -vec + B*Vstar
    q = -vec + B*Vstar;
    
    %using Yuval Tassa's Newton-based LCP solver, download from http://www.mathworks.com/matlabcentral/fileexchange/20952
    z0 = V_stacked-Vstar; lbnd = zeros(I*J*Nz,1); ubnd = Inf*ones(I*J*Nz,1); % set initial value z0 where we guess v(n+1)=v(n) and lower/upper bound for solution z
    z = LCP(B,q,lbnd,ubnd,z0,0); % last argument is display of iteration data
   
    LCP_error = max(abs(z.*(B*z + q))); % check the accuracy of the LCP solution (should all be 0)
    if LCP_error > 10^(-6)
        disp('LCP not solved')
        %break
    end
    
    V_stacked = z+Vstar; %calculate value function
    
    
    V(:,:,1) = reshape(V_stacked(1:I*J),I,J); % this part is again the same as in the without adjustment case
    V(:,:,2) = reshape(V_stacked(I*J+1:I*J*Nz),I,J);
    
    
    Vchange = V - v;
    v = V;

    dist(n) = max(max(max(abs(Vchange))));
    disp(['With adjustment, adjustment, iter ', int2str(n), ', dist ' num2str(dist(n)) ]);
    if dist(n)<tol
        disp(['Value Function Converged, Iteration = ' int2str(n)]);
        break
    end
end



plot(dist(n-20:n))

ss = zzz + Rb.*bbb - c; % solution for savings in liquid asset b

adj = V_stacked==Vstar; %INDICATOR FOR ADJUSTING
adj = (abs(V_stacked - Vstar)<10^(-6)); % allow for some tolerance in defining solution for adjustment decision

adjRegion1 = reshape(adj(1:I*J),I,J).*ones(I,J); % multiplication by ones(I,J) transforms logicals into floats
adjRegion2 = reshape(adj(I*J+1:end),I,J).*ones(I,J);
adjRegion(:,:,1)=adjRegion1; adjRegion(:,:,2)=adjRegion2; % this reshapes the adj. decision into (I, J, Nz)

%AJUSTMENT TARGETS
bAdj1 = lininterp1(x',bAdjAux(1,:)',sumGrid); % this defines the solution for (b',a') conditional on adjusting on the original grid (b, a) from the grid on x
bAdj2 = lininterp1(x',bAdjAux(2,:)',sumGrid);
aAdj1 = lininterp1(x',aAdjAux(1,:)',sumGrid);
aAdj2 = lininterp1(x',aAdjAux(2,:)',sumGrid);

bAdj(:,:,1) = reshape(bAdj1,I,J);
bAdj(:,:,2) = reshape(bAdj2,I,J);
aAdj(:,:,1) = reshape(aAdj1,I,J);
aAdj(:,:,2) = reshape(aAdj2,I,J);

%% KF EQUATION

if max(abs(sum(A, 2))) > 1e-10
    warning('Some rows of A do not sum to 0')
end

% get list of points in the adjustment region
pointslist = (1:L)';
adjpoints = pointslist(adj);
noadjpoints = pointslist(~adj);

% BUILD M MATRIX
% first, put 1's on diagonal for points not in adjustment region
M = sparse(noadjpoints, noadjpoints, ones(length(noadjpoints),1), L, L);
% next, take care of points in adjustment region
pointsmat = reshape(pointslist, I, J, Nz);
count_idwpoints = 0; % count number of times inverse distance weights are used to allocate mass to be rebalanced
for mcount = 1:length(adjpoints)
    m = adjpoints(mcount);
    if m <= I*J % get (b',a') from point (b,a)
        bprime = bAdj1(m);
        aprime = aAdj1(m);
        zind = 1;
    else
        bprime = bAdj2(m - I*J);
        aprime = aAdj2(m - I*J);
        zind = 2;
    end
    bprime_left = discretize(bprime, b); % this gives the index of the left edge for liquid asset
    aprime_left = discretize(aprime, a); % same for illiquid asset
    bprime_right = bprime_left + 1;
    aprime_right = aprime_left + 1;

    % map from grid indexes to points
    point11 = pointsmat(bprime_left, aprime_left, zind);
    point12 = pointsmat(bprime_left, aprime_right, zind);
    point21 = pointsmat(bprime_right, aprime_left, zind);
    point22 = pointsmat(bprime_right, aprime_right, zind);
    neighpoints = [point11 point12;
                    point21 point22];
    % check if each point is not an adjustment point
    pointscheck = reshape(ismember(neighpoints(:), adjpoints), 2, 2); % 0 = not an adjustment point, so keep it
    pointscheck = 1 - pointscheck;
    if sum(sum(pointscheck)) == 4 && BilinW == 1 % use weights proportional to area opposite the point
        totarea = (b(bprime_right) - b(bprime_left)) * (a(aprime_right) - a(aprime_left)); % this is just da * db with equispaced grid 
        weights = totarea^(-1) * [b(bprime_right) - bprime; bprime - b(bprime_left)] * [a(aprime_right) - aprime aprime - a(aprime_left)];
    else
        count_idwpoints = count_idwpoints + 1; % use inverse distance weights
        neigh_bvals = [ones(1,2) * b(bprime_left);
                     ones(1,2) * b(bprime_right)];
        neigh_avals = [ones(2,1) * a(aprime_left) ones(2,1) * a(aprime_right)];
        bprime_mat = ones(2,2) * bprime;
        aprime_mat = ones(2,2) * aprime;
        dist_points = sqrt((neigh_bvals - bprime_mat).^2 + (neigh_avals - aprime_mat).^2); % this is the Euclidean distance
        inverse_dist = (1 ./ dist_points) .* pointscheck; % kill the points which are adjustment points
        if sum(sum(isinf(inverse_dist))) == 1
            weights = isinf(inverse_dist) .* ones(2,2);
        else
            weights = inverse_dist / sum(sum(inverse_dist));
        end
    end
    M = M + sparse(m * ones(4, 1), neighpoints(:), weights(:), L, L);
end

display(['Number of points where IDWs were used = ', num2str(count_idwpoints)])

AMT = transpose(A*M);

% Solve for stationary distribution (time-iteration)
Deltat = 25;
g_stacked = ones(L,1)./(L*da*db); % initial distribution
for n = 1:maxit
    g_nplus12 = transpose(M)*g_stacked; % apply M to ensure g has no mass in adjustment region
    g_nplus1 = (speye(L) - Deltat*AMT)\g_nplus12;   % implicit method
    if (sum(abs(g_nplus1-g_stacked))) < tol
        disp(['Distribution Converged, Iteration = ' int2str(n)]);
        break
    end
    g_stacked = g_nplus1;
end

if max(abs(AMT * g_stacked)) > 1e-6
    warning('Solution to KFE not accurate')
end
if sum(g_stacked < -1e-6) ~= 0
    warning('Solution to KFE has negative values')
end

g(:,:,1) = reshape(g_stacked(1:I*J),I,J);
g(:,:,2) = reshape(g_stacked(I*J+1:I*J*2),I,J);

%{
% Solve for stationary distribution (one-shot method)

D = sparse(adjpoints, adjpoints, ones(length(adjpoints),1), L, L);
AtildeT = D + AMT;

% Fix one value so matrix isn't singular:
vec = zeros(L,1);
rng(10,'twister');
iFix = noadjpoints(randi(length(noadjpoints),1));
vec(iFix)=.01;
AtildeT(iFix,:) = [zeros(1,iFix-1),1,zeros(1,L-iFix)];

% Solve system:
g_stacked = AtildeT\vec;
g_sum = g_stacked'*ones(L,1)*da*db;
g_stacked = g_stacked./g_sum;
if max(abs(AMT * g_stacked)) > 1e-6
    warning('Solution to KFE not accurate')
end
if sum(g_stacked < -1e-6) ~= 0
    warning('Solution to KFE has negative values')
end

g(:,:,1) = reshape(g_stacked(1:I*J),I,J);
g(:,:,2) = reshape(g_stacked(I*J+1:I*J*2),I,J);
%}

%%%%%%%%%%%
% FIGURES %
%%%%%%%%%%%

figure
set(gcf,'PaperPosition',[0 0 15 5])
subplot(1,2,1)
surf(b,a,adjRegion1')
view([0 90])
xlabel('Liquid Wealth (b)','FontSize',14, 'interpreter','latex')
ylabel('Illiquid Wealth (a)','FontSize',14, 'interpreter','latex')
title('Low Productivity')
shading flat
grid off
box on

subplot(1,2,2)
surf(b,a,adjRegion2')
view([0 90])
xlabel('Liquid Wealth (b)','FontSize',14, 'interpreter','latex')
ylabel('Illiquid Wealth (a)','FontSize',14, 'interpreter','latex')
title('High Productivity')
shading flat
grid off
box on
print -depsc adjRegion.eps

figure
set(gcf,'PaperPosition',[0 0 15 5])
subplot(1,2,1)
surf(b,a,adjRegion1')
view([0 90])
xlabel('Liquid Wealth (b)','FontSize',14, 'interpreter','latex')
ylabel('Illiquid Wealth (a)','FontSize',14, 'interpreter','latex')
title('Low Productivity')
shading flat
grid off
box on

subplot(1,2,2)
surf(b,a,adjRegion2')
view([0 90])
xlabel('Liquid Wealth (b)','FontSize',14, 'interpreter','latex')
ylabel('Illiquid Wealth (a)','FontSize',14, 'interpreter','latex')
title('High Productivity')
shading flat
grid off
box on

ss_plot = ss;
ss_plot(adjRegion == 1)=NaN;

figure
subplot(1,2,1)
plot(b,ss_plot(:,:,1)')
xlabel('Liquid Wealth (b)','FontSize',14, 'interpreter','latex')
title('Liquid Savings, Low Income')

subplot(1,2,2)
plot(b,ss_plot(:,:,2)')
xlabel('Liquid Wealth (b)','FontSize',14, 'interpreter','latex')
title('Liquid Savings, High Income')

c_plot = c;
c_plot(adjRegion == 1)=NaN;

figure
subplot(1,2,1)
plot(b,c_plot(:,:,1)')
xlabel('Liquid Wealth (b)','FontSize',14, 'interpreter','latex')
title('Consumption, Low Income')

subplot(1,2,2)
plot(b,c_plot(:,:,2)')
xlabel('Liquid Wealth (b)','FontSize',14, 'interpreter','latex')
title('Consumption, High Income')

bAdj_plot = bAdj;
bAdj_plot(adjRegion ~= 1)=NaN;
% bAdj_plot = bAdj_plot ./ bbb - 1;

aAdj_plot = aAdj;
aAdj_plot(adjRegion ~= 1)=NaN;
% aAdj_plot = aAdj_plot ./ aaa - 1;

% adjustment targets
figure
set(gcf,'PaperPosition',[0 0 15 5])
subplot(1,2,1)
surf(b,a,bAdj_plot(:,:,1)')
xlabel('Liquid Wealth (b)','FontSize',14, 'interpreter','latex')
title('Liquid Adj. Target, Low Income')

subplot(1,2,2)
surf(b,a,bAdj_plot(:,:,2)')
xlabel('Liquid Wealth (b)','FontSize',14, 'interpreter','latex')
title('Liquid Adj. Target, High Income')
print -depsc adjTarget_b.eps

figure
set(gcf,'PaperPosition',[0 0 15 5])
subplot(1,2,1)
surf(b,a,aAdj_plot(:,:,1)')
xlabel('Liquid Wealth (b)','FontSize',14, 'interpreter','latex')
title('Illiquid Adj. Target, Low Income')

subplot(1,2,2)
surf(b,a,aAdj_plot(:,:,2)')
xlabel('Liquid Wealth (b)','FontSize',14, 'interpreter','latex')
title('Illiquid Adj. Target, High Income')
print -depsc adjTarget_a.eps

% KF equation
figure
set(gcf,'PaperPosition',[0 0 15 5])
subplot(1,2,1)
surf(b,a,g(:,:,1)')
xlabel('Liquid Wealth (b)','FontSize',14, 'interpreter','latex')
ylabel('Illiquid Wealth (a)','FontSize',14, 'interpreter','latex')
title('Stationary Density, Low Income')
subplot(1,2,2)
surf(b,a,g(:,:,2)')
xlabel('Liquid Wealth (b)','FontSize',14, 'interpreter','latex')
ylabel('Illiquid Wealth (a)','FontSize',14, 'interpreter','latex')
title('Stationary Density, High Income')
print -depsc distribution.eps
