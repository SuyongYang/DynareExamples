function results = KFE_SS(parameters,results)
% This codes computes the KFE
I    = parameters.I;                % numer of points
J    = parameters.J;                % numer of points
da   = results.da;         
dz   = results.dz;
AT   = (results.A)';
b = zeros(I*J,1);

%need to fix one value, otherwise matrix is singular

i_fix   = 1;
b(i_fix)= -1;


%Solve linear system
gg = AT\b;

gg = max(gg,0); % To avoid numerical errors 
g_sum = gg'*ones(I*J,1)*da*dz;
gg = gg./g_sum;

results.g = reshape(gg,I,J);

