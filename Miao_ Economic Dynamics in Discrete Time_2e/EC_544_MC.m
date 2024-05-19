%% EC 544
%% Markov Chain Examples

close all
clear all

%%Finite State Example 1
P = zeros(6,6);
P(1,1) = 0.5;
P(1,2) = 0.5;
P(2,3) = 1;
P(3,1) = 1/3;
P(3,4) = 1/3;
P(3,5) =1/3;
P(4,4) = 1/2;
P(4,5) =1/2;
P(5,6) = 1;
P(6,5) = 1;
mc = dtmc(P);
figure;
 graphplot(mc,'LabelEdges',true,'ColorNodes',true)
print -depsc ex1.eps
isreducible(mc)
asymptotics(mc)
classify(mc)
isergodic(mc)
rng(1); % For reproducibility
numSteps = 20;
x0= zeros(1,6);
x0(1)=2;
X = simulate(mc,numSteps,'X0',x0)
figure;
simplot(mc,X);


%% Limiting distribution example 3
P =[0 0 3/4 1/4;
    0 0 1/4 3/4;
    3/4 1/4 0 0;
    1/4 3/4 0 0];

mc = dtmc(P);
figure;
 graphplot(mc,'LabelEdges',true,'ColorNodes',true)
print -depsc ex3.eps

%% Stationary distribution
P = [.9 .075 .025 ; .8 .15 .05; .25 .25 .5];
mc = dtmc(P,'StateNames',["Bull" "Bear" "Stagnant"]);
figure;
graphplot(mc,'LabelEdges',true,'ColorNodes',true)
title(" Market MC1")
stat_distribution_1 = asymptotics(mc)
stat_distribution_2 = P^100
[eig_verifiy, eigs] = eig(P');
stat_distribution_3 = eig_verifiy/sum(eig_verifiy)


P = [0 .5 .5 ; .5 .5 0; 0 0 1];
mc = dtmc(P,'StateNames',["Bull" "Bear" "Stagnant"]);
figure;
graphplot(mc,'LabelEdges',true,'ColorNodes',true)
title(" Market MC2")



%% Periodicity 
periodic_mat = [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 0 0 0 ];
periodic_MC = dtmc(periodic_mat)
[bins,ClassStates,ClassRecurrence,ClassPeriod] = classify(periodic_MC)

figure
graphplot(periodic_MC,'ColorEdges',true,'ColorNodes',true)
title('Periodic MC with d = 4')

aperiodic_mat =  [.4 .6; .6 0.4]; 
aperiodic_MC = dtmc(aperiodic_mat)
[bins_2,ClassStates_2,ClassRecurrence_2,ClassPeriod_2] = classify(aperiodic_MC)
figure
graphplot(aperiodic_MC,'ColorEdges',true,'ColorNodes',true)
title('Aperiodic MC')

%% Many states
numStates = 8;
Zeros = 50;
stateNames = strcat(repmat("Regime ",1,8),string(1:8));

rng(1) % For reproducibility
mc = mcmix(8,'Zeros',Zeros,'StateNames',stateNames);
figure();
graphplot(mc,'ColorNodes',true);
title('Complicated MC ')

