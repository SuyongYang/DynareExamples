%% neal2.m%%  amended version, tjs Sept 13, 2001%% corrections to match value functions in the text (Chao had used old ones%% before we changed the model.%% %% this version generates graphs neal1.new.eps and neal2new.eps for the second%% printing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%			Neal's model of career choice
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %
%  Neal's worker chooses career, job(theta,epsilon) pair given the distribution 
%  of career and job, which is uniform in this program. He can either reject 
%  current pair in hand, or reject current job and keep career, or retain his
%  existing pair

clear
diary neal.out
disp('model of carrer choice');
disp('');
%
%  set parameter values
%
%   form job and career grid
%   
% maxjob = 1.0;                     % maximum value of job grid     maxjob=5.0;  incjob=.1;  
%incjob = 0.1;                     % size of job grid increments
njob = round(maxjob/incjob);      % number of grid points
%epsilon = [0.1:0.1:1.0]';         % the states of jobepsilon=[incjob:incjob:maxjob]';
maxcar = 5.0;                     % maximum value of career grid   
inccar = 0.1;                     % size of career grid increments
ncar = round(maxcar/inccar);      % number of grid points
theta = [inccar:inccar:5.0]';           % the states of career

n=njob;
% uniform distributionsfprob = (1/n)*ones(n,1);           % probability of states of jobs   gprob =(1/n)*ones(n,1);           % probability of states of careers beta=0.95;
   %
   %  initialize some variables
   %
   v       = zeros(n);      
   decis   = zeros(1,n*n); % written into row vector for future use
   iter    = 0;
   test    = 10;
   
   %  iterate on Bellman's equation and get the decision 
   %  rules and the value function at the optimum         
   %
   while test ~= 0;
      
      %  value function v will be n by n matrix, its (i,j) element denotes the maximum
      %  value given the ith career and jth job in hand
      
      v1 = (theta*ones(1,n)+ones(n,1)*epsilon')/(1-beta); % value of not drawing 
      v2 = theta*ones(1,n)+(ones(n,1)*epsilon'+beta*v)*fprob*ones(1,n); % value of drawing a new job
      v3 = gprob'*(theta*ones(1,n)+ones(n,1)*epsilon'+beta*v)*fprob*ones(n,n); % value of drawing both job and career
      
      v11 = reshape(v1,[1 n*n]); %
      v22 = reshape(v2,[1 n*n]); % reshape values into 1 by n*n vector, first by 
      v33 = reshape(v3,[1 n*n]); % careers, then by jobs
      
      vstack = [v11;v22;v33]; % stack value functions for comparation across column
      
      [tvstack,tdecis] = max(vstack); 
      
      
      tv = reshape(tvstack,[n n]);
      test=max(any(tdecis-decis));
      v=tv;
      decis=tdecis;
      iter=iter+1;
      
   end;
   
   % decis=1 means value of not drawing is highest
   % decis=2 means value of redrawing a job is highest
   % decis=3 means value of redrawing both is highest

decis = reshape(decis,[n n]);
 
disp('v')
v
disp('decis')
decis

[x,y] = meshgrid(theta,epsilon);
colormap(hsv)
surf(x,y,v')
title('Graph of Career Choice')
xlabel('career choice (\theta)')
ylabel('job choice (\epsilon)')
zlabel('v(\theta,\epsilon)')print -depsc neal1new
figure(2)contourf(theta,epsilon,decis')%caxis([-5 5])%colormap(cool)colormap(gray)text(2,2.5,'new life')text(4.25,2.5,'new job')xlabel('\theta'), ylabel('\epsilon')print -depsc neal2newreturn
figure
[x,y] = meshgrid(theta,epsilon);
colormap(cool)
surf(x,y,decis')
title('Graph of Career Choice')
xlabel('career choice')
ylabel('job choice')
zlabel('decision rule of career-job pair')

diary off