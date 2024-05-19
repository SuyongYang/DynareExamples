% Moll lecture practice 
% consumption and capital diagram

% dot lambda = lambda(rho + delta -F'(k))
% dot k = F(k) - delta*k - (lambda/alpha)^(1/(alpha-1))
% dot lambda = u''(c) dot c
% u''(c)*dot c = u'(c)(rho + delta - F'(k))
% u' = alpha*c^(alpha-1)
% u''= alpha*(alpha-1)*c^(alpha-2)
clear all

%calibration

alpha = 0.3; %preference param
Alpha = 0.34;
dt = 0.01;
%ss lam, k, Alpha
 
rho = 0.05; 
delta = 0.08;
%kei = 2:.1:5;
%con = .5:.1:4;

%lamm=double.empty(length(kei),length(lambda),0);
%kk=double.empty(length(kei),length(lambda),0);

SSk = (Alpha / (rho+delta))^(1/(1-Alpha));
SSc = SSk^Alpha - delta*SSk;
%SSk = 4.29, SSc=1.29

kei = [1.0 2.0 3 4.2 SSk 4.4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
con = [.1 .3 .5 .7 1 1.2 SSc 1.4 2 3 4 5 6 7 8 9 10 11];

for j = 1:length(kei) 

     for jj = 1:length(con) %:0.1:0.2
        k(1)=kei(j);
        c(1)=con(jj);
        for i=1:10000
            k(i+1)=k(i) + (k(i)^(Alpha) - delta*k(i) - c(i))*dt;
            c(i+1)=c(i)+(c(i)/(alpha-1)*(rho + delta-  Alpha*k(i)^(Alpha-1)))*dt;
            % u' = alpha*c^(alpha-1)
            % u''= alpha*(alpha-1)*c^(alpha-2)
            if k(i+1)>0 & c(i+1)>0
                continue
            else
                break
            end

        end
        %hold off
        
        plot(k,c)
        consumption=c;
        kk=k;
        clear c k
        hold on

        
    end
end
hold off