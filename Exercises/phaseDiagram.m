% Phase Diagram from Moll lectures _ self practice

% dot lambda = lambda(rho + delta -F'(k))
% dot k = F(k) - delta*k - (lambda/alpha)^(1/(alpha-1))
clear all

%calibration

alpha = 0.3; %preference param
Alpha = 0.34;
dt = 0.01;
%ss lam, k, Alpha

rho = 0.05;
delta = 0.08;
kei = 0.2:0.1:0.4;
lambda = 1:0.01:2;
%lamm=double.empty(length(kei),length(lambda),0);
%kk=double.empty(length(kei),length(lambda),0);

for j = 1:length(kei)

     for jj = 1:length(lambda) %:0.1:0.2
        k(1)=kei(j);
        lam(1)=lambda(jj);
        for i=1:5000
            k(i+1)=k(i) + (k(i)^(Alpha) - delta*k(i) - (lam(i)/alpha)^(1/(alpha-1)))*dt;
            lam(i+1)=lam(i) +(lam(i)*(rho + delta-  Alpha*k(i)^(Alpha-1)))*dt;
            if k(i+1)>0 & lam(i+1)>0
                continue
            else
                break
            end

        end
        %hold off
        
        plot(lam,k)
        lamm=lam;
        kk=k;
        clear lam k
        hold on

        
    end
end
hold off
% 
% for i =0:1
%     for j=1:10
%         
%         plot(i+1:i+10)
%         hold on
%     end
% end
%  plot(1:10)
%  hold on 
%  plot(11:20)
%  hold off
