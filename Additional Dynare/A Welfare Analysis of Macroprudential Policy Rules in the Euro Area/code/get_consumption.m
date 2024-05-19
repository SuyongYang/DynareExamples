function [rez] = get_consumption(param)
% find the asymetric level of ss consumption
% asymetry is driven by relative country size

Y			= param(1);
gy			= param(2);
alpha_C_h	= param(3);
alpha_C_f	= param(4);
alpha_I_h	= param(5);
alpha_I_f	= param(6);
I			= param(7);
n			= param(8);
C			= param(9);

func = @(myx)([	Y*(1-gy) - (1-alpha_C_h)*myx(1) - (1-alpha_I_h)*I - ((1-n)/n)*( alpha_C_f*myx(2) + alpha_I_f*I ),
				Y*(1-gy) - (1-alpha_C_f)*myx(2) - (1-alpha_I_f)*I - (n/(1-n))*( alpha_C_h*myx(1) + alpha_I_h*I ) ,		
				]);
myx0=[C,C]; 
options=optimset('display','off','MaxFunEvals',3000,'MaxIter',3000); 
rez=fsolve(func,myx0,options);

end