% Equivalente de segunda ordem com atraso

s=tf('s'); G=1/((1+s)*(1+0.5*s)*(1+0.25*s));
%G = exp(-s)/(s+1)^4;
% figure(1);
% step(feedback(G,1));

% 1ª ordem
T = 0;
taui = [1, 0.5, 0.25];
theta0 = 0;
K = 1;
tau1 = taui(1) + taui(2)/2;
theta = theta0 + taui(2)/2 + sum(taui(3:end)) + sum(T);
G1a = K*exp(-theta*s)/(tau1*s + 1);

% Calculando os parâmetros do controlador PI de acordo com o SIMC
Kp1 = (1/K)*(tau1/(tau_c+theta));
Ti1 = min(tau1, 4*(tau_c + theta));
Ki1 = Kp1/Ti1;

Kpi = Kp1 + Ki1/s;
H1 = feedback(G*Kpi,1);
% figure(2);
% step(H1);

% 2ª ordem
T = 0;
taui = [1, 0.5, 0.25];
theta0 = 0;
K = 1;
tau1 = taui(1);
tau2 = taui(2) + taui(3)/2;
theta = theta0 + taui(3)/2 + sum(taui(4:end)) + sum(T);
G2a = K*exp(-theta*s)/((tau1*s + 1)*(tau2*s + 1));

% Calculando os parâmetros do controlador PID de acordo com o SIMC
Kp2 = (1/K)*(tau1/(tau_c+theta));
Ti2 = min(tau1, 4*(tau_c + theta));
Td2 = tau2;
Ki2 = Kp2/Ti2;
Kd2 = Kp2*Td2;

Kpid2 = Kp2 + Ki2/s + Kd2*s/(1+0.01*s);
H2 = feedback(G*Kpid2,1);
figure(3);
step(H2);
hold on;

% Criterios antes da otimização
MSb=norm(feedback(1,G*Kpid2),inf)
MTb=norm(feedback(G*Kpid2,1),inf)
Jvb=norm(feedback(G/s,Kpid2),inf)
Jub=norm(feedback(Kpid2, G),inf)

% 
MS_max=1.7;
MT_max=1.3;
Ju_max=250;
x = [Kp2, Ki2, Kd2];

options = optimset('Algorithm','active-set');
x=fmincon(@(x) objfun(x,s,G),x,[],[],[],[],...
[], [], @(x)confun(x,s,G,MS_max,MT_max,Ju_max), options);

Kp = x(1); Ki = x(2); Kd = x(3);
Kpid = Kp + Ki/s + Kd*s/(1+0.01*s);
H = feedback(G*Kpid,1);

%Critérios após a otimização
Jv=norm(feedback(G/s,Kpid),inf)
Ju=norm(feedback(Kpid, G),inf)
MS=norm(feedback(1,G*Kpid),inf)
MT=norm(feedback(G*Kpid,1),inf)
step(H);
legend('Não otimizado', 'Otimizado');






