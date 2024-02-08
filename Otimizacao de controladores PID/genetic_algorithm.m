%% Projeto de controladores PID por meio de Genetic Algorithm (GA)

s=tf('s'); 
% G=1/(1+s)^3;                        % G1
G=exp(-0.3*s)/((1+s)*(1+0.5*s));   % G2


% Obtendo as aproximaçãoes de Half-rule
T = 0;
% taui = [1, 1, 1];     % G1
taui = [1, 0.5, 0];   % G2
theta0 = 0.3;
K = 1;
tau1 = taui(1);
tau2 = taui(2) + taui(3)/2;
theta = theta0 + taui(3)/2 + sum(taui(4:end)) + sum(T);
G2a = K*exp(-theta*s)/((tau1*s + 1)*(tau2*s + 1));

% Calculando os parâmetros do controlador PID de acordo com o SIMC
tau_c = theta;
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
MSb=norm(feedback(1,pade(G)*Kpid2),inf)
MTb=norm(feedback(pade(G)*Kpid2,1),inf)
Jvb=norm(feedback(pade(G)/s,Kpid2),inf)
Jub=norm(feedback(Kpid2, pade(G)),inf)

% Critérios de restrição
MS_max=1.6;
MT_max=1.05;
%Jv_max=0.66;
Ju_max=100;
% Gm_Max = 3;
% Pm_Max = 45;

% x = [Kp2, Ki2, Kd2];

% Otimização do controlador por meio do algoritmo genético
lb = [2, 0.5, 0.1];
ub = [30, 10, 10];
options = optimoptions('ga', 'display', 'iter');
x = ga(@(x) objfun(x,s,G),3,[],[],[],[],...
lb,ub,@(x)confun(x,s,G,MS_max,MT_max,Ju_max), options);

Kp = x(1); Ki = x(2); Kd = x(3);
Kpid = Kp + Ki/s + Kd*s/(1+0.01*s);
H = feedback(G*Kpid,1);

% Critérios após a otimização
Jv=norm(feedback(pade(G)/s,Kpid),inf)
Ju=norm(feedback(Kpid, pade(G)),inf)
MS=norm(feedback(1,pade(G)*Kpid),inf)
MT=norm(feedback(pade(G)*Kpid,1),inf)
step(H);
stepinfo(H)
title('Planta G_1');
legend('Controlador SIMC', 'Controlador otimizado', 'location', 'best','FontSize', 10);
grid on;
