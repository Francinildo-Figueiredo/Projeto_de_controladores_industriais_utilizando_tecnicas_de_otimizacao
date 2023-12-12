%% Otimização de controladores por meio de um algoritmo bioinspirado em um enxame de partículas
clc;close all; clear;
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
Kp1 = (1/K)*(tau1/(tau_c+theta));
Ti1 = min(tau1, 4*(tau_c + theta));
Td1 = tau2;
Ki1 = Kp1/Ti1;
Kd1 = Kp1*Td1;

Kpid1 = Kp1 + Ki1/s + Kd1*s/(1+0.01*s);
H1 = feedback(G*Kpid1,1);
figure(3);
step(H1);
hold on;

% Criterios antes da otimização
MSb=norm(feedback(1,pade(G)*Kpid1),inf)
MTb=norm(feedback(pade(G)*Kpid1,1),inf)
Jvb=norm(feedback(pade(G)/s,Kpid1),inf)
Jub=norm(feedback(Kpid1, pade(G)),inf)
[Gmb,Pmb] = margin(G*Kpid1)

% Critérios de restrição
MS_max=1.5;
MT_max=1.01;
%Jv_max=0.66;
Ju_max=100;
Gm_Max = 3;
Pm_Max = 65;

x1 = [Kp1, Ki1, Kd1];

% Otimização do controlador por meio do algoritmo genético
lb = [1, eps, eps];
ub = [10, 10, 10];
x2 = PSO_non_linear_constraint(@(x) objfun(x,s,G), @(x)confunC(x,s,G,...
MS_max,MT_max,Ju_max, Gm_Max, Pm_Max),@(x)[],lb,ub,3,600,20);

Kp2 = x2(1); Ki2 = x2(2); Kd2 = x2(3);
Kpid2 = Kp2 + Ki2/s + Kd2*s/(1+0.01*s);
H2 = feedback(G*Kpid2,1);

% Critérios após a otimização
Jv=norm(feedback(pade(G)/s,Kpid2),inf)
Ju=norm(feedback(Kpid2, pade(G)),inf)
MS=norm(feedback(1,pade(G)*Kpid2),inf)
MT=norm(feedback(pade(G)*Kpid2,1),inf)
[Gm,Pm] = margin(G*Kpid2)

step(H2);
% title('Planta G_1');
legend('SIMC', 'PSO', 'best','FontSize', 10);
grid on;