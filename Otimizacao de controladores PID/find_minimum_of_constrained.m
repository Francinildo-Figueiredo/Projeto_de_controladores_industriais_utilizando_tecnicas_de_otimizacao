%% Otimização de controladores por meio da função 

clc   %clears the command window
clear    % clears the previous work space
close all    % closes the privous graphical objects(figures)

s=tf('s'); 
% G=1/(1+s)^3;                     % G1
G=exp(-0.3*s)/((1+s)*(1+0.5*s));   % G2


% Obtendo as aproximaçãoes de Half-rule
T = 0;
% taui = [1, 1, 1];   % G1
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

t = 0:0.1:10;
y1 = step(H1, t);
ISE1 = sum((1-y1).^2)
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
MS_max=1.7;
MT_max=1.05;
Ju_max=100;
% Gm_Max = 3;
% Pm_Max = 45;

x1 = [Kp1, Ki1, Kd1];
% x1 = [0,0,0];

% Otimização do controlador por meio da fmincon
lb = [0, 0, 0];
ub = [30, 10, 10];
options = optimoptions('fmincon', 'display', 'iter');
x2=fmincon(@(x) objfun(x,s,G),x1,[],[],[],[],...
lb, ub, @(x)confun(x,s,G,MS_max,MT_max,Ju_max), options);

Kp2 = x2(1); Ki2 = x2(2); Kd2 = x2(3);
Kpid2 = Kp2 + Ki2/s + Kd2*s/(1+0.01*s);
H2 = feedback(G*Kpid2,1);

% Critérios após a otimização
Jv_fmincon=norm(feedback(pade(G)/s,Kpid2),inf)
Ju_fmincon=norm(feedback(Kpid2, pade(G)),inf)
MS_fmincon=norm(feedback(1,pade(G)*Kpid2),inf)
MT_fmincon=norm(feedback(pade(G)*Kpid2,1),inf)
[Gm_fmincon,Pm_fmincon] = margin(G*Kpid2)

t = 0:0.1:10;
y2 = step(H2, t);
ISE2 = sum((1-y2).^2)
step(H2);
% title('Planta G_1');
legend('SIMC', 'fmincon', 'best','FontSize', 10);
grid on;
x2

