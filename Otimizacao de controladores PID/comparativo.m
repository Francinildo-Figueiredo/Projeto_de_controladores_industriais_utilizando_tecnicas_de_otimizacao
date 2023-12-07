%% Análise comparativa entre os algoritmos de otimização 

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
Pm_Max = 45;

x1 = [Kp1, Ki1, Kd1];

% Otimização do controlador por meio do algoritmo genético
lb = [0, 0, 0];
ub = [30, 10, 10];
options = optimoptions('ga', 'display', 'iter');
x1 = ga(@(x1) objfun(x1,s,G),3,[],[],[],[],...
lb,ub,@(x1)confun(x1,s,G,MS_max,MT_max,Ju_max, Gm_Max, Pm_Max), options);

Kp2 = x1(1); Ki2 = x1(2); Kd2 = x1(3);
Kpid2 = Kp2 + Ki2/s + Kd2*s/(1+0.01*s);
H2 = feedback(G*Kpid2,1);

% Critérios após a otimização
Jv_ga=norm(feedback(pade(G)/s,Kpid2),inf)
Ju_ga=norm(feedback(Kpid2, pade(G)),inf)
MS_ga=norm(feedback(1,pade(G)*Kpid2),inf)
MT_ga=norm(feedback(pade(G)*Kpid2,1),inf)
[Gm_ga,Pm_ga] = margin(G*Kpid2)

% Otimização dos parâmetros do controlador por meio da fmincon
lb = [0, 0, 0];
ub = [30, 10, 10];
x2 = [Kp1, Ki1, Kd1];
options = optimoptions('fmincon', 'display', 'iter');
x2=fmincon(@(x2) objfun(x2,s,G),x2,[],[],[],[],...
lb, ub, @(x2)confun(x2,s,G,MS_max,MT_max,Ju_max, Gm_Max, Pm_Max), options);

Kp3 = x2(1); Ki3 = x2(2); Kd3 = x2(3);
Kpid3 = Kp3 + Ki3/s + Kd3*s/(1+0.01*s);
H3 = feedback(G*Kpid3,1);

% Critérios após a otimização
Jv_fmincon=norm(feedback(pade(G)/s,Kpid3),inf)
Ju_fmincon=norm(feedback(Kpid3, pade(G)),inf)
MS_fmincon=norm(feedback(1,pade(G)*Kpid3),inf)
MT_fmincon=norm(feedback(pade(G)*Kpid3,1),inf)
[Gm_fmincon,Pm_fmincon] = margin(G*Kpid3)

step(H2); hold on;
step(H3);
% title('Planta G_1');
legend('SIMC', 'GA', 'fmincon', 'best','FontSize', 10);
grid on;
