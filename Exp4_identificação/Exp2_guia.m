% Experimento sobre técnicas de sincronia para controladores PID
%%

% Resposta ao degrau em malha aberta
s = tf('s');
G = exp(-s)/(s+1)^4;

figure(1);
step(G);
title('Resposta ao degrau em malha aberta');
grid on;

% Parâmetros do controlador PI

% Parâmetros obtidos pela regra de aproximação de Half-rule
L = 7/2; K = 1; T = 3/2;

% Método Cohen-Coon
Kp1 = (T/(K*L))*(0.9 + L/(12*T));
Ti1 = (L*(30*T + 3*L))/(9*T + 20*L);

Cpi1 = pid(Kp1, Kp1/Ti1);
H1 =feedback(Cpi1*G,1);

% figure(2);
% step(H1);
%title('Controlador PI pelo método de Cohen-Coon');
% grid on;
% stepinfo(H1)

% Método de SIMC com \tau_c = \theta
Kp2 = T/(2*K*L);
Ti2 = min(T, 8*L);

Cpi2 = pid(Kp2, Kp2/Ti2);
H2 =feedback(Cpi2*G,1);

% figure(3);
% step(H2);
%title('Controlador PI pelo método de SIMC com \tau_c = \theta');
% grid on;
% stepinfo(H2)

% Parâmetros do controlador PID

% Método Cohen-Coon
Kp3 = (T/(K*L))*(4/3 + L/(4*T));
Ti3 = (L*(32*T + 6*L))/(13*T + 8*L);
Td3 = (4*L*T)/(11*T + 2*L);

Cpid1 = pid(Kp3, Kp3/Ti3, Kp3*Td3);
H3 =feedback(Cpid1*G,1);

figure(4);
step(H3);
title('Controlador PID pelo método de Cohen-Coon');
grid on;
stepinfo(H3)

% Parâmetros obtidos pela regra de aproximação de Half-rule de 2ªordem
K = 1; tau1 = 1; tau2 = 3/2; theta = 5/2;

% Método de SIMC com \tau_c = \theta
Kp4 = tau1/(2*K*theta);
Ti4 = min(tau1, 8*theta);
Td4 = tau2;

Cpid2 = pid(Kp4, Kp4/Ti4, Kp4*Td4);
H4 =feedback(Cpid2*G,1);

figure(5);
step(H4);
title('Controlador PID pelo método de SIMC com \tau_c = \theta');
grid on;
stepinfo(H4)
