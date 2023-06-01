% Técnicas de Sintonia para Controladores PID - CSTR
%%

% Ponto de operação
Fs   = 1;     % Vazão volumétrica (m^3/h)
Cafs = 10;    % Concentração de A no fluxo de alimentação ((kg*mol)/m^3)
Tfs  = 298;   % Temperatura de alimentação (K)
Tjs  = 298;   % Temperatura da capa (K)
Ts   = 311.2; % Temperatura do reator (K)
Cas  = 8.564; % Concentração de A no reator ((kg*mol)/m^3)

RefCa = 9; dCafs = 0; dTfs = 0; dTjs = 0;

% Função de transferência identificada que relaciona a entrada F e a saída
% Ca
G = tf(3.153, [2.186, 1]);
G.ioDelay = 0.106;

% Método Cohen-Coon
L = 0.106; K = 3.153; T = 2.186;

% Controlador PI
Kp1 = (T/(K*L))*(0.9 + L/(12*T));
Ti1 = (L*(30*T + 3*L))/(9*T + 20*L);
Ki1 = Kp1/Ti1;

%Cpi1 = pid(Kp1, Kp1/Ti1);

% Controlador PID
Kp2 = (T/(K*L))*(4/3 + L/(4*T));
Ti2 = (L*(32*T + 6*L))/(13*T + 8*L);
Td2 = (4*L*T)/(11*T + 2*L);
Ki2 = Kp2/Ti2;
Kd2 = Kp2*Td2;

%Cpid = pid(Kp2, Kp2/Ti2, Kp2*Td2);

% Método SIMC para \tau_c = \theta
K = 3.153; theta = 0.106; tau1 = 2.186; tau_c = theta;

Kp3 = (1/K)*(tau1/(tau_c+theta));
Ti3 = min(tau1, 4*(tau_c + theta));
Ki3 = Kp3/Ti3;
Kpid3 = Kp3 + Ki3/s;

figure(1);
out1 = sim('CSTR_NL_SIMC', 'ReturnWorkspaceOutputs', 'on');
plot(out1.CaNL);
grid on;
hold on;

% Criterios antes da otimização
MSb=norm(feedback(1,pade(G)*Kpid3),inf)
MTb=norm(feedback(pade(G)*Kpid3,1),inf)
Jvb=norm(feedback(pade(G)/s,Kpid3),inf)
Jub=norm(feedback(Kpid3, pade(G)),inf)

MS_max=1.7;
MT_max=1.3;
%Jv_max=0.66;
Ju_max=20;
x = [Kp3, Ki3, 0];

options = optimset('Algorithm','active-set');
x=fmincon(@(x) objfun(x,s,G),x,[],[],[],[],...
[], [], @(x)confun(x,s,G,MS_max,MT_max,Ju_max), options);

Kp = x(1); Ki = x(2); Kd = x(3);
Kpid = Kp + Ki/s + Kd*s/(1+0.01*s);

%Critérios após a otimização
Jv=norm(feedback(pade(G)/s,Kpid),inf)
Ju=norm(feedback(Kpid, pade(G)),inf)
MS=norm(feedback(1,pade(G)*Kpid),inf)
MT=norm(feedback(pade(G)*Kpid,1),inf)
out1 = sim('CSTR_NL_Otimizado', 'ReturnWorkspaceOutputs', 'on');
plot(out1.CaNL);
grid on;
title('Planta G_{CSTR}');
legend('Controlador SIMC', 'Controlador otimizado', 'location', 'best','FontSize', 10);
xlabel('t (h)');
ylabel('Ca (kgmol/m^3)');