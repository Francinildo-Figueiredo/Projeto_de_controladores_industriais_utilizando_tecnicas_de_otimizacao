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
G11 = tf(3.153, [2.186, 1]);
G11.ioDelay = 0.106;

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

% Controlador PI
%Cpi2 = pid(Kp3, Kp3/Ti3);
