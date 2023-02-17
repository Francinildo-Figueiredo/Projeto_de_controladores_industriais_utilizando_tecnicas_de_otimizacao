% Definição dos parâmetros do reator
V     = 1;         % Volume do reator (m^3)
dE    = 11843;     % Energia de ativação (J/mol)
rhocp = 500;       % Densidade*Capacidade térmica (J/(K*m^3))
R     = 1.987;     % Constante de gás ideal (L/(mol*K))
UA    = 150;       % Coef. geral transf. calor*área para troca de calor (J/(Kh))
k0    = 9703*3600; % Fator pré-exponencial (h^-1)
dH    = -5960;     % Calor de reação (J/mol)

% Ponto de operação
Fs   = 1;     % Vazão volumétrica (m^3/h)
Cafs = 10;    % Concentração de A no fluxo de alimentação ((kg*mol)/m^3)
Tfs  = 298;   % Temperatura de alimentação (K)
Tjs  = 298;   % Temperatura da capa (K)
Ts   = 311.2; % Temperatura do reator (K)
Cas  = 8.564; % Concentração de A no reator ((kg*mol)/m^3)

% Modelando as equações dinâmicas na forma de variáveis de estado

% Definindo a expressão para a taxa de reação por unidade de volume
r = k0*exp(-dE/(R*Ts))*Cas;

% Da equação de balanço no componente A, obtemos
dCadt = (Fs/V)*(Cafs - Cas) - r;

% Da equação de balanço de energia, obtemos
dTdt  = (Fs/V)*(Tfs - Ts) + (-dH/rhocp)*r - (UA/(V*rhocp))*(Ts - Tjs);

% Obtendo o modelo linearizado em espaço de estados
drdCas = k0*exp(-dE/(R*Ts));
drdTs   = dE/(R*Ts^2)*r;

A11 = -Fs/V - drdCas;
A12 = -drdTs;
A21 = (-dH/rhocp)*drdCas;
A22 = -Fs/V + (-dH/rhocp)*drdTs - UA/(V*rhocp);
A = [A11, A12; A21, A22];

Bu1 = (Cafs - Cas)/V;
Bu2 = (Tfs - Ts)/V;
Bu = [Bu1; Bu2];

Bd11 = Fs/V;
Bd12 = 0;
Bd13 = 0;
Bd21 = 0;
Bd22 = Fs/V;
Bd23 = UA/(V*rhocp);
Bd = [Bd11, Bd12, Bd13; Bd21, Bd22, Bd23];

B = horzcat(Bu, Bd);
C = [1, 0];
D = [0, 0, 0, 0];

G = ss(A, B, C, D);

% Convertendo os modelos do tempo contínuo para discreto


