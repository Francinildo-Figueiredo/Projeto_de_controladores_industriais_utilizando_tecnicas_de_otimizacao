%% Experimento 2: Linearização e Discretização de Sistemas (STH)
%%
% Definição dos parâmetros da simulação
V      = 283.168e-3; % Volume do tanque [m^3]
Vj     = 28.319e-3;  % Volume da jaqueta [m^3]
rhocp  = 2.165e3;    % Densidade*Cap. calorífica do fluido no tanque [W/(ºC*m^3)]
rhocpj = 2.165e3;    % Densidade*Cap. calorífica do fluido na jaqueta [W/(ºC*m^3)]
UA     = 3.065;      % Coef. de transf. de calor*Área da superfície [W/(ºC*min)]

% Ponto de operação
Fs   = 4.720e-4; % Vazão no tanque [m^3/s]
Fjs  = 7.079e-4; % Vazão na jaqueta [m^3/s]
Tis  = 10;       % Temperatura entrando no tanque [ºC]
Tjis = 93.33;    % Temperatura entrando na jaqueta [ºC]
Tjs  = 65.56;    % Temperatura saindo da jaqueta [ºC]
Ts   = 51.67;    % Temperatura saindo do tanque [ºC]

%% Equações
%%
% Equação diferencial que descreve ocomportamento dinâmico da temperatura:
dTdt  = Fs/V*(Tis - Ts) + (UA*(Tjs - Ts))/(V*rhocp);

% Equação diferencial para a temperatura do fluido na jaqueta:
dTjdt = Fjs/Vj*(Tjis - Tjs) - (UA*(Tjs - Ts))/(Vj*rhocpj);

%% Obtendo o modelo linearizado em espaço de estados
%%
A11 = -Fs/V - UA/(V*rhocp);
A12 = UA/(V*rhocp);
A21 = UA/(Vj*rhocpj);
A22 = -Fjs/Vj - UA/(Vj*rhocpj);
A = [A11, A12; A21, A22];

Bu1 = 0; Bu2 = (Tjis - Tjs)/Vj;
Bu = [Bu1; Bu2];

Bd11 = (Tis - Ts)/V; Bd12 = Fs/V; Bd13 = 0;
Bd21 = 0; Bd22 = 0; Bd23 = Fjs/Vj;
Bd = [Bd11, Bd12, Bd13; Bd21, Bd22, Bd23];

B = horzcat(Bu, Bd);
C = [1, 0];
D = [0, 0, 0, 0];

Gss = ss(A, B, C, D);

%% Obtendo o modelo no domínio de Laplace
%%
s   = tf('s');
I   = eye(2,2); 
G   = C*(s*I - A)^(-1)*Bu;
Gd1 = C*(s*I - A)^(-1)*B(:,1);
Gd2 = C*(s*I - A)^(-1)*B(:,2);
Gd3 = C*(s*I - A)^(-1)*B(:,3);

%% Comparando as simulações do modelo linear e não linear, para uma pequena
%% variação na entrada
%%
figure(1);
out1 = sim('STH_NL', 'ReturnWorkspaceOutputs', 'on');
plot(out1.T_NL);
hold on;
out2 = sim('STH_Linear', 'ReturnWorkspaceOutputs', 'on');
plot(out2.T_Linear, 'k');
hold on;
out3 = sim('STH_LinearTF', 'ReturnWorkspaceOutputs', 'on');
plot(out3.T_LinearTF, '--');
legend('Não linear', 'Linear', 'Linear TF', 'Location', 'Best');
title('Comparando as saídas das três simulações realizadas');
xlabel('t (s)');
ylabel('T (ºC)');
grid on;

%% Convertendo os modelos do tempo contínuo para discreto Ta = 50s, 150s, 300s
%% 
Gdisc1 = c2d(Gss, 50, 'zoh');
Gdisc2 = c2d(Gss, 150, 'zoh');
Gdisc3 = c2d(Gss, 300, 'zoh');

figure(2);
out1 = sim('STH_Discreto', 'ReturnWorkspaceOutputs', 'on');
plot(out1.T_Discreto);
legend('Contínuo', 'Ta = 50 s', 'Ta = 150 s', 'Ta = 300 s', 'Location', 'Best');
title('Saídas do sistema contínuo e discretizado');
xlabel('t (s)');
ylabel('T (ºC)');
grid on;
