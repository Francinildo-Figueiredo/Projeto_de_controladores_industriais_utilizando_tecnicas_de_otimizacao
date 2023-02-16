%%
%% Exercício 1
%%
% Função de transferência
G = tf([1, 5], [1, 4.5, 7, 4.5, 1]);
G.ioDelay=1; % Aplicando um atraso de 1

% Controlador PID
s = tf('s');
C = 0.25 + 0.04*(1/s) + 0.4*s;

% Modelo de zeros, polos e ganhos
G1 = zpk(G);

% Verificando a estabilidade do sistema
polos = pole(G);
zeros = zero(G);
pzmap(G);

% Obtendo a resposta ao degrau
step(G);
info = stepinfo(G);

% Função de tranferência do controlador PID baseado nos parâmetros Kp, Ki e
% Kd
Kp = 0.25; Ki = 0.04; Kd = 0.4;
C = pid(Kp, Ki, Kd);

% Definindo a função de ganho de malha
L = G*C;

% Análise gáfica da função de transferência G(s)
rlocus(G);
grid on
bode(G);
nyquist(G);

% Obtendo as margens de granho e fase
[gm, pm, wg,wp] = margin(G);

