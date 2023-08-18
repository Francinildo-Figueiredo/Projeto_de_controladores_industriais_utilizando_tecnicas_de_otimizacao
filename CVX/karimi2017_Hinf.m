%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projeto do controlador PI centralizado usando formulacao H_inf do artigo:
% A data-driven approach to robust control of multivariable systems by
% convex optimization - equacao 15
%
% Controlador; K = X Y^-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all

% Definicao da Matriz de Funcoes de Transferencia - Sistema TITO - Wood-Berry
num11 = [12.8]; den11 = [16.7 1]; td11 = 1; g11 = tf(num11,den11,'iodelay',td11);
num12 = [-18.9]; den12 = [21 1]; td12 = 3; g12 = tf(num12,den12,'iodelay',td12);
num21 = [6.6]; den21 = [10.9 1]; td21 = 7; g21 = tf(num21,den21,'iodelay',td21);
num22 = [-19.4]; den22 = [14.4 1]; td22 = 3; g22 = tf(num22,den22,'iodelay',td22);
G = [g11 g12; g21 g22];
n = 2; % dimensao
I = eye(n); 

% controlador inicial
% boyd2016 - MIMO PID tuning via iterated LMI restriction
ki = 0.1*I/dcgain(G);
c11 = pid(0,ki(1,1)); c12 = pid(0,ki(1,2));
c21 = pid(0,ki(2,1)); c22 = pid(0,ki(2,2));
% dhanyaram2015 - MIMO PID tuning via iterated LMI restriction
% c11 = pid(0.3140,0.0471); c12 = pid(-0.3058,-0.04587);
% c21 = pid(0.1068,0.01602); c22 = pid(-0.2072,-0.03108);

% Dados projeto
N = 300; %numero de pontos
w = logspace(-3,3,N); %faixa de frequencia
ak = 0.1; %parametro de W1 - ponderacao da sensibilidade
W2 = I; %ponderacao da entrada
bd = 0.2;
Qmax = 0.7381;
Smax = 1.2;

% Resposta em frequencia do processo
Gw = freqresp(G,w);

% resolucao do problema de otimizacao convexa usando cvx
cvx_begin sdp %inicializa o cvx
    variables gamma2 xp(n,n) xi(n,n) %variaveis 
    %gamma2: objetivo
    %xp e xi: ganhos proporcional e integral do controlador
    dual variables S{N} %variavel dual, usada para construir 
    %as restricoes no intervalo de frequencia
    minimize(gamma2) %funcao objetivo
    subject to %restricoes
        for k = 1 : N %construcao das restrições no intervalo de frequencia
            W1 = 1/(Smax)*((1i*w(k)+bd)/(1i*w(k)+0.001*bd)); %ponderacao da sensibilidade
            %W2 = 100/(Qmax)*((1i*w(k)+bd)/(1i*w(k)+100*bd));
            W2 = 0;
            Y = [(1j*w(k)) 0; 0 (1j*w(k))]; %Y e Yc, como está sendo utilizado o controlador PI: Y = Yc = [s 0; 0 s]
            X = [xp(1,1)*1j*w(k)+xi(1,1) xp(1,2)*1j*w(k)+xi(1,2); xp(2,1)*1j*w(k)+xi(2,1) xp(2,2)*1j*w(k)+xi(2,2)]; %X 
            Z = Y+Gw(:,:,k)*X; %P da equacao 15 do artigo
            Xc = [c11.Kp*1j*w(k)+c11.ki c12.Kp*1j*w(k)+c12.ki; c21.Kp*1j*w(k)+c21.ki c22.Kp*1j*w(k)+c22.ki]; %Xc
            Zc = Y+Gw(:,:,k)*Xc; %Pc da equacao 15 do artigo 
            S{k} : [Z'*Zc+Zc'*Z-Zc'*Zc, (W1*Y)',  (W2*X)';...
                   (W1*Y),              gamma2*I,  zeros(2);...
                   (W2*X),              zeros(2), gamma2*I] == hermitian_semidefinite(3*n); %restricao equacao 15 do artigo
        end       
cvx_end %finaliza o cvx

% Controlador projetado
c11 = pid(xp(1,1),xi(1,1)); c12 = pid(xp(1,2),xi(1,2));
c21 = pid(xp(2,1),xi(2,1)); c22 = pid(xp(2,2),xi(2,2));
C = [c11 c12; c21 c22];

% Simulacao da malha fechada obtida
T = (I+G*C)^-1*G*C; % malha fechada
S = (I+G*C)^-1;
Q = C*S;
Msopt = max(max(sigma(S,w)));
Mtopt = max(max(sigma(T,w)));
Mqopt = max(max(sigma(Q,w)));
r = [ones(200,1) zeros(200,1); zeros(200,1) zeros(200,1); zeros(200,1) ones(200,1)]; %sinal referencia
t = 0:1:599;
lsim(T,r,t) %simulacao


