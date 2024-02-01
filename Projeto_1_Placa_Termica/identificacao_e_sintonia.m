%% Dados obtidos
dobtidos = load('exp18092023174722.mat');
plot(dobtidos.pv2(1:end-1));
title('Medições de temperatura obtidas na PV2');
xlim([0, 1900]);
xlabel('t(s)');
ylabel('T(°C)');

%% Aquecimento da Malha 1
pv1_aq = dobtidos.pv1(173:434);
mv1_aq = dobtidos.mv1(173:434);
t = 0:1:length(mv1_aq)-1;

% Função de transferência entre a MV1 e PV1 no aquecimento
[K11_aq, T11_aq, L11_aq] = parametrosFOPTD(pv1_aq - pv1_aq(1), mv1_aq(1) - 20, t(end) - t(end-1));
G11_aq = tf(K11_aq, [T11_aq 1], 'iodelay', L11_aq);
t = 0:1:length(pv1_aq)-1;
Sis_id = lsim(G11_aq, mv1_aq - 20, t);

figure(1);
plot(t,pv1_aq-pv1_aq(1), t, Sis_id);
title('MV1 e PV1 no aquecimento');

% Função de transferência entre a MV1 e PV2 no aquecimento
pv2_aq = dobtidos.pv2(173:434);

[K21_aq, T21_aq, L21_aq] = parametrosFOPTD(pv2_aq - pv2_aq(1), mv1_aq(1) - 20, t(end) - t(end-1));
G21_aq = tf(K21_aq, [T21_aq 1], 'iodelay', L21_aq);
t = 0:1:length(pv2_aq)-1;
Sis_id = lsim(G21_aq, mv1_aq - 20, t);

figure(2);
plot(t,pv2_aq-pv2_aq(1), t, Sis_id);
title('MV1 e PV2 no aquecimento');

%% Resfriamento da Malha 1

pv1_resf = dobtidos.pv1(434:875);
mv1_resf = dobtidos.mv1(434:875);

% Função de transferência entre a MV1 e PV1 no resfriamento
[K11_resf, T11_resf, L11_resf] = parametrosFOPTD(pv1_resf - pv1_resf(1),-10, t(end) - t(end-1));
G11_resf = tf(K11_resf, [T11_resf 1], 'iodelay', L11_resf);
t = 0:1:length(pv1_resf)-1;
Sis_id = lsim(G11_resf, -10*ones(size(t)), t);

figure(3);
plot(t,pv1_resf - pv1_resf(1), t, Sis_id);
title('MV1 e PV1 no resfriamento');

% Função de transferência entre a MV1 e PV2 no resfriamento
pv2_resf = dobtidos.pv2(434:875);

[K21_resf, T21_resf, L21_resf] = parametrosFOPTD(pv2_resf - pv2_resf(1), -10, t(end) - t(end-1));
G21_resf = tf(K21_resf, [T21_resf 1], 'iodelay', L21_resf);
t = 0:1:length(pv2_resf)-1;
Sis_id = lsim(G21_resf, -10*ones(size(t)), t);

figure(4);
plot(t,pv2_resf - pv2_resf(1), t, Sis_id);
title('MV1 e PV2 no resfriamento');

%% Aquecimento da Malha 2

pv1_aq = dobtidos.pv1(1445:1658);
mv2_aq = dobtidos.mv2(1445:1658);
t = 1:1:length(dobtidos.pv1(1:end-1));

% Função de transferência entre a MV2 e PV1 no aquecimento
[K12_aq, T12_aq, L12_aq] = parametrosFOPTD(pv1_aq - pv1_aq(1), mv2_aq(1), t(end) - t(end-1));
G12_aq = tf(K12_aq, [T12_aq 1], 'iodelay', L12_aq);
t = 0:1:length(pv1_aq)-1;
Sis_id = lsim(G12_aq, mv2_aq, t);

figure(5);
plot(t,pv1_aq-pv1_aq(1), t, Sis_id);
title('MV2 e PV1 no aquecimento');

% Função de transferência entre a MV2 e PV2 no aquecimento
pv2_aq = dobtidos.pv2(1445:1658);

[K22_aq, T22_aq, L22_aq] = parametrosFOPTD(pv2_aq - pv2_aq(1), mv2_aq(1), t(end) - t(end-1));
G22_aq = tf(K22_aq, [T22_aq 1], 'iodelay', L22_aq);
t = 0:1:length(pv2_aq)-1;
Sis_id = lsim(G22_aq, mv2_aq, t);

figure(6);
plot(t,pv2_aq-pv2_aq(1), t, Sis_id);
title('MV2 e PV2 no aquecimento');

%% Resfriamento da Malha 2

pv1_resf = dobtidos.pv1(1659:end);
mv2_resf = dobtidos.mv2(1659:end);

% Função de transferência entre a MV2 e PV1 no resfriamento
[K12_resf, T12_resf, L12_resf] = parametrosFOPTD(pv1_resf - pv1_resf(1),-10, t(end) - t(end-1));
G12_resf = tf(K12_resf, [T12_resf 1], 'iodelay', L12_resf);
t = 0:1:length(pv1_resf)-1;
Sis_id = lsim(G12_resf, -10*ones(size(t)), t);

figure(7);
plot(t,pv1_resf - pv1_resf(1), t, Sis_id);
title('MV2 e PV1 no resfriamento');

% Função de transferência entre a MV2 e PV2 no resfriamento
pv2_resf = dobtidos.pv2(1659:end);

[K22_resf, T22_resf, L22_resf] = parametrosFOPTD(pv2_resf - pv2_resf(1), -10, t(end) - t(end-1));
G22_resf = tf(K22_resf, [T22_resf 1], 'iodelay', L22_resf);
t = 0:1:length(pv2_resf)-1;
Sis_id = lsim(G22_resf, -10*ones(size(t)), t);

figure(8);
plot(t,pv2_resf - pv2_resf(1), t, Sis_id);
title('MV2 e PV2 no resfriamento');

%% Modelos médios

pv1 = [pv1_aq, pv1_resf + (pv1_aq(end)-pv1_resf(1))];

% Modelo médio entre a MV1 e PV1
K11_med = (K11_aq + K11_resf)/2;
T11_med = (T11_aq + T11_resf)/2;
L11_med = (L11_aq + L11_resf)/2;

G11_med = tf(K11_med, [T11_med 1], 'iodelay', L11_med);
t = 0:1:length(mv1_aq)-1;
Sis_id_aq = lsim(G11_med, mv1_aq - 20, t);
t = 0:1:length(pv1_resf)-1;
Sis_id_resf = lsim(G11_med, -10*ones(size(t)), t);

Sis_id = [Sis_id_aq', (Sis_id_resf-Sis_id_resf(end))'];
t = 0:1:length(pv1)-1;

figure(9);
plot(t(end-441:end),pv1(end-441:end)-pv1(1), t(end-441:end), Sis_id(end-441:end));
title('MV1 e PV1 modelo médio');

% Modelo médio entre a MV1 e PV2
K21_med = (K21_aq + K21_resf)/2;
T21_med = (T21_aq + T21_resf)/2;
L21_med = (L21_aq + L21_resf)/2;

G21_med = tf(K21_med, [T21_med 1], 'iodelay', L21_med);

% Modelo médio entre a MV2 e PV1
K12_med = (K12_aq + K12_resf)/2;
T12_med = (T12_aq + T12_resf)/2;
L12_med = (L12_aq + L12_resf)/2;

G12_med = tf(K12_med, [T12_med 1], 'iodelay', L12_med);

% Modelo médio entre a MV2 e PV2
K22_med = (K22_aq + K22_resf)/2;
T22_med = (T22_aq + T22_resf)/2;
L22_med = (L22_aq + L22_resf)/2;

G22_med = tf(K22_med, [T22_med 1], 'iodelay', L22_med);

% figure(9);
% plot(t,pv2-pv2(1), t, Sis_id);
% title('MV2 e PV2 modelo médio');

%% Projetando os controladores PI ulizando o método SIMC, com \tau_c = \theta
%% Com base nos parâmetros dos modelos médios de G11 e G22

% PI para o modelo médio de G11
Kp11 = T11_med/(2*K11_med*L11_med);
Ti11 = min(T11_med, 8*L11_med);

Cpi11 = pid(Kp11, Kp11/Ti11);
H11 =feedback(Cpi11*G11_med,1);
step(H11);
hold on;
stepinfo(H11);
t = 0:1:250;
Sis_id_resf = lsim(H11, -1*ones(size(t)), t);
plot(t+250, Sis_id_resf);

% PI para o modelo médio de G22
Kp22 = T22_med/(2*K22_med*L22_med);
Ti22 = min(T22_med, 8*L22_med);

Cpi22 = pid(Kp22, Kp22/Ti22);
H22 =feedback(Cpi22*G22_med,1);
step(H22);
hold on;
stepinfo(H22);

% Controle em malha fechada da PV1
load('exp18092023174722.mat')
plot(pv1(2:end)-pv1(2)); hold on; plot(sp1(2:end)-sp1(2));
xlim([0, 405]);
title('Controle em malha fechada da variável de processo 1');
legend('Temperatura de saída','Temperatura de referência', 'Location', 'best')
xlabel('t(s)');
ylabel('T(°C)');

% Controle em malha fechada da PV2
load('exp18092023174722.mat')
plot(pv2(440:end)-pv2(440)); hold on; plot(sp2(440:end)-sp2(440));
xlim([0, 518]);
title('Controle em malha fechada da variável de processo 2');
legend('Temperatura de saída','Temperatura de referência', 'Location', 'best')
xlabel('t(s)');
ylabel('T(°C)');


%% Otimização dos controlador PI da malha 1

s = tf('s');

% Criterios antes da otimização
MSb=norm(feedback(1,pade(G11_med)*Cpi11),inf)
MTb=norm(feedback(pade(G11_med)*Cpi11,1),inf)
Jvb=norm(feedback(pade(G11_med)/s,Cpi11),inf)
Jub=norm(feedback(Cpi11, pade(G11_med)),inf)

% Critérios de restrição
MS_max=1.5;
MT_max=1;
% Jv_max=inf;
Ju_max=100;
Ki11 = Kp11/Ti11;
x = [Kp11, Kp11/Ti11, 0];
% x = [5 0.1 0.2];
% x = [0,0,0];

% Otimização dos parâmetro do controlador PI
options = optimoptions('fmincon', 'display', 'iter');
x=fmincon(@(x) objfun(x,s,G11_med),x,[],[],[],[],...
[], [], @(x)confun(x,s,G11_med,MS_max,MT_max,Ju_max), options);

Kp = x(1); Ki = x(2); Kd = x(3);
Kpid1 = Kp + Ki/s + Kd*s/(1+0.01*s);
H = feedback(G11_med*Kpid1,1);

% Critérios após a otimização
MS=norm(feedback(1,pade(G11_med)*Kpid1),inf)
MT=norm(feedback(pade(G11_med)*Kpid1,1),inf)
Jv=norm(feedback(pade(G11_med)/s,Kpid1),inf)
Ju=norm(feedback(Kpid1, pade(G11_med)),inf)

t = 0:1:100;
y11 = step(H11, t);
e_abs11 = abs(1-y11);
IAE11 = sum(e_abs11)
step(H11, t);
title('Resposta ao degrau para a malha 1');
hold on;
y = step(H,t);
e_abs = abs(1-y);
IAE = sum(e_abs)
step(H,t);
legend('SIMC', 'fmincon', 'Location','southeast');
Kp
Ti = Kp/Ki
Td = Kp/Kd

%% Otimização dos controlador PI da malha 2

s = tf('s');

% Criterios antes da otimização
MSb=norm(feedback(1,pade(G22_med*Cpi22)),inf)
MTb=norm(feedback(pade(G22_med)*Cpi22,1),inf)
Jvb=norm(feedback(pade(G22_med)/s,Cpi22),inf)
Jub=norm(feedback(Cpi22, pade(G22_med)),inf)

% Critérios de restrição
MS_max=1.6;
MT_max=1;
% Jv_max=inf;
Ju_max=100;
Ki22 = Kp22/Ti22;
% x = [Kp22, Ki22, 0];
x = [4, 0.1, 0.2];

% Otimização dos parâmetro do controlador PI
options = optimset('Algorithm','active-set');
x=fmincon(@(x) objfun_J(x),x,[],[],[],[],...
[], [], @(x)confun(x,s,G22_med,MS_max,MT_max,Ju_max), options);

Kp = x(1); Ki = x(2); Kd = x(3);
Kpid2 = Kp + Ki/s + Kd*s/(1+0.01*s);
H = feedback(G22_med*Kpid2,1);

% Critérios após a otimização
MS=norm(feedback(1,pade(G22_med)*Kpid2),inf)
MT=norm(feedback(pade(G22_med)*Kpid2,1),inf)
Jv=norm(feedback(pade(G22_med)/s,Kpid2),inf)
Ju=norm(feedback(Kpid2, pade(G22_med)),inf)

t = 0:1:150;
y22 = step(H22, t);
e_abs22 = abs(1-y22);
IAE22 = sum(e_abs22)
step(H22, t);
title('Resposta ao degrau para a malha 2');
hold on;
y = step(H,t);
e_abs = abs(1-y);
IAE = sum(e_abs)
step(H,t);
legend('SIMC', 'Otimizado', 'Location','southeast');
Kp
Ti = Kp/Ki
Td = Kp/Kd

%% Otimização com controle multivariável

G = [G11_med, G12_med; G21_med, G22_med];
C = [Cpi11, 0; 0, Cpi22];
I = eye(2);

% Criterios antes da otimização
MSb=norm(feedback(I,pade(G)*C),inf)
MTb=norm(feedback(pade(G)*C,I),inf)
Jvb=norm(feedback(pade(G)/s,C),inf)
Jub=norm(feedback(C, pade(G)),inf)

% Critérios de restrição
MS_max=1.6;
MT_max=1.3;
Ju_max=100;
% x1 = [Kp11, Ki11, 0, Kp22, Ki22, 0];
x1 = [3,1,0.2,5,1,0.3];

% Otimização do controlador por meio da fmincon
lb = [0,0,0,0,0,0];
ub = [10,10,5,10,10,5];
options = optimoptions('fmincon', 'display', 'iter');
x2=fmincon(@(x) objfun(x,s,G),x1,[],[],[],[],...
lb, ub, @(x)confun(x,s,G,MS_max,MT_max,Ju_max), options);

Kpid1 = x2(1) + x2(2)/s; + x2(3)*s/(1+0.01*s);
Kpid2 = x2(4) + x2(5)/s; + x2(6)*s/(1+0.01*s);
Kpid  = [Kpid1, 0; 0, Kpid2];

% Critérios após a otimização
MS=norm(feedback(I,pade(G)*Kpid),inf)
MT=norm(feedback(pade(G)*Kpid,I),inf)
Jv=norm(feedback(pade(G)/s,Kpid),inf)
Ju=norm(feedback(Kpid, pade(G)),inf)
