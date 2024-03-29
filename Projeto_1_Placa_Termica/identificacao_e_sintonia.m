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
Kp11_SIMC = T11_med/(2*K11_med*L11_med);
Ti11_SIMC = min(T11_med, 8*L11_med);

Cpi11 = pid(Kp11_SIMC, Kp11_SIMC/Ti11_SIMC);
H11_SIMC =feedback(Cpi11*G11_med,1);
step(H11_SIMC);
hold on;
stepinfo(H11_SIMC);
t = 0:1:250;
Sis_id_resf = lsim(H11_SIMC, -1*ones(size(t)), t);
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


%% Otimização do controlador PI da malha 1 através do método SQP

s = tf('s');

% Criterios antes da otimização
MSb=norm(feedback(1,pade(G11_med)*Cpi11),inf)
MTb=norm(feedback(pade(G11_med)*Cpi11,1),inf)
Jvb=norm(feedback(pade(G11_med)/s,Cpi11),inf)
Jub=norm(feedback(Cpi11, pade(G11_med)),inf)

% Critérios de restrição
MS_max=1.5;
MT_max=1.1;
% Jv_max=inf;
Ju_max=100;
Ki11_SIMC = Kp11_SIMC/Ti11_SIMC;
x11_SIMC = [Kp11_SIMC, Ki11_SIMC];
% x = [5 0.1 0.2];
% x = [0,0,0];

% Otimização do controlador por meio da fmincon
lb = [0, 0.085];
ub = [30, 10];
options = optimoptions('fmincon', 'display', 'iter');
x11_SQP=fmincon(@(x) objfun(x,s,G11_med),x11_SIMC,[],[],[],[],...
lb, ub, @(x)confun(x,s,G11_med,MS_max,MT_max,Ju_max), options);

Kp11 = x11_SQP(1); Ki11 = x11_SQP(2);
Kpi11_SQP = Kp11 + Ki11/s;
H11_SQP = feedback(G11_med*Kpi11_SQP,1);

% Critérios após a otimização
Jv11_SQP=norm(feedback(pade(G11_med)/s,Kpi11_SQP),inf)
Ju11_SQP=norm(feedback(Kpi11_SQP, pade(G11_med)),inf)
MS11_SQP=norm(feedback(1,pade(G11_med)*Kpi11_SQP),inf)
MT11_SQP=norm(feedback(pade(G11_med)*Kpi11_SQP,1),inf)
[Gm11_SQP,Pm11_SQP] = margin(G11_med*Kpi11_SQP)
%%
% Otimização do controlador por meio do algoritmo genético
lb = [0, 0.085];
ub = [6, 0.1];
options = optimoptions('ga', 'display', 'iter');
x11_GA = ga(@(x) objfun(x,s,G11_med),2,[],[],[],[],...
lb,ub,@(x)confun(x,s,G11_med,MS_max,MT_max,Ju_max), options);

Kp11 = x11_GA(1); Ki11 = x11_GA(2);
Kpi11_GA = Kp11 + Ki11/s;
H11_GA = feedback(G11_med*Kpi11_GA,1);

% Critérios após a otimização
Jv11_GA=norm(feedback(pade(G11_med)/s,Kpi11_GA),inf)
Ju11_GA=norm(feedback(Kpi11_GA, pade(G11_med)),inf)
MS11_GA=norm(feedback(1,pade(G11_med)*Kpi11_GA),inf)
MT11_GA=norm(feedback(pade(G11_med)*Kpi11_GA,1),inf)

t = 0:1:100;
y11_SIMC = step(H11_SIMC, t);
IAE11_SIMC = sum(abs(1-y11_SIMC))
step(H11_SIMC, t);
hold on;
y11_SQP = step(H11_SQP,t);
IAE11_SQP = sum(abs(1-y11_SQP))
step(H11_SQP,t);
hold on;
y11_GA = step(H11_GA,t);
IAE11_GA = sum(abs(1-y11_GA))
step(H11_GA,t);
legend('SIMC', 'SQP', 'GA', 'Location','southeast');
title('Resposta ao degrau para a malha 1');
Ti11_SQP = x11_SQP(1)/x11_SQP(2);
Ti11_GA  = x11_GA(1)/x11_GA(2);

%% Resultados obtidos no experimento para cada método de sintonia

dobt_SIMC = load('exp03022024125543.mat');
dobt_SQP  = load('exp03022024143419.mat');
dobt_GA   = load('exp03022024152205.mat');

% Controle realizado pelo método de otimização GA
subplot(2,1,1);
plot(dobt_GA.pv1(154:576), 'LineWidth', 2);hold on;
plot(dobt_GA.pv2(154:576), '-.', 'LineWidth', 1.5);hold on;
plot(dobt_GA.sp1(154:576), '--', 'LineWidth', 1.5);
xlim([0, 430]);
ylim([30, 50]);
xlabel('t(s)', 'fontsize', 12);
ylabel('T(°C)', 'fontsize', 12);
legend('PV1', 'PV2', 'SP1','location','southeast', 'fontsize', 10);

subplot(2,1,2);
plot(dobt_GA.mv1(154:576), 'LineWidth', 2);hold on;
plot(dobt_GA.mv2(154:576), '-.', 'LineWidth', 1.5);
xlim([0, 430]);
ylim([0, 50]);
xlabel('t(s)', 'fontsize', 12);
ylabel('MV(%)', 'fontsize', 12);
legend('MV1', 'MV2','location','best', 'fontsize', 10);

% Métricas
t = 0:1:length(dobt_GA.sp1(154:548))-1;
IAE_GA = sum(abs(dobt_GA.sp1(154:548)-dobt_GA.pv1(154:548)))
ISE_GA = sum((dobt_GA.sp1(154:548)-dobt_GA.pv1(154:548)).^2)
ITAE_GA = sum(t.*abs(dobt_GA.sp1(154:548)-dobt_GA.pv1(154:548)))

%% Controle realizado pelo método de otimização SQP
subplot(2,1,1);
plot(dobt_SQP.pv1(132:526), 'LineWidth', 2);hold on;
plot(dobt_SQP.pv2(132:526), '-.', 'LineWidth', 1.5);hold on;
plot(dobt_SQP.sp1(132:526), '--', 'LineWidth', 1.5);
% xlim([0, 430]);
ylim([30, 50]);
xlabel('t(s)', 'fontsize', 12);
ylabel('T(°C)', 'fontsize', 12);
legend('PV1', 'PV2', 'SP1','location','southeast', 'fontsize', 10);

subplot(2,1,2);
plot(dobt_SQP.mv1(132:526), 'LineWidth', 2);hold on;
plot(dobt_SQP.mv2(132:526), '-.', 'LineWidth', 1.5);
% xlim([0, 430]);
ylim([0, 50]);
xlabel('t(s)', 'fontsize', 12);
ylabel('MV(%)', 'fontsize', 12);
legend('MV1', 'MV2','location','best', 'fontsize', 10);

% Métricas
t = 0:1:length(dobt_SQP.sp1(132:526))-1;
IAE_SQP = sum(abs(dobt_SQP.sp1(132:526)-dobt_SQP.pv1(132:526)))
ISE_SQP = sum((dobt_SQP.sp1(132:526)-dobt_SQP.pv1(132:526)).^2)
ITAE_SQP = sum(t.*abs(dobt_SQP.sp1(132:526)-dobt_SQP.pv1(132:526)))

