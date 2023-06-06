dobtidos = load('exp23052023115431.mat');

%% Aquecimento da Malha 1
pv1_aq = dobtidos.pv1(275:667);
mv1_aq = dobtidos.mv1(275:667);
t = 1:1:length(dobtidos.pv1(1:end-1));

% Função de transferência entre a MV1 e PV1 no aquecimento
[K11_aq, T11_aq, L11_aq] = parametrosFOPTD(pv1_aq - pv1_aq(1), mv1_aq(1) - 20, t(end) - t(end-1));
G11_aq = tf(K11_aq, [T11_aq 1], 'iodelay', L11_aq);
t = 0:1:length(pv1_aq)-1;
Sis_id = lsim(G11_aq, mv1_aq - 20, t);

figure(1);
plot(t,pv1_aq-pv1_aq(1), t, Sis_id);
title('MV1 e PV1 no aquecimento');

% Função de transferência entre a MV1 e PV2 no aquecimento
pv2_aq = dobtidos.pv2(275:667);

[K21_aq, T21_aq, L21_aq] = parametrosFOPTD(pv2_aq - pv2_aq(1), mv1_aq(1) - 20, t(end) - t(end-1));
G21_aq = tf(K21_aq, [T21_aq 1], 'iodelay', L21_aq);
t = 0:1:length(pv2_aq)-1;
Sis_id = lsim(G21_aq, mv1_aq - 20, t);

figure(2);
plot(t,pv2_aq-pv2_aq(1), t, Sis_id);
title('MV1 e PV2 no aquecimento');

%% Resfriamento da Malha 1

pv1_resf = dobtidos.pv1(668:end-7);
mv1_resf = dobtidos.mv1(668:end-7);

% Função de transferência entre a MV1 e PV1 no resfriamento
[K11_resf, T11_resf, L11_resf] = parametrosFOPTD(pv1_resf - pv1_resf(1),-10, t(end) - t(end-1));
G11_resf = tf(K11_resf, [T11_resf 1], 'iodelay', L11_resf);
t = 0:1:length(pv1_resf)-1;
Sis_id = lsim(G11_resf, -10*ones(size(t)), t);

figure(3);
plot(t,pv1_resf - pv1_resf(1), t, Sis_id);
title('MV1 e PV1 no resfriamento');

% Função de transferência entre a MV1 e PV2 no resfriamento
pv2_resf = dobtidos.pv2(668:end-7);

[K21_resf, T21_resf, L21_resf] = parametrosFOPTD(pv2_resf - pv2_resf(1), -10, t(end) - t(end-1));
G21_resf = tf(K21_resf, [T21_resf 1], 'iodelay', L21_resf);
t = 0:1:length(pv2_resf)-1;
Sis_id = lsim(G21_resf, -10*ones(size(t)), t);

figure(4);
plot(t,pv2_resf - pv2_resf(1), t, Sis_id);
title('MV1 e PV2 no resfriamento');

%% Aquecimento da Malha 2

pv1_aq = dobtidos.pv1(275:667);
mv2_aq = dobtidos.mv1(275:667);
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
pv2_aq = dobtidos.pv2(275:667);

[K22_aq, T22_aq, L22_aq] = parametrosFOPTD(pv2_aq - pv2_aq(1), mv2_aq(1), t(end) - t(end-1));
G22_aq = tf(K22_aq, [T22_aq 1], 'iodelay', L22_aq);
t = 0:1:length(pv2_aq)-1;
Sis_id = lsim(G22_aq, mv2_aq, t);

figure(6);
plot(t,pv2_aq-pv2_aq(1), t, Sis_id);
title('MV2 e PV2 no aquecimento');

%% Resfriamento da Malha 2

pv1_resf = dobtidos.pv1(668:end-7);
mv2_resf = dobtidos.mv2(668:end-7);

% Função de transferência entre a MV2 e PV1 no resfriamento
[K12_resf, T12_resf, L12_resf] = parametrosFOPTD(pv1_resf - pv1_resf(1),-10, t(end) - t(end-1));
G12_resf = tf(K12_resf, [T12_resf 1], 'iodelay', L12_resf);
t = 0:1:length(pv1_resf)-1;
Sis_id = lsim(G12_resf, -10*ones(size(t)), t);

figure(7);
plot(t,pv1_resf - pv1_resf(1), t, Sis_id);
title('MV2 e PV1 no resfriamento');

% Função de transferência entre a MV2 e PV2 no resfriamento
pv2_resf = dobtidos.pv2(668:end-7);

[K22_resf, T22_resf, L22_resf] = parametrosFOPTD(pv2_resf - pv2_resf(1), -10, t(end) - t(end-1));
G22_resf = tf(K22_resf, [T22_resf 1], 'iodelay', L22_resf);
t = 0:1:length(pv2_resf)-1;
Sis_id = lsim(G22_resf, -10*ones(size(t)), t);

figure(8);
plot(t,pv2_resf - pv2_resf(1), t, Sis_id);
title('MV2 e PV2 no resfriamento');

%% Modelos médios

pv1 = [pv1_aq, pv1_resf];

% Modelo médio entre a MV1 e PV1
K11_med = (K11_aq + K11_resf)/2;
T11_med = (T11_aq + T11_resf)/2;
L11_med = (L11_aq + L11_resf)/2;

G11_med = tf(K11_med, [T11_med 1], 'iodelay', L11_med);
t = 0:1:length(pv1_aq)-1;
Sis_id_aq = lsim(G11_med, mv1_aq - 20, t);
t = 0:1:length(pv1_resf)-1;
Sis_id_resf = lsim(G11_med, -10*ones(size(t)), t);

Sis_id = [Sis_id_aq', (Sis_id_resf-Sis_id_resf(end))'];
t = 0:1:length(pv1)-1;

figure(9);
plot(t,pv1-pv1(1), t, Sis_id);
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

%% Projetando os controladores PI ulizando o método SIMC, com \tau_c = \theta
%% Com base nos parâmetros dos modelos médios de G11 e G22

% PI para o modelo médio de G11
Kp11 = T11_med/(2*K11_med*L11_med);
Ti11 = min(T11_med, 8*L11_med);

Cpi11 = pid(Kp11, Kp11/Ti11);
H11 =feedback(Cpi11*G11_med,1);

% PI para o modelo médio de G22
Kp22 = T22_med/(2*K22_med*L22_med);
Ti22 = min(T22_med, 8*L22_med);

Cpi22 = pid(Kp22, Kp22/Ti22);
H22 =feedback(Cpi22*G22_med,1);
