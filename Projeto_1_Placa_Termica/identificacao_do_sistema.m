dobtidos = load('exp23052023115431.mat');

%% Aquecimento da Malha 1
pv1_aq = dobtidos.pv1(275:667);
mv1_aq = dobtidos.mv1(275:667);
t = 1:1:length(dobtidos.pv1(1:end-1));

% Função de transferência entre a MV1 e PV1 no aquecimento
[K11, T11, L11] = parametrosFOPTD(pv1_aq - pv1_aq(1), mv1_aq(1) - 20, t(end) - t(end-1));
G11 = tf(K11, [T11 1], 'iodelay', L11);
t = 0:1:length(pv1_aq)-1;
Sis_id = lsim(G11, mv1_aq - 20, t);

figure(1);
plot(t,pv1_aq-pv1_aq(1), t, Sis_id);

% Função de transferência entre a MV1 e PV2 no aquecimento
pv2_aq = dobtidos.pv2(275:667);

[K21, T21, L21] = parametrosFOPTD(pv2_aq - pv2_aq(1), mv1_aq(1) - 20, t(end) - t(end-1));
G21 = tf(K21, [T21 1], 'iodelay', L21);
t = 0:1:length(pv2_aq)-1;
Sis_id = lsim(G21, mv1_aq - 20, t);

figure(2);
plot(t,pv2_aq-pv2_aq(1), t, Sis_id);

%% Resfriamento da Malha 1

pv1_resf = dobtidos.pv1(668:end-7);
mv1_resf = dobtidos.mv1(668:end-7);

% Função de transferência entre a MV1 e PV1 no resfriamento
[K11, T11, L11] = parametrosFOPTD(pv1_resf - pv1_resf(1),-10, t(end) - t(end-1));
G11_resf = tf(K11, [T11 1], 'iodelay', L11);
t = 0:1:length(pv1_resf)-1;
Sis_id = lsim(G11_resf, -10*ones(size(t)), t);

figure(3);
plot(t,pv1_resf - pv1_resf(1), t, Sis_id);

% Função de transferência entre a MV1 e PV2 no resfriamento
pv2_resf = dobtidos.pv2(668:end-7);

[K21, T21, L21] = parametrosFOPTD(pv2_resf - pv2_resf(1), -10, t(end) - t(end-1));
G21_resf = tf(K21, [T21 1], 'iodelay', L21);
t = 0:1:length(pv2_resf)-1;
Sis_id = lsim(G21_resf, -10*ones(size(t)), t);

figure(4);
plot(t,pv2_resf - pv2_resf(1), t, Sis_id);

%% Aquecimento da Malha 2

pv1_aq = dobtidos.pv1(275:667);
mv2_aq = dobtidos.mv1(275:667);
t = 1:1:length(dobtidos.pv1(1:end-1));

% Função de transferência entre a MV2 e PV1 no aquecimento
[K12, T12, L12] = parametrosFOPTD(pv1_aq - pv1_aq(1), mv2_aq(1), t(end) - t(end-1));
G12 = tf(K12, [T12 1], 'iodelay', L12);
t = 0:1:length(pv1_aq)-1;
Sis_id = lsim(G12, mv2_aq, t);

figure(5);
plot(t,pv1_aq-pv1_aq(1), t, Sis_id);

% Função de transferência entre a MV2 e PV2 no aquecimento
pv2_aq = dobtidos.pv2(275:667);

[K22, T22, L22] = parametrosFOPTD(pv2_aq - pv2_aq(1), mv2_aq(1), t(end) - t(end-1));
G22 = tf(K22, [T22 1], 'iodelay', L22);
t = 0:1:length(pv2_aq)-1;
Sis_id = lsim(G22, mv2_aq, t);

figure(6);
plot(t,pv2_aq-pv2_aq(1), t, Sis_id);

%% Resfriamento da Malha 2

pv1_resf = dobtidos.pv1(668:end-7);
mv2_resf = dobtidos.mv2(668:end-7);

% Função de transferência entre a MV2 e PV1 no resfriamento
[K12, T12, L12] = parametrosFOPTD(pv1_resf - pv1_resf(1),-10, t(end) - t(end-1));
G12_resf = tf(K12, [T12 1], 'iodelay', L12);
t = 0:1:length(pv1_resf)-1;
Sis_id = lsim(G12_resf, -10*ones(size(t)), t);

figure(7);
plot(t,pv1_resf - pv1_resf(1), t, Sis_id);

% Função de transferência entre a MV2 e PV2 no resfriamento
pv2_resf = dobtidos.pv2(668:end-7);

[K22, T22, L22] = parametrosFOPTD(pv2_resf - pv2_resf(1), -10, t(end) - t(end-1));
G22_resf = tf(K22, [T22 1], 'iodelay', L22);
t = 0:1:length(pv2_resf)-1;
Sis_id = lsim(G22_resf, -10*ones(size(t)), t);

figure(8);
plot(t,pv2_resf - pv2_resf(1), t, Sis_id);