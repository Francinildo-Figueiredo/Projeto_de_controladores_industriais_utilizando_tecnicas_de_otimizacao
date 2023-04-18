%% Experimento 2: Identificação de Sistemas (CSRT)
%%

% Ponto de operação
Fs   = 1;     % Vazão volumétrica (m^3/h)
Cafs = 10;    % Concentração de A no fluxo de alimentação ((kg*mol)/m^3)
Tfs  = 298;   % Temperatura de alimentação (K)
Tjs  = 298;   % Temperatura da capa (K)
Ts   = 311.2; % Temperatura do reator (K)
Cas  = 8.564; % Concentração de A no reator ((kg*mol)/m^3)

% Identificação da função de transferência G11(s): saída Ca, entrada F
dFj = 0.01; dCaf = 0; dTf = 0; dTj = 0;
out1 = sim('CSRT_NL', 'ReturnWorkspaceOutputs', 'on');
t = out1.get('tout');
Ca_F = out1.get('CaNL').Data;
Fstep = out1.get('Fstep').Data;
[K11, T11, L11] = parametrosFOPTD(Ca_F - Ca_F(1), Fstep(1) - Fs, t(end) - t(end-1));
G11 = tf(K11, [T11 1], 'iodelay', L11);
Ca_Fsim = lsim(G11, Fstep - Fs, t);
emqG11 = mean((Ca_F - (Ca_Fsim+Ca_F(1))).^2);
disp(['Erro médio ao quadrado de G11: ' num2str(emqG11)]);
figure(1);
plot(t,Ca_F, '--', t,Ca_Fsim+Ca_F(1), '-');
title('Teste do degrau e simulação de G11 identificado para dF = 0.01 m^3/h');
grid on;
legend('Modelo não linear', 'Modelo G11 identificado', 'Location', 'Best');
xlabel('t (h)');
ylabel('Ca (kgmol/m^3)');

% Identificação da função de transferência G12(s): saída Ca, perturbação Caf
dFj = 0; dCaf = 0.007; dTf = 0; dTj = 0;
out1 = sim('CSRT_NL', 'ReturnWorkspaceOutputs', 'on');
t = out1.get('tout');
Ca_Caf = out1.get('CaNL').Data;
Cafstep = out1.get('Cafstep').Data;
[K12, T12, L12] = parametrosFOPTD(Ca_Caf - Ca_Caf(1), Cafstep(1) - Cafs, t(end) - t(end-1));
G12 = tf(K12, [T12 1], 'iodelay', L12);
Ca_Cafsim = lsim(G12, Cafstep - Cafs, t);
emqG12 = mean((Ca_Caf - (Ca_Cafsim+Ca_Caf(1))).^2);
disp(['Erro médio ao quadrado de G12: ' num2str(emqG12)]);
figure(2);
plot(t,Ca_Caf, '--', t,Ca_Cafsim+Ca_Caf(1), '-');
title('Teste do degrau e simulação de G12 identificado para dCaf = 0.007 kgmol/m^3');
grid on;
legend('Modelo não linear', 'Modelo G12 identificado', 'Location', 'Best');
xlabel('t (h)');
ylabel('Ca (kgmol/m^3)');

% Identificação da função de transferência G13(s): saída Ca, perturbação Tf
dFj = 0; dCaf = 0; dTf = 1; dTj = 0;
out1 = sim('CSRT_NL', 'ReturnWorkspaceOutputs', 'on');
t = out1.get('tout');
Ca_Tf = out1.get('CaNL').Data;
Tfstep = out1.get('Tfstep').Data;
[K13, T13, L13] = parametrosFOPTD(Ca_Tf - Ca_Tf(1), Tfstep(1) - Tfs, t(end) - t(end-1));
G13 = tf(K13, [T13 1], 'iodelay', L13);
Ca_Tfsim = lsim(G13, Tfstep - Tfs, t);
emqG13 = mean((Ca_Tf - (Ca_Tfsim+Ca_Tf(1))).^2);
disp(['Erro médio ao quadrado de G13: ' num2str(emqG13)]);
figure(3);
plot(t,Ca_Tf, '--', t,Ca_Tfsim+Ca_Tf(1), '-');
title('Teste do degrau e simulação de G13 identificado para dTf = 1 K');
grid on;
legend('Modelo não linear', 'Modelo G13 identificado', 'Location', 'Best');
xlabel('t (h)');
ylabel('Ca (kgmol/m^3)');

% Identificação da função de transferência G14(s): saída Ca, perturbação Tj
dFj = 0; dCaf = 0; dTf = 0; dTj = 1;
out1 = sim('CSRT_NL', 'ReturnWorkspaceOutputs', 'on');
t = out1.get('tout');
Ca_Tj = out1.get('CaNL').Data;
Tjstep = out1.get('Tjstep').Data;
[K14, T14, L14] = parametrosFOPTD(Ca_Tj - Ca_Tj(1), Tjstep(1) - Tjs, t(end) - t(end-1));
G14 = tf(K14, [T14 1], 'iodelay', L14);
Ca_Tjsim = lsim(G14, Tjstep - Tjs, t);
emqG14 = mean((Ca_Tj - (Ca_Tjsim+Ca_Tj(1))).^2);
disp(['Erro médio ao quadrado de G14: ' num2str(emqG14)]);
figure(4);
plot(t,Ca_Tj, '--', t,Ca_Tjsim+Ca_Tj(1), '-');
title('Teste do degrau e simulação de G14 identificado para dTj = 1 K');
grid on;
legend('Modelo não linear', 'Modelo G14 identificado', 'Location', 'Best');
xlabel('t (h)');
ylabel('Ca (kgmol/m^3)');