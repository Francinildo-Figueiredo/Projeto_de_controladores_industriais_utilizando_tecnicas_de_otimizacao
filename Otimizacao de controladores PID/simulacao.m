%% Simulação dos sistemas de controle otimizados

%% Exemplo 1

s=tf('s'); 
G=1/(1+s)^3;                     % G1

% Obtendo as aproximaçãoes de Half-rule
T = 0;
taui = [1, 1, 1];
theta0 = 0;
K = 1;
tau1 = taui(1);
tau2 = taui(2) + taui(3)/2;
theta = theta0 + taui(3)/2 + sum(taui(4:end)) + sum(T);
G2a = K*exp(-theta*s)/((tau1*s + 1)*(tau2*s + 1));

% Calculando os parâmetros do controlador PID de acordo com o SIMC
tau_c = theta;
Kp = (1/K)*(tau1/(tau_c+theta));
Ti = min(tau1, 4*(tau_c + theta));
Td = tau2;
Ki = Kp/Ti;
Kd = Kp*Td;

out1 = sim('rejeicao_a_disturbios', 'ReturnWorkspaceOutputs', 'on');
y_SIMC = out1.simout.Data(:,1);
disturb = out1.simout.Data(:,2);
t = out1.tout;

% Controlador obtido por meio do SQP

Kp = 2.0159; Ki = 0.7581; Kd = 0.9416;

out2 = sim('rejeicao_a_disturbios', 'ReturnWorkspaceOutputs', 'on');
y_SQP = out2.simout.Data(:,1);

% Controlador obtido por meio do GA

Kp = 2.0198; Ki = 0.7784; Kd = 0.9798;

out3 = sim('rejeicao_a_disturbios', 'ReturnWorkspaceOutputs', 'on');
y_GA = out3.simout.Data(:,1);

IAE_SIMC = sum(abs(1-y_SIMC))
IAE_SQP = sum(abs(1-y_SQP))
IAE_GA = sum(abs(1-y_GA))
ISE_SIMC = sum((1-y_SIMC).^2)
ISE_SQP = sum((1-y_SQP).^2)
ISE_GA = sum((1-y_GA).^2)
ITAE_SIMC = sum(t.*abs(1-y_SIMC))
ITAE_SQP = sum(t.*abs(1-y_SQP))
ITAE_GA = sum(t.*abs(1-y_GA))

plot(t,y_SIMC,t,y_SQP,'r',t,y_GA,'k',t,disturb,'linewidth',1);
ylim([0 1.4]);
yline(1,'--','color',[0.8500 0.3250 0.0980],'linewidth',1.5);
legend('SIMC','SQP','GA','Disturbance','Setpoint','Location','best','FontSize',10);
ylabel('Amplitude');
xlabel('Time(s)');

%% Exemplo 2

s=tf('s'); 
G=exp(-0.3*s)/((1+s)*(1+0.5*s));   % G2

% Obtendo as aproximaçãoes de Half-rule
T = 0;
taui = [1, 0.5, 0];
theta0 = 0.3;
K = 1;
tau1 = taui(1);
tau2 = taui(2) + taui(3)/2;
theta = theta0 + taui(3)/2 + sum(taui(4:end)) + sum(T);
G2a = K*exp(-theta*s)/((tau1*s + 1)*(tau2*s + 1));

% Calculando os parâmetros do controlador PID de acordo com o SIMC
tau_c = theta;
Kp = (1/K)*(tau1/(tau_c+theta));
Ti = min(tau1, 4*(tau_c + theta));
Td = tau2;
Ki = Kp/Ti;
Kd = Kp*Td;

out1 = sim('rejeicao_a_disturbios', 'ReturnWorkspaceOutputs', 'on');
y_SIMC = out1.simout.Data(:,1);
disturb = out1.simout.Data(:,2);
t = out1.tout;

% Controlador obtido por meio do SQP

Kp = 2.3247; Ki = 1.6661; Kd = 0.8332;

out2 = sim('rejeicao_a_disturbios', 'ReturnWorkspaceOutputs', 'on');
y_SQP = out2.simout.Data(:,1);

% Controlador obtido por meio do GA

Kp = 2.3306; Ki = 1.6733; Kd = 0.8298;

out3 = sim('rejeicao_a_disturbios', 'ReturnWorkspaceOutputs', 'on');
y_GA = out3.simout.Data(:,1);

IAE_SIMC = sum(abs(1-y_SIMC))
IAE_SQP = sum(abs(1-y_SQP))
IAE_GA = sum(abs(1-y_GA))
ISE_SIMC = sum((1-y_SIMC).^2)
ISE_SQP = sum((1-y_SQP).^2)
ISE_GA = sum((1-y_GA).^2)
ITAE_SIMC = sum(t.*abs(1-y_SIMC))
ITAE_SQP = sum(t.*abs(1-y_SQP))
ITAE_GA = sum(t.*abs(1-y_GA))

plot(t,y_SIMC,t,y_SQP,'r',t,y_GA,'k',t,disturb,'linewidth',1);
ylim([0 1.4]);
yline(1,'--','color',[0.8500 0.3250 0.0980],'linewidth',1.5);
legend('SIMC','SQP','GA','Disturbance','Setpoint','Location','best','FontSize',10);
ylabel('Amplitude');
xlabel('Time(s)');

