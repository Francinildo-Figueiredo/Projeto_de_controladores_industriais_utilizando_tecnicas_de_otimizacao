s=tf('s'); G=1/((1+s)*(1+0.5*s)*(1+0.25*s));
%G=1/(1+s)^3;
MS_max=1.7; % Margem de estabilidade máxima
MT_max=1.3; % Margem de amortecimento máxima
Kinf=15;    % Ganho em alta frequência
min_x=[3 0.4 0.6];      % Limites superiores e inferiores
max_x=[5 0.9 1];        % de Ki, tau e zeta
x=(min_x+max_x)/2;      % Valor inicial = valor médio
options = optimset('Algorithm','active-set');
x=fmincon(@(x) objfun(x,s,G,Kinf),x,[],[],[],[],...
min_x, max_x, @(x)confun(x,G,MS_max,MT_max,Kinf), options);
Ki=x(1); tau=x(2); zeta=x(3); beta=Kinf/(Ki*tau);
K=tf(Ki*[tau^2 2*zeta*tau 1], [tau/beta 1 0]);