%% Exemplo 1fdf: Wood-Berry distillation

s = tf('s');
G11 = exp(-s)*12.8/(16.7*s+1);
G12 = exp(-3*s)*(-18.9)/(21*s+1);
G21 = exp(-7*s)*6.6/(10.9*s+1);
G22 = exp(-3*s)*(-19.4)/(14.4*s+1);

G = [G11 G12; G21 G22];

Kp1 = 1; Ki1 = 0.01; Kd1 = 0.1;
rho1 = [Kp1 Ki1 Kd1];
phi = [1 1/s s/(0.1*s+1)];
C1 = rho1*phi';

Kp2 = 1; Ki2 = 0.01; Kd2 = 0.1;
rho2 = [Kp2 Ki2 Kd2];
phi = [1 1/s s/(0.1*s+1)];
C2 = rho2*phi';

C = [C1 0; 0 C2];
for j = 1:size(G,2)
    Gw(j) = freqresp(G(j,j));
    Gw(j) = squeeze(Gw);
end
