Kp = 1e-2:1e-1:5;
Ki = 1e-2:1e-1:5;
Jv = zeros(length(Kp), length(Ki));

s = tf('s');
%G=1/((1+s)*(1+0.5*s)*(1+0.25*s));
G = 1/(1+s)^3;

MS_max=1.8;
MT_max=1.3;
%Jv_max=0.66;
Ju_max=20;

for i = 1:1:length(Kp)
    for j = 1:1:length(Ki)
        x0 = [Kp(i), Ki(j), 0];
        K=Kp(i) + Ki(j)/s; %+ Kd*s/(1+0.01*s); % Controlador PID
        stab=norm(feedback(pade(G)*K,1)); % Norma 2 da função de sensibilidade complementar
        if stab<inf 
            options = optimset('Algorithm','active-set');
            x = fmincon(@(x0) objfun(x0,s,G),x0,[],[],[],[],...
            [], [], @(x0)confun(x0,s,G,MS_max,MT_max,Ju_max), options);
            Jv(i,j) = objfun([x(1), x(2), x(3)] , s, G);
        else
            Jv(i,j) = NaN; 
        end
        disp(i);
    end
end

%[X,Y] = meshgrid(Kp,Ki);
%Jv = Jv/max(max(Jv));
%Jv = mag2db(Jv);
surf(Kp, Ki, Jv);
%zlim([0 1]);
xlabel('Kp');
ylabel('Ki');
zlabel('Jv');
grid on;

    