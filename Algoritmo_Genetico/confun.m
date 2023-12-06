function [C,Ceq]=confun(x,s,G,MS_max,MT_max, Ju_max, Gm_Max, Pm_Max)
    Kp=x(1); Ki=x(2); Kd=x(3);
    K=Kp + Ki/s + Kd*s/(1+0.01*s); % Controlador PID
    MS=norm(feedback(1,pade(G)*K),inf);                 % Margem de estabilidade
    MT=norm(feedback(pade(G)*K,1),inf);                 %
    Ju=norm(feedback(K, pade(G)),inf);
    [Gm,Pm] = margin(G*K);
%     Jv=norm(feedback(pade(G)/s,K),inf);
    stab=norm(feedback(pade(G)*K,1)); % Norma 2 da função de sensibilidade complementar
    if stab<inf, C=[MS-MS_max; MT-MT_max; Ju-Ju_max; Gm-Gm_Max; Pm-Pm_Max];
    else C=[1; 1; 1; 1; 1]; end
    Ceq=[];
return