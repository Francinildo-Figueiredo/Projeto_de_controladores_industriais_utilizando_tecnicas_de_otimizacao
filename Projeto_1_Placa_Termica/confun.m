function [C,Ceq]=confun(x,s,G,MS_max,MT_max,Ju_max)
    K1 = x(1) + x(2)/s; + x(3)*s/(1+0.01*s);
    K2 = x(4) + x(5)/s; + x(6)*s/(1+0.01*s);
    K  = [K1, 0; 0, K2];
    I = eye(2);
    MS=norm(feedback(I,pade(G)*K),inf);             
    MT=norm(feedback(pade(G)*K,I),inf);  
    Ju=norm(feedback(K, pade(G)),inf);
    stab=norm(feedback(pade(G)*K,I)); % Norma 2 da função de sensibilidade complementar
    if stab<inf, C=[MS-MS_max; MT-MT_max; Ju-Ju_max];
    else C=[1; 1; 1]; end
    Ceq=[];
return