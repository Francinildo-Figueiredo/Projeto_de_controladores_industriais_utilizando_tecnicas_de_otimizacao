function Jv=objfun(x,s,G)

    K1 = x(1) + x(2)/s; + x(3)*s/(1+0.01*s);
    K2 = x(4) + x(5)/s; + x(6)*s/(1+0.01*s);
    K  = [K1, 0; 0, K2];
    Jv=norm(feedback(pade(G)/s,K),inf);
return 