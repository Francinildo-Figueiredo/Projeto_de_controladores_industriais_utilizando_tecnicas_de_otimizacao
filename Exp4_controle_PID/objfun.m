function Jv=objfun(x,s,G)
    Kp=x(1); Ki=x(2); Kd=x(3);
    Fd = Kp + Ki/s + Kd*s/(1+0.01*s);
    Jv=norm(feedback(pade(G)/s,Fd),inf);
    %Ju=norm(feedback(Fd, pade(G)),inf);
return 