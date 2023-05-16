wc = 0.1;
beta = 4;
zeta = 1;
x = 0:0.1:4;
K_fase = @(x) -90 + acos((1-x^2)/sqrt((1-x^2)^2+(2*zeta*x)^2))...
         - atan(x/beta);
plot(x, K_fase(x))
     