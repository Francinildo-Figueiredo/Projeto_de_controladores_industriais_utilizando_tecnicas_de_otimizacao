Greal = tf(2, [1 2 1], 'iodelay', 1);
h = 5;
DeltaT = 0.2;
t = 0:DeltaT:10;
u = h*ones(size(t));
randn('seed', 314);
y = lsim(Greal, u, t) + 0.2*randn(length(t), 1);
[G0, T1, L] = parametrosFOPTD(y, h, DeltaT);
Gid = tf(G0, [T1 1], 'iodelay', L);
yy = lsim(Gid, u, t);
plot(t,y,t,yy);
legend('y', 'yy');