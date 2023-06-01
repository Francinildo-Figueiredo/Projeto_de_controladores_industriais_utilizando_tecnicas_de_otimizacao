dobtidos = load('exp23052023115431.mat');
pv1_aq = dobtidos.pv1(275:667);
mv1_aq = dobtidos.mv1(275:667);
duty_op = dobtidos.mv1(250);
t = 1:1:length(dobtidos.pv1(1:end-1));

[K, T, L] = parametrosFOPTD(pv1_aq - pv1_aq(1), mv1_aq(1) - duty_op, t(end) - t(end-1));
G = tf(K, [T 1], 'iodelay', L);
t = 0:1:length(pv1_aq)-1;
Sis_id = lsim(G, mv1_aq - duty_op, t);

figure(1);
plot(t,pv1_aq, t, Sis_id+pv1_aq(1));

pv1_resf = dobtidos.pv1(668:end-7);
mv1_resf = dobtidos.mv1(668:end-7);


