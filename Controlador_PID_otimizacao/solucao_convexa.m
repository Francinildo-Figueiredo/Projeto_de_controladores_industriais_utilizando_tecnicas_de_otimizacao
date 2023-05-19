Kp = 0:0.1:15;
Ki = 0:0.1:15;
Jv = zeros(1, 2*length(Kp));

s = tf('s');
G=1/((1+s)*(1+0.5*s)*(1+0.25*s));

k = 1;
for i = 1:1:length(Kp)
    for j = 1:1:length(Ki)
        Jv(k) = objfun([Kp(i), Ki(j), 0] , s, G);
        k = k + 1;
        disp(k);
    end
end


    