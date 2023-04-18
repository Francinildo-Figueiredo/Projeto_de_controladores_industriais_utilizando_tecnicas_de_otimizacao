function [G0, T1, L] = parametrosFOPTD(y,h,DeltaT)
theta = zeros(1, 3);
R = zeros(size(theta,2));
f = zeros(length(theta),1);

N = length(y);
for k = 1:N
    phi = [h*k*DeltaT, -h, -y(k)]';
    R = R + phi*phi';
    A = sum(y(1:k))*DeltaT;
    f = f + phi*A;
end
theta = (R\f)';
G0 = theta(1); T1 = theta(3); L = theta(2)/theta(1);

