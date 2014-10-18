function y = beam_diff2(t,y,K,omega)
n = length(y)/2;
X = y(1:n);
V = y(n+1:end);% this is dx/dt

% move to a system of 1-D
% M*x_tt+K*x = F
% y = x_t
% y_t = x_tt
% new sys:
% x_t = y
% y_t = -M^(-1)*K*x + M^(-1)*F


F = forcingFunctionVector(n,omega,t);

dX = V;
dV = -K*X+F;

y = [dX;dV];

end
