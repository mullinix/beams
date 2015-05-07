function Ydot = beam_ode_perturb(t,Y,K,G,P,Sigma,C,X2inv,k,epsilon,Omega)
% function Ydot = beam_ode_perturb(t,Y,K,G,P,Sigma,C,X2inv,k,epsilon,Omega)
% A perturbed version of beam_ode.
x=Y(1:end/2);
y=Y(end/2+1:end);

Om = calc_omega(t,k,epsilon,Omega);
Om_dot = Om(2);
Om = Om(1);

X0 = (K+Om^2*(Sigma-P)+Om_dot*G);
X1 = 2*G*Om+C;

xdot = y;
ydot = -X2inv*(X1*y+X0*x);

Ydot = [xdot;ydot];

end