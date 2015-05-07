function Ydot = beam_ode(~,Y,X0,X1,X2inv)
% function Ydot = beam_ode(~,Y,X0,X1,X2inv)
% an ode file for the beam problem without forcing or nonhomogeneous terms

x=Y(1:end/2);
y=Y(end/2+1:end);

xdot = y;
ydot = -X2inv*(X1*y+X0*x);

Ydot = [xdot;ydot];

end