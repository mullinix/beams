function dXdt = coupled_diff(~,X,A,B)

Y1 = X(1:end/2);
Y2 = X(end/2+1:end);

dY1 = Y2;
dY2 = A*Y2+B*Y1;

dXdt = [dY1;dY2];

end
