function [k,m] = febeam1(el,xi,len,area,rho,ipt)
% element stiffness matrix, eqn 8.1.11
c = el*xi/len^3;
k = c*[12,      6*len,      -12,        6*len;...
       6*len,   4*len^2,    -6*len,     2*len^2;...
       -12,     -6*len,     12,         -6*len;...
       6*len,   2*len^2,    -6*len,     4*len^2];
if ipt==1
    mm=rho*area*len/420;
    m = mm*[156,    22*len,     54,     -13*len;...
            22*len, 4*len^2,    13*len, -3*len^2;...
            54,     13*len,     156,    -22*len;...
            -13*len,-3*len^2,   -22*len,4*len^2];
elseif ipt==2
    mass=rho*area*len;
    m = diag([mass/2,0,mass/2,0]);
else
    mass=rho*area*len;
    m = mass.*diag([0.5,len^2/78,0.5,len^2/78]);
end
end