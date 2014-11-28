% This file generates the 'p'  and 'm' matrices from symbolic computations.
%% create symbolics. declare real.
syms a x xe xe1 L
a = sym(a,'real');
x = sym(x,'real');
xe = sym(xe,'real');
xe1 = sym(xe1,'real');
L = xe1-xe;
%% create N vectors.
Nv1 = [i1-3*x^2/L^2+2*x^3/L^3; 
           x-2*x^2/L+x^3/L^2; 
         3*x^2/L^2-2*x^3/L^3; 
              -x^2/L+x^3/L^2;];
Nw1 = Nv1;
Ns1 = [1-x/L; x/L];
Ns = [Ns1(1); 0; 0; 0; 0; Ns1(2); 0; 0; 0; 0;];
Nv = [0; Nv1(1); Nv1(2); 0; 0; 0; Nv1(3); Nv1(4); 0; 0;];
Nw = [0; 0; 0; Nw1(1); Nw1(2); 0; 0; 0; Nw1(3); Nw1(4);];
% %% find first derivatives.
% Nvx = diff(Nv,x);
% Nwx = diff(Nw,x);
%% compute outer products.
NNs = Ns*Ns';
NNv = Nv*Nv';
NNw = Nw*Nw';
% %% construct the polynomial from the integrand.
% polly = a*(L-x)+0.5*(L^2-x^2);
%% sum outer products.
NspNv = NNs+NNv;
NspNvpNw = NspNv + NNw;
%% define integrands
p_intgrand = NspNv;
m_intgrand = NspNvpNw;
%% solve definite integral
p_I = int(p_intgrand,x,xe,xe1);
m_I = int(m_intgrand,x,xe,xe1);
p_I = simplify(p_I);
m_I = simplify(m_I);
%% write matrix to file.
p_fh = fopen('p_syms.dat','w');
m_fh = fopen('m_syms.dat','w');
[m,n] = size(p_I);
for i=1:m
    for j=1:n
        fprintf(p_fh,'P(%d,%d): %s\n',i,j,char(p_I(i,j)));
        fprintf(m_fh,'M(%d,%d): %s\n',i,j,char(m_I(i,j)));
    end
end
fclose(p_fh);
fclose(m_fh);