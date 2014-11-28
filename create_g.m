% This file generates the 'g' matrix from symbolic computations.
%% create symbolics. declare real.
syms a x xe xe1 L
a = sym(a,'real');
x = sym(x,'real');
xe = sym(xe,'real');
xe1 = sym(xe1,'real');
L = sym(L,'real');
% L = xe1-xe;

%% create N vectors.
Nv1 = [1-3*x^2/L^2+2*x^3/L^3; x-2*x^2/L+x^3/L^2; 3*x^2/L^2-2*x^3/L^3; -x^2/L+x^3/L^2;];
Ns1 = [1-x/L; x/L];
Ns = [Ns1(1); 0; 0; 0; 0; Ns1(2); 0; 0; 0; 0;];
Nv = [0; Nv1(1); Nv1(2); 0; 0; 0; Nv1(3); Nv1(4); 0; 0;];

%% compute outer products.
NvNs = Nv*Ns';
NsNv = Ns*Nv';

%% define integrand
intgrand = NvNs-NsNv;

%% solve definite integral
I = int(intgrand,x,0,L);
I = collect(simplify((60/L)*I),L);

%% write matrix to file.
fh = fopen('g_syms.dat','w');
[m,n] = size(I);
for i=1:m
    for j=1:n
        fprintf(fh,'G(%d,%d): %s\n',i,j,char(I(i,j)));
    end
end
fclose(fh);