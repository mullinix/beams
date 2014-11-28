% This file generates the 'sigma' matrix from symbolic computations.
%% create symbolics. declare real.
syms a x xe xe1 L l
a = sym(a,'real');
x = sym(x,'real');
xe = sym(xe,'real');
xe1 = sym(xe1,'real');
L = sym(L,'real');
l = sym(l,'real');
% l = xe1-xe;
%% create N vectors.
Nv1 = [1-3*x^2/l^2+2*x^3/l^3; x-2*x^2/l+x^3/l^2; 3*x^2/l^2-2*x^3/l^3; -x^2/l+x^3/l^2;];
Nw1 = Nv1;
Nv = [0; Nv1(1); Nv1(2); 0; 0; 0; Nv1(3); Nv1(4); 0; 0;];
Nw = [0; 0; 0; Nw1(1); Nw1(2); 0; 0; 0; Nw1(3); Nw1(4);];
%% find first derivatives.
Nvx = diff(Nv,x);
Nwx = diff(Nw,x);
%% compute outer products.
NNvx = Nvx*Nvx';
NNwx = Nwx*Nwx';
%% construct the polynomial from the integrand.
polly = a*(L-(x+xe))+0.5*(L^2-(x+xe)^2);
%% sum outer products.
NvpNw = NNvx+NNwx;
%% define integrand
intgrand = polly*NvpNw;
%% solve definite integral
aa=0;
bb=l;
I = int(intgrand,x,aa,bb);
I = collect(simplify(420*l*I),l);
%% write matrix to file.
fh = fopen('sigma_syms.dat','w');
[m,n] = size(I);
for i=1:m
    for j=1:n
        fprintf(fh,'S(%d,%d): %s\n',i,j,char(I(i,j)));
    end
end
fclose(fh);