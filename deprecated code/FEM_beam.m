function sys = fem_beam(L)
% Generates a FEM model of a 3D cantilever Timoshenko beam
% Replace arbitrarily many of the preset parameters by inputs to the file!
% 
% Please cite this model as:
% H. Panzer, J. Hubele et al.:
% Generating a Parametric Finite Element Model
%   of a 3D Cantilever Timoshenko Beam Using Matlab
% Technical reports on Automatic Control, vol. TRAC-4,
% Institute of Automatic Control, Technische Universitat Munchen, 2009
%
% Usage example: sys = fem_beam(1); % for beam of length L=1m
%
% Authors:  Jorg Hubele (joerg.hubele@mytum.de)
%           Heiko Panzer (panzer@tum.de)
% Last Change: 20 Nov 2009
%
% Visit our webpage for more code downloads: www.rt.mw.tum.de

% ========== parameters ==========
% (input) L:total beam length (x) [m]
N = 20; %number of elements [1]
t = .01; %t: beam thickness (y) [m]
h = .15; %h: beam height (z) [m]
rho = 7850; %density of steel [kg/m^3]
m = L*t*h*rho; %m: total beam mass [kg]
E = 210e9; %E: Young's modulus of steel [N/m^2]
nu = 3/10; %nu: Poisson's ratio
d1 = 8e-6; d2 = 8; %d1, d2: dampening ratio (D = d1*K + d2*M)

% ========== physical quantities ==========
G = E/2/(1+nu); % G: Shear modulus [N/m^2]
l = L/N; % l: beam element length
A = t*h ; % beam area [m^2]
ASy = 5/6*A ;
ASz = 5/6*A ; % effective area of shear
Iy = 1/12*h^3*t;
Iz = 1/12*t^3*h; % second moments of area [m^4]
Ip = 1/12*t*h*(h^2+t^2); % polar moment of inertia [m^4]
It = min([h t])^3*max([h t])/7; % torsion constant [m^4]
Py = 12*E*Iz/(G*ASy*l^2) ; Pz = 12*E*Iy/(G*ASz*l^2) ; % Phi
% ========== element mass and stiffness matrices ==========
% element mass matrix
M11 = zeros(6,6);
M11(1,1) = 1/3;
M11(2,2) = 13/35 + 6*Iz/(5*A*l^2);
M11(3,3) = 13/35 + 6*Iy/(5*A*l^2);
M11(4,4) = Ip/(3*A);
M11(5,5) = l^2/105 + 2*Iy/(15*A);
M11(6,6) = l^2/105 + 2*Iz/(15*A);
M11(6,2) = 11*l/210 + Iz/(10*A*l);
M11(2,6) = M11(6,2) ;
M11(5,3) = -11*l/210 - Iy/(10*A*l);
M11(3,5) = M11(5,3) ;
M22 = -M11 + 2*diag(diag(M11));
M21 = zeros(6,6);
M21(1,1) = 1/6;
M21(2,2) = 9/70 - 6*Iz/(5*A*l^2);
M21(3,3) = 9/70 - 6*Iy/(5*A*l^2);
M21(4,4) = Ip/(6*A);
M21(5,5) = -l^2/140 - Iy/(30*A);
M21(6,6) = -l^2/140 - Iz/(30*A);
M21(6,2) = -13*l/420 + Iz/(10*A*l);
M21(2,6) = -M21(6,2);
M21(5,3) = 13*l/420 - Iy/(10*A*l);
M21(3,5) = -M21(5,3);
Me = m/N*[M11, M21'; M21, M22];

% element stiffness matrix
K11 = zeros(6,6);
K11(1,1) = E*A/l ;
K11(2,2) = 12*E*Iz/(l^3*(1+Py)) ;
K11(3,3) = 12*E*Iy/(l^3*(1+Pz)) ;
K11(4,4) = G*It/l ;
K11(5,5) = (4+Pz)*E*Iy/(l*(1+Pz)) ;
K11(6,6) = (4+Py)*E*Iz/(l*(1+Py)) ;
K11(2,6) = 6*E*Iz/(l^2*(1+Py)) ;
K11(6,2) = K11(2,6) ;
K11(3,5) = -6*E*Iy/(l^2*(1+Pz)) ;
K11(5,3) = K11(3,5) ;
K22 = -K11 + 2*diag(diag(K11));
K21 = K11 - 2*diag(diag(K11));
K21(5,5) = (2-Pz)*E*Iy/(l*(1+Pz)) ;
K21(6,6) =(2-Py)*E*Iz/(l*(1+Py)) ;
K21(2,6) = -K21(6,2);
K21(3,5) = -K21(5,3);
Ke = [K11, K21'; K21, K22];

% ========== assembly ==========
% global mass and stiffness matrices: N*(6+1) dofs
M = zeros(6*(N+1), 6*(N+1));
K = zeros(6*(N+1), 6*(N+1));
% assembly
for i = 1:N
    a=1+6*(i-1):6*(i+1);
    K(a,a) = K(a,a) + Ke;
    M(a,a) = M(a,a) + Me;
end
% boundary condition
K(1:6,:)=[]; K(:,1:6)=[];
M(1:6,:)=[]; M(:,1:6)=[];
N = N-1; % JAM - experimenting with BCs
% dampening matrix 
D = zeros(6*(N+1));%d1*K + d2*M ;
% global load and output vectors
B=zeros(6*(N+1),1); C=zeros(1,6*(N+1));
B((N+1)*6-3,1)=-1; % input: beam tip, -z
C(1,(N+1)*6-3)=1; % output: beam tip, +z
% ========== second order --> first order ==========C(1,N*6-3)=1;
F = eye(6*(N+1));
E = [F zeros(6*(N+1),6*(N+1)); zeros(6*(N+1),6*(N+1)) M];
A = [zeros(6*(N+1),6*(N+1)) F; -K -D] ;
B = [zeros(6*(N+1),size(B,2)); B];
C = [C, zeros(size(C,1), 6*(N+1))];
sys = dss(A,B,C,0,E);
end