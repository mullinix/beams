function Ydot = beam_ode_forced(t,Y,varargin)
% function Ydot = beam_ode_forced(t,Y,varargin)
% Ode function file.
%
% Varargin: 13 elements
%
% X2inv, mode,  epsilon, Omega,  F_idx, 
% K,     Sigma, P,       G,      C,    
% F,     s_pts, v_full 

%% set parameters from input
X2Inv=varargin{1}; % inverted mass matrix
mode=varargin{2}; % eigenfrequency
epsilon=varargin{3}; % forcing strength
Omega=varargin{4}; % rotation frequency
F_idx=varargin{5}; % where force is applied
% system matrices
K=varargin{6}; 
Sigma=varargin{7};
P=varargin{8};
G=varargin{9};
C=varargin{10};
% inertial force vector
F=varargin{11};
% where inertial forces are applied
s_idx=varargin{12};
v_idx=varargin{13};

%% grab ICs
x=Y(1:end/2);
y=Y(end/2+1:end);

%% calculate Omega
Om = calc_omega(t,epsilon,Omega);
Om_dot = Om(2);
Om = Om(1);

%% set custom forces
Pf = epsilon*sin(mode*t).*F_idx;%epsilon;%
% Pf=1e3*F_idx*(t<1e-6);
% F=F*(t<50);

%% set inertial forces, add custom forces
F=F.*(Om^2).*s_idx+F.*Om_dot.*v_idx+Pf;%+F_idx.*Om;

%% calculate matrices for this omega 
X0 = (K+Om^2*(Sigma-P)+Om_dot*G);
X1 = 2*G*Om+C;

%% set system values
xdot = y;
ydot = X2Inv*(F-X1*y-X0*x);

%% set output
Ydot = [xdot;ydot];

end