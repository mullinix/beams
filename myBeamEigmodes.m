function [evals,shape,omega_sqr] = myBeamEigmodes(N)
L = 1; % length of beam [m] (0.45)
% N = 5000; % number of elements [1]

rho = 2700; % density of steel [kg/m^3] (steel: 7850, Alum: 2700)
% m = L*w*h*rho; % m: total beam mass [kg]
E = 7.2e10; % E: Young's modulus of steel [N/m^2] (steel: 210e9, alum: 7.2e10)
% nu = 3/10; % nu: Poisson's ratio (steel: 0.3, alum: 0.35)

% rods
d = 0.01; % d: beam diameter [m] (0.01)
A = pi*(d/2)^2; % A: beam cross-sectional area [m] 
I = pi/4*(d/2)^4; % I: moment of inertia  

% rectangular beam
% w = 0.02; % w: beam width (lateral) [m] (0.02)
% h = 0.003; % h: beam height (vertical) [m] (0.003)
% A = w*h; % cross section area [m^2], rectangular beam
% I = (w*h^3)/12; % moment of inertia, rectangular beam

l = L/N; % equidistant node length

% element mass matrix
M = (rho*A*l/420)*[156,22*l,54,-13*l;...
					22*l,4*l*l,13*l,-3*l*l;...
					54,13*l,156,-22*l;...
					-13*l,-3*l*l,-22*l,4*l*l];

% element stiffness matrix
K = ((E*I)/(l^3))*[12,6*l,-12,6*l;...
				    6*l,4*l*l,-6*l,2*l*l;...
					-12,-6*l,12,-6*l;...
					6*l,2*l*l,-6*l,4*l*l];

% assembly block
% initialize global matrices
MM = zeros(2*(N+1));
KK = zeros(2*(N+1));
% loop to assemble
for i=1:N
	j = 2*i-1;
	k = j+3;
    %[j:k]
	MM(j:k,j:k) = MM(j:k,j:k) + M;
	KK(j:k,j:k) = KK(j:k,j:k) + K;
end

% apply BC's - fixed -> free; alternatively, set top left 2x2 to identity
KK([1,2],[1,2]) = eye(2);%[]; KK(:,[1,2]) = []; 
MM([1,2],[1,2]) = eye(2);%[]; MM(:,[1,2]) = []; 

% fprintf(1,'Mass Matrix: ');
% MM
% fprintf(1, 'Stiffness Matrix: ');
% KK

% generalized eigenvalue problem
% [M](w^2)v +[K]v=0

% move to a system of 1-D
% M*x_tt+K*x = F
% y = x_t
% y_t = x_tt
% new sys:
% x_t = y
% y_t = -M^(-1)*K*x + F
[shape,omega_sqr] = eig(KK,MM);
P = length(shape);
omega_sqr = diag(omega_sqr);
evals = real(sqrt(omega_sqr))./(2*pi);

for idx=1:P
    shape(5:2:end,idx) = shape(5:2:end,idx)./max(abs(shape(5:2:end,idx)));
end

num_modes = 5;
cc=lines(num_modes);% cool, lines, copper
legend_text = cell(num_modes,1);
for idx=1:num_modes
    plot(linspace(0,L,P/2)',shape(1:2:end,idx+2),'color',cc(idx,:),'LineWidth',3);
    hold on;
    legend_text{idx} = sprintf('Mode %d',idx);
end
set(gca,'FontSize',15);
legend(legend_text);
grid on
xlabel('X: Profile of Beam');
ylabel('Y: Vertical Displacement');
title('Eigenshapes of Beam Displacement');

end