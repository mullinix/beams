clear;
%% load values
% tic;
% load('w-bending_unclamped_disconnected.mat');
% tload = toc;
% fprintf(1,'Time to load mat: %.5es\n',tload);
timoshenko2D(2,1);
linearization = 3;

%% build vars
nn = num_beams;
angles = (0:nn-1).*(2*pi/nn);
a_locs = [cos(angles)',sin(angles)']*props.a;
locs_perm = [a_locs(2:end,:);a_locs(1,:)];
locs_diff = locs_perm-a_locs;
a_lens = locs_diff(:,1)'*locs_diff(:,1)+locs_diff(:,2)'*locs_diff(:,2);
a_lens = sqrt(a_lens);
l = props.a*2*pi/3; % only for axial coupling
k=props.E*props.I./a_lens;
beta=1e4;
k_round = beta*props.E*props.I./props.L;
kk = [k(1),-k(1);-k(1),k(1)];
kk_round = [ k_round,-k_round;
            -k_round, k_round];
% num_elts = 1;
num_nodes = num_elts+1;
node_dofs = 5;
bc_nodes = (0:num_beams-1)*num_nodes+1;
free_dofs = bc_nodes*node_dofs;
adjusted_dofs = free_dofs-(1:nn).*(node_dofs-1);
% spy(K)
% pause
%% modify stiffness matrix
for i=1:nn-1
    gdofs = adjusted_dofs(i:i+1);
    K(gdofs,gdofs) = K(gdofs,gdofs)+kk_round; %#ok<SAGROW>
end
gdofs = adjusted_dofs([nn,1]);
K(gdofs,gdofs) = K(gdofs,gdofs)+kk_round;
% K(1,:)=[]; K(:,1)=[];
% M(1,:)=[]; M(:,1)=[];
% G(1,:)=[]; G(:,1)=[];
% P(1,:)=[]; P(:,1)=[];
% Sigma(1,:)=[]; Sigma(:,1)=[];
%% set parameters for system solution
delta = props.delta;
% gamma = props.gamma;
gamma = 1;
Omega = gamma/T;
% Omega = props.Omega;
T = props.T;
T_squared = props.T_squared;

%% set matrices for solution
X0 = (K+Omega^2*(Sigma-P));
X1 = 2*G*Omega;
X2 = M;
%% adimensionalize
X1 = X1./T;
X2 = X2./T_squared;
%% solve normal modes using linearization 3
% spy(K)
tic;
Z = zeros(size(K));
if(linearization==1)
    % std linearization
    [shape,omega] = polyeig(-X0,-X1,X2);
    evals = omega;
elseif(linearization==3)
    % L3 linearization
    Z = zeros(size(K));
    A = [Z,-X0; X2, Z];
    B = [X2, X1; Z, X2];
    [shape,omega] = eig(A,B);
    evals = diag(omega);
    omega = evals;
elseif(linearization==4)
    % L4 linearization
    Z = zeros(size(K));
    A = [X0,Z; X1, X0];
    B = [Z, -X0; X2, Z];
    [shape,omega] = eig(A,B);
    evals = diag(omega);
    omega = evals;
end
omega = abs(omega);
% omega = omega/T;
omega = sort(omega,'ascend');

freqs = omega./(2*pi);
[freqs,ix] = sort(abs(freqs),'ascend');
omega = omega(ix);
shape = shape(:,ix);
tf = toc;
err = sum(abs(1-(abs(real(evals))+abs(imag(evals)))./abs(evals)));
display(err);
fprintf(1,'eigenvalue solve time: %.5e\n',tf);
