clear;
%% load values
% tic;
% load('w-bending_unclamped_disconnected.mat');
% tload = toc;
% fprintf(1,'Time to load mat: %.5es\n',tload);
% diary on;
nelts=20;
% bad=0;
% while ~bad 
timoshenko2D(nelts);
num_beams = 3;
linearization = 1;
node_dofs = 5;
num_bcs = 4;
gamma = 50;% turning rate: 0=0rpm, 1=900rpm, 10=9krpm, 50=50krpm
beta=1e3;

%% remove truncation error
% K_max = max(max(abs(K)));
% K_low = 10^(real(floor(log10(K_max)))-15);
% K(abs(K)<K_low)=0;
% G_max = max(max(abs(G)));
% G_low = 10^(real(floor(log10(G_max)))-15);
% G(abs(G)<G_low)=0;
% M_max = max(max(abs(M)));
% M_low = 10^(real(floor(log10(M_max)))-15);
% M(abs(M)<M_low)=0;
% P_max = max(max(abs(P)));
% P_low = 10^(real(floor(log10(P_max)))-15);
% P(abs(P)<P_low)=0;
% S_max = max(max(abs(Sigma)));
% S_low = 10^(real(floor(log10(S_max)))-15);
% Sigma(abs(Sigma)<S_low)=0;

%% build multi beams
KK = zeros(size(K)*num_beams); MM=KK; GG=KK; PP=KK; SS=KK;
for b_idx = 1:num_beams
    start_idx = (b_idx-1)*(num_nodes*node_dofs-num_bcs)+1;
    end_idx = b_idx*(num_nodes*node_dofs-num_bcs);
    KK(start_idx:end_idx,start_idx:end_idx)=K;
    GG(start_idx:end_idx,start_idx:end_idx)=G;
    MM(start_idx:end_idx,start_idx:end_idx)=M;
    PP(start_idx:end_idx,start_idx:end_idx)=P;
    SS(start_idx:end_idx,start_idx:end_idx)=Sigma;
end

K=KK; G=GG; M=MM; P=PP; Sigma=SS;
clear KK GG MM PP SS;

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

k_round = beta*props.E*props.I./props.L;
kk = [ k,-k;
      -k, k];
kk_round = [ k_round,-k_round;
            -k_round, k_round];
% num_elts = 1;
num_nodes = num_elts+1;
node_dofs = 5;
[r,c] = size(K);
len=r/3;

bc_nodes = (0:num_beams-1)*num_nodes+1;
free_dofs = bc_nodes*node_dofs;
first_dofs = free_dofs-(1:nn).*(node_dofs-1);
adjusted_dofs = [[first_dofs(1:end-1)',first_dofs(2:end)'];
                 [first_dofs(1),       first_dofs(end)   ]];
% spy(K)
% pause
%% modify stiffness matrix
for i=1:nn
    gdofs = adjusted_dofs(i,:);
    K(gdofs,gdofs) = K(gdofs,gdofs)+kk_round; 
end
% gdofs = [1,adjusted_dofs(nn)];
% K(gdofs,gdofs) = K(gdofs,gdofs)+kk;
% K(1,:)=[]; K(:,1)=[];
% M(1,:)=[]; M(:,1)=[];
% G(1,:)=[]; G(:,1)=[];
% P(1,:)=[]; P(:,1)=[];
% Sigma(1,:)=[]; Sigma(:,1)=[];
%% Apply permutation
% II = eye(len);
% ZZ = zeros(len);
% Rho = [ZZ,ZZ,II;
%        II,ZZ,ZZ;
%        ZZ,II,ZZ];% S3 permutation matrix
% % R2 = Rho^2;
% K = K*Rho;
% spy(K);
% pause;
%% set parameters for system solution
delta = props.delta;
T = props.T;
Omega = gamma/T;
T_squared = props.T_squared;

%% set matrices for solution
X0 = (K+Omega^2*(Sigma-P));
X1 = 2*G*Omega*1i;
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
    a=norm(X2);
    b=norm(X1);
    c=norm(X0);
    alpha = sqrt(c/a);
    beta = 2/(c+alpha*b);
    X0 = X0*beta;
    X1 = X1*alpha*beta;
    X2 = X2*alpha^2*beta;
    [shape,omega] = polyeig(-X0,-X1,X2);
    evals = omega;
    omega = omega*(alpha);
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
[omega,ix] = sort(omega,'ascend');
shape = shape(:,ix);
evals = evals(ix);
freqs = abs(evals);
% savename = sprintf('%d-elts.mat',nelts);
% save(savename,'omega','shape');
tf = toc;
err = sum(abs(1-(abs(real(evals))+abs(imag(evals)))./abs(evals)));
display(err);
fprintf(1,'eigenvalue solve time: %.5e\n',tf);
% if(isnan(err))
%     bad=1;
% end
% nelts=nelts+1;
% end
% diary off;
