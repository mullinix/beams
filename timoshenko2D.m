function timoshenko2D(num_elts,num_beams)
%function err = timoshenko2D(params)

%% System parameters
% L = 0.1524; % length of beam [m] (0.45)
% N = 2; % number of elements [1]
base = pwd;
%% Write custom parameters to file for assembly
% Valid parameters (case sensitive), (default setting), description:
%     N,            (100),      # of nodes
%     L,            (0.1524),   length of rod
%     g,            (9.81),     accel. due to gravity
%     rho,          (800),      density
%     E,            (2.07e11),  Young's modulus
%     G,            (E/2.6),    Shear modulus
%     mu,           (0.833),    "shear coefficient"
%     w,            (0.0254),   width of beam at base
%     h,            (0.0046),   height of beam at base
%     twist_angle,  (45),       twist angle at tip wrt base
%     breadth_taper,(2.56),     taper ratio, width
%     depth_taper,  (2.29),     taper ratio, height
%     Omega,        (100),      rotation speed, rps
%     offset        (0),        dist. from center of rot. to base of beam
% Parser format: "N = 100" on each line, N = from list above, 100 = value
% to set.
% Whitespace will be trimmed

params = [0.1,0,3,0];
% params = [delta,gamma,linearization,print];

linearization = params(3);

rho = 7820;
E = 210e9;
depth = 0.25;
breadth = 0.25;
A = depth*breadth;
I = breadth*(depth^3)/12;
alpha = 70;
L = alpha*sqrt(I/A);
T_squared = (rho*A*L^4)/(E*I);
T = sqrt(T_squared);


delta = params(1);
gamma = params(2);
a = delta*L;%0.1;%
Omega = gamma/T;%20*(2*pi);%
props.rho = rho;
props.E = E;
props.depth = depth;
props.breadth = breadth;
props.A = A;
props.I = I;
props.L = L;
props.T_squared = T_squared;
props.T = T;
props.delta = delta;
props.gamma = gamma;
props.a = a;
props.Omega = Omega;
% display(props);
% 
% delete('*.dat');
% delete('*.txt');

%% build material
fmat =fopen('mat.txt','w');
fprintf(fmat,'E = %.15e\n',E);
fprintf(fmat,'rho = %.15e\n',rho);
fprintf(fmat,'geometry_file = elements.txt\n');
fclose(fmat);

%% build geometry
% num_elts = 70;
% num_beams = 3;
num_nodes = num_elts+1;
fgeo = fopen('elements.txt','w');
for j=0:num_beams-1
    for idx=j*(num_nodes)+1:(j+1)*num_elts+j
        fprintf(fgeo,'%d\t%.15e\t%.15e\t%d,%d\n',idx-j,depth,breadth,idx,idx+1);
    end
end
%     for idx=1:num_elts
%         fprintf(fgeo,'%d\t%.15e\t%.15e\t%d,%d\n',idx,depth,breadth,idx,idx+1);
%     end
%     for idx=(num_nodes)+1:2*num_elts+1
%         fprintf(fgeo,'%d\t%.15e\t%.15e\t%d,%d\n',idx-1,depth,breadth,idx,idx+1);
%     end
%     for idx=2*(num_nodes)+1:3*num_elts+2
%         fprintf(fgeo,'%d\t%.15e\t%.15e\t%d,%d\n',idx-2,depth,breadth,idx,idx+1);
%     end
% fprintf(fgeo,'%d\t%.15e\t%.15e\t%d,%d\n',idx,depth,breadth,1,num_elts+1);
% fprintf(fgeo,'%d\t%.15e\t%.15e\t%d,%d\n',idx,depth,breadth,num_elts+1,2*num_elts+1);
% fprintf(fgeo,'%d\t%.15e\t%.15e\t%d,%d\n',idx,depth,breadth,2*num_elts+1,1);
fclose(fgeo);

%% build nodes
dx = L/num_elts;
nodes = (0:num_elts)';
fnodes = fopen('node_locations.txt','w');
for j=0:num_beams-1
    bar = [(nodes*dx+a).*cos(j*2*pi/num_beams),(nodes*dx+a).*sin(j*2*pi/num_beams)];
    for idx = 1:num_nodes
        fprintf(fnodes,'%d\t%.15e,%.15e,%.15e\n',j*num_nodes+idx,bar(idx,1),bar(idx,2),0);
    end
end
% bar_1 = [nodes*dx+a,zeros(num_elts+1,1)];
% bar_2 = [(nodes*dx+a).*cos(2*pi/3),(nodes*dx+a).*sin(2*pi/3)];
% bar_3 = [(nodes*dx+a).*cos(4*pi/3),(nodes*dx+a).*sin(4*pi/3)];
% for idx = 1:num_nodes
%     fprintf(fnodes,'%d\t%.15e,%.15e,%.15e\n',idx,bar_1(idx,1),bar_1(idx,2),0);
% end
% for idx = 1:num_nodes
%     fprintf(fnodes,'%d\t%.15e,%.15e,%.15e\n',num_nodes+idx,bar_2(idx,1),bar_2(idx,2),0);
% end
% for idx = 1:num_nodes
%     fprintf(fnodes,'%d\t%.15e,%.15e,%.15e\n',2*num_nodes+idx,bar_3(idx,1),bar_3(idx,2),0);
% end
fclose(fnodes);

%% build loads
floads = fopen('loads.txt','w');
fprintf(floads,'1\t0\t0,-1,0\n');
fclose(floads);

%% build bcs
fbcs = fopen('bcs.txt','w');
fprintf(fbcs,'1\t1,1,1,1,0\n');
for i=1:num_beams-1
    fprintf(fbcs,'%d\t1,1,1,1,0\n',i*num_nodes+1);
end
fclose(fbcs);

%% write inputs file

fh = fopen('inputs.txt','w');
%fprintf(fh,'N=%d\n',N);
%fprintf(fh,'breadth_taper=%.15e\n',1);
%fprintf(fh,'depth_taper=%.15e\n',1);
%fprintf(fh,'h=%.15e\n',0.0254);
%fprintf(fh,'twist_angle=%.15e\n',0);
fprintf(fh,'composite_number = %d\n',1);
fprintf(fh,'material_files = %s\n','mat.txt');
fprintf(fh,'loads_file = %s\n','loads.txt');
fprintf(fh,'nodes_file = %s\n','node_locations.txt');
fprintf(fh,'bcs_file = %s\n','bcs.txt');
fprintf(fh,'displacement_dofs = %s\n','1,2,2');
fprintf(fh,'total_length = %.15e\n',L);
fprintf(fh,'a = %.15e\n',a);
fprintf(fh,'Omega = %.15e\n',Omega);
fclose(fh);
%% Calculate stiffness and mass matrices
tic;
system([base,'/analyzer']);
tf = toc;
fprintf(1,'<matlab> c execution time: %.5e\n',tf);
%% Load stiffness and mass matrices
tic;
K = load('K.dat');
M = load('M.dat');
G = load('G.dat');
Sigma = load('Sigma.dat');
P = load('P.dat');
tf = toc;
fprintf(1,'data load time: %.5e\n',tf);

% K(1:8,:) = []; K(:,1:8) = [];
% M(1:8,:) = []; M(:,1:8) = [];

% %% Establish polynomial eigenvalue coefficient matrices
% X0 = (K+Omega^2*(Sigma-P));
% X1 = 2*G*Omega;
% X2 = M;
% 
% %% Adimensionalize
% X1 = X1./T;
% X2 = X2./T_squared;
% 
% %% Apply rescaling
% % % Original scales of matrices
% % X2_min = real(floor(log10(min(abs(X2(X2~=0)))))); 
% % X2_max = real(floor(log10(max(abs(X2(X2~=0))))));
% % fprintf(1,'Pre-Scale of X2: [10^%d->10^%d]\n',X2_min,X2_max);
% % X1_min = real(floor(log10(min(abs(X1(X1~=0)))))); 
% % X1_max = real(floor(log10(max(abs(X1(X1~=0))))));
% % fprintf(1,'Pre-Scale of X1: [10^%d->10^%d]\n',X1_min,X1_max);
% % X0_min = real(floor(log10(min(abs(X0(X0~=0)))))); 
% % X0_max = real(floor(log10(max(abs(X0(X0~=0))))));
% % fprintf(1,'Pre-Scale of X0: [10^%d->10^%d]\n',X0_min,X0_max);
% % 
% % % space_rescaling_factor = 10^(-X1_min*0);%1e-4;% alpha
% % time_rescaling_factor = 10^((ceil((X1_max-X1_min)/2)+X1_min-X2_min));%1e17;% beta
% % time_rescaling_factor2 = 10^((ceil((X0_max-X0_min)/2)+X0_min-X1_min));%1e17;% beta
% % 
% % X1 = X1.*time_rescaling_factor;
% % X2 = X2.*time_rescaling_factor;
% % KK=sparse(KK); MM = sparse(MM);
% 
% % % New scales
% % X2_min = real(floor(log10(min(abs(X2(X2~=0)))))); 
% % X2_max = real(floor(log10(max(abs(X2(X2~=0))))));
% % fprintf(1,'New-Scale of X2: [10^%d->10^%d]\n',X2_min,X2_max);
% % X1_min = real(floor(log10(min(abs(X1(X1~=0)))))); 
% % X1_max = real(floor(log10(max(abs(X1(X1~=0))))));
% % fprintf(1,'New-Scale of X1: [10^%d->10^%d]\n',X1_min,X1_max);
% % X0_min = real(floor(log10(min(abs(X0(X0~=0)))))); 
% % X0_max = real(floor(log10(max(abs(X0(X0~=0))))));
% % fprintf(1,'New-Scale of X0: [10^%d->10^%d]\n',X0_min,X0_max);
% 
% %% Solve the generalized eigenvalue problem
% % u(v,t) = v*sin(omega*t)
% % [M](w^2)v -[K]v=0
% 
% tic;
% 
% if(linearization==1)
%     % std linearization
%     [shape,omega] = polyeig(-X0,-X1,X2);
%     evals = omega;
% elseif(linearization==3)
%     % L3 linearization
%     Z = zeros(size(K));
%     A = [Z,-X0; X2, Z];
%     B = [X2, X1; Z, X2];
%     [shape,omega] = eig(A,B);
%     evals = diag(omega);
%     omega = evals;
% elseif(linearization==4)
%     % L4 linearization
%     Z = zeros(size(K));
%     A = [X0,Z; X1, X0];
%     B = [Z, -X0; X2, Z];
%     [shape,omega] = eig(A,B);
%     evals = diag(omega);
%     omega = evals;
% end
% 
% omega = abs(omega);
% % omega = omega/T;
% omega = sort(omega,'ascend');
% 
% freqs = omega./(2*pi);
% [freqs,ix] = sort(abs(freqs),'ascend');
% omega = omega(ix);
% tf = toc;
% fprintf(1,'eigenvalue solve time: %.5e\n',tf);
% 
% if(params(4)==0)
%     assignin('base','omega',omega);
%     assignin('base','evals',evals);
%     assignin('base','freqs',freqs);
%     assignin('base','shape',shape);
    assignin('base','M',M);
    assignin('base','K',K);
    assignin('base','G',G);
    assignin('base','P',P);
    assignin('base','Sigma',Sigma);
%     assignin('base','X0',X0);
%     assignin('base','X1',X1);
%     assignin('base','X2',X2);
    assignin('base','T',T);
    assignin('base','props',props);
    assignin('base','num_beams',num_beams);
    assignin('base','num_elts',num_elts);
% else
%     fout = fopen(sprintf('omega-d%d-g%d-lin%d.txt',delta,gamma,linearization),'w');
%     fprintf(fout,'%1.6e\n',omega(1:30));
%     fclose(fout);
% end
% 
% % assignin('base','time_rescaling_factor',time_rescaling_factor);
% % assignin('base','space_rescaling_factor',space_rescaling_factor);
% 
% % mode_to_try = 1;%3*round(length(omega_sqr)/4);
% 
% %% calculate eigenvalue error
% % result = 0.0;
% % % omega_sqr = omega;
% % % omega_sqr = omega_sqr(ix);
% % for mode_to_try=1:length(evals)
% %     tmp_result = (X0 + evals(mode_to_try)*X1 + (evals(mode_to_try)^2)*X2)*shape(:,mode_to_try);
% %     %M*omega_sqr(mode_to_try)*shape(:,mode_to_try)-K*shape(:,mode_to_try);
% %     result = result + tmp_result'*tmp_result;
% % end
% % % result = result/(num_elts*dofs)^2;
% % fprintf(1,'Eigenvalue error (SSE): %.5e\n',result);
% 
% %% calculate solution accuracy
% % plot(real(evals),imag(evals),'*');
% err = sum(abs(1-(abs(real(evals))+abs(imag(evals)))./abs(evals)));
% display(err);
end