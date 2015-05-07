% Rcurves
%%%%%%%%%%%%%%%%%%%
% plots the anticipated shape of the fra curves. not terribly meaningful;
% another representation of plot_gamma_vs_omega.

%% get matrices
% set to 1 to reassemble the matrices
rebuild = 0;

% if reassembly is required, run code, else load matrices
if(rebuild == 1)
    build_coupling_stiffness;
else
load matrices;
end

X2 = M;

%% setup forced response params
% set perturbation strength
epsilon = 1e0;%5e20;
% create damping matrix - combination of mass and stiffness
s1 = 1e2;
s2 = s1;
C = s1*M+s2*K;

%% get curves
phyles = strsplit(ls([pwd,'/data/omega_evals_sweep/datastruct-*.mat']),'[\n]+');
% ignore "empty" strings from the split.
phyles=phyles(~strcmp(phyles,''));
% instantiate vecror of gamma values
g = zeros(length(phyles),1);
% initialize Y for plotting
Y=zeros(length(phyles),48);
% for each file...
for idx=1:length(phyles)
    % load the data
    load(phyles{idx}); 
    % here we grab the value of gamma from the file name
    % separate name by the dash
    gam=strsplit(phyles{idx},'[-]+');
    % the end has the value of gamma and the extension
    gam=gam{end};
    % remove the extension
    gam=gam(1:end-4);
    %convert to a number
    gam=str2double(gam);
    g(idx)=gam;
    Y(idx,:) = data2save.omega;
end
[g,idx]=sort(g);
Y = Y(idx,:);
gg=linspace(min(g),max(g),1e4);
YY=zeros(length(gg),length(Y(1,:)));
for ii=1:length(Y(1,:))
    YY(:,ii)=spline(g,Y(:,ii),gg);
end

%% ODE params
% initial displacement is random,small; init. vel. is zero.
% Y0 = zeros(length(X2)*2,1);%[(rand(length(X2),1)-0.5)*1e-12,zeros(length(X2),1)];%
% F_idx = zeros(length(X2)/3,1);%double(grab_elts(length(X2)/3,length(X2),1,3))';
% F_idx(3)=1;
% F_idx = repmat(F_idx,3,1);

% % set time vector
% numpts = 5e2;
% % num_cycles = 2;
% num_beams = props.num_beams;
% t_final = 1e3;%;(2*pi*num_cycles)/freq;%
% t = linspace(0,t_final,numpts);%[0,t_final];%

% find the inverse only once
% X2inv = inv(X2);

% ode_opts = odeset('Stats','on');%,'Refine',10);
kin=1:10;
R = zeros(length(kin),length(gg));
% a=norm(X2);
% 


%% looping params
% Omega is set in a loop. Set the factor k.   
plot_gamma_vs_omega;
dpoints_t = cell(length(YY(1,:)),length(kin));
dpoints_y = cell(length(YY(1,:)),length(kin));
figure(1); hold on;
for k_idx=1:length(kin)
k = kin(k_idx);
    for g_idx=1:length(gg)
%         Omega=g(g_idx)/props.T;
% %         X0 = (K+Omega^2*(Sigma-P));
%         X1 = 2*G*Omega+C;
%         b=norm(X1);%c=norm(X0);
%         Rtemp=epsilon./sqrt(abs(k*g(g_idx)-Y(g_idx,:))+Y(g_idx,:));
        Rtemp = 1./(epsilon+abs(YY(g_idx,:).^2-(k*gg(g_idx)).^2));
        R(k_idx,g_idx)=sum(Rtemp);
    end
%     R(k_idx,:) =  R(k_idx,:)./max( R(k_idx,:) );
    for w_idx=1:length(YY(1,:))
        [dpoints_t{w_idx,k},dpoints_y{w_idx,k}]=dataset_intersection(gg,YY(:,w_idx)',gg,k*gg);
        plot(dpoints_t{w_idx,k},dpoints_y{w_idx,k},'r*');
    end
end
figure(2);clf;
plot(gg,R);