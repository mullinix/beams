function forced_response_loop(kin,gamma_in)
% function forced_response_loop(kin,gamma_in)
% loop over k and gamma values and run forced_response.m


%% forced response loop
%% start clean
% clear;

% %% get matrices
% % set to 1 to reassemble the matrices
% rebuild = 0;
% 
% % if reassembly is required, run code, else load matrices
% if(rebuild == 1)
%     build_coupling_stiffness;
% else
load matrices;
% end

X2 = M;

%% setup forced response params
% set perturbation strength
epsilon = 1e-2;
% create damping matrix - combination of mass and stiffness
s1 = 1e-1;
s2 = 1e-1;
C = s1*M+s2*K;

%% ODE params
% initial displacement is random,small; init. vel. is zero.
Y0 = [(rand(length(X2),1)-0.5)*1e-12;zeros(length(X2),1)];

% set time vector
numpts = 1e3;
% num_cycles = 2;
num_beams = props.num_beams;
t_final = 1e3;%;(2*pi*num_cycles)/freq;%
t = linspace(0,t_final,numpts);%[0,t_final];%

% find the inverse only once
X2inv = inv(X2);

%% looping params
% Omega is set in a loop. Set the factor k.
for k_idx=1:length(kin)
k = kin(k_idx);
    for g_idx=1:length(gamma_in)
        Omega=gamma_in(g_idx)/props.T;

        tic;
        %% ode solution
        [t,Y] = ode23t(@(t,Y)beam_ode_perturb(t,Y,K,G,P,Sigma,C,X2inv,k,epsilon,Omega),t,Y0);
        tf=toc;
        fprintf(1,'k=%d, gamma=%.1f, time to solve ODE: %.1fs\n',k,gamma_in(g_idx),tf);

        %% store max response values, SV and W
        % TODO: find max(Y(SV)) and max(Y(W)) to store for this Omega, k
        spts = grab_elts(5,length(Y(1,:))/props.num_beams/2,1,2);
        spts = repmat(spts,1,3);
        v_pts = grab_elts(5,length(Y(1,:))/props.num_beams/2,1,3);
        v_pts = repmat(v_pts,1,3);
        w_pts = grab_elts(5,length(Y(1,:))/props.num_beams/2,1,5);
        w_pts = repmat(w_pts,1,3);
        v_full = grab_elts(5,length(Y(1,:))/props.num_beams/2,2,3);
        v_full = repmat(v_full,1,3);
        w_full = grab_elts(5,length(Y(1,:))/props.num_beams/2,2,5);
        w_full(1) = (1==0);
        w_full = repmat(w_full,1,3);

        numpts = length(t);

        x=linspace(0,props.L,props.nelts+1);
        x=x(end-1:end);

        displacement.U=zeros(numpts,props.num_beams);
        displacement.V=zeros(numpts,props.num_beams);
        displacement.W=zeros(numpts,props.num_beams);
        Si = Y(:,spts);
        Vi = Y(:,v_pts); 
        Wi = Y(:,w_pts);
        Vfull = Y(:,v_full);
        Wfull = Y(:,w_full);
        tic;
        len=size(Vi,2);
        len2=size(Vfull,2);
        for j=1:num_beams
            plot_pts = len/num_beams*j;
            plot_pts2 = len2/num_beams*j-3:len2/num_beams*j;
            H_vw = hvw_analysis(x,Vfull(:,plot_pts2),Wfull(:,plot_pts2));
            Ui = Si(:,plot_pts)-H_vw;
            u=real(Ui);           
            v=real(Vi(:,plot_pts)); 
            w=real(Wi(:,plot_pts)); 
            displacement.U(:,j)=u;
            displacement.V(:,j)=v;
            displacement.W(:,j)=w;
        end
        tf=toc;
        fprintf(1,'k=%d, gamma=%.1f, nested loop time: %.1fs\n',k,gamma_in(g_idx),tf);    
        fname = sprintf('./data/forced_response_sweep/forced_response_data-k-%d-gamma-%.1f.mat',k,gamma_in(g_idx));
        save(fname,'displacement');
    end
end


end
