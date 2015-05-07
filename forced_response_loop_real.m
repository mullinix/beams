function forced_response_loop_real(kin,gamma_in)
% function forced_response_loop_real(kin,gamma_in)
% A function to loop over k values and rotation speeds and force the system
% at the hub. Not used in the final product.

%% forced response loop
%% start clean
% clear;
fprintf(1,'start time: %s\n',datestr(now));
%% get matrices
% set to 1 to reassemble the matrices
rebuild = 0;

% if reassembly is required, run code, else load matrices
if(rebuild == 1)
    build_coupling_stiffness;
else
load matrices; load modes;
end

X2 = M./props.T./props.T;


%% setup forced response params
% set perturbation strength
epsilon = 1e1;
% create damping matrix - combination of mass and stiffness
s1 = 1e-5;
s2 = s1;
C = s1*M+s2*K;
G=G./props.T; %#ok<NODEF>
C=C./props.T;
num_beams=props.num_beams;

mode=modes.vals(3);

%% ODE params
% initial displacement is random,small; init. vel. is zero.
F_idx = zeros(length(X2)/num_beams,1);%double(grab_elts(length(X2)/3,length(X2),1,3))';
F_idx = repmat(F_idx,num_beams,1);F_idx(length(F_idx)/3+5)=1;
% s_idx=grab_elts(5,length(X2)/num_beams,1,2);
% v_idx=grab_elts(5,length(X2)/num_beams,2,3);
% w_idx=grab_elts(5,length(X2)/num_beams,2,3);
% s_idx = repmat(s_idx,1,num_beams)';
% v_idx = repmat(v_idx,1,num_beams)';
% 
% s_idx = repmat(s_idx,1,num_beams)';

%chunk_size, array_length, group_size, idx_start
s_pts = grab_elts(5,length(X2)/num_beams,1,2);
s_pts = repmat(s_pts,1,num_beams)';
v_pts = grab_elts(5,length(X2)/num_beams,1,3);
v_pts = repmat(v_pts,1,num_beams)';
w_pts = grab_elts(5,length(X2)/num_beams,1,5);
w_pts = repmat(w_pts,1,num_beams)';
v_full = grab_elts(5,length(X2)/num_beams,2,3);
v_full = repmat(v_full,1,num_beams)';
w_full = grab_elts(5,length(X2)/num_beams,2,5);
w_full(1) = (1==0);
w_full = repmat(w_full,1,num_beams)';

Y0 = zeros(length(X2),1);
% Y0(s_idx|v_idx)=shape_store;
Y0 = [Y0;zeros(length(X2),1)];%[eye(length(X2),1)*1e-12,zeros(length(X2),1)];


% find the inverse only once
X2inv = inv(X2);
% MM=eye(size(X2)*2);
% MM(end/2+1:end,end/2+1:end)=X2;

ode_opts = odeset('Stats','on');%,'Refine',10);%,'Mass',MM);

%% looping params
% Omega is set in a loop. Set the factor k.
for k_idx=1:length(kin)
k = kin(k_idx);
    for g_idx=1:length(gamma_in)
        Omega=gamma_in(g_idx);
        % set time vector
        numpts = 1e3;
        num_cycles = 10;
        t_final = num_cycles*(2*pi/Omega);
        t = linspace(0,t_final,numpts);%[0,t_final];%
        tic;
        %% ode solution
        [t,Y] = ode23t(@(T,Y)beam_ode_forced(T,Y,X2inv,mode,epsilon,...
                                             Omega,F_idx,K,Sigma,...
                                             P,G,C,F,s_pts,v_full,...
                                             w_pts),t,Y0,ode_opts);
        tf=toc;
        fprintf(1,'freq=%.1f, gamma=%.1f, time to solve ODE: %.1fs\n',mode,gamma_in(g_idx),tf);

        %% store max response values, SV and W
        % TODO: find max(Y(SV)) and max(Y(W)) to store for this Omega, k
        numpts = length(t);
        % need relative node locations for Hvw solutions
        x=linspace(0,props.L,props.nelts+1);
        x=x(2:end);
        % initialize storage
        Si = Y(:,s_pts);
        Vi = Y(:,v_pts); 
        Wi = Y(:,w_pts);
        Vfull = Y(:,v_full);
        Wfull = Y(:,w_full);
        displacement.U=zeros(numpts,size(Si,2));
        displacement.V=zeros(numpts,size(Si,2));
        displacement.W=zeros(numpts,size(Si,2));

        tic;
        len=size(Vi,2);
        len2=size(Vfull,2);
        for j=1:num_beams
            plot_pts = len/num_beams*(j-1)+1:len/num_beams*j;
            plot_pts2 = len2/num_beams*(j-1)+1:len2/num_beams*j;
            H_vw = hvw_analysis(x,Vfull(:,plot_pts2),Wfull(:,plot_pts2));
            Ui = Si(:,plot_pts)-H_vw;
            u=real(Ui);           
            v=real(Vi(:,plot_pts)); 
            w=real(Wi(:,plot_pts)); 
            displacement.U(:,plot_pts)=u;
            displacement.V(:,plot_pts)=v;
            displacement.W(:,plot_pts)=w;
        end
        displacement.time=t;
        displacement.props=props;
        displacement.modes=modes;
        tf=toc;
        fprintf(1,'freq=%.1f, gamma=%.1f, nested loop time: %.1fs\n',mode,gamma_in(g_idx),tf);    
        fname = sprintf('./data/forced_response_sweep/forced_response_data_real-freq-%.1f-gamma-%.1f.mat',mode,gamma_in(g_idx));
        save(fname,'displacement');
    end
end

fprintf(1,'finish time: %s\n',datestr(now));

end
