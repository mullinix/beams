%% start clean
clear;

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
epsilon = 1e-1;

% Omega is set in a loop. Set the factor k.
k = 1;

% create damping matrix - combination of mass and stiffness
s1 = 1e-3;
s2 = 1e-3;
C = s1*M+s2*K;

%% ODE params
% initial displacement is random,small; init. vel. is zero.
Y0 = [(rand(length(X2),1)-0.5)*1e-12;zeros(length(X2),1)];

% get values from data runs
% TODO: get the frequency to determine the length of a run and Omega value
Omega=25/props.T;%props.Omega;

% set time vector
numpts = 1e4;
num_cycles = 2;
num_beams = props.num_beams;
t_final = 1e2;%;(2*pi*num_cycles)/freq;%
t = linspace(0,t_final,numpts);%[0,t_final];%

% find the inverse only once
X2inv = inv(X2);

tic;
%% ode solution
[t,Y] = ode23t(@(t,Y)beam_ode_perturb(t,Y,K,G,P,Sigma,C,X2inv,k,epsilon,Omega),t,Y0);
tf=toc;
fprintf(1,'time to solve ODE: %.1fs\n',tf);

%% store max response values, SV and W
% TODO: find max(Y(SV)) and max(Y(W)) to store for this Omega, k
spts = grab_elts(5,length(Y(1,:))/props.num_beams/2,1,2);
spts = repmat(spts,1,3);
s_max = max(max(abs(Y(:,spts))));
v_pts = grab_elts(5,length(Y(1,:))/props.num_beams/2,1,3);
v_pts = repmat(v_pts,1,3);
v_max = max(max(abs(Y(:,v_pts))));
w_pts = grab_elts(5,length(Y(1,:))/props.num_beams/2,1,5);
w_pts = repmat(w_pts,1,3);
w_max = max(max(abs(Y(:,w_pts))));
v_full = grab_elts(5,length(Y(1,:))/props.num_beams/2,2,3);
v_full = repmat(v_full,1,3);
w_full = grab_elts(5,length(Y(1,:))/props.num_beams/2,2,5);
w_full(1) = (1==0);
w_full = repmat(w_full,1,3);

numpts = length(t);

x=linspace(0,props.L,props.nelts+1);
x=x(2:end);

U=zeros(numpts,props.num_beams);
V=zeros(numpts,props.num_beams);
W=zeros(numpts,props.num_beams);
tic;
for i=1:length(t)
    Si = Y(i,spts);  
    Vi = Y(i,v_pts); 
    Wi = Y(i,w_pts); 
    Vfull = Y(i,v_full);
    Wfull = Y(i,w_full);
    len=length(Vi);
    len2=length(Vfull);
    for j=1:num_beams
        plot_pts = (len/num_beams*(j-1)+1):len/num_beams*j;
        plot_pts2 = (len2/num_beams*(j-1)+1):len2/num_beams*j;
        H_vw = hvw(x,Vfull(plot_pts2),Wfull(plot_pts2));
        Ui = Si(plot_pts)-H_vw;
        u=real(Ui);           
        v=real(Vi(plot_pts)); 
        w=real(Wi(plot_pts)); 
        U(i,j)=u(end);
        V(i,j)=v(end);
        W(i,j)=w(end);
    end
end
tf=toc;
fprintf(1,'nested loop time: %.1fs\n',tf);
