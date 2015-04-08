function myBeamForced()

%% Setup mass and stiffness element matrices
L = 0.45; % length of beam [m] (0.45)
N = 10; % number of elements [1]

rho = 7850; % density of steel [kg/m^3] (steel: 7850, Alum: 2700)
% m = L*w*h*rho; % m: total beam mass [kg]
E = 21e10; % E: Young's modulus of steel [N/m^2] (steel: 21e10, alum: 7.2e10)
% nu = 3/10; % nu: Poisson's ratio (steel: 0.3, alum: 0.35)

% % rods
% d = 0.01; % d: beam diameter [m] (0.01)
% A = pi*(d/2)^2; % A: beam cross-sectional area [m] 
% I = pi/4*(d/2)^4; % I: moment of inertia  

% % rectangular beam
w = 0.02; % w: beam width (lateral) [m] (0.02)
h = 0.003; % h: beam height (vertical) [m] (0.003)
A = w*h; % cross section area [m^2], rectangular beam
I = (w*h^3)/12; % moment of inertia, rectangular beam

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

%% Assembly block
% initialize global matrices
MM = zeros(2*(N+1));
KK = zeros(2*(N+1));
% loop to assemble
for i=1:N;%
	j = 2*i-1;
	k = j+3;
	MM(j:k,j:k) = MM(j:k,j:k) + M;
	KK(j:k,j:k) = KK(j:k,j:k) + K;
end

spy(MM);
pause;

%% Apply rescaling
% Original scales of matrices
MM_min = real(floor(log10(min(abs(MM(MM~=0)))))); 
MM_max = real(floor(log10(max(abs(MM(MM~=0))))));
fprintf(1,'Pre-Scale of M: [10^%d->10^%d]\n',MM_min,MM_max);
KK_min = real(floor(log10(min(abs(KK(KK~=0)))))); 
KK_max = real(floor(log10(max(abs(KK(KK~=0))))));
fprintf(1,'Pre-Scale of K: [10^%d->10^%d]\n',KK_min,KK_max);

space_rescaling_factor = 10^(-KK_min*0);%1e-4;% alpha
time_rescaling_factor = 10^(KK_min/2-MM_min);%1e17;% beta

KK = KK.*space_rescaling_factor;
MM = MM.*time_rescaling_factor.*space_rescaling_factor;
% KK=sparse(KK); MM = sparse(MM);

% New scales
MM_min = min(abs(MM(MM~=0))); MM_max = max(abs(MM(MM~=0)));
fprintf(1,'Scale of M: [10^%d->10^%d]\n',real(floor(log10(MM_min))),real(floor(log10(MM_max))));
KK_min = min(abs(KK(KK~=0))); KK_max = max(abs(KK(KK~=0)));
fprintf(1,'Scale of K: [10^%d->10^%d]\n',real(floor(log10(KK_min))),real(floor(log10(KK_max))));

%% Apply BC's - fixed -> free; set top left 2x2 to identity
KK([1,2],:) = []; KK(:,[1,2]) = [];
MM([1,2],:) = []; MM(:,[1,2]) = [];
N=N-1;
% KK([1,2],[1,2]) = eye(2); 
% MM([1,2],[1,2]) = eye(2); 



%% Solve the generalized eigenvalue problem
% u(v,t) = v*sin(omega*t)
% [M](w^2)v +[K]v=0

[shape,omega_sqr] = eig(KK,MM);
omega_sqr = diag(omega_sqr);
omega = sqrt(omega_sqr);
evals = real(omega)./(2*pi);
[evals,ix] = sort(evals,'ascend');
shape = shape(:,ix);

% shape_omega = shape*diag(sin(omega));
%% Verify eigenvalue solution
mode_to_try = 1;
omega_sqr = omega_sqr(ix);
result = MM*omega_sqr(mode_to_try)*shape(:,mode_to_try)-KK*shape(:,mode_to_try);
result = result'*result;
fprintf(1,'Eigenvalue error (SSE): %e\n',result);

%% Setup integration, integrate in time
omega = evals(mode_to_try)*(2*pi);% corresponds to Hz -> radians
t0 = 0;
numCycles=2;
tf = numCycles*(1/omega)*2*pi;% end after four full cycles
numpts = 1e3;
thyme = linspace(t0,tf,numpts);
X = shape(1:2:end,mode_to_try);
y0 = [shape(:,mode_to_try)./max(abs(X));zeros(2*N+2,1)];%
X = X./max(abs(X));
% MM_inv = inv(MM);

% [t,y] = ND_integrator('beam_diff',thyme,y0,KK,MM_inv,omega,'bd4');

% MM_big = [eye(2*N+2),zeros(2*N+2);zeros(2*N+2),MM_inv];
MM_big = [eye(2*N+2),zeros(2*N+2);zeros(2*N+2),MM];
% options = odeset('Mass',MM_big,'RelTol',1e-6);
options = odeset('RelTol',1e-8,'Mass',MM_big);
[t,y] = ode23t(@(t,y) beam_diff2(t,y,KK,omega), thyme, y0,options);

[~,cols] = size(y);
y_out_length = cols;
BAR = linspace(0,L,y_out_length/4); %#ok<*NASGU>

%% Calculate SSE, print to command window
err_vec = y(end,1:2:cols/2-1)'-X;
err_vec = err_vec'*err_vec/numCycles;
fprintf(1,'Error/cycle wrt IC: %1.5e\n',err_vec);
% 
skip_pts = 2*numpts/1e3;

%% Show the 3D representation of one (last) period. show few points relative
% to numpts. this would be 1/4th of points at numpts=1000, or 250pts of fs,
% here we show 1/16th of points since we are taking the last 1/4 points.
% figure(1);
% [Y,T] = meshgrid(linspace(0,1,cols/4),t(1:skip_pts:(numpts/4)+1));
% surf(T,Y,y(end-(numpts/4):skip_pts:end,1:2:cols/2));%,'EdgeColor','None');
% shading interp;
% xlabel('Time, (s)');
% ylabel('Rod (m)');
% zlabel('Displacement, (arb)');
% title('Final cycle of rod vibration');

%% Show the wiggling bar, this wiggles for 4 periods. show few points relative
% % to numpts. this would be 1/4th of points at numpts=1000, or 250pts of fs.
wiggleBar = 1;
if(wiggleBar)
    figure(2); %#ok<*UNRCH>
    plot_idx = skip_pts;
    % % instantiate plotting var
    z = y(plot_idx,1:2:cols/2);
    % % instantiate axes plots
    h1 = plot(BAR,X);
    hold on;
    h2 = plot(BAR,z,'g');
    % % set data sources for plots
    set(h1,'YDataSource','X');
    set(h1,'XDataSource','BAR');
    set(h2,'YDataSource','z');
    set(h2,'XDataSource','BAR');
    % % create legend
    legend({'eig soln';'perturbed';});
    % % Static-ize legend to improve performance
    set(gca,'LegendColorbarListeners',[]); 
    setappdata(gca,'LegendColorbarManualSpace',1);
    setappdata(gca,'LegendColorbarReclaimSpace',1);
    % % set title
    title(sprintf('idx: %d t: %e',plot_idx,t(plot_idx)));
    % % static axes to improve performance
    xlim([0,L]); ylim([-1.25,1.25]);
    xlabel('Rod (m)');
    ylabel('Displacement (arb)')
    drawnow;
    % % loop over data and draw. this is far more efficient than drawnow, and we
    % % can actually control the display speed with skip_pts. pause is required
    % % to get animation.
    for plot_idx = 2*skip_pts:skip_pts:length(t)
        z = y(plot_idx,1:2:cols/2)';
        set(h2,'XData',BAR,'YData',z);
        title(sprintf('idx: %d t: %e',plot_idx,t(plot_idx)));
        pause(0.01);
    end
end

%% Show drift errors
if(y(1,cols/2-1)>0)
    sign=1;
else
    sign=-1;
end

wave = sign.*cos(omega*t);

err = y(:,cols/2-1)-wave;
SSE = err'*err/numCycles;
fprintf(1,'Error/cycle wrt COS: %.5e\n',SSE);
assignin('base','evals',evals*sqrt(time_rescaling_factor));
% 
% figure(3);
% 
% % [t_peaks,peeks] = findPeaks([t',y(:,cols/2-1)],[],[],1); %#ok<NASGU>
% % 
% % periods = t_peaks-t_peaks(1);
% % period_err = periods./periods(2); %#ok<NASGU>
% hold on;
% plot(t,wave,'g');
% xlabel('time (s)');
% ylabel('displacement (arb)')
% legend({'Beam tip displacement';'Peaks found';'Wave'});
% hold off;
end