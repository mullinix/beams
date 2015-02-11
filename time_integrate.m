%% Setup integration, integrate in time
mode_to_try=1;
omega = evals(mode_to_try)*(2*pi);% corresponds to Hz -> radians
t0 = 0;
numCycles=2;
tf = numCycles*(1/omega)*2*pi;% end after four full cycles
numpts = 1e3;
thyme = linspace(t0,tf,numpts);
X = shape(2*end/3+3:5:end,end-mode_to_try+1);
N = length(X);
y0 = shape(2*end/3+3:5:end,end-mode_to_try+1)./max(abs(X));
Y0 = zeros(3*N*2,1);

for idx=1:num_beams
    idx_range = (1:num_dofs)*(idx*num_dofs-1)+1
    Y0(idx_range) = y0;
end

A = -X2\X1;
B = -X2\X0;
% options = odeset('RelTol',1e-8,'Mass',MM_big);
[t,y] = ode23t(@(t,y) coupled_diff(t,y,A,B), thyme, y0);%,options);