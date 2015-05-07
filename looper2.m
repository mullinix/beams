% looper2
%%%%%%%%%%%%%%%
% this is the improved loop_omega that uses eigs instead of polyeig.

tic;
load matrices;
gamma=0.0;
delta = props.delta;
T = props.T;
Omega = gamma/T;
T_squared = props.T_squared;

%% set matrices for solution
X0 =(K+Omega^2*(Sigma-P));
X1 = 2*G*Omega*1i;
X2 = M;
%% adimensionalize
X1 = X1./T;
X2 = X2./T_squared;

% a=norm(X2);
% b=norm(X1);
% c=norm(X0);
% alpha = sqrt(c/a);
% beta = 2/(c+alpha*b);
% X0 = X0*beta;
% X1 = X1*alpha*beta;
% X2 = X2*alpha^2*beta;
% Z = zeros(size(K));
% A =[Z,-X0; X2, Z];
% B =[X2, X1; Z, X2];
mx_Cell={-X0,-X1,X2};
n = length(mx_Cell{1,1});
p=length(mx_Cell)-1;
A = eye(n*p);
A(1:n,1:n) = mx_Cell{1,1};
if p == 0
    B = eye(n);
    p = 1;
else
    B = diag(ones(n*(p-1),1),-n);
    j = 1:n;
    for k = 1:p
        B(1:n,j) = - mx_Cell{1,k+1};
        j = j+n;
    end
end
%[shape,omega]= polyeig(-X0,-X1,X2);
[shape,omega]= eigs(A,B);
omega = diag(omega);
% omega = evals;
evals = omega;
% omega = omega*(alpha);
omega = abs(omega);
[omega,ix]= sort(omega,'ascend');
shape = shape(:,ix);
evals = evals(ix);
freqs = abs(evals);
%  savename = sprintf('%d-elts.mat',nelts);
% save(savename,'omega','shape');
tf = toc;
err = sum(abs(1-(abs(real(evals))+abs(imag(evals)))./abs(evals)));
display(err);
fprintf(1,'eigenvalue solve time: %.5e\n',tf);