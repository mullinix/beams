clear all; close all;

N=100;
% rng(1);
M1 = (rand(N)+1)*1e2;% make the matrix from 10^0 -> 10^-5
M1=M1*M1';% matrix is now symmetric
% rng(2);
M2 = (rand(N)+1)*1e6;
M2=M2*M2';
ind = rand(N)>0.5;
M = M1.*ind+M2.*(1-ind);

% rng(3);
K1 = (rand(N)+1)*1e-5;
K1=K1*K1';
% rng(4);
K2 = (rand(N)+1)*1e+8;
K2=K2*K2';
ind = rand(N)>0.5;
K = K1.*ind+K2.*(1-ind);

%K = ones(N);

MM_min = real(floor(log10(min(abs(M(M~=0)))))); 
MM_max = real(floor(log10(max(abs(M(M~=0))))));
fprintf(1,'Pre-Scale of M: [10^%d->10^%d]\n',MM_min,MM_max);
KK_min = real(floor(log10(min(abs(K(K~=0)))))); 
KK_max = real(floor(log10(max(abs(K(K~=0))))));
fprintf(1,'Pre-Scale of K: [10^%d->10^%d]\n',KK_min,KK_max);
% K=K+K';
% M=M+M';

KK_mid = (KK_max-KK_min)/2+KK_min
rescale_factor = 10^(KK_mid-MM_min)

M = M.*rescale_factor;

[shape,evals] = eig(K,M);
evals=diag(evals);

SSE1=zeros(1,N);

for ii=1:N
 result = M*evals(ii)*shape(:,ii)-K*shape(:,ii);
 SSE1(ii) = result'*result;
end

plot(1:N,SSE1,'bo-')


[shape,evals] = eig(K,M);
evals=diag(evals);
