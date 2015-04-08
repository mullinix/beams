mode = 1;
% evec_test = shape(:,mode);
% len = length(evals)/6;
% II = eye(len);
% ZZ = zeros(len);
% Rho = [ZZ,ZZ,II;
%        II,ZZ,ZZ;
%        ZZ,II,ZZ];% S3 permutation matrix
% perm = 2;
% spy(Rho^perm)
% evec_perm1 = (Rho^(perm))*evec_test;
% res = (X2*evals(mode)^2-X1*evals(mode)-X0)*evec_perm1;
% norm(res)