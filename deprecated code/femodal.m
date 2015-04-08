function [Omega,Phi,ModF]=femodal(M,K,F)
% function [Omega,Phi,ModF]=femodal(M,K,F)
%     Inputs:
%         M = Mass Matrix
%         K = Stiffness Matrix
%         F = Forcing function
%     Output:
%         Omega = Natural freqs in ascending order
%         Phi = Modal matrix - each column corresponding to the eigenvector
%         ModF = Modal of input F matrix
% From: "FEM Using Matlab" by Kwon and Bang, CRC press, 1997, p.292

% book checks for size, then doesn't do anything with it...?
% [a,b] = size(M);
% [c,d] = size(F);
[V,D] = eig(K,M);
[~,k] = sort(diag(D));
V = V(:,k);
Factor = diag(V'*M*V);
Phi = V/(sqrt(diag(Factor)));
Omega = diag(sqrt(Phi'*K*Phi));% book has Vnorm, but Vnorm:=Phi
ModF = Phi'*F;% book has Vnorm, but Vnorm:=Phi
end