rng(1);
if(~exist('N','var'))
    N=10;
end
X0 = rand(N);X0 = X0+X0';X0=X0*1e8;
X1 = rand(N);X1 = X1+X1';X1=X1*1e8;
X2 = rand(N);X2 = X2+X2';X2=X2*1e8;

[Ve,Xe] = eig(X0);Xe=diag(Xe);
[Veg,Xeg] = eig(X0,-X2);Xeg=diag(Xeg);
[Vp,Xp] = polyeig(X0,X1,X2);


sum_eg = 0;
sum_p = 0;
for i=1:N
    err = (X0+X2*Xeg(i))*Veg(:,i);
    sum_eg = sum_eg + err'*err;
    lam = Xp(i);
    err = (X0+X1*0*lam+X2*lam^2)*Vp(:,i);
    sum_p = sum_p + err'*err;
    lam = Xp(i+N);
    err = (X0+X1*lam+X2*lam^2)*Vp(:,i+N);
    sum_p = sum_p + err'*err;
end
% for i=1:2*N
%     lam = Xp(i);
%     err = (X0+X1*0*lam+X2*lam^2)*Vp(:,i);
%     sum_p = sum_p + err'*err;
% end

fprintf(1,'Gnrl Eval Err: %1.3e\n',sum_eg);
fprintf(1,'Poly Eval Err: %1.3e\n',sum_p);

% Xp=Xp.^2;
% 
% plot(real(Xeg),imag(Xeg),'o',real(Xp),imag(Xp),'r+');hold off;
% legend({'General Evals'; 'Poly Evals'});
