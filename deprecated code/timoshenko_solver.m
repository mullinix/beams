function timoshenko_solver()

params = [1,1];

[x,f] = fsolve(@timoshenko2D,params);

fh = fopen('timo_solver_results.txt','w');
fprintf(fh,'error: %1.5e\n',f);
fprintf(fh,'b,rho/L: %1.5e\n',x);

end