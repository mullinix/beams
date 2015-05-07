% loop_omega
%%%%%%%%%%%%%%%%%%%%%%
% loop over omega rotation speed, grab the first 50 eigenvalues, save the
% results to disk. 
%% TODO:
% Replace the full eigenvalue solver (polyeigs) with the "eigs" command so
% that only 50 values are solved, not the entire system. Will be much
% faster!
%% run the code
global gamma;
for gamma=0:0.1:50
    build_coupling_stiffness
    data2save.omega = (omega(1:48));
    data2save.evals = (evals(1:48));
    data2save.freqs = (freqs(1:48));
    savename = sprintf('datastruct-%.1f.mat',gamma);
    save(savename,'data2save');
end