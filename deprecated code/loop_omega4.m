global gamma;
for gamma=46:0.1:50
    build_coupling_stiffness
    data2save.omega = (omega(1:48));
    data2save.evals = (evals(1:48));
    data2save.freqs = (freqs(1:48));
    savename = sprintf('datastruct-%.1f.mat',gamma);
    save(savename,'data2save');
end