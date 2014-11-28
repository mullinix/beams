function linearization_error()
delta = [0,1,5];
gamma = 0:10;
lins = [1,3,4];
Errs = zeros(size(lins));

for i=1:length(lins)
    for j=1:length(delta)
        for k=1:length(gamma)
            Errs(i) = Errs(i)+timoshenko2D([delta(j);gamma(k);lins(i);1]);
        end
    end
    fprintf(1,'Error for linearization %d: %1.5e\n',lins(i),Errs(i));
end

end