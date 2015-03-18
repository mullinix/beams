err = zeros(length(SVevals),1);
for i=1:length(SVevals)
    res = (SVX2*SVevals(i)^2-SVX1*SVevals(i)-SVX0)*SVshape(:,i);
    err(i) = norm(res);
%     if err(i)>1e-8
%         fprintf(1,'Error: Eigen value/shape %d is not a solution!\n',i);
%         fprintf(1,'Residual error: %.5e\n', err(i));
%         break;
%     end
end
fprintf(1,'Max error: %.5e\n',max(err));
large.logic = err>1e-8;
if(sum(large.logic))
    large.err = err(large.logic);
    idx = (1:length(evals))';
    large.idx = idx(large.logic);
    display([large.err,large.idx]);
end
    