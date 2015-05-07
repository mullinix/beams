function [out,idx] = uniquetol(A,tol)
% function [out,idx] = uniquetol(A,tol)
% Returns unique values in a dataset (unique to within a tolerance)
tolinv=1/tol;
out = unique(round(A.*tolinv)./tolinv);
idx=zeros(length(out),1);
for i=1:length(out)
    idx(i)=find(abs(A-out(i))<tol,1,'first');
end
out = A(idx);
    
end