function[w,R] = newtRahp(fn,df,w,varargin)
%fparams,dfparams)
if(nargin>3 && ~isempty(varargin(1)))
    myFunc=@(x) fn(x,varargin{1});
else
    myFunc=@(x) fn(x);
end
if(nargin>4 && ~isempty(varargin(2)))
    mydFunc=@(x) df(x,varargin{2});
else
    mydFunc=@(x) df(x);
end
fw=myFunc(w);
dR=mydFunc(w);
R=fw;
i=0;
maxiter=1e4;
tol=1e-8;
while(abs(R)>tol)
    i=i+1;
    if(abs(dR)<1e-15)
        dR=1e-15;
    end
    dw=dR\(-R);
    w=w+dw;
    fw=myFunc(w);
    dR=mydFunc(w);
    R=fw;
    if(i>maxiter)
        error('Exceeded maxiter=%d!',maxiter);
    end
    fprintf(1,'iter: %d, R: %.2e\n',i,R);
end

end



