function[x0,R] = newtRaph(fn,df,x0,varargin)
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
fx=myFunc(x0);
dR=mydFunc(x0);
R=fx;
iter=0;
maxiter=1e4;
tol=1e-8;
while(norm(R)>tol)
    iter=iter+1;
    if(abs(dR)<1e-15)
        dR=1e-15;
    end
    dx=dR\(-R);
    x0=x0+dx;
    fx=myFunc(x0);
    dR=mydFunc(x0);
    R=fx;
    if(iter>maxiter)
        error('Exceeded maxiter=%d!',maxiter);
    end
%     fprintf(1,'iter: %d, R: %.2e\n',iter,R);
end

end



