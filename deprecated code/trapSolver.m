function Y = trapSolver(fn,x0,t,varargin)
params=varargin;
Y=zeros(length(t),length(x0));
Y(1,:)=feval(fn,t(1),x0,params{:});

for idx=2:length(t)
    dt=t(idx)-t(idx-1);
    x0=newtRaph(@trapRes,@dtrapRes,x0,{fn,x0,t(idx-1),dt,params{:}},{dt});
%     x0=fsolve(@(x)trapRes(x,{fn,x0,t(idx-1),dt,params{:}}),x0);
    Y(idx,:)=feval(fn,t(idx),x0,params{:});
end

end