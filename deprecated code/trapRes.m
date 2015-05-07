function R = trapRes(w,params)
fn=params{1};% function handle
vn=params{2};% previous x value
tn=params{3};% previous time value
dt=params{4};% time step
params=params(5:end);
a=fn(tn+dt,w,params{:});
b=fn(tn,vn,params{:});
R = w-vn-(dt/2)*(a+b);

end