function dR = dtrapRes(fw,params)
    dt=params{1};
    n=length(fw);
    dR = (eye(n)+dt/2)\eye(n);
end