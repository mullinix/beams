function retval = hvw_analysis(x,v,w)
% function retval = hvw_analysis(x,v,w)
% evaluates and returns h_v + h_w for the purpose of extracting u from the
% stretch coordinate s.

retval = zeros(size(v,1),length(x));
for xidx=1:length(x)-1
    el = x(xidx+1)-x(xidx);
    idx=2*(xidx-1)+1:2*(xidx+1);
    y1=v(:,idx)';
    y2=w(:,idx)';
    A = [1, 0,    0,      0;
         0, 1,    0,      0;
         1,el, el^2,   el^3;
         0, 1, 2*el, 3*el^3];
    c1=A\y1;
    c2=A\y2;
    for i=1:length(retval(:,1))
        cprime1 = polyder(c1);
        cprime2 = polyder(c2);
        p1 = polyint(conv(cprime1,cprime1));
        p2 = polyint(conv(cprime2,cprime2));
        retval(i,xidx) = 0.5*(polyval(p1,el)+polyval(p2,el));
    end
end
end