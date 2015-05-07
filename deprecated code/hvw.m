function retval = hvw(x,v,w)
retval = zeros(1,length(x));
for i=1:length(x)-1
el = x(i+1)-x(i);
idxi=2*(i-1)+1;
idx=idxi:idxi+3;
y1=v(idx)';
y2=w(idx)';
A = [1, 0,    0,      0;
     0, 1,    0,      0;
     1,el, el^2,   el^3;
     0, 1, 2*el, 3*el^3];
c1=A\y1;
c2=A\y2;
cprime1 = polyder(c1);
cprime2 = polyder(c2);
p1 = polyint(conv(cprime1,cprime1));
p2 = polyint(conv(cprime2,cprime2));
retval(i) = 0.5*(polyval(p1,el)+polyval(p2,el));
end
end