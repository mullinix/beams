function [t,y] = dataset_intersection(t1,x1,t2,x2)
% function [t,y] = dataset_intersection(t1,x1,t2,x2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the intersection of two datasets x1(t1), x2(t2).
% Uses linear interpolation.
%
% We need t1==t2 OR min(t1)==min(t2) AND max(t1)==max(t2)!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% James Mullinix, 4/23/2015

if(t1~=t2) % if not the same indep var, spline!
    if(min(t1)~=min(t2)||max(t1)~=max(t2))
        error('Cannot spline! Independent variable bounds must agree!');
    end
    l1=length(t1); l2=length(t2);
    if(l1~=l2)
        splen=max([l1,l2])*2;
        T=linspace(min(t1),max(t1),splen);
        x1 = spline(t1,x1,T);
        x2 = spline(t2,x2,T);
    end
else % or its the same, and set T.
    T=t1;
end
% find the sign of the difference function. -1 or 1. then take the
% derivative of this function; it's discrete, so we will either have a
% slope that is zero, or we will have crossed zero.
idx = (diff(sign(x1-x2))~=0);
% two issues: diff has one less index, and we need to get the points before
% and after the intersection for the calculation.
% 1 is before intersection, 2 is after
idx1 = [idx,1==0];
idx2 = [1==0,idx];
% find difference of independent var
dT = T(idx2)-T(idx1);
% find slope of line
m1 = (x1(idx2)-x1(idx1))./dT;
m2 = (x2(idx2)-x2(idx1))./dT;
% find intercept of line
b1 = x1(idx1)-m1.*T(idx1);
b2 = x2(idx1)-m2.*T(idx1);
% get all of the independent var intersection pts
t = (b2-b1)./(m1-m2);
% calculate dependent var intersection pts
y = m1.*t+b1;

end