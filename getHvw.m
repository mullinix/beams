function Hret = getHvw(x0,x1)
%% Instantiate return variable
Hret = zeros(4,1);
%% Create inperpolatory shape function
el = x1-x0;
H = [1,0,-3/el^2,2/el^3;
     0,1,-2/el,1/el^2;
     0,0,3/el^2,-2/el^3;
	 0,0,-1/el,1/el^2];
%% Solve stretch value addition
for i=1:4
    Hp = polyder(H(i,:)); % dv/dtheta
    Hpi = polyint(conv(Hp,Hp));% integral of (dv/dtheta)^2
    Hret(i) = polyval(Hpi,x1)-polyval(Hpi,0);% Solution of definite integral
end

end