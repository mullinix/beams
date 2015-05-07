function omega_out = calc_omega(t,epsilon,Omega)
% function omega_out = calc_omega(t,epsilon,Omega)
% calculates the Omega function for time integration.

%% perturbations of rotation frequency. add k to input.
% omega_out  = [Omega.*(1+epsilon.*sin(k.*Omega.*t));
%               k.*Omega.^2.*epsilon.*cos((k.*Omega.*t))];

%% reproduce paper: smooth acceleration/deceleration of omega
% omega_out=zeros(2,length(t));
if(t<4*pi/Omega)
    omega_out  = [Omega*sin(Omega/8*t);
                  Omega^2/8*cos(Omega/8*t)];
% elseif(t<40)
%     omega_out=[Omega;0];
else
    omega_out = [Omega; 0];%[50-t(t>=40)+(5/pi)*sin(pi/5*t(t>=40));
                    %-1+cos(pi/5*t(t>=40))];
end

%% Constant rotation
% omega_out = [Omega; 0];
          
end