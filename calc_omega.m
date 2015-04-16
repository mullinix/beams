function omega_out = calc_omega(t,k,epsilon,Omega)

omega_out  = [Omega.*(1+epsilon.*cos(k.*Omega.*t));
              -k.*Omega.^2.*epsilon.*sin((k.*Omega.*t))];
end