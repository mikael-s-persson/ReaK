function ddt = compute_derivative_travel_time(beta, p_deltas_0, norm_delta, coefs, dt_amax)

N = length(beta);

[p_deltas_1,A] = compute_p_deltas(beta,coefs,dt_amax);
    
c = (norm_delta - p_deltas_1(:,1) - p_deltas_1(:,2) + p_deltas_0(1,1) + p_deltas_0(1,2));


p_deltas_dot_1 = compute_dp_deltas(beta,coefs,dt_amax,A);
    
c_dot = -p_deltas_dot_1(:,1) - p_deltas_dot_1(:,2);
    
ramp_time = zeros(N,2);
for I=1:N
    if (A(I,1) < dt_amax)
        ramp_time(I,1) = coefs(5) * (coefs(5) * beta(I) - coefs(1) * coefs(2)) / dt_amax;
    else
        ramp_time(I,1) = coefs(5) * (coefs(5) * beta(I) - coefs(1) * coefs(2)) / A(I,1);
    end
    if(A(I,2) < dt_amax)
        ramp_time(I,2) = coefs(5) * (coefs(5) * beta(I) - coefs(3) * coefs(4)) / dt_amax;
    else
        ramp_time(I,2) = coefs(5) * (coefs(5) * beta(I) - coefs(3) * coefs(4)) / A(I,2);
    end
end
    
ddt = ramp_time(:,1) + ramp_time(:,2) - abs(c) ./ (beta .* beta) + (sign(c) .* c_dot) ./ beta;




