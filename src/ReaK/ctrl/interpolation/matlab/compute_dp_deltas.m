function dp_deltas = compute_dp_deltas(beta, coefs, dt_amax,A)


N = length(beta);
dp_deltas = zeros(N,2);

for I=1:N
    if (A(I,1) < dt_amax)
        dp_deltas(I,1) = beta(I) * coefs(5) * (beta(I) * coefs(5) - coefs(1) * coefs(2)) / dt_amax + coefs(1) * coefs(1) * (1 - coefs(2) * coefs(2)) / dt_amax + 0.5 * dt_amax;
    else
        dp_deltas(I,1) = beta(I) * coefs(5) * (beta(I) * coefs(5) - coefs(1) * coefs(2)) / A(I,1) + 0.5 * ((1.0 + dt_amax * dt_amax / (A(I,1) * A(I,1))) * coefs(1) * coefs(1) * (1 - coefs(2) * coefs(2)) / A(I,1) + dt_amax);
    end
    if(A(I,2) < dt_amax)
        dp_deltas(I,2) = beta(I) * coefs(5) * (beta(I) * coefs(5) - coefs(3) * coefs(4)) / dt_amax + coefs(3) * coefs(3) * (1 - coefs(4) * coefs(4)) / dt_amax + 0.5 * dt_amax;
    else
        dp_deltas(I,2) = beta(I) * coefs(5) * (beta(I) * coefs(5) - coefs(3) * coefs(4)) / A(I,2) + 0.5 * ((1.0 + dt_amax * dt_amax / (A(I,2) * A(I,2))) * coefs(3) * coefs(3) * (1 - coefs(4) * coefs(4)) / A(I,2) + dt_amax);
    end
end







