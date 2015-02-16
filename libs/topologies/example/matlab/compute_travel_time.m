function dt = compute_travel_time(beta, p_deltas_0, norm_delta, coefs, dt_amax)


[p_deltas_1,A] = compute_p_deltas(beta,coefs,dt_amax);

N = length(beta);

for I=1:N
    if (A(I,1) < dt_amax)
        A(I,1) = sqrt(4.0 * A(I,1) * dt_amax);
    else
        A(I,1) = A(I,1) + dt_amax;
    end
    if(A(I,2) < dt_amax)
        A(I,2) = sqrt(4.0 * A(I,2) * dt_amax);
    else
        A(I,2) = A(I,2) + dt_amax;
    end
end
    
c = (norm_delta - p_deltas_1(:,1) - p_deltas_1(:,2) + p_deltas_0(1,1) + p_deltas_0(1,2));
    
dt = A(:,1) + A(:,2) + abs(c) ./ beta;





