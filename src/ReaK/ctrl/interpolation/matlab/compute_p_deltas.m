function [p_deltas,A] = compute_p_deltas(beta, coefs, dt_amax)

p_deltas = zeros(length(beta),2);
A = zeros(length(beta),2);

A(:,1) = sqrt( coefs(1) * coefs(1) - (beta * coefs(5)) .* (2.0 * coefs(1) * coefs(2) - beta * coefs(5)) );
t_0 = coefs(1) * coefs(2) / coefs(5);
A(:,2) = sqrt( coefs(3) * coefs(3) - (beta * coefs(5)) .* (2.0 * coefs(3) * coefs(4) - beta * coefs(5)) );
t_1 = coefs(3) * coefs(4) / coefs(5);
    
for I=1:length(beta)
    if(A(I,1) < dt_amax)
        p_deltas(I,1) = 0.5 * sqrt( dt_amax * A(I,1) ) * ( t_0 + 3.0 * beta(I) );
    else
        p_deltas(I,1) = 0.5 * ((dt_amax + A(I,1)) * (beta(I) + t_0) + dt_amax * dt_amax / A(I,1) * (beta(I) - t_0));
    end

    if(A(I,2) < dt_amax)
        p_deltas(I,2) = 0.5 * sqrt( dt_amax * A(I,2) ) * (t_1 + 3.0 * beta(I));
    else
        p_deltas(I,2) = 0.5 * ((dt_amax + A(I,2)) * (beta(I) + t_1) + dt_amax * dt_amax / A(I,2) * (beta(I) - t_1));
    end
end





