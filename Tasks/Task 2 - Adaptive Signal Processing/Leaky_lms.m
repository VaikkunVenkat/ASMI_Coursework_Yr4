function [error , weights] = Leaky_lms(x,order,mu,gamma)
    N = length(x);
    xHat = zeros(N, 1);
    weights = zeros(order, N+1);
    error = zeros(N, 1);
    for i = order+1:N
        xHat(i) = weights(:, i)' * x(i-1:-1:i-order);
        error(i) = x(i) - xHat(i);
        weights(:, i+1) = (1-mu*gamma)*weights(:, i) + mu * error(i) * x(i-1:-1:i-order);
    end
end
