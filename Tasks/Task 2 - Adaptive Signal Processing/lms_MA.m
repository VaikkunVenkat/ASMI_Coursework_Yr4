function [error , xHat,weights] = lms_MA(x,input,order,mu,gamma)
    N = length(x);
    xHat = zeros(N, 1);
    weights = zeros(order, N);
    error = zeros(N, 1);
    for i = order+1:N
        input_past = input(i:-1:i-order+1);
        xHat(i) = weights(:, i)' * input_past;
        error(i) = x(i) - xHat(i);
        weights(:, i+1) = (1-mu*gamma)*weights(:, i) + mu * error(i) * input_past;
    end
    weights = weights(:,2:end);
end

