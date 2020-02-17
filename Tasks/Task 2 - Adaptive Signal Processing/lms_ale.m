function [error,predicted_x,weights] = lms_ale(s_realisation,mu,delay,M)
    N = length(s_realisation);
    predicted_x = zeros(N, 1);
    weights = zeros(M, N);
    error = zeros(N, 1);
    for i = delay+M:N
        u = s_realisation(i-delay:-1:i-delay-M+1);
        predicted_x(i) = weights(:, i)' * u;
        error(i) = s_realisation(i) - predicted_x(i);
        weights(:, i+1) = weights(:, i) + mu * error(i) * u;
    end
    weights = weights(:,2:end);
end

