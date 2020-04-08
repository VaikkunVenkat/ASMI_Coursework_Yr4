function [w_n,y_hat,e_n] = lms4_dynamical_bias(y,mu,order,a,w_n)
    N=length(y);
    e_n = zeros(1,N);
    y_hat = zeros(1,N);
    for i = order+1:N
        inputPrevious = [1;y(i-1:-1:i-order)];
        y_hat(1,i) = a*tanh(w_n(:,i)' * inputPrevious);
        e_n(1,i) = y(i) - y_hat(1,i);
        w_n(:,i+1) = w_n(:,i) + mu*e_n(1,i)* inputPrevious;
    end
    w_n = w_n(:,2:end);
end