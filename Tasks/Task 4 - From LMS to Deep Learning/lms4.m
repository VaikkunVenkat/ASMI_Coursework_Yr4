function [w_n,y_hat,e_n] = lms4(y,mu,order,dynamic,a)
%lms4 implements 4th order linear prediction using input data in y
    N=length(y);
    w_n = zeros(order,N);
    e_n = zeros(1,N);
    y_hat = zeros(1,N);
    for i = order+1:N
        inputPrevious = y(i-1:-1:i-order);
        if(dynamic == "dynamic")
            y_hat(1,i) = a*tanh(w_n(:,i)' * inputPrevious);
        elseif(dynamic == "normal")
            y_hat(1,i) = w_n(:,i)' * inputPrevious;
        end
        e_n(1,i) = y(i) - y_hat(1,i);
        w_n(:,i+1) = w_n(:,i) + mu*e_n(1,i)* inputPrevious;
    end
    w_n = w_n(:,2:end);
end
