function [w_init,WeightPath,Error,Rp] = pre_train(y,mu,order,a)
    N=length(y);nEpochs = 100;
    w_n = zeros(order+1,N);
    Error = zeros(nEpochs,1);
    WeightPath = zeros(order+1,nEpochs);
    Rp = zeros(nEpochs,1);
    for j = 1:nEpochs
        e_n = zeros(1,N);
        y_hat = zeros(1,N);
        for i = order+1:N
            inputPrevious = [1;y(i-1:-1:i-order)];
            y_hat(1,i) = a*tanh(w_n(:,i)' * inputPrevious);
            e_n(1,i) = y(i) - y_hat(1,i);
            w_n(:,i+1) = w_n(:,i) + mu*e_n(1,i)* inputPrevious;
        end
    Error(j,1) = mean(e_n(1,:));
    WeightPath(:,j) = w_n(:,end);
    Rp(j,1) = 10*log10(var(y_hat)/var(e_n));
    w_n = w_n(:,2:end);
    w_temp = w_n(:,end);
    w_n = zeros(order+1,N);
    w_n(:,order+1) = w_temp;
    end
    w_init = w_temp;
end