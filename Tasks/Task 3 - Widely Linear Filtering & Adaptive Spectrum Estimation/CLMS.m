function [h_n,y_hat,e_n] = CLMS(y,x,mu,numCoeffs)
    N = length(x);
    y_hat = complex(zeros(1, N));
    h_n = complex(zeros(numCoeffs, N));
    e_n = complex(zeros(1, N));
    for i = numCoeffs+1:N
        x_n = x(i:-1:i-numCoeffs+1);
        y_hat(i) = h_n(:, i)' * x_n;
        e_n(i) = y(i) - y_hat(i);
        h_n(:,i+1) = h_n(:, i) + mu * conj(e_n(i)) * x_n;
    end
    h_n = h_n(:,2:end);
end

