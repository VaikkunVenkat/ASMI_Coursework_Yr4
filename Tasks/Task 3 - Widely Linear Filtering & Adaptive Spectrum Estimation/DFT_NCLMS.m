function [w_n, yHat, error] = DFT_NCLMS(y,x_n,mu)
    N = length(x_n);
    w_n = complex(zeros(N,N));
    error = complex(zeros(N,1));
    yHat = complex(zeros(N,1));
    for i=1:N
        x_i = x_n(:,i);
        yHat(i,1) = w_n(:,i)' * x_i;
        error(i,1) = y(i,1)-yHat(i,1);
        w_n(:,i+1) = w_n(:,i) + (mu/(x_i'*x_i) + 0.002)*conj(error(i))*x_i;
    end
    w_n = w_n(:,2:end);
end
