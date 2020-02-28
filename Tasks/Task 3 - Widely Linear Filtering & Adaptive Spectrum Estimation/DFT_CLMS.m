function [w_n,y_hat,e_n] = DFT_CLMS(y,x_n,mu,gamma)
%Implementation of the DFT-CLMS algorithm
%filter weights adapted using CLMS algorithm, and will converge to DFT
%solution if learning rate = 1
    N=length(y);
    y_hat = complex(zeros(N,1));
    w_n = complex(zeros(N,N));
    e_n = complex(zeros(N,1));
    for i =1:N
        x_i = x_n(:,i);
        y_hat(i,1) = w_n(:,i)' * x_i;
        e_n(i,1) = y(i,1) - y_hat(i,1);
        w_n(:,i+1) = (1-mu*gamma)*w_n(:,i) + mu * conj(e_n(i,1)) * x_i;
    end
    w_n = w_n(:,2:end);
end

