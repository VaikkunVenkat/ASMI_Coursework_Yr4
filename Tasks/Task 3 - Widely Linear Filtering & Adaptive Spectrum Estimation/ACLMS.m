function [h_n,g_n,xHat,e_n] = ACLMS(y,x,mu,numCoeffs)
    N = length(x);
    h_n = complex(zeros(numCoeffs,N));
    g_n = complex(zeros(numCoeffs,N));
    e_n = complex(zeros(1,N));
    xHat = complex(zeros(1,N));
    for i=numCoeffs+1:N
        x_n = x(i:-1:i-numCoeffs+1);
        xHat(i) = h_n(:,i)'*x_n + g_n(:,i)'*conj(x_n);
        e_n(i) = y(i)-xHat(i);
        h_n(:,i+1) = h_n(:,i) + mu*conj(e_n(i))*x_n;
        g_n(:,i+1) = g_n(:,i) + mu*conj(e_n(i))*conj(x_n);
    end
    h_n = h_n(:,2:end);
    g_n = g_n(:,2:end);
end

