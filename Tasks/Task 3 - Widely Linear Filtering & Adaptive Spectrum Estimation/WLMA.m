function y = WLMA(b1,b2,x)    
    N = length(x);
    y = complex(zeros(1,N));
    for i=2:N
        y(i) = x(i) + b1*x(i-1) + b2*conj(x(i-1));
    end
end

