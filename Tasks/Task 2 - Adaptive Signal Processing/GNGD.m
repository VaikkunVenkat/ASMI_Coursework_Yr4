function [xhat, error, Weight, etaAdaptive] = GNGD(x, input, order, rho,mu)
    N = length(x);
    Weight = zeros(order,N);
    error = zeros(1,N);
    xhat = zeros(1,N);
    etaAdaptive = zeros(1,N);
    epsilon = ones(1,N);
    for i=order+1:N
        xPrevious = input(i:-1:i-order+1);
        xhat(i) = Weight(:,i)'*xPrevious';
        error(i) = x(i)-xhat(i);
        Weight(:,i+1) = Weight(:,i) + mu/((xPrevious'*xPrevious) + epsilon(i))*error(i)*xPrevious';
        %Weight(:,i+1) = Weight(:,i) + etaAdaptive(i)*error(i)*xPrevious';
        %etaAdaptive(i+1) = mu/((xPrevious'*xPrevious) + epsilon(i));
        Numerator = epsilon(i)*epsilon(i-1)*xPrevious'*xPrevious;
        Denominator = (epsilon(i-1) + (xPrevious'*xPrevious)).^2;
        epsilon(i+1) = epsilon(i) - rho*mu*(Numerator/Denominator);
    end
    Weight = Weight(:,2:end);
end
