function [xhat, error, Weight, muAdaptive] = GASSlms(x, input, order, rho,alpha ,type)
    N = length(x);
    Weight = zeros(order,N);
    error = zeros(1,N);
    xhat = zeros(1,N);
    psi = zeros(order,N);
    muAdaptive = zeros(1,N);
    for i=order+1:N
        xPrevious = input(i:-1:i-order+1);
        xhat(i) = Weight(:,i)'*xPrevious';
        error(i) = x(i)-xhat(i);
        Weight(:,i+1) = Weight(:,i) + muAdaptive(i)*error(i)*xPrevious';
        muAdaptive(i+1) = muAdaptive(i) + rho*error(i)*xPrevious*psi(:,i);
        switch type
            case 'Benveniste'
                psi(:,i+1) = (eye(order)-muAdaptive(i)*xPrevious'*xPrevious)*psi(:,i) + error(i)*xPrevious';
            case 'Farhang'
                psi(:,i+1) = alpha*psi(:,i) + error(i)*xPrevious';
            case 'Mathews'
                psi(:,i+1) = error(i)*xPrevious';
        end
    end
    Weight = Weight(:,2:end);
end
