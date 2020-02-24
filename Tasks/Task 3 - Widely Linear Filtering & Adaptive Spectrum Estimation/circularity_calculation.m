function [rho_z,nu] = circularity_calculation(x) % Output [circularity_quotient, circularity_coefficient]
    MeanSignal = mean(x,2);
    xR = real(MeanSignal); xI = imag(MeanSignal);
    Cov = cov(xR,xI);
    PseudoCovariance = Cov(1,1)-Cov(2,2) + 2*1i*Cov(1,2);
    Covariance = mean(abs(MeanSignal).^2);
    rho_z = PseudoCovariance / Covariance ;
    nu = abs(PseudoCovariance)/Covariance;
end

