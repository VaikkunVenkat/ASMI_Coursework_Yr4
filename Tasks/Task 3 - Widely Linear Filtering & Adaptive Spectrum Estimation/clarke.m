function ClarkeVoltage = clarke(V_amp,V_Phi,fsig,fsamp,nsamp,phi)
%Determine the general Complex-valued representation of the clarke
%transform of the three phases
A = sqrt(6)/6 * (V_amp(1) + V_amp(2)*exp(1i * (V_Phi(1))) + V_amp(3)*exp(1i * V_Phi(2)));
B = sqrt(6)/6 * (V_amp(1) + V_amp(2)*exp(-1i * (V_Phi(1) + 2*pi/3)) + V_amp(3)*exp(-1i * (V_Phi(2) - 2*pi/3)));
ClarkeVoltage = A*exp(1i*(2*pi*(fsig/fsamp)*1:nsamp + phi)) + B*exp(-1i*(2*pi*(fsig/fsamp)*1:nsamp + phi));
end

