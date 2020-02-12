% Define spatial signal
dx   = 1e-3;
Fsx  = 1/dx; % sampling frequency
X    = 0:dx:2;
L    = length(X);     % = 2001
Lfft = 2^nextpow2(L); % = 2048
dfx  = 1/(Lfft*dx);
fsig = 10;
Z   = sin(2*pi*fsig*X);
Z_dirac = [1, zeros(1,L-1)];

% Compute PSD from FFT
fourier_image  = fftshift(fft(Z,Lfft)/Lfft);
PSDfromFFT     = abs(fourier_image).^2;
freq           = -Fsx/2:dfx:Fsx/2-dfx; % Frequency axis
% Compute PSD from Autocorrelation
Rxx = xcorr(Z,'biased');
lag = (-(L-1):1:(L-1))*dx;
PSDfromRxx = abs(fftshift(fft(Rxx,Lfft)/Lfft));
% Parseval sums
fprintf('\nParseval sums : \n\n');
fprintf('Spatial sum = %f\n', sum(X(:).^2));
fprintf('Frequency sum (from FFT) = %f\n', sum(PSDfromFFT(:))*Lfft);
fprintf('Frequency sum (from ACF) = %f\n', sum(PSDfromRxx(:))*Lfft);
% Plots
figure;
subplot(3,1,1); plot(X,Z,'-b'); title('Signal sin(2\pi10n)'); xlabel('Time[s]');ylabel('Amplitude') ; grid on;
subplot(3,1,2); plot(lag,Rxx,'-r'); title('Autocorrelation of sin(2\pi10n)');xlabel('lag');ylabel('r(k)') ; grid on;
subplot(3,1,3);
plot(freq,PSDfromFFT,'r'); title('PSD based on FFT'); hold on;
plot(freq, PSDfromRxx, 'b'); title('PSD based on autocorrelation');
legend('PSD based on FFT','PSD based on autocorrelation'); hold off;
xlabel('Frequency [Hz]');ylabel('|P(\omega)|');grid on;