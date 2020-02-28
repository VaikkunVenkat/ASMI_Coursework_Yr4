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
% sine
fourier_image  = fftshift(fft(Z,Lfft)/Lfft);
PSDfromFFT     = abs(fourier_image).^2;
freq           = -Fsx/2:dfx:Fsx/2-dfx; % Frequency axis
% dirac
fourier_image_dirac = fftshift(fft(Z_dirac,Lfft)/Lfft);
PSDfromFFT_dirac = abs(fourier_image_dirac).^2;

% Compute PSD from Autocorrelation
% Sine
Rxx = xcorr(Z,'biased');
lag = (-(L-1):1:(L-1))*dx;
PSDfromRxx = abs(fftshift(fft(Rxx,Lfft)/Lfft));

%Dirac
Rxx_dirac = xcorr(Z_dirac , 'biased');
PSDfromRXX_dirac = abs(fftshift(fft(Rxx_dirac,Lfft)/Lfft));

% Parseval sums
fprintf('\nParseval sums : \n\n');
fprintf('Spatial sum = %f\n', sum(X(:).^2));
fprintf('Frequency sum (from FFT) = %f\n', sum(PSDfromFFT(:))*Lfft);
fprintf('Frequency sum (from ACF) = %f\n', sum(PSDfromRxx(:))*Lfft);
% Plots
figure(1);
plot(lag,Rxx,'-r'); title('Autocorrelation of Sinusoid');xlabel('lag');ylabel('r(k)') ; grid on;grid minor;
set(gca,'FontSize',18);
figure(2);
plot(freq,10*log10(PSDfromFFT),'r'); hold on;
plot(freq, 10*log10(PSDfromRxx), 'b'); title('PSD based on FFT and Autocorrelation of Sinusoid');
legend('PSD based on FFT','PSD based on autocorrelation'); hold off;
xlabel('Frequency [Hz]');ylabel('Power [dB]');grid on; grid minor;
set(gca,'FontSize',18);
figure(3);
plot(lag,Rxx_dirac,'-r'); title('Autocorrelation of Kronecker Delta');xlabel('lag');ylabel('r_{\delta}(k)') ; 
grid on; grid minor;
set(gca,'FontSize',18);
figure(4);
plot(freq,10*log10(PSDfromFFT_dirac),'r'); title('PSD based on FFT and Autocorrelation'); hold on;
plot(freq, 10*log10(PSDfromRXX_dirac), 'b');
legend('PSD based on FFT','PSD based on autocorrelation'); hold off;
xlabel('Frequency [Hz]');ylabel('|P(\omega)| [dB]');grid on; grid minor;
set(gca,'FontSize',18);