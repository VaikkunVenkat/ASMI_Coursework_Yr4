% Generate Time-Changing Frequency and Time
N = 1500;
f_n = frequency_n();
phi = integrateFrequency(f_n);
figure;
plot(1:1500,f_n,'b','LineWidth',2); hold on;
xlabel('Time Index [n]');ylabel('f(n)');
title('Frequency of signal over time');grid on;grid minor;set(gca,'FontSize',18);
figure;
plot(1:1500,phi,'r','LineWidth',2); hold on;
xlabel('Time Index [n]');ylabel('\phi(n)');
title('Phase of signal over time');grid on;grid minor;set(gca,'FontSize',18);



% Generate FM Signal
noiseVariance = 0.05;
fs = 5000;
eta = wgn(1500,1,10*log10(noiseVariance),'complex');
signal = exp(1i * 2*pi/fs .* phi);
y = signal + eta;

% determine AR(1) coefficient and static power spectrum
[a,e] = aryule(y,1);
[h,w] = freqz(1,a,N);
TheoreticalPSD = abs(h).^2;
FFT_y = fft(y);
PSDy = (1/(2*pi*N)) * abs(FFT_y).^2;
freq = 0:(2*pi)/N:2*pi-(2*pi)/N;
EmpPSD = abs(fftshift(FFT_y)).^2;
figure;
plot(w/pi,10*log10(TheoreticalPSD),'b','LineWidth',2); hold on;
plot(freq/pi,10*log10(PSDy),'r','LineWidth',2);
xlabel('Normalised Frequency (\times\pi rad/sample)');ylabel('Power/Frequency (dB/rad/sample)');
grid on; grid minor; title('Theoretical AR(1) and Empirical PSD of FM signal');
set(gca,'FontSize',18);legend('Theoretical AR(1)','Empirical');
%%
% Determine Empirical and AR(1) power spectra in each of the three regions
% Linear Section
y_flat = y(1:500);
y_linear = y(501:1000);
y_quadratic = y(1001:1500);
y_signals = [y_flat , y_linear , y_quadratic];
N_section = 500;
for i =1:3
    y_section = y_signals(:,i);
    [a,e] = aryule(y_section,1);
    [h,w] = freqz(1,a,N_section);
    TheoreticalPSD = abs(h).^2;
    FFT_y = fft(y_section);
    PSDy = (1/(2*pi*N_section)) * abs(FFT_y).^2;
    freq = 0:(2*pi)/N_section:2*pi-(2*pi)/N_section;
    EmpPSD = abs(fftshift(FFT_y)).^2;
    figure;
    plot(w/pi,10*log10(TheoreticalPSD),'b','LineWidth',2); hold on;
    plot(freq/pi,10*log10(PSDy),'r','LineWidth',2);
    xlabel('Normalised Frequency (\times\pi rad/sample)');ylabel('PSD [dB]');
    grid on; grid minor; title('AR(1) and Empiral PSD FM signal estimate, section ' + string(i));
    set(gca,'FontSize',18);legend('Theoretical AR(1)','Empirical');
end    
%%
% CLMS based spectrum estimation.
L = 1024;
a1_hat = zeros(1,N); H = zeros(L,N);mu=0.6;order=1;
x = [0;y(1:end-1)];
[a1_hat(1,:),~,~] = CLMS(y,x,mu,order);
for n = 1:N
    [h ,w]= freqz(1 , [1, -conj(a1_hat(1,n))], L);
    H(:, n) = abs(h).^2;
end
medianH = 50*median(median(H));
H(H > medianH) = medianH;
figure;
surf(1:N, (w/pi)*1000, H, 'LineStyle','none');
colormap spring;colorbar;
view(2);
xlabel('Time Index [n]');
ylabel('f(n) [Hz]');
set(gca, 'Fontsize', 18);
title('Time-Frequency Estimation of noisy FM signal, \mu = ' + string(mu));
grid on;grid minor;

