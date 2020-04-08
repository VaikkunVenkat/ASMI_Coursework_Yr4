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

% Determine circularity distribution
figure;
plot(y,'b*');xlabel('y_{r}(n)');ylabel('y_{i}(n)');title('Complex Distribution of FM Signal');
grid on; grid minor; set(gca,'FontSize',18);
%%
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
surf(1:N, (w/pi)*(fs/2), H, 'LineStyle','none');
colormap winter;colorbar;
view(2);ylim([0 1000]);
xlabel('Time Index [n]');
ylabel('f(n) [Hz]');
set(gca, 'Fontsize', 18);
title('Time-Frequency Estimation of noisy FM signal, \mu = ' + string(mu));
grid on;grid minor;
%%
mu = [0.05,0.3,0.6];num_mu=length(mu);
a1_hat = zeros(num_mu,N);error_n = zeros(num_mu,N);


for i =1:length(mu)
    [a1_hat(i,:),~,error_n(i,:)] = CLMS(y,x,mu(i),order);
end

figure;plot(abs(a1_hat(3,:)),'r','LineWidth',2); hold on; plot(abs(a1_hat(2,:)),'g','LineWidth',2);
plot(abs(a1_hat(1,:)),'b','LineWidth',2)
hold off;
xlabel('Time Index [n]');ylabel('$$|\hat{a}_1(n)|$$','Interpreter','Latex');title('Trajectory of a_1(n) on FM Signal using CLMS');
grid on; grid minor;set(gca,'FontSize',18); legend('\mu=0.6','\mu=0.3','\mu=0.05');

figure;
plot(10*log10(abs(error_n(3,:)).^2),'r','LineWidth',2); hold on; plot(10*log10(abs(error_n(2,:)).^2),'g','LineWidth',2);plot(10*log10(abs(error_n(1,:).^2)),'b','LineWidth',2);
hold off;
xlabel('Time Index [n]');ylabel('Error Power [dB]');title('Learning curves for TF estimation of FM signal using CLMS');
grid on; grid minor;set(gca,'FontSize',18); legend('\mu=0.6','\mu=0.3','\mu=0.05');

figure;plot(conj(a1_hat(3,:)),'r*'); hold on; plot(conj(a1_hat(2,:)),'g*');
plot(conj(a1_hat(1,:)),'b*')
hold off;xlabel('Re(a_1(n))');ylabel('Im(a_1(n))');title('Trajectory of a_1(n) on FM Signal using CLMS');
grid on; grid minor;set(gca,'FontSize',18); legend('\mu=0.6','\mu=0.3','\mu=0.05');

