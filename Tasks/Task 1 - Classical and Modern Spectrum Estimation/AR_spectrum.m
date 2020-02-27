N=1000;
noise = randn(N,1);
a = [2.76 -3.81 2.65 -0.92];
a = [1 -a];
x = filter(1,a,noise);
x = x(500:end);noise = noise(500:end); N_length = length(x);
figure(1);
hold on;plot(x,'r','LineWidth',2);plot(noise,'k','LineWidth',2);xlabel('Time Index');ylabel('Amplitude');title('x[n] = 2.76x[n-1]-3.81x[n-2]+2.65x[n-3]-0.92x[n-4] + w[n]');
legend('AR(4) Process','N(0,1) noise')
grid on; grid minor;
set(gca,'FontSize',18)

figure(2);
[H,F] = freqz(1,a,[],1);
orders = [2 5 20];

for i = 1:length(orders)
    [Pxx,F_yulear] = pyulear(x,orders(i),1024,1); hold on;
    plot(F_yulear,10*log10(Pxx),'LineWidth',2);
end
plot(F,20*log10(abs(H)),'LineWidth',2); xlabel('Normalized Frequency (Hz)'); ylabel('PSD (dB/Hz)');title('Actual and Estimated PSD of AR(4) process, N = ' + string(N));
grid on ; grid minor; legend('AR(2)','AR(5)','AR(18)','Theoretical');
set(gca,'FontSize',18)

figure(3);
[arcoefs,E,K] = aryule(x,20);
pacf=-K;
stem(pacf,'filled','k');
xlabel('lag');ylabel('Partial ACF'); title('Partial Autocorrelation Sequence');xlim([1 20]);
uconf = 0.2;
lconf = -uconf; hold on;
plot([1 20],[1 1]'*[lconf uconf],'r--');grid on; grid minor;
legend('PACF of AR(4) process')
set(gca,'FontSize',18)



%%
figure(3);
xacf = xcorr(x,'unbiased'); % Autocorrelation of AR(4) process
EmpPSD = abs(fftshift(fft(xacf)));
freq = linspace(-0.5,0.5,2*N_length);
plot(freq, 10*log10(EmpPSD),'LineWidth',2); xlabel('Normalized Frequency (Hz)');ylabel('PSD (dB/Hz)');title('Empircal PSD Based on Data');
grid on; grid minor; set(gca, 'FontSize',18);legend('AR(4) Empirical PSD');



