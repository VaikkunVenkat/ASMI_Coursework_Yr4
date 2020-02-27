%% Periodogram based estimation to sunspots timeseries
load sunspot.dat
ZurichNumber=sunspot(:,2);
Year = sunspot(:,1);
N = length(ZurichNumber);
zeromeanDetrendZurich = detrend(ZurichNumber - mean(ZurichNumber));
ZurichLog = log(ZurichNumber+eps) - mean(log(ZurichNumber+eps));
figure(1);
plot(Year,ZurichNumber ,'g','LineWidth',2);
hold on;
plot(Year,zeromeanDetrendZurich,'b','LineWidth',2);
plot(Year, ZurichLog, 'r','LineWidth',2);
legend('original','zero-mean-detrend' , 'log');
xlabel('Year');
ylabel('Zurich Number'); title('Sunspot Plots');
grid on; grid minor;
set(gca,'FontSize',18)

[pOriginal,fOriginal] = periodogram(ZurichNumber,hanning(N),[],1);
[pZeroMeanDetrend,fZeroMean] = periodogram(zeromeanDetrendZurich , hanning(N),[],1);
[pLog,fLog] = periodogram(ZurichLog , hanning(N) , [] , 1);

figure(2);
plot(fOriginal,10*log10(pOriginal),'g','LineWidth',2);hold on;
plot(fZeroMean,10*log10(pZeroMeanDetrend),'b','LineWidth',1);
plot(fLog,10*log10(pLog) ,'r','LineWidth',2);
legend('original','zero-mean-detrend','log');
xlabel('Cycles/Year')
ylabel('dB / (Cycles/Year)')
title('Periodogram of Relative Sunspot Number Data')
set(gca,'FontSize',18); grid on; grid minor;

%% Periodogram Based estimation of EEG
% 96000 samples, fs = 1200Hz, duration = 80s. Determine X, the rate at
% which the patient was flashed with the visual stimulus using spectral
% estimation
load Data/EEG_Data/EEG_Data_Assignment1
dt=1/fs; N=length(POz); samplesPerHertz = 5; Nfft = fs*samplesPerHertz; nOverlap = 0;
POz_zeroMean = POz - mean(POz);
t = 0:dt:(N-1)*dt;
figure(1);
plot(t,POz);xlabel('Time[s]');ylabel('EEG [V]');title('SSVEP');grid on; grid minor;set(gca,'FontSize',18);
figure(2);
[pxx,f] = periodogram(POz_zeroMean,hanning(N),Nfft ,fs);
plot(f,10*log10(pxx), 'g' , 'LineWidth',2);grid on;grid minor;
legend('Standard');title('EEG Periodogram with standard method');
xlabel('Frequency (Hz)');
ylabel('PSD (dB)');
xlim([10 50]);
ylim([-150 -100]);
set(gca,'FontSize',18);

%Bartlett
windows = [10 5 1];
PeriodogramWelch = []; frequencyWelch = [];
for i = 1:length(windows)
    windowSize = windows(i) * fs;
    [PeriodogramWelch(:,i),frequencyWelch(:,i)] = pwelch(POz_zeroMean, bartlett(windowSize), nOverlap, Nfft, fs);
end

welch = PeriodogramWelch(:,1);
welch = welch(1:(end-60));
variance_standard = (var(pxx));variance_average = var(welch);
figure(3);
plot(frequencyWelch(:,1),10*log10(PeriodogramWelch(:,1)), 'g' , 'LineWidth',2);
hold on;
plot(frequencyWelch(:,2),10*log10(PeriodogramWelch(:,2)), 'b' , 'LineWidth',2);
plot(frequencyWelch(:,3),10*log10(PeriodogramWelch(:,3)), 'r' , 'LineWidth',2);
legend('\Deltat = 10s', '\Deltat = 5s', '\Deltat = 1s' )
grid on;grid minor;title('EEG Periodograms with Bartlett Method');
xlabel('Frequency (Hz)');
ylabel('PSD (dB)');
xlim([10 50]);
ylim([-150 -100]);
set(gca,'FontSize',18);
hold off;
