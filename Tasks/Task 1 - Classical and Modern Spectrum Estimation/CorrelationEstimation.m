%% Correlogram part a)

%generate different signals
fsamp = 25600;
fsig = 1000; % sampling frequency
nsamp = 128;
Lfft = 2^nextpow2(nsamp);
dfx  = 1/(Lfft*tsamp);
WGN = randn(nsamp,1);    %noise
noisySine   = sinegen(fsamp,fsig,nsamp) + WGN'; %Noise Sine
windowLength = 10;
filterCoeff = ones(1, windowLength)/windowLength;
filteredNoise = filter(filterCoeff,1 , WGN);    %Filtered WGN
%plot different signals
figure(1);
plot(t,noisySine,'LineWidth',2); hold on; 
plot(t,WGN,'LineWidth',2);
plot(t,filteredNoise,'Linewidth',2); legend('noisy sine','noise','filtered noise'); hold off;
grid on; grid minor;
xlabel('Time[s]');ylabel('Amplitude');title('Plot of noise, filtered noise, and noisy sine');
set(gca,'FontSize',18);

signal = {WGN, filteredNoise , noisySine};
label = ["White Gaussian nNise", "Filtered White Gaussian Noise", "Noisy Sinusoid"];
numSignals = length(signal);
acfUnbiased = cell(numSignals,1);
acfBiased = cell(numSignals,1);
psdUnbiased = cell(numSignals,1);
psdBiased = cell(numSignals,1);

for i =1:numSignals
   acfUnbiased{i} = xcorr(signal{i},'unbiased');
   acfBiased{i} = xcorr(signal{i}, 'biased');
   psdUnbiased{i} = abs(fftshift(fft(acfUnbiased{i},Lfft)/Lfft));
   psdBiased{i} = abs(fftshift(fft(acfBiased{i},Lfft)/Lfft));
end

lag = (-(nsamp-1):1:(nsamp-1))*tsamp;

figure(2)
for i = 1:numSignals
   subplot(3,1,i);
   plot(lag,acfUnbiased{i},'r','LineWidth',2);hold on;plot(lag,acfBiased{i},'b','LineWidth',2);
   xlabel('lag [k]');ylabel('$$\hat{r}(k)$$','Interpreter','LaTex');title('ACF estimate of ' + label(i));
   grid on; grid minor; legend('Unbiased','Biased');set(gca,'FontSize',18);
end

freq = linspace(-fsamp/2,fsamp/2 , nsamp ); % Frequency axis
figure(3)
for i = 1:numSignals
   subplot(3,1,i);
   plot(freq,psdUnbiased{i},'r','LineWidth',2);hold on;plot(freq,psdBiased{i},'b','LineWidth',2);
   xlabel('frequency');ylabel('$$\hat{P}_{c}(\omega)$$','Interpreter','LaTex');title('Correlogram of ' + label(i));
   grid on; grid minor; legend('Unbiased','Biased');set(gca,'FontSize',18);
end

%%
% Use your code from the previous section (only the biased ACF estimator) to generate the PSD estimate of several [5]
% realisations of a 2 sinewaves corrupted with noise and plot them
fsig1 = 100; fsig2 = 200; numRealisations = 100; fsamp = 500; nsamp=128;
randomProcess = [];BiasedACF_realisation = [];BiasedPSD_realisation = [];
for i = 1:numRealisations
    randomProcess(i,:) = sinegen(fsamp,fsig1,nsamp) + sinegen(fsamp,fsig2,nsamp) + randn(1,nsamp) ;
    BiasedACF_realisation(i,:) = xcorr(randomProcess(i,:),'unbiased');
    BiasedPSD_realisation(i,:) = abs(fftshift(fft(BiasedACF_realisation(i,:),Lfft)/Lfft));
end
Ensemble_average = mean(randomProcess);AveragePSD = mean(BiasedPSD_realisation);stdPSD = std(BiasedPSD_realisation);

figure(1);
for i = 1:numRealisations
    plot(0 : 1/fsamp : (nsamp-1)*1/fsamp , randomProcess(i,:) , 'b');hold on;
end
plot(0 : 1/fsamp : (nsamp-1)*1/fsamp , Ensemble_average , 'r','LineWidth',2);
xlabel('Time[s]');ylabel('Amplitude');title('Overlay of 100 realisations of random process');grid on; grid minor;
set(gca,'FontSize',18);
    

figure(2);subplot(2,1,1);freq = linspace(-fsamp/2,fsamp/2 , nsamp );
for i = 1:numRealisations
    plot(freq,BiasedPSD_realisation(i,:),'b','LineWidth',1); hold on;
end
plot(freq,AveragePSD,'r','LineWidth',2); 
grid on; grid minor;xlabel('Frequency [Hz]');ylabel('$$\hat{P}_{c}(\omega)$$','Interpreter','LaTex');title('PSD of 100 realisations of sin(2\pi100/25600n) + sin(2\pi200/25600n) + w(n)');
set(gca,'FontSize',18);
subplot(2,1,2);
plot(freq,stdPSD,'b','LineWidth',2); 
grid on; grid minor;xlabel('Frequency [Hz]');ylabel('$$\hat{\sigma}_{C}(\omega)$$','Interpreter','LaTex');title('Standard Deviation of Correlogram Estimates of Noisy Sinusoids');
set(gca,'FontSize',18);

%% dB representation
figure(3); subplot(2,1,1);
for i = 1:numRealisations
    plot(freq,10*log10(BiasedPSD_realisation(i,:)),'k','LineWidth',1); hold on;
end
plot(freq,10*log10(AveragePSD),'m','LineWidth',2); 
grid on; grid minor;xlabel('Frequency [Hz]');ylabel('$$\hat{P}_{c}(\omega), [dB]$$','Interpreter','LaTex');title('PSD (in DB) of 100 realisations of sin(2\pi100/25600n) + sin(2\pi200/25600n) + w(n)');
set(gca,'FontSize',18);
subplot(2,1,2);
plot(freq,stdPSD,'g','LineWidth',2); 
grid on; grid minor;xlabel('Frequency [Hz]');ylabel('$$\hat{\sigma}_{C}(\omega)$$','Interpreter','LaTex');title('Standard Deviation (in dB) of Correlogram Estimates of Noisy Sinusoids');
set(gca,'FontSize',18);
%% Periodogram with Frequency Resolution
clear all; clc; close all;
numPoints = [20 30 50];colors = ['r','b','g'];
for i = 1:length(numPoints)
    n = 0:numPoints(i);
    noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
    x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n) + noise;
    subplot(length(numPoints),1,i);
    [pxx,w] = periodogram(x,hamming(length(x)));
    plot(w/(2*pi) , 10*log10(pxx),colors(i),'LineWidth',1);
    xlabel('Normalized Frequency (x \pi rad/sample)');
    ylabel('P(\omega) [dB]'); grid on; grid minor;
    title('Periodogram with N = ' + string(numPoints(i)) + ' Points');
    set(gca,'FontSize',18);
end

%% MUSIC algorithm
numPoints = [20 30 50];numRealisations = 50;N_points = length(numPoints);
f1 = 0.3 ; f2 = 0.32;
noisePower = 0.2;
realisationPSD_lengths = cell(N_points,1);meanPSD_lengths = cell(N_points,1);
stdPSD_lengths = cell(N_points,1);
for i = 1:N_points
    n = 0:numPoints(i);
    clean = exp(1j*2*pi*f1*n)+exp(1j*2*pi*f2*n);
    realisationPSD = cell(numRealisations,1);
    PSDOnlyRealisation = [];
    for j = 1:numRealisations
        noise = noisePower/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
        realisation = clean + noise;
        [X,R] = corrmtx(realisation,14,'mod');
        [S,F] = pmusic(R,2,[],1,'corr');
        realisationPSD{j} = [S F];
        PSDOnlyRealisation(j,:) = S;
    end
    realisationPSD_lengths{i} = realisationPSD;
    meanPSD_lengths{i} = mean(PSDOnlyRealisation);
    stdPSD_lengths{i} = std(PSDOnlyRealisation);
    
end

figure(1)
for i = 1:N_points
    subplot(N_points,1,i);
    for j = 1:numRealisations
        frequency = realisationPSD_lengths{i}{j}(:,2);
        PSD = realisationPSD_lengths{i}{j}(:,1);
        plot(frequency , PSD , 'b' , 'LineWidth' , 1); hold on;
    end
    plot(frequency, meanPSD_lengths{i}, 'r','LineWidth',2);
    set(gca,'xlim',[0.25 0.40]); grid on; grid minor; xlabel('Hz');ylabel('P(\omega)');
    title('PseudoSpectrum with N = ' + string(numPoints(i)) + ' Points');set(gca,'FontSize',18);
    hold off;
end

figure(2)
for i = 1:N_points
    subplot(N_points,1,i);
    plot(frequency, stdPSD_lengths{i}, 'r','LineWidth',2);
    set(gca,'xlim',[0.25 0.40]); grid on; grid minor; xlabel('Hz');ylabel('\sigma_{P}(\omega)');
    title('Standard Deviation of Pseudospectrum with N = ' + string(numPoints(i)) + ' Points');set(gca,'FontSize',18);
    
end
        
    

