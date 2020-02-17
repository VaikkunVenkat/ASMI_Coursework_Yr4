%% Implement adaptive line enhancer.
N = 1000;numRealisations=100;
b= [1 0 0.5]; a=1; 
v = wgn(N,numRealisations,0);
eta = filter(b,a,v);
w0 = 0.01 * pi;
x = sin(w0 * (1:N));
s = x' + eta;
mu=0.01;
delay = [1,2,3,4,5];
M=3; % change delay and order of predictor
error = zeros(N,numRealisations,length(delay));
xHat=zeros(N,numRealisations,length(delay));
Weights_realisations = zeros(M,N,numRealisations,length(delay));
MSPE = zeros(1,length(delay));
for d = 1:length(delay)
    for i =1:numRealisations
        s_realisation = s(:,i);
        [error(:,i,d) , xHat(:,i,d), Weights_realisations(:,:,i,d)] = lms_ale(s_realisation,mu,delay(d),M);
    end
    MSPE(d) = mean(mean((repmat(x', 1, numRealisations)-squeeze(xHat(:,:,d))).^2));
end

% Plot the Weight Trajectories
meanWeights = squeeze(mean(Weights_realisations,3));
figure(1);chosenDelay = delay(3);
for m=1:M
    plot(meanWeights(m,:,chosenDelay),'LineWidth',2);hold on;
end
xlabel('Iteration Number'); ylabel('Weight Values');title('Temporal Trajectory of Weights,\mu=0.01, M = ' + string(M) + ', \Delta = ' + string(chosenDelay));
grid on; grid minor; set(gca,'FontSize',18);legend('w1','w2','w3');

% Plot the noise-corrupted sine, the de-noised sine, and the true sine
% Manually change the delay to plot
figure(2);
plot(mean(s,2),'LineWidth',2);hold on;plot(mean(xHat(:,:,chosenDelay),2),'LineWidth',2);
plot(x','LineWidth',2);grid on; grid minor;legend('Noise-Corrupted','De-Noised','Clean');
xlabel('Time Index [n]');ylabel('Amplitude');title('Noise-Corrupted,De-Noised and Clean Sinewaves for M = ' + string(M) + ', \Delta = ' + string(chosenDelay));
set(gca,'FontSize',18);

figure;
plot(s,'r','LineWidth',2);hold on;plot(xHat(:,:,chosenDelay),'b','LineWidth',2);
plot(x','g','LineWidth',2);grid on; grid minor;
xlabel('Time Index [n]');ylabel('Amplitude');title('100 realisations of corrupted,de-noised and clean sine, M = ' + string(M) + ', \Delta = ' + string(chosenDelay));
set(gca,'FontSize',18);
%% Investigate effect of delay and order on MPSE

delay = 1:25;
M = [5,10,15,20];
error = zeros(N,numRealisations,length(delay));
xHat=zeros(N,numRealisations,length(delay));
MSPE = zeros(length(M),length(delay));
for m = 1:length(M)
    chosenOrder = M(m);
    Weights_realisations = zeros(chosenOrder,N,numRealisations,length(delay));
    for d = 1:length(delay)
        chosenDelay = delay(d);
        for i =1:numRealisations
            s_realisation = s(:,i);
            [error(:,i,d) , xHat(:,i,d), Weights_realisations(:,:,i,d)] = lms_ale(s_realisation,mu,chosenDelay,chosenOrder);
        end
        MSPE(m,d) = mean(mean((repmat(x', 1, numRealisations)-squeeze(xHat(:,:,d))).^2));
    end
end

figure(4);
for i =1:length(M)
    plot(MSPE(i,:),'LineWidth',2);hold on;
end

hold off; grid on; grid minor;xlabel('Delay \Delta');ylabel('MSPE');
title('MSPE across \Delta \in (1,25) and M \in {5,10,15,20}');
legend('M=5','M=10','M=15','M=20');set(gca,'FontSize',18);

%% Minimum occurs at delta = 3. Determine the best model order selection by 
%determining which one has the smallest MSPE.

delay = 3;
M=1:20;
error = zeros(N,numRealisations);
xHat=zeros(N,numRealisations);
MSPE = zeros(length(M),1);
for m = 1:length(M)
    chosenOrder = M(m);
    Weights_realisations = zeros(chosenOrder,N,numRealisations);
    for i =1:numRealisations
        s_realisation = s(:,i);
        [error(:,i) , xHat(:,i), Weights_realisations(:,:,i)] = lms_ale(s_realisation,mu,delay,chosenOrder);
    end
    MSPE(m) = mean(mean((repmat(x', 1, numRealisations)-squeeze(xHat(:,:))).^2));
end

figure(5);
plot(MSPE,'LineWidth',2);xlabel('Model Order, M');ylabel('MSPE');
[minimum_MSPE ,I] = min(MSPE); hold on;
plot([1 I],[minimum_MSPE minimum_MSPE],'r--');
plot([I I],[0 minimum_MSPE],'r--');
title('MSPE across model orders for \Delta=3');
grid on; grid minor;set(gca,'FontSize',18);

%% Implement Adaptive Noise Cancellation
M=5;delay=3; gamma=0; % Optimal / no leakage
nuHat = zeros(N,numRealisations);
xhat = zeros(N,numRealisations);
Weights_realisations = zeros(M,N,numRealisations);
for i =1:numRealisations
    s_realisation = s(:,i); eta_realisation = eta(:,i);
   [xhat(:,i), nuHat(:, i), Weights_realisations(:, :, i)] = lms_MA(s_realisation,eta_realisation, M,mu,gamma); 
end
MSPE_ANC = mean(mean((repmat(x', 1, numRealisations)-xhat(:,:)).^2));

% Compute MSPE for ALE and ANC and compare
Models = 1:20;MSPE_ANC = zeros(length(Models),1);MSPE_ALE = zeros(length(Models),1);
for m =1:length(Models)
    chosenModel = Models(m);
    Weights_realisations_ALE = zeros(chosenModel,N,numRealisations);
    Weights_realisations_ANC = zeros(chosenModel,N,numRealisations);
    for i =1:numRealisations
        s_realisation = s(:,i); eta_realisation = eta(:,i);
       [xhat(:,i), nuHat(:, i), Weights_realisations_ANC(:, :, i)] = lms_MA(s_realisation,eta_realisation, chosenModel,mu,gamma); 
       [error(:,i) , xHat(:,i), Weights_realisations_ALE(:,:,i)] = lms_ale(s_realisation,mu,delay,chosenModel);
    end
    MSPE_ANC(m) = mean(mean((repmat(x', 1, numRealisations)-xhat(:,:)).^2));
    MSPE_ALE(m) = mean(mean((repmat(x', 1, numRealisations)-squeeze(xHat(:,:))).^2));
end


% Plot ANC
figure(6);
plot(s,'r','LineWidth',2);
hold on;plot(xhat,'b','LineWidth',2);
plot(x','g','LineWidth',2);grid on; grid minor;
xlabel('Time Index [n]');ylabel('Amplitude');title('100 realisations of corrupted,de-noised and clean sine with ANC, M = ' + string(M));
set(gca,'FontSize',18);

figure(7);
plot(mean(s,2),'r','LineWidth',2);
hold on;plot(mean(xhat,2),'b','LineWidth',2);
plot(x','g','LineWidth',2);grid on; grid minor;
xlabel('Time Index [n]');ylabel('Ampltitude');title('Averaged Noisy, De-Noised and Clean Sinewaves with ANC, M = ' + string(M));
set(gca,'FontSize',18);legend('Noise-Corrupted','De-Noised with ANC','Clean');

figure(8);
plot(MSPE_ANC,'m','LineWidth',2);
hold on; plot(MSPE_ALE, 'k','LineWidth',2);
xlabel('Model Order, M');ylabel('MSPE');title('MSPE for different model orders with ANC and ALE, \Delta=3');
grid on; grid minor; set(gca,'FontSize',18);legend('ANC' , 'ALE');

%% Remove 50Hz from mains of EEG experiment
load Data/EEG_Data/EEG_Data_Assignment2
%%
dt=1/fs; N=length(POz); samplesPerHertz = 5; Nfft = fs*samplesPerHertz; nOverlap = 0;
POz_zeroMean = POz - mean(POz);
t = 0:dt:(N-1)*dt;

% Plot original Time Domain SSVEP
% figure(9);
% plot(t,POz_zeroMean,'r');xlabel('Time [s]');ylabel('SSVEP [mV]');title('SSVEP from EEG recording');grid on; grid minor; set(gca,'FontSize',18);
windows = 5;
windowSize = windows * fs;
[PeriodogramWelch_original,frequencyWelch_original] = pwelch(POz_zeroMean, bartlett(windowSize), nOverlap, Nfft, fs);



% Plot original Power Spectral Density Estimate
% figure(10);
% plot(frequencyWelch(:,1),10*log10(PeriodogramWelch(:,1)),'g','LineWidth',2);
% grid on;grid minor;title('EEG Welch Periodogram, \Deltat = ' + string(windows(1)));
% xlabel('Frequency (Hz)');
% ylabel('PSD (dB)');
% ylim([-150 -100]);
% set(gca,'FontSize',18);
% hold off;

% Plot original Spectrogram
K = 2^12;
L = 16384;
overlap = 0.75;
figure
spectrogram(POz_zeroMean, hamming(K), floor(overlap*K),L, fs,'yaxis');
colormap bone
ylim([0 100])
xlabel('Time (mins)')
ylabel('Frequency (Hz)')
set(gca, 'Fontsize', 22)
title(['EEG Original Spectrogram'])

% Remove 50Hz component through adaptive noise cancellation
F_mains = 50;
M= 20 %[5 10 15 20];    %Experiment with Model Order

mu = [0.005 0.01 0.02 0.04];  %Experiment with step size;
gamma = 0; %No Leaky

% Generate reference noisy sinusoid input
Reference_input = sin(2*pi*(F_mains/fs)*(1:N))' + wgn(N,1,0);

% Plot reference input
% figure(12);
% plot(Reference_input,'k','LineWidth',2);xlabel('Time Index');ylabel('Amplitude');
% grid on; grid minor; legend('Noise Reference Input');
% title('Noise Reference Input');set(gca,'FontSize',18);

% Perform ANC using LMS algorithm applied to POz data with different step
% sizes and model orders
frequencyWelch = []; PeriodogramWelch = [];
for i = 1:length(mu)
    weights_POz = zeros(M,N); %Change to M(i) if changing model orders
    [xHat, RefInputHat, weights_POz(:,:)] = lms_MA(POz_zeroMean,Reference_input, M,mu(i),gamma);
    figure
    spectrogram(xHat, hamming(K), floor(overlap*K),L, fs,'yaxis');
    ylim([0 100])
    xlabel('Time (mins)')
    ylabel('Frequency (Hz)')
    set(gca, 'Fontsize', 22)
    title(['EEG Spectrogram post ANC, \mu = ' + string(mu(i)) + ', M = ' + string(M)])
    windowSize = 5 * fs;
    [PeriodogramWelch(:,i),frequencyWelch(:,i)] = pwelch(xHat, bartlett(windowSize), nOverlap, Nfft, fs);
end

figure(13);
plot(frequencyWelch_original,10*log10(PeriodogramWelch_original),'LineWidth',2); hold on;
for i = 1:length(mu)
    plot(frequencyWelch(:,i),10*log10(PeriodogramWelch(:,i)),'LineWidth',2);
end
hold off;
xlabel('Frequency [Hz]');ylabel('Power [dB]');title('ANC Signal Power for Different Step-Sizes of constant Model Order M = '+ string(M));
grid on; grid minor;legend('Original','\mu=0.005','\mu=0.01','\mu=0.02','\mu=0.04');set(gca,'FontSize',18);


