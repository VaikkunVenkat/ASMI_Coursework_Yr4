N = 1000; % data length
realisations = 100; % number of realisations

% filter
b = [1 0 0.5];
a = 1;
order = length(b);

% params
mu = 0.01;
delta = [1 2 3 4];
M = 5;

xhat = zeros(length(delta),N,realisations);
error = zeros(length(delta),N,realisations);
A = zeros(length(delta),M,N,realisations);
s = zeros(realisations,N);
for d=1:length(delta)
    
    x = sin(0.01*pi*(1:N));
    for i=1:realisations
        v = normrnd(0, sqrt(1), 1, N);
        eta = filter(b, a, v);
        s(i,:) = x + eta;
        
        [xhat(d,:,i), error(d, :, i), A(d, :, :, i)] = lmsALE(s(i,:), mu, M, delta(d));
    end
    MSPE(d) = mean(mean((repmat(x', 1, realisations)-squeeze(xhat(d,:,:))).^2));
end

figure
plot(s', 'r','linewidth',2);
hold on
plot(squeeze(xhat(3,:,:)), 'b','linewidth',2); % vary for different plots
plot(x', 'k','linewidth',2);
xlabel('n')
ylabel('Signal Amplitude')
set(gca, 'Fontsize', 22)
title('MSPE=0.3077', 'Fontsize', 35)

figure
plot(mean(xhat(3,:,:),3),'linewidth',2)
xlabel('n')
ylabel('Denoised signal')
set(gca, 'Fontsize', 22)
title('Ensemble average Delta=4', 'Fontsize', 35)

% b 

N = 1000; % data length
realisations = 100; % number of realisations

% filter
b = [1 0 0.5];
a = 1;
order = length(b);

% params
mu = 0.01;
delta = 1:25;
M = [5 10 15 20];


MSPE = zeros(length(M),length(delta));
for m=1:length(M)
    xhat = zeros(length(delta),N,realisations);
    error = zeros(length(delta),N,realisations);
    A = zeros(length(delta),M(m),N,realisations);
    s = zeros(realisations,N);
    for d=1:length(delta)

        x = sin(0.01*pi*(1:N));
        for i=1:realisations
            v = normrnd(0, sqrt(1), 1, N);
            eta = filter(b, a, v);
            s(i,:) = x + eta;

            [xhat(d,:,i), error(d, :, i), A(d, :, :, i)] = lmsALE(s(i,:), mu, M(m), delta(d));
        end
        MSPE(m,d) = mean(mean((repmat(x', 1, realisations)-squeeze(xhat(d,:,:))).^2));
    end
end

figure
plot(MSPE','linewidth',2)
xlabel('Delta')
ylabel('MSPE')
set(gca, 'Fontsize', 22)
legend('M=5','M=10','M=15','M=20')
title('MSPE vs delta and M', 'Fontsize', 35)

figure
plot(M,MSPE(:,3),'linewidth',2)
xlabel('order')
ylabel('MSPE')
set(gca, 'Fontsize', 22)
title('MSPE vs order', 'Fontsize', 35)

% c
%%
b = [1 0 0.5];
a = 1;
order = length(b);

% params
mu = 0.01;
delta = [3];
M = 5;

xhat = zeros(N,realisations);
error = zeros(N,realisations);
A = zeros(M,N,realisations);
s = zeros(realisations,N);

    
x = sin(0.01*pi*(1:N));
for i=1:realisations
    v = normrnd(0, sqrt(1), 1, N);
    eta = filter(b, a, v);
    s(i,:) = x + eta;

    [xhat(:,i), error(:, i), A(:, :, i)] = lmsG(s(i,:), eta, mu, M); %x, in, mu, order)
end
MSPE = mean(mean((repmat(x', 1, realisations)-xhat(:,:)).^2));


figure
plot(s', 'r','linewidth',2);
hold on
plot(error(:,:), 'b','linewidth',2); % vary for different plots
plot(x', 'k','linewidth',2);
xlabel('n')
ylabel('Signal Amplitude')
set(gca, 'Fontsize', 22)
title('ANC Delta=3 M=5', 'Fontsize', 35)

%%
% d

load Coursework/Data/EEG_Data/EEG_Data_Assignment2

x = detrend(POz);

N = length(x);
% signal generation
f = 50/fs;
var = 1;
sine = sin(2*pi*f*(1:N));
w = normrnd(0, sqrt(var), 1, N);
eps = sine + w;

% Spectogram
K = 2^12;
L = 16384;
overlap = 0.75;
figure
spectrogram(x, hamming(K), floor(overlap*K),L, fs,'yaxis');
ylim([0 100])
xlabel('Time (min)')
ylabel('Frequency (Hz)')
set(gca, 'Fontsize', 22)
title(['Original signal'])

% varying mu and M : comment/uncomment accordingly
% mu = [0.001, 0.005, 0.01, 0.03];
% M = 20;
mu = 0.01;
M = [5 10 20 30];

l = length(M); % change mu-M accordingly

% change init according to what tested
xhat = zeros(N,l);
error = zeros(N,l);
% PSD = zeros(L,l);
for i=1:l
    
%     A = zeros(M(i),N);
%     [xhat(:,i), error(:, i), A(:, :)] = lmsG(x, eps, mu, M(i)); %x, in, mu, order)

    A = zeros(M(i),N);
    [xhat(:,i), error(:, i), A(:, :)] = lmsG(x, eps, mu, M(i)); %x, in, mu, order)
    
    figure
    spectrogram(squeeze(error(:,i)), hamming(K), floor(overlap*K),L, fs,'yaxis');
    ylim([0 100])
    xlabel('Time (min)')
    ylabel('Frequency (Hz)')
    set(gca, 'Fontsize', 22)
    title(['M=', num2str(M(i)), ' \mu=', num2str(mu)])
    
    xHat = reshape(error(:,i)-mean(x), 5*fs, []);
    PSD(:,i) = mean(periodogram(xHat,hamming(length(xHat)),5*fs,fs), 2);
%     PSD(:,i) = mean(periodogram(xHat,hamming(length(xHat)), [], L, 'centered'), 2);   
end

xr = reshape(x, 5*fs, []);
% PSDOrig = mean(periodogram(xr, [], L, 'centered'), 2);
PSDOrig = mean(periodogram(xr,hamming(length(xHat)),5*fs,fs), 2);
fp = (fs/2)*linspace(0, 1, length(PSDOrig));

figure
plot(fp,10*log10(PSDOrig), 'linewidth', 2);
hold on
plot(fp,10*log10(PSD(:,1)), 'linewidth', 2);
plot(fp,10*log10(PSD(:,2)), 'linewidth', 2);
plot(fp,10*log10(PSD(:,3)), 'linewidth', 2);
plot(fp,10*log10(PSD(:,4)), 'linewidth', 2);
xlim([0 80])
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
set(gca, 'Fontsize', 22)
legend('original','M=5','M=10','M=20','M=30')
title(['Periodogram for varying M'])