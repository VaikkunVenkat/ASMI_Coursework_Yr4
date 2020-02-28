N = 1500;fs = 5000;
f_n = frequency_n();
phi = integrateFrequency(f_n);
freqRange = (0:N-1)*fs/N ;
noiseVariance = 0.05;
fs = 5000;
eta = wgn(1500,1,10*log10(noiseVariance),'complex');
signal = exp(1i * 2*pi/fs .* phi);
y = signal + eta;

%Circularity plot:
plot(y,'r*');xlabel('y_{r}(n)');ylabel('y_{i}(n)');title('Circularity Distribution of y(n)');
grid on; grid minor; set(gca,'FontSize',18); 



alpha = exp(1i * 2*pi/N);
input = (1/N)*exp(1i*2*pi/N * (0:N-1)'*(0:N-1));
weights = zeros(N,N);
error = zeros(1,N);
prediction = zeros(1,N);
mu=1;gamma=0.05;

[weights, prediction,error] = DFT_CLMS(y, input, mu,gamma);


figure;
surf(1:N, freqRange, abs(weights), 'LineStyle','none');
ylim([0 600]);
colormap winter;colorbar;
view(2);
xlabel('Time Index [n]');
ylabel('f(n) [Hz]');
set(gca, 'Fontsize', 18);
title('DFT-CLMS Estimation of noisy FM signal, \mu = ' + string(mu) + ', \gamma = ' + string(gamma));
grid on;grid minor;

[weights, prediction,error] = DFT_NCLMS(y,input,mu);


figure;
surf(1:N, freqRange, abs(weights), 'LineStyle','none');
ylim([0 600]);
colormap winter;colorbar;
view(2);
xlabel('Time Index [n]');
ylabel('f(n) [Hz]');
set(gca, 'Fontsize', 18);
title('DFT-NCLMS Estimation of noisy FM signal');
grid on;grid minor;

%% Perform DFT CLMS on the POz data
clc;clear all; close all;
load('Data/EEG_Data/EEG_Data_Assignment1');
%%
N=1200;
segmentPOz = POz(200:200+1200-1);
freqRange = (0:N-1)*fs/N ;

input = (1/N)*exp(1i*2*pi/N * (0:N-1)'*(0:N-1));
mu=1;gamma=0;

[weights, prediction,error] = DFT_CLMS(segmentPOz, input, mu,gamma);


figure;
surf(1:N, freqRange, abs(weights), 'LineStyle','none');
ylim([0 600]);
colormap winter;colorbar;
view(2);
xlabel('Time Index [n]');
ylabel('f(n) [Hz]');
set(gca, 'Fontsize', 18);
title('DFT-CLMS Estimation of noisy FM signal, \mu = ' + string(mu) + ', \gamma = ' + string(gamma));
grid on;grid minor;



