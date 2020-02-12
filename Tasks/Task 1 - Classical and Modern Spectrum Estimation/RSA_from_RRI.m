load Trials-RRI/Trial1RRI
load Trials-RRI/Trial2RRI
load Trials-RRI/Trial3RRI

Ts1 = 1/fsRRI ; Ts2 = 1/fsRRI2; Ts3 = 1/fsRRI3; Fs = fsRRI ; % sampling frequency for all RRI are same
n1=length(xRRI);n2=length(xRRI2);n3=length(xRRI3);  numRRI = 3;
t1=0:Ts1:(n1-1)*Ts1 ; t2=0:Ts2:(n2-1)*Ts2 ; t3=0:Ts3:(n3-1)*Ts3 ;
windows = [50 150];

figure(1);
subplot(3,1,1);
plot(t1,xRRI,'b','LineWidth',2);xlabel('Time [s]');ylabel('RRI [mV]');title('RRI of Trial 1 (Normal)'); grid on; grid minor; legend('Normal'); set(gca,'FontSize',18)
subplot(3,1,2);
plot(t2,xRRI2,'r','LineWidth',2);xlabel('Time [s]');ylabel('RRI [mV]');title('RRI of Trial 2 (Fast)'); grid on; grid minor; legend('Fast'); set(gca,'FontSize',18)
subplot(3,1,3);
plot(t3,xRRI3,'g','LineWidth',2);xlabel('Time [s]');ylabel('RRI [mV]');title('RRI of Trial 3 (Slow)'); grid on; grid minor; legend('Slow'); set(gca,'FontSize',18)
speed = ['Normal' 'Fast' 'Slow']
RRIsignals = {detrend(xRRI-mean(xRRI)), detrend(xRRI2 - mean(xRRI2)) , detrend(xRRI3 - mean(xRRI3))};
PSDStandard = cell(numRRI,1);
nSamples = [n1 n2 n3]; nFft=2048; colors = ['b','g','r'];
figure(2);
for i = 1:numRRI
    [PSDStandard{i}, f] = periodogram(RRIsignals{i}, hamming(nSamples(i)));
    subplot(numRRI,1,i);
    plot(f,10*log10(PSDStandard{i}),colors(i),'LineWidth',2);
    xlabel('Frequency [Hz]');ylabel('P(\omega) [dB]');title('Periodogram of Trial ' + string(i));
    grid on; grid minor;set(gca,'FontSize',18);
    
end
figure(3);
AvgPSD = {numRRI , length(windows)};
for i = 1:numRRI
    subplot(numRRI,1,i)
    for j = 1:length(windows)
        % number of samples of window
        lengthWindow = windows(j) * Fs;
        [AvgPSD{i, j}, fWelch] = pwelch(RRIsignals{i}, hamming(lengthWindow), 0, nFft, Fs);
        plot(fWelch, 10*log10(AvgPSD{i,j}),'LineWidth',2); hold on;
        legend('\deltat = ' + string(windows(j)));
    end
    xlabel('Frequency [Hz]');ylabel('P(\omega) [dB]');grid on;grid minor;title('Averaged Periodogram for Trial ' + string(i));
    set(gca , 'FontSize' , 18);legend('\deltat = ' + string(50) , '\deltat = ' + string(150));
end

% AR estimation of Power Spectrum

ModelOrders = 2:3:12 ; numOrders = length(ModelOrders);
var_AR_estimation = [] ; F = cell(numRRI,1);
figure(4)
for i = 1: numRRI
    subplot(numRRI,1,i);
    [PSDStandard{i}, f] = periodogram(RRIsignals{i}, hamming(nSamples(i)),nFft,Fs);
    plot(f,10*log10(PSDStandard{i}),colors(i),'LineWidth',1); hold on;
    xlabel('Frequency [Hz]');ylabel('P(\omega) [dB]');title('Periodogram of Trial ' + string(i));
    grid on; grid minor;set(gca,'FontSize',18);
    for j = 1: numOrders
        [Pxx,F_yulear] = pyulear(RRIsignals{i},ModelOrders(j),nFft,Fs);
        plot(F_yulear,10*log10(Pxx),'LineWidth',2);
    end
    legend('Original Periodogram','AR(2)','AR(5)','AR(8)','AR(11)')
end

    

