% Generate circular white gaussian noise of zero-mean unit-variance
numRealisations=100;N=1000;
x = wgn(N,numRealisations,0,'complex');

[circ_quot_x , circ_coeff_x] = circularity_calculation(x);

% initial Complex Distribution of x:
figure(1);
plot(x,'b.','LineWidth',2);xlabel('Real');ylabel('Imaginary');title('Circular White Gaussian Noise, \mu = 0');
legend('\rho_z = ' + string(round(circ_quot_x,2)) + ' \eta = ' + string(round(circ_coeff_x,2)))
grid on; grid minor;set(gca,'FontSize',18);

% Perform Widely Linear MA(1) Filtering
b1 = 1.5 + 1i; b2 = 2.5 - 0.5*1i;order=1;
y = complex(zeros(N,numRealisations));
for i = 1:numRealisations
    y(:,i) = WLMA(b1,b2,x(:,i));
end
[circ_quot_y , circ_coeff_y] = circularity_calculation(y);
figure(2);
plot(y,'r.','LineWidth',2);xlabel('Real');ylabel('Imaginary');title('WLMA(1) Process');
legend('\rho_z = ' + string(round(circ_quot_y,2)) + ' \eta = ' + string(round(circ_coeff_y,2)))
grid on; grid minor;set(gca,'FontSize',18);

% Perform system identification

%Initialise CLMS Variables
mu=0.005;numCoeffs = 2*order;
weights_CLMS = complex(zeros(numCoeffs,N,numRealisations));
prediction_CLMS = complex(zeros(N,numRealisations));
error_CLMS = complex(zeros(N,numRealisations));

%Initialise ACLMS Variables
weights_h = complex(zeros(numCoeffs,N,numRealisations));
weights_g = complex(zeros(numCoeffs,N,numRealisations));
prediction_ACLMS = complex(zeros(N,numRealisations));
error_ACLMS = complex(zeros(N,numRealisations));

for i = 1:numRealisations
    [weights_CLMS(:,:,i),prediction_CLMS(:,i),error_CLMS(:,i)] = CLMS(y(:,i),x(:,i),mu,numCoeffs);
    [weights_h(:,:,i), weights_g(:,:,i),prediction_ACLMS(:,i),error_ACLMS(:,i)] = ACLMS(y(:,i), x(:,i), mu, numCoeffs);
end

% plot convergence of CLMS weights
MeanWeights = conj(mean(weights_CLMS,3));
b2 = MeanWeights(2,:);
b2R = real(b2);b2I = imag(b2); 
figure(3);
plot(b2R,'g','LineWidth',2); hold on; plot(b2I,'g','LineWidth',2);
xlabel('Iternation Number');ylabel('Coefficient Value');title('Weight Trajectories of WLMA(1) Process using CLMS - \mu=' + string(mu))
legend('Re(b_1)','Im(b_1)');
grid on; grid minor; set(gca,'FontSize',18); 

% plot error for CLMS
MeanError = mean(error_CLMS,2);
SquaredErrordB = 10*log10(abs(MeanError).^2);
figure(4);
plot(SquaredErrordB,'r','LineWidth',2);xlabel('Iteration Number');ylabel('Error Power [dB]');
title('Learning Curve');grid on; grid minor;
set(gca,'FontSize',18);

% plot convergence of ACLMS weights
MeanWeights_h = conj(mean(weights_h,3));
MeanWeights_g = conj(mean(weights_g,3));
b1 = MeanWeights_h(2,:); b1R = real(b1);b1I = imag(b1);
b2 = MeanWeights_g(2,:); b2R = real(b2);b2I = imag(b2);
figure(5);
plot(b1R,'g','LineWidth',2); hold on;plot(b1I,'g','LineWidth',2);
plot(b2R,'b','LineWidth',2); hold on; plot(b2I,'b','LineWidth',2);
xlabel('Iternation Number');ylabel('Coefficient Value');title('Weight Trajectories of WLMA(1) Process using ACLMS - \mu=' + string(mu))
legend('Re(b_1)','Im(b_1)','Re(b_2)','Im(b_2)');
grid on; grid minor; set(gca,'FontSize',18); 


% Plot error curve of ACLMS
MeanError = mean(error_ACLMS,2);
SquaredErrordB = 10*log10(abs(MeanError).^2);
figure(6);
plot(SquaredErrordB,'r','LineWidth',2);xlabel('Iteration Number');ylabel('Error Power [dB]');
title('WLMA(1) Process Learning Curve - ACLMS');grid on; grid minor;
set(gca,'FontSize',18);


%% Wind Data

%% Import the Wind Data
clc; clear all; close all;
numSpeeds = 3;
windComplexSignals = cell(numSpeeds,1);
names = ["high-wind","medium-wind","low-wind"];
for i = 1:numSpeeds
    windSpeed = names(i);
    load('Data/wind-dataset/'+windSpeed)
    windComplexSignals{i} = complex(v_east,v_north);
end

%%

% Plot circularity distribution for each wind speed
colors = ["b.","r.","g."];
for i = 1:numSpeeds
    windSignal = windComplexSignals{i};
    [circ_quotient, circ_coeff] = circularity_calculation(windSignal);
    figure(i);
    plot(windSignal,colors(i),'LineWidth',2);
    xlabel('v_{east}');ylabel('v_{north}');title(names(i)+ " complex distribution");
    legend('\rho_z = ' + string(round(circ_quotient,2)) + ', \eta = ' + string(round(circ_coeff,2)))
    grid on; grid minor;
    set(gca,'FontSize',18);
end

windHigh =  windComplexSignals{3};
windMedium = windComplexSignals{2};
windLow = windComplexSignals{1};

% One Step Ahead Prediction. delay input by one.

%High Wind speed
filterOrders = [1 5 15]; numOrders = length(filterOrders); N = length(windHigh);
ErrorHigh = complex(zeros(N,numOrders));
ErrorHighACLMS = complex(zeros(N,numOrders));
mu = 0.01;
weightCLMS = cell(numOrders,1); predictionCLMS = complex(zeros(N,numOrders));
weight_h_ACLMS = cell(numOrders,1); predictionACLMS = complex(zeros(N,numOrders));
weight_g_ACLMS = cell(numOrders,1);
for k=1:numOrders
    weights_h = complex(zeros(filterOrders(k),N));
    weights_h_ACLMS = complex(zeros(filterOrders(k),N));
    weights_g_ACLMS = complex(zeros(filterOrders(k),N));
    delaywindHigh = [0; windHigh(1:end-1)];
    [weights_h, predictionCLMS(:,k),ErrorHigh(:,k)] = CLMS(windHigh, delaywindHigh, mu, k);
    [weights_h_ACLMS,weights_g_ACLMS , predictionACLMS(:,k) ,ErrorHighACLMS(:,k)] = ACLMS(windHigh, delaywindHigh, mu, k);
    weightCLMS{k} = weights_h;
    weight_h_ACLMS{k} = weights_h_ACLMS;
    weight_g_ACLMS{k} = weights_g_ACLMS;
end

figure;
for k =1:numOrders
    plot(10*log10(abs(ErrorHigh(:,k)).^2),'LineWidth',1); hold on;
end
xlabel('Iteration Number');ylabel('Error Power [dB]');title('High Wind Speed - CLMS Learning Curves');
grid on; grid minor; set(gca,'FontSize',18);legend('M=1','M=5','M=15')

figure;
for k =1:numOrders
    plot(10*log10(abs(ErrorHighACLMS(:,k)).^2),'LineWidth',1); hold on;
end
xlabel('Iteration Number');ylabel('Error Power [dB]');title('High Wind Speed - CLMS Learning Curves');
grid on; grid minor; set(gca,'FontSize',18);legend('M=1','M=5','M=15')

for i = 1:numOrders
    windSignal = windHigh;
    [~, circ_coeff_wind] = circularity_calculation(windSignal);
    [~ , circ_coeff_prediction] = circularity_calculation(predictionCLMS(:,i));
    varPrediction = var(predictionCLMS(:,i)); 
    figure;
    plot(windSignal,'.'); hold on;
    plot(predictionCLMS(:,i),'.');
    xlabel('v_{east}');ylabel('v_{north}');title("One-step update of high-speed wind, M = " + string(filterOrders(i)) + '\sigma_{z}^{2} = ' + string(round(varPrediction,3)));
    legend('\eta_{original} = ' + string(round(circ_coeff_wind,2)),'\eta_{update} = ' + string(round(circ_coeff_prediction,2)));
    grid on; grid minor;
    set(gca,'FontSize',18);
end

% CONTINUE SAME FOR LOW AND MEDIUM WIND SPEED

%%
clear all; clc; close all;





