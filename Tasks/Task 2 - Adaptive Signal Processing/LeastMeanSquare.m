%% Implement LMS for second order AR predictor
N=1000;numRealisations = 100;
b=1; a = [0.1 0.8];order=length(a); a = [1 -a];variance = 0.25; 
nu = sqrt(variance)*randn(N,numRealisations);
x=filter(b,a,nu);

mu = [0.01 0.05];gamma=0;
error_realisations = zeros(N,numRealisations);
Weights_realisations = zeros(2,N+1,numRealisations);
ErrorVals = cell(length(mu),1);
WeightVals = cell(length(mu),1);
for l = 1:length(mu)
    for k = 1:numRealisations
        x_realisation = x(:,k);
        [E,W] = lms(x_realisation,order,mu(l) , gamma);
        error_realisations(:,k) = E;
        Weights_realisations(:,:,k) = W;
    end
    ErrorVals{l} = error_realisations;
    WeightVals{l} = Weights_realisations;
end

figure(1);
errorRealisation_mu1 = ErrorVals{1}(:,1);
errorRealisation_mu2 = ErrorVals{2}(:,1);
plot(10*log10(errorRealisation_mu1.^2 + eps),'g','LineWidth',2); hold on;
plot(10*log10(errorRealisation_mu2.^2 + eps),'b','LineWidth',2); hold on;
xlabel('Iteration Number');ylabel('Error Power [dB]');title('Error Power of LMS for AR(2) process');
grid on; grid minor; legend('\mu=0.01','\mu=0.05');set(gca,'FontSize',18);

figure(2);
Realisations_mu1 = ErrorVals{1};
Realisations_mu2 = ErrorVals{2};
plot(10*log10(mean(Realisations_mu1.^2 +eps,2)),'g','LineWidth',2);hold on;
plot(10*log10(mean(Realisations_mu2.^2,2)),'b','LineWidth',2);
xlabel('Iteration number');ylabel('Error Power [dB]');title('Learning Curve of LMS for AR(2) process');
grid on; grid minor; legend('\mu=0.01','\mu=0.05');set(gca,'FontSize',18);

figure(3);
plot(WeightVals{2}(1,:,1));hold on;plot(WeightVals{2}(2,:,1));

%%

MSE_mu1 = mean(Realisations_mu1(600:end,:).^2 , 1);
MSE_mu2 = mean(Realisations_mu2(600:end,:).^2 , 1);

EMSE_mu1 = mean(MSE_mu1 - variance,2);
EMSE_mu2 = mean(MSE_mu2 - variance,2);

misadjustment_mu1 = EMSE_mu1 / variance; % misadjustmentmu1 = -0.9988
misadjustment_mu2 = EMSE_mu2 / variance; % misdjustmentmu2 = -1 (perfect)

Correlation_matrix = [0.9259 0.4630 ; 0.4630 0.9259];
M_lms_mu1 = (mu(1)* trace(Correlation_matrix))/2;
M_lms_mu2 = (mu(2)* trace(Correlation_matrix))/2;

%%
MeanWeights_mu1 = mean(WeightVals{1},3);
MeanWeights_mu2 = mean(WeightVals{2},3);
figure(4);
plot(MeanWeights_mu1(1,:),'k','LineWidth',2);hold on;plot(MeanWeights_mu1(2,:),'k','LineWidth',2);
plot(MeanWeights_mu2(1,:),'m','LineWidth',2);plot(MeanWeights_mu2(2,:),'m','LineWidth',2);
plot([1 1000],[0.1 0.1],'r--');plot([1 1000],[0.8 0.8],'r--')
xlabel('Iteration Number');ylabel('w[n]');title('Average of Leaky LMS Trajectory of AR(2) Coefficients Over 100 Realisations  \gamma = ' + string(gamma));
grid on; grid minor;legend('\mu = 0.01, a_1=0.1','\mu = 0.01, a_2=0.8','\mu=0.05, a_1=0.1','\mu=0.05, a_2=0.8');
set(gca,'FontSize',18);


%Take final coefficient values and compare to truth
% mu1
a_1_est_mu1 = mean(MeanWeights_mu1(1,800:end)); error_a1_mu1 = abs(a_1_est_mu1 - 0.1);
a_2_est_mu1 = mean(MeanWeights_mu1(2,800:end)); error_a2_mu1 = abs(a_2_est_mu1 - 0.8);
a_1_est_mu2 = mean(MeanWeights_mu2(1,800:end)); error_a1_mu2 = abs(a_1_est_mu2 - 0.1);
a_2_est_mu2 = mean(MeanWeights_mu2(2,800:end)); error_a2_mu2 = abs(a_2_est_mu2 - 0.8);





