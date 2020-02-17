N=1000;
Variance = 0.5;numRealisations = 100;order=1;
nu = sqrt(Variance)*randn(N,numRealisations);
b = [1 0.9];a=1;
x = filter(b,a,nu);
%%
mu = [0.01 0.1];
error = zeros(N,numRealisations);
Weights= zeros(N,numRealisations);
ErrorVals = cell(length(mu),1);
WeightVals = cell(length(mu),1);
AlgorithmNames = ["Benveniste","Farhang","Mathews"]; numAlgorithms=3;
WeightsTrajectories = cell(numAlgorithms,1);
Errortrajectories = cell(numAlgorithms,1);
Weights_standardLMS = cell(length(mu),1);
Error_standardLMS = cell(length(mu),1);
rho=0.005;alpha=0.8;
for n = 1:length(AlgorithmNames)
    chosenAlgorithm = AlgorithmNames(n);
    for k = 1:numRealisations
        nu_realisation = nu(:,k);
        [Prediction,error(:,k),Weights(:,k),Mu] = GASSlms(x(:,k),nu_realisation,order,rho,alpha,chosenAlgorithm);
    end
    WeightsTrajectories{n} = Weights;
    Errortrajectories{n} = error;
end
error = zeros(N,numRealisations);
Weights= zeros(N,numRealisations);
xHat = zeros(N,numRealisations);
for n = 1:length(mu)
    chosenMu = mu(n);
    for k =1:numRealisations
        nu_realisation = nu(:,k);
        [error(:,k),xHat(:,k),Weights(:,k)] = lms_MA(x(:,k),nu_realisation,order,chosenMu,0);
    end
    Weights_standardLMS{n} = Weights;
    Error_standardLMS{n} = error;
end

% Weight Error Computation
AverageWeightsBen = mean(WeightsTrajectories{1},2);
AverageWeightsFarhang = mean(WeightsTrajectories{2},2);
AverageWeightsMathews = mean(WeightsTrajectories{3},2);
AverageWeight_standardMu1 = mean(Weights_standardLMS{1},2);
AverageWeight_standardMu2 = mean(Weights_standardLMS{2},2);
figure(1)
plot(b(1)-AverageWeightsBen,'r','LineWidth',2); hold on;
plot(b(1)-AverageWeightsFarhang,'g','LineWidth',2);plot(b(1)-AverageWeightsMathews,'b','LineWidth',2);
plot(b(1)-AverageWeight_standardMu1,'k','LineWidth',2);plot(b(1)-AverageWeight_standardMu2,'m','LineWidth',2)
grid on; grid minor; xlabel('Iteration Number');ylabel('0.9 - w(n)');title('Weight Error Curves for different GASS Algorithms');
legend('Benveniste','Farhang-Ang','Mathews','Standard LMS - \mu=0.01','Standard LMS - \mu=0.1');set(gca,'FontSize',18);

% MSE 
figure(2);
Realisations_mu1 = Error_standardLMS{1};
Realisations_mu2 = Error_standardLMS{2};
Realisations_BenViste = Errortrajectories{1};
Realisations_Farhang = Errortrajectories{2};
Realisations_Mathews = Errortrajectories{3};
plot(10*log10(mean(Realisations_BenViste.^2 ,2)),'r','LineWidth',1); hold on;
plot(10*log10(mean(Realisations_Farhang.^2,2)) ,'g','LineWidth',1);plot(10*log10(mean(Realisations_Mathews.^2 ,2)),'b','LineWidth',2);
plot(10*log10(mean(Realisations_mu1.^2,2)) ,'k','LineWidth',1);
plot(10*log10(mean(Realisations_mu2.^2,2)) ,'m','LineWidth',1);
xlabel('Iteration number');ylabel('Error Power [dB]');title('Learning Curves of LMS for MA(1) process');
grid on; grid minor; legend('Benveniste','Farhang-Ang','Mathews','Standard LMS - \mu=0.01','Standard LMS - \mu=0.1');set(gca,'FontSize',18);

figure(3);
plot(mean(xHat,2));

%% 

error = zeros(N,numRealisations);
Weights= zeros(N,numRealisations);
ErrorVals = cell(length(mu),1);
WeightVals = cell(length(mu),1);
AlgorithmName = "Benveniste";
WeightsTrajectories = cell(2,1);
Errortrajectories = cell(2,1);
rho=0.004;alpha=0.8;mu=0.5;

for k = 1:numRealisations
    nu_realisation = nu(:,k);
    [Prediction,error(:,k),Weights(:,k),Mu] = GASSlms(x(:,k),nu_realisation,order,rho,alpha,AlgorithmName);
end

WeightsTrajectories{1} = Weights;
Errortrajectories{1} = error;
for k = 1:numRealisations
    nu_realisation = nu(:,k);
    [Prediction,error(:,k),Weights(:,k),Mu] = GNGD(x(:,k),nu_realisation,order,rho,mu);
end

WeightsTrajectories{2} = Weights;
Errortrajectories{2} = error;
AverageWeightsBen = mean(WeightsTrajectories{1},2);
AverageWeightsGNDG = mean(WeightsTrajectories{2},2);
figure(3);
plot(b(1)-AverageWeightsBen,'r','LineWidth',2); hold on;
plot(b(1)-AverageWeightsGNDG,'g','LineWidth',2);
grid on; grid minor; xlabel('Iteration Number');ylabel('0.9 - w(n)');title('Weight Error Curves for different Benveniste and GNGD LMS Algorithms \mu = ' + string(mu) + ', \rho = ' + string(rho));
legend('Benveniste','GNGD');set(gca,'FontSize',18);

figure(4);
ErrorBenveniste = Errortrajectories{1};
ErrorGNDG = Errortrajectories{2};
plot(10*log10(mean(ErrorBenveniste.^2 ,2)),'r','LineWidth',1); hold on;
plot(10*log10(mean(ErrorGNDG.^2,2)) ,'g','LineWidth',1);
xlabel('Iteration number');ylabel('Error Power [dB]');title('Learning Curves of Benveniste and GNGD LMS for MA(1) process \mu = ' + string(mu) + ', \rho = ' + string(rho)');
grid on; grid minor; legend('Benveniste','GNGD');set(gca,'FontSize',18);


