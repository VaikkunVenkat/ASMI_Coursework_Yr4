load /PCAPCR
% Determine Singular Values of X and X_noise 
S_X= svd(X);
S_Xnoise = svd(Xnoise);
rank_X = rank(X);
rank_Xnoise = rank(Xnoise);
figure(1);
stem(S_X, 'b','filled');hold on; stem(S_Xnoise, 'r','filled'); hold off;
xlabel('Singular Value Index');ylabel('Singular Value');grid on; grid minor;
legend('Singular Values of X','Singular Values of Xnoise')
title('Plot of Singular Values of X and Xnoise');set(gca,'FontSize',18);

error = S_X - S_Xnoise;
Squared_error = error.^2;

figure(2);
stem(Squared_error,'g','LineWidth',2);xlabel('Singular Value Index');ylabel('e(n)^2');title('Squared Error between Singular Values of X and Xnoise');
grid on; grid minor;set(gca,'FontSize',18);

%%
% Construct low-rank approximation of matrix Xnoise
[U_Xnoise , S_Xnoise , V_Xnoise] = svd(Xnoise);
k = rank_X;
Xnoise_tilde = U_Xnoise(:,1:k)*S_Xnoise(1:k,1:k)*V_Xnoise(:,1:k)';
S_XnoiseTilde = svd(Xnoise_tilde);

Norm_X_XNoise = norm(X - Xnoise,'fro'); Norm_X_XnoiseTilde = norm(X-Xnoise_tilde , 'fro');


error_tilde = S_X - S_XnoiseTilde; squaredErrorTilde = error_tilde.^2;
figure(3);
stem(Squared_error,'g','LineWidth',2); hold on; stem(squaredErrorTilde,'k','LineWidth',2); xlabel('Singular Value Index');ylabel('e(n)^2');title('Squared Error of Singular Values of Noisy and De-Noised Matrices');
grid on; grid minor; legend('Noisy X', 'De-Noised X');set(gca,'FontSize',18);

% Compute Error for different ranks and show that rank 3 is the most
% optimal rank:
E_X_XnoiseTilde = []; E_XnoiseTilde_Xnoise = []
for i =1:rank_Xnoise
    Xnoise_tilde = U_Xnoise(:,1:i)*S_Xnoise(1:i,1:i)*V_Xnoise(:,1:i)';
    E_X_XnoiseTilde(i) = norm(X - Xnoise_tilde,'fro');
    E_XnoiseTilde_Xnoise(i) = norm(Xnoise_tilde - Xnoise , 'fro');
end

figure(4);
plot(1:rank_Xnoise , E_X_XnoiseTilde, 'r', 'LineWidth',2);
hold on; plot(1:rank_Xnoise , E_XnoiseTilde_Xnoise, 'k', 'LineWidth',2)
xlabel('Rank'); ylabel('Frobenius Norm of Error');title('Errors Between True, Noisy and Reconstructed Input Signal Matrices');
grid on; grid minor; legend('Difference between true and de-noised X', 'Difference between de-noised and noisy X');set(gca,'FontSize',18);

%% Compute OLS and PCR solutions of B and find errors

Bhat_OLS = (Xnoise' * Xnoise)\ (Xnoise' * Y);   %OLS solution to Bhat
[U_Xnoise , S_Xnoise , V_Xnoise] = svd(Xnoise); %SVD of Xnoise
Xnoise_tilde = U_Xnoise(:,1:k)*S_Xnoise(1:k,1:k)*V_Xnoise(:,1:k)';
V1_r = V_Xnoise(:,1:k); Sigma1_r = S_Xnoise(1:k,1:k); U1_r = U_Xnoise(:,1:k); 
Bhat_PCR = V1_r * (Sigma1_r\(U1_r'*Y));
%csvwrite('B_hat.csv',[Bhat_OLS zeros(10,1) Bhat_PCR],0,0); Write the
%Estimated coefficient Matrices to a csv file

Yhat_OLS = Xnoise * Bhat_OLS ; Yhat_PCR = Xnoise * Bhat_PCR;
error_OLS = Y - Yhat_OLS ; error_PCR = Y - Yhat_PCR;
errornorm_OLS = norm(error_OLS,'fro'); errornorm_PCR = norm(error_PCR, 'fro');
figure(1);
subplot(2,2,1);
plot(mean(error_OLS,2),'r','LineWidth',2);grid on; grid minor;xlabel('Time Index [n]');ylabel('Error');
title('Error Plot between True and OLS-Estimated-Values');set(gca,'FontSize',18);
subplot(2,2,2);
plot(mean(error_PCR,2),'k','LineWidth',2);grid on; grid minor;xlabel('Time Index [n]');ylabel('Error');
title('Error Plot between True and PCR-Estimated-Values');set(gca,'FontSize',18);
subplot(2,2,3);
histogram(mean(error_OLS,2),'facecolor','r');grid on; grid minor;xlabel('Error');ylabel('Frequency');
title('Hisogram of OLS Error');legend('\mu = ' + string(mean(mean(error_OLS,2))) + ', \sigma = ' + string(std(mean(error_OLS,2)))) ;  set(gca,'FontSize',18);
subplot(2,2,4);
histogram(mean(error_PCR,2),'facecolor','k');grid on; grid minor;xlabel('Error');ylabel('Frequency');
title('Hisogram of PCR Error');legend('\mu = ' + string(mean(mean(error_PCR,2))) + ', \sigma = ' + string(std(mean(error_PCR,2)))); set(gca,'FontSize',18);
%%
[U_Xtest , S_Xtest , V_Xtest] = svd(Xtest);
k = 3;
Xtest_tilde = U_Xtest(:,1:k)*S_Xtest(1:k,1:k)*V_Xtest(:,1:k)';
Yhat_OLS = Xtest * Bhat_OLS ; Yhat_PCR = Xtest_tilde * Bhat_PCR;
error_OLS = Y - Yhat_OLS ; error_PCR = Y - Yhat_PCR;
mseOLS = sum(sum((error_OLS).^2))/numel(Yhat_OLS); msePCR = sum(sum((error_PCR).^2))/numel(Yhat_PCR);  
figure(1);
subplot(2,2,1);
plot(mean(error_OLS,2),'b','LineWidth',2);grid on; grid minor;xlabel('Time Index [n]');ylabel('Error');
title('Error Plot between True Test and OLS-Estimated-Values');set(gca,'FontSize',18);
subplot(2,2,2);
plot(mean(error_PCR,2),'g','LineWidth',2);grid on; grid minor;xlabel('Time Index [n]');ylabel('Error');
title('Error Plot between True Test and PCR-Estimated-Values');set(gca,'FontSize',18);
subplot(2,2,3);
histogram(mean(error_OLS,2),'facecolor','b');grid on; grid minor;xlabel('Error');ylabel('Frequency');
title('Hisogram of Test OLS Error');legend('\mu = ' + string(mean(mean(error_OLS,2))) + ', \sigma = ' + string(std(mean(error_OLS,2)))) ;  set(gca,'FontSize',18);
subplot(2,2,4);
histogram(mean(error_PCR,2),'facecolor','g');grid on; grid minor;xlabel('Error');ylabel('Frequency');
title('Hisogram of Test PCR Error');legend('\mu = ' + string(mean(mean(error_PCR,2))) + ', \sigma = ' + string(std(mean(error_PCR,2)))); set(gca,'FontSize',18);

%% regval

[Yhat_OLS_reg , YOLS] = regval(Bhat_OLS);
ERROR = YOLS - Yhat_OLS_reg;
MSE_OLS = sum(sum((YOLS-Yhat_OLS_reg).^2))/numel(YOLS);


[Yhat_PCR_reg , YPCR] = regval(Bhat_PCR);
ERROR = YPCR - Yhat_PCR_reg;
MSE_PCR = sum(sum((YPCR-Yhat_PCR_reg).^2))/numel(YPCR);

