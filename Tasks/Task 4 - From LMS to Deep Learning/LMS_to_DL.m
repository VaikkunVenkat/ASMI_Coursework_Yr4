%% LOAD TIME SERIES DATA
load('Data/time-series')
%% Question 4
mu=1e-5;a=81;order=4;N=length(y);w_n=zeros(order+1,N);dynamic="normal";
[weightLMS,yHatLMS,errorLMS] = lms4(y,mu,order,dynamic,a);
[weight,yHat,error] = lms4_dynamical_bias(y,mu,order,a,w_n);
figure;
plot(1:N,y,'r','LineWidth',2);hold on;plot(1:N,yHat,'b','LineWidth',2);
xlabel('Time Index [n]');ylabel('Amplitude');title('Dynamical neuron with bias predicted and original signal');
grid on; grid minor;set(gca,'FontSize',18);legend('Original','Dynamical Neuron with bias');
figure;
plot(1:N,y,'r','LineWidth',2);hold on;plot(1:N, yHatLMS, 'b','LineWidth',2);
xlabel('Time Index [n]');ylabel('Amplitude');title('simple neuron predicted and original signal');
grid on; grid minor;set(gca,'FontSize',18);legend('Original','Simple Neuron');
figure;
plot(weight','LineWidth',2);hold on;
grid on; grid minor;xlabel('Time Index [n]');ylabel('w(n)');title('Weight trajectory of dynamical neuron with bias');
set(gca,'FontSize',18);
figure;
plot(weightLMS','LineWidth',2);hold on;
grid on; grid minor;xlabel('Time Index [n]');ylabel('w(n)');title('Weight trajectory of simple neuron with LMS');
set(gca,'FontSize',18);
figure;
plot(10*log10(error.^2),'k','LineWidth',2);hold on;plot(10*log10(errorLMS.^2),'m','LineWidth',2);
grid on; grid minor;xlabel('Time Index [n]');ylabel('Error Power');title('Learning curve of dynamical neuron with bias and simple neuron with LMS');
set(gca,'FontSize',18);legend('Dynamical Neuron with Bias','Simple Neuron');

%% QUESTION 5
[w_init,WeightPath,Error,Rp] = pre_train(y,mu,order,a);
w_n = zeros(order+1,N);w_n(:,order+1) = w_init;
[weights,yHat_pretrain,e_n] = lms4_dynamical_bias(y,mu,order,a,w_n);
figure;
plot(1:N,y,'r','LineWidth',2);hold on;plot(1:N, yHat_pretrain, 'b','LineWidth',2);
xlabel('Time Index [n]');ylabel('Amplitude');title('Pre-trained dynamical neuron with bias predicted and original signals');
grid on; grid minor;set(gca,'FontSize',18);legend('pre-trained dynamical neuron with bias','Original');
figure;
plot(1:N,10*log10(e_n.^2),'k','LineWidth',2);hold on;plot(1:N,10*log10(error.^2),'m','LineWidth',2);
xlabel('Time Index [n]');ylabel('Error Power');title('Learning Curves of normal and pre-trained dynamical neuron with bias');
grid on; grid minor;set(gca,'FontSize',18);legend('With pre-training','Without pre-training');