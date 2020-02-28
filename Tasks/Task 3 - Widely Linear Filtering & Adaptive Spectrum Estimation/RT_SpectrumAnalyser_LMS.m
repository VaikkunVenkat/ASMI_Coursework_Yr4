N = 1500;
f_n = frequency_n();
phi = integrateFrequency(f_n);
noiseVariance = 0.05;
fs = 5000;
eta = wgn(1500,1,10*log10(noiseVariance),'complex');
signal = exp(1i * 2*pi/fs .* phi);
y = signal + eta;

w = zeros(K, N);
e = zeros(1, N);
input = (1/K)*exp(1i*2*pi*(0:N-1)'*(0:K-1)/K);

[w, e] = DFT_CLMS(y, input, mu, gamma);