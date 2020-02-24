function y = cosinegen(fsamp,fsig,nsamp,phi,delta) % phi incorporates the phase-shift by 2pi/3 for phases b and c
% Generate a cosine waveform of a certain signal frequency and phase (and phase distortion), given the
% sampling frequency and the number of samples
tsamp = 1/fsamp;
t = 0:tsamp:(nsamp-1)*tsamp;
y = cos((2*pi*fsig*t) + phi + delta);
end

