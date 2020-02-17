load RAW_NEW
%1 : 0.001
% Trial 1: 0 - 240s
% Trial 2: 260 - 490s
% Trial 3: 520 - 750;
DataChannel = data(:,2);
xECG1 = DataChannel(1:240000);
xECG2 = DataChannel(260000: 490000);
xECG3 = DataChannel(520000: 750000);

[xRRI1_Vaikkun,fsRRI1Vaikkun]=ECG_to_RRI(xECG1,fs);
[xRRI2_Vaikkun, fsRRI2Vaikkun] = ECG_to_RRI(xECG2,fs);
[xRRI_Vaikkun, fsRRI3Vaikkun] = ECG_to_RRI(xECG3,fs);


