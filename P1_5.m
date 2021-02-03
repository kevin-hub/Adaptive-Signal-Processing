%% 1.5 Real World Signals: Respiratory Sinus Arrhythmia from RR-Intervals
clear all; close all; clc
% EXTRACTING RRI FROM ECG (ONLY NEEDS TO BE DONE ONE TIME)

%iAmp_import_v40('ECG_Data/RAW_1.BIN'); %convert data to .mat
%iAmp_import_v40('ECG_Data/RAW_2.BIN');

% raw1 = load('ECG_Data/RAW_1.mat');
% raw2 = load('ECG_data/RAW_2.mat');
% 
% data = raw1.data;
% fsfs = raw1.fs;
% data_mean = data - mean(data);
% data_detrend = detrend(data_mean);
% plot(data_detrend);
% 
% s1 = 924;  s2 = 246245; %splits of three trials for raw1
% s3 = 324667; s4 = 504322; 
% s5 = 544900; s6 = 779359;

% s1 = 1560;  s2 = 250950; %splits of three trials for raw2
% s3 = 260160; s4 = 497777; 
% s5 = 507860; s6 = 748015;

% ecg_normal = data_detrend(s1:s2);
% ecg_fast = data_detrend(s3:s4);
% ecg_slow = data_detrend(s5:s6);
% 
% [rri, fs] = ECG_to_RRI(ecg_normal,fsfs);
% [rri_fast, fsrrii] = ECG_to_RRI(ecg_fast,fs);
% [rri_slow, fsrriii] = ECG_to_RRI(ecg_slow,fs);

%% a)
normal = load('ECG_Data/rri_normal_noanom2');
fast = load('ECG_Data/rri_fast_noanom2');
slow = load('ECG_Data/rri_slow_noanom2');
load ECG_Data/ECG_Data.mat

fs = normal.fs;
rri = {xRRI1 xRRI2 xRRI3};

nfft = 1024;
noverlap = 0;
window_t = [50 150];
window = window_t * fs;

figure(1)
hold on
for i = 1:length(rri)
    [psdStnd{i}, f] = periodogram(rri{i},rectwin(length(rri{i})), nfft, fs);
    psdW50{i} = pwelch(rri{i}, window(1), noverlap, nfft, fs);
    psdW150{i} = pwelch(rri{i}, window(2), noverlap, nfft, fs);
    
    subplot(3,1,i)
    hold on
    plot(f,mag2db(psdStnd{i}), f,mag2db(psdW50{i}), f,mag2db(psdW150{i}))
end

subplot(3,1,1)
title('PSD of RRI Under Normal Breathing')
lgd = legend('Standard','W = 50','W = 150');


subplot(3,1,2)
title('PSD of RRI Under Fast Breathing')
ylabel('PSD (dB)')

subplot(3,1,3)
title('PSD of RRI Under Slow Breathing')
xlabel('Frequency (Hz)')

for p = 2:50
    [coeffEst{p}, varEst(p)] = aryule(xRRI1,p);
    h = freqz(sqrt(varEst(p)),coeffEst{p},length(psdStnd{1}));
    psdEst{p} = abs(h).^2;
    %subplot(3,1,1)
    %plot(f, mag2db(psdEst{p}),'linewidth',1.1)
end

figure(2)
subplot(3,1,1)
hold on
plot(f,mag2db(psdStnd{1}),'k')
plot(f, mag2db(psdEst{2}), f, mag2db(psdEst{20}), f, mag2db(psdEst{50}),'linewidth',1.1)
title('PSD Estimate of RRI Under Normal Breathing')
lgd = legend('Standard','p = 2','p = 20','p = 50');

for p = 2:50
    [coeffEst{p}, varEst(p)] = aryule(xRRI2,p);
    h = freqz(sqrt(varEst(p)),coeffEst{p},length(psdStnd{2}));
    psdEst{p} = abs(h).^2;
end

subplot(3,1,2)
hold on
plot(f,mag2db(psdStnd{2}),'k')
plot(f, mag2db(psdEst{2}), f, mag2db(psdEst{20}), f, mag2db(psdEst{50}),'linewidth',1.1)
title('PSD Estimate of RRI Under Fast Breathing')
ylabel('PSD (dB)')

for p = 2:50
    [coeffEst{p}, varEst(p)] = aryule(xRRI3,p);
    h = freqz(sqrt(varEst(p)),coeffEst{p},length(psdStnd{3}));
    psdEst{p} = abs(h).^2;
end

subplot(3,1,3)
hold on
plot(f,mag2db(psdStnd{3}),'k')
plot(f, mag2db(psdEst{2}), f, mag2db(psdEst{20}), f, mag2db(psdEst{50}),'linewidth',1.1)
title('PSD Estimate of RRI Under Slow Breathing')
xlabel('Frequency (Hz)')
	