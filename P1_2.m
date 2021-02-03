%% 1.2 Periodogram-based Methods Applied to Real World Data
clear all; close all; clc
load sunspot.dat
load EEG_Data/EEG_Data_Assignment1.mat
%% a)
xaxis = sunspot(:,1);
sunspot = sunspot(:,2);
sun_mean = sunspot - mean(sunspot);
sun_detrend = detrend(sunspot);
sun_mean_detrend = detrend(sun_mean);
sun_log = log(sunspot);
sun_log(isinf(sun_log)) = min(sun_log(~isinf(sun_log)));
sun_log_mean = sun_log - mean(sunspot);

[pxx,w] = periodogram(sunspot);
[pxx_mean_detrend,w1] = periodogram(sun_mean_detrend);
[pxx_log_mean, w2] = periodogram(sun_log_mean,hamming(length(xaxis)));

figure(1)
hold on
plot(xaxis,sunspot, xaxis,sun_mean, xaxis,sun_detrend, xaxis,sun_mean_detrend, xaxis,sun_log_mean, 'linewidth',1.1)
legend('Original','Mean removal','Detrend','Mean-detrend','Log-mean')
title('Sunspot Time Series')
xlabel('Year')
ylabel('Wolf Number')

figure(2)
hold on
plot(w/pi,mag2db(pxx), w1/pi,mag2db(pxx_mean_detrend), w2/pi,mag2db(pxx_log_mean), 'linewidth',1.1)
legend('Original','Mean-detrend','Log-mean')
title('Periodogram of Sunspot Time Series')
xlabel('Normalised Frequency (\pi rad/sample)')
ylabel('PSD (dB)')
ylim([-100 100])

%% b) The basis for brain computer interface (BCI)

ts = 1/fs;
nfft = 5 * fs;
[pxx_poz, f_poz] = periodogram(POz,[],nfft,fs);

figure(3)
plot(f_poz, mag2db(pxx_poz))
title('Periodogram of EEG Signal')
xlabel('Frequency (Hz)')
ylabel('PSD (dB)')
xlim([0 60])

window_t = [1 5 10];


for i = 1:length(window_t)
    window_len = window_t(i) * fs;
    pxx_avg{i} = pwelch(POz,window_len,0,nfft,fs);

end
figure(4)
hold on
plot(f_poz,mag2db(pxx_avg{3}),'linewidth',1.1)
hold on
plot(f_poz,mag2db(pxx_avg{2}),'linewidth',1.1)
hold on
plot(f_poz,mag2db(pxx_avg{1}),'linewidth',1.1)
legend('W = 10', 'W = 5', 'W = 1')
title('Averaged Periodogram of EEG Signal')
xlabel('Frequency (Hz)')
ylabel('PSD (dB)')
xlim([0 60])
