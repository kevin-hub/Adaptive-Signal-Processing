%% 2.3 Adaptive Noise Cancellation
clear all; close all; clc;

%% a)
fs = 1;
nSamples = 1000;
t = (0:nSamples-1)/fs;
nRps = 100;
nTrans = 0;
maxDelay = 10;
order = 5;

A = 1;
w0 = 0.01 * pi;
x = A * sin(w0 * t);

coeff = [0 0.5];
var = 1;

MAmdl = arima('MA', coeff, 'Constant', 0, 'Variance', var);
[MAsig,innovation] = simulate(MAmdl,nSamples,'NumPaths',nRps);

cn = MAsig';
wn = innovation';

for j = 1:maxDelay
    for i = 1:nRps
        s = cn(i,:) + x;
        [u] = pp(s,order,j);
        
        [y(i,:),~,~] = lms(u,s,0.01,0);
        
        e(i,:) = (x-y(i,:)).^2;
    end
    y_ALE(j,:) = mean(y);
    mse(j,:) = mean(e.^2);
    mpse(j) = mean(mse(j,:));
end

figure(1)
plot(1:maxDelay,10*log10(mpse))
title('Error of ALE (M = 5)')
ylabel('MPSE (dB)')
xlabel('\Delta (samples)')

%% b)
minDelay = 3;
maxDelay = 25;
delay = minDelay:maxDelay;
M = linspace(5,20,4);

for j = 1:length(delay)
    for i = 1:length(M)
        for k = 1:nRps
            s_M = cn(k,:) + x;
            [u_M] = pp(s_M,M(i),delay(j));

            [y_M(k,:),~,~] = lms(u_M,s_M,0.01,0);
            
            e_M(k,:) = (x-y_M(k,:)).^2;
        end
        mpse_M(j,i) = mean(mean(e_M));
    end
end

figure(2)
hold on
plot(delay, 10*log10(mpse_M(:,1)))
plot(delay, 10*log10(mpse_M(:,2)))
plot(delay, 10*log10(mpse_M(:,3)))
plot(delay, 10*log10(mpse_M(:,4)))
xlim([minDelay maxDelay])
title('Error of ALE with Various Filter Orders')
xlabel('\Delta(samples)')
ylabel('MSPE (dB)')
legend('M=5','M=10','M=15','M=20','Location','northwest')

%% c)
alpha = 0.9;
beta = 0.05;
maxDelay = 10;


for i = 1:nRps
    s_ANC = cn(i,:) + x;

    sn = alpha * cn(i,:) + beta;
    [snpp] = pp(sn,order,0);

    [n(i,:),~,~] = lms(snpp,s_ANC,0.005,0);

    y_ANC(i,:) = s_ANC - n(i,:);

    e_ANC(i,:) = (x-y_ANC(i,:)).^2;
end
mpse_ANC = mean(mean(e_ANC.^2));

figure(3)
hold on
plot(s_ANC,'k','LineWidth',1.1)
plot(mean(y_ANC),'b','LineWidth',1.1)
plot(y_ALE(3,:),'c','LineWidth',1.1)
plot(x,'r','LineWidth',1.1)
title('Prediction Signals of Adaptive Filtering Configurations')
xlabel('Iteration (samples)')
ylabel('Amplitude')
legend('Noise-corrupted','ANC(M=5)','ALE(M=5,\Delta=3)','Clean')

%% d)
load EEG_Data/EEG_Data_Assignment2.mat

nSamples = length(POz);
nfft = 8192;
t = (0:nSamples-1)/fs;
order = 5:5:20;
M = length(order);
step = [1e-2, 1e-3, 1e-4];
u = length(step);

A = 1;
f0 = 50;
sine = sin(2*pi*50*t);
noise = randn(1, nSamples);
ref = sine + noise;

for i = 1:M
    [refpp] = pp(ref,order(i),0);
    for j = 1:u
        [nPOz,~,~] = lms(refpp,POz',step(j),0);
        yPOz{i,j} = POz' - nPOz;
        ePOz{i,j} = (POz'-yPOz{i,j}).^2;
        mspePOz{i,j} = mean(ePOz{i,j});
    end
end

figure(4);
spectrogram(POz,2^12,round(0.5*2^12),2^13,fs,'yaxis');
ylim([0 60]);
title('Spectrogram of POz');
% specteograms by ALE with different order and step size
for i = 1:M
    figure;
    for j = 1:u
        subplot(u, 1, j);
        spectrogram(yPOz{i,j}, 2^12,round(0.5*2^12),2^13,fs,'yaxis');
        ylim([0 60]);
        title(['Spectrogram of ANC Signal: M = ', num2str(order(i)), sprintf(', \\mu = '), num2str(step(j))]);
    end
end

order_opt = 10;
step_opt = 1e-3;

M_opt = find(order == order_opt);
u_opt = find(step == step_opt);

[psdPOz, fPOz] = periodogram(POz,[],nfft,fs);
psd_opt = periodogram(yPOz{M_opt,u_opt},[],nfft,fs);

figure;
plot(fPOz,mag2db(psdPOz), fPOz,mag2db(psd_opt))
title('Periodogram of Original and ANC EEG Signal')
xlabel('Frequency (Hz)')
ylabel('PSD (dB)')
legend('Original','ANC')
xlim([0 60])


