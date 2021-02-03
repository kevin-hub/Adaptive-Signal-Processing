%% 1.3 Correlation Estimation
clear all; close all; clc
%% a)
fs = 1;
N = 1024;
lags = -N+1:N-1;
w = lags./(2*N) * fs;
t = ((0:(N-1))/fs).';
sinef = [0.2 0.4];

WGN = wgn(N,1,0);
filtWGN = filter([1 1],1,WGN);
noisySine = sin(2*pi*sinef(1)*t) + sin(2*pi*sinef(2)*t) + WGN;

signals = {WGN, filtWGN, noisySine};

for i = 1:length(signals)
    r_bias_i = xcorr(signals{i},'biased');
    r_unbias_i = xcorr(signals{i},'unbiased');
    
    acf_bias{i} = r_bias_i;
    acf_unbias{i} = r_unbias_i;
    
    psd_bias_i = real(fftshift(fft(ifftshift(r_bias_i))));
    psd_unbias_i = real(fftshift(fft(ifftshift(r_unbias_i))));
    
    psd_bias{i} = psd_bias_i;
    psd_unbias{i} = psd_unbias_i;
end

%correlogram/autocorrelation plot
figure(1)
subplot(3,1,1)
plot(lags,acf_unbias{1}, lags,acf_bias{1})
title('Correlogram of White Gaussian Noise')
xlim([-N+1 N+1])
lgd = legend('Unbiased','Biased');
subplot(3,1,2)
plot(lags,acf_unbias{2}, lags,acf_bias{2})
title('Correlogram of Filtered White Gaussian Noise')
xlim([-N+1 N+1])
ylabel('Autocorrelation')
subplot(3,1,3)
plot(lags,acf_unbias{3}, lags,acf_bias{3})
title('Correlogram of Sinusoidal Signals with Noise')
xlim([-N+1 N+1])
xlabel('Lags/k')

%psd (correlogram-based estimation) plot
figure(2)
subplot(3,1,1)
plot(w,psd_unbias{1}, w,psd_bias{1})
title('PSD of White Gaussian Noise')
%xlim([-N+1 N+1])
lgd = legend('Unbiased','Biased');

subplot(3,1,2)
plot(w,psd_unbias{2}, w,psd_bias{2})
title('PSD of Filtered White Gaussian Noise')
%xlim([-N+1 N+1])
ylabel('PSD')
subplot(3,1,3)
plot(w,psd_unbias{3}, w,psd_bias{3})
title('PSD of Sinusoidal Signals with Noise')
%xlim([-N+1 N+1])
xlabel('Normalised Frequency (\pi rad/sample)')

%% b)
nRls = 100;

figure(3)
subplot(2,1,1)
hold on
for i = 1:nRls
    WGN = wgn(N,1,0);
    noisySine = sin(2*pi*sinef(1)*t) + sin(2*pi*sinef(2)*t) + WGN;
    
    acfRps{i} = xcorr(noisySine,'biased');
    psdRps{i} = real(fftshift(fft(ifftshift(acfRps{i}))));
    
    plot(w,psdRps{i},'c')
end
psdRps_mean = mean(cell2mat(psdRps),2);
plot(w,psdRps_mean,'b')
title('PSD by Biased Estimators of Different Realisations and Mean')
xlabel('Normalised Frequency (\pi rad/sample)')
ylabel('PSD')

%% c)
subplot(2,1,2)
psdRps_std = std(cell2mat(psdRps).');
plot(w,psdRps_std,'r')
title('Standard Deviation of the PSD Estimates')
xlabel('Normalised Frequency (\pi rad/sample)')
ylabel('PSD')


figure(4)
subplot(2,1,1)
hold on

for i = 1:nRls
    plot(w,mag2db(psdRps{i}),'c')
end

plot(w,mag2db(psdRps_mean),'b')
title('PSD by Biased Estimators of Different Realisations and Mean')
xlabel('Normalised Frequency (\pi rad/sample)')
ylabel('PSD (dB)')

subplot(2,1,2)
psdRps_std = std(mag2db(cell2mat(psdRps).'));
plot(w,psdRps_std,'r')
title('Standard Deviation of the PSD Estimates')
xlabel('Normalised Frequency (\pi rad/sample)')
ylabel('PSD (dB)')

%% d)
fs = 1;
N = [30 50];
nfft = 128;

figure(5)
hold on
for i = 1:length(N)
    n = 0:N(i)-1;
    noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n))); 
    x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+ noise;
    [pxx,f] = periodogram(x,rectwin(length(x)),nfft,fs);
    plot(f,mag2db(pxx),'linewidth',1.1)
end
title('Periodogram of Complex Exponentials')
xlabel('Normalised Frequency (\pi rad/sample)')
ylabel('PSD (dB)')
legend('N=30','N=50')

%% e)

N = [20 30 40];
for i = 1:length(N)
    for k = 1:nRls
        n = 0:N(i)-1;
        noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n))); 
        x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+ noise;

        [X,R] = corrmtx(x,14,'mod');
        [S{k},F] = pmusic(R,2,[],1,'corr');
    end
    MUpsd{i} = S;
    MUmean{i} = mean(cell2mat(S).');
    MUstd{i} = std(cell2mat(S).');
end



figure(6) 
subplot(3,2,1)
hold on
for i = 1:nRls
    plot(F,MUpsd{1}{i},'c')
end
plot(F,MUmean{1},'b'); 
set(gca,'xlim',[0.25 0.40]);
grid on;
title('MUSIC PSD Estimates and Mean: N=20')
subplot(3,2,2)
plot(F,MUstd{1},'r')
set(gca,'xlim',[0.25 0.40]);
grid on;
title('Standard Deviation of MUSIC PSD Estimates: N=20')

subplot(3,2,3)
hold on
for i = 1:nRls
    plot(F,MUpsd{2}{i},'c')
end
plot(F,MUmean{2},'b'); 
set(gca,'xlim',[0.25 0.40]);
grid on;
title('                                                         N=30')
ylabel('Pseudospectrum');
subplot(3,2,4)
plot(F,MUstd{2},'r')
set(gca,'xlim',[0.25 0.40]);
grid on;
title('                                                         N=30')

subplot(3,2,5)
hold on
for i = 1:nRls
    plot(F,MUpsd{3}{i},'c')
end
plot(F,MUmean{3},'b'); 
set(gca,'xlim',[0.25 0.40]);
grid on;
title('                                                         N=40')
xlabel('Hz'); 
subplot(3,2,6)
plot(F,MUstd{3},'r')
set(gca,'xlim',[0.25 0.40]);
grid on;
title('                                                         N=40')
xlabel('Hz'); 



