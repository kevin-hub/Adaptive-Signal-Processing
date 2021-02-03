%% 3.2 Adaptive AR Model Based Time-Frequency Estimation
clear; close all; clc;

%% a)
fs = 1500;
N = 1500;
n = 1:N;

f = 100+(n-500)/2;
f(1:500) = 100;
f(1001:1500) = 100 + ((n(1001:1500)-1000)/25).^2;

var = 0.05;
noise = wgn(N,1,pow2db(var),'complex');

phase = cumsum(f);
modulation = exp(1j*(2*pi/fs*phase));
y =  exp(1j*(2*pi/fs*phase)) + noise';

figure(1)
subplot(2,1,1)
plot(f)
title('Frequency Variation the FM Signal')
xlabel('Sample')
ylabel('Frequency (Hz)')
subplot(2,1,2)
plot(angle(modulation))
title('Phase Variation the FM Signal')
xlabel('Sample')
ylabel('Angle (rad)')

coeff = aryule(y,1);
[h,w] = freqz(1,coeff,N,fs);
psd = abs(h).^2;

figure(2)
plot(w,pow2db(psd))
title('Complete Estimate of the FM Signal using AR(1)')
xlabel('Frequency (Hz)')
ylabel('PSD (dB)')
% 
% coeff1 = aryule(y(1:500),1);
% coeff2 = aryule(y(501:1000),1);
% coeff3 = aryule(y(1001:1500),1);
% 
% h1 = freqz(1,coeff1,N,fs);
% h2 = freqz(1,coeff2,N,fs);
% h3 = freqz(1,coeff3,N,fs);
% 
% psd1 = abs(h1).^2;
% psd2 = abs(h2).^2;
% psd3 = abs(h3).^2;

% plot(w,pow2db(psd1))
% plot(w,pow2db(psd2))
% plot(w,pow2db(psd3))

M = [1 2 5 10];
figure(3)
for i = 1:length(M)
    coeff1 = aryule(y(1:500),M(i));
    coeff2 = aryule(y(501:1000),M(i));
    coeff3 = aryule(y(1001:1500),M(i));

    h1 = freqz(1,coeff1,N,fs);
    h2 = freqz(1,coeff2,N,fs);
    h3 = freqz(1,coeff3,N,fs);

    psd1 = abs(h1).^2;
    psd2 = abs(h2).^2;
    psd3 = abs(h3).^2;
    
    subplot(3,1,1)
    hold on
    plot(w,pow2db(psd1))
    subplot(3,1,2)
    hold on
    plot(w,pow2db(psd2))
    subplot(3,1,3)
    hold on
    plot(w,pow2db(psd3))
    
end
subplot(3,1,1)
title('First Segment Estimate of the FM Signal using AR(1)')
legend('M=1','M=2','M=5','M=10')
subplot(3,1,2)
title('Second Segment Estimate of the FM Signal using AR(1)')
ylabel('PSD (dB)')
subplot(3,1,3)
title('Third Segment Estimate of the FM Signal using AR(1)')
xlabel('Sample')

%% b)
ypp = pp(y,1,1);
u = [1, 0.1, 0.01, 0.001];

for i = 1:length(u)
    [~, ~, w_h] = clms(ypp,y,u(i),0);
    for n = 1:N
        [h,w] = freqz(1, [1 -conj(w_h(n))], 1024,fs);
        H(:,n) = abs(h).^2;
    end
    % remove outliers
    medianH = 50 * median(median(H));
    H(H > medianH) = medianH;
    
    figure;
    surf(1:N,w,H,'EdgeColor','none');
    view(2);
    c = colorbar;
    c.Label.String = 'PSD (dB)';
    ylim([0 500])
    xlabel('Sample')
    ylabel('Frequency (Hz)')
end
figure(4)
title('CLMS-AR(1) Spectrogram \mu = 1')
figure(5)
title('CLMS-AR(1) Spectrogram \mu = 0.1')
figure(6)
title('CLMS-AR(1) Spectrogram \mu = 0.01')
figure(7)
title('CLMS-AR(1) Spectrogram \mu = 0.001')