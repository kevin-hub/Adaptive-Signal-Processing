%% 3.3 A Real Time Spectrum Analyser Using Least Mean Square
clear; close all; clc;

%% c)
fs = 1500;
N = 1500;
n = 1:N;
W = 1500;
w = (0:W-1);

f = 100+(n-500)/2;
f(1:500) = 100;
f(1001:1500) = 100 + ((n(1001:1500)-1000)/25).^2;


var = 0.05;
noise = wgn(N,1,pow2db(var),'complex');

phase = cumsum(f);
modulation = exp(1j*(2*pi/fs*phase));
y =  exp(1j*(2*pi/fs*phase)) + noise';

x = 1/W * exp(n' * 1j * 2 * pi * w / W)';

gamma = [0 0.01 0.1 0.5];
for i = 1:length(gamma)
    [~, ~, h] = clms(x,y',1,gamma(i));
    
    H = abs(h).^2;
    
    medianH = 50 * median(median(H));
    H(H > medianH) = medianH;
    
    figure;
    surf(n,w,H,'EdgeColor','none');
    view(2);
    c = colorbar;
    c.Label.String = 'PSD (dB)';
    xlabel('Sample')
    ylabel('Frequency (Hz)')
    
end
figure(1)
title('DFT-CLMS Spectrogram \gamma = 0')
figure(2)
title('DFT-CLMS Spectrogram \gamma = 0.01')
figure(3)
title('DFT-CLMS Spectrogram \gamma = 0.1')
figure(4)
title('DFT-CLMS Spectrogram \gamma = 0.5')
%% d)
load 'EEG_Data/EEG_Data_Assignment1.mat'

W = 1200;
w = 0:W-1;
a = 1000;

POz_seg = POz(a:a+W-1);

N = length(POz_seg);
n = 1:N;

x = 1/W * exp(n' * 1j * 2 * pi * w / W)';

gamma = [0 0.001];
for i = 1:length(gamma)
    [~, ~, h] = clms(x,POz_seg',1,gamma(i));
    
    H = abs(h).^2;
    
    medianH = 1000 * median(median(H));
    H(H > medianH) = medianH;
    
    figure;
    surf(n,w,H,'EdgeColor','none');
    view(2);
    c = colorbar;
    c.Label.String = 'PSD (dB)';
    xlabel('Sample')
    ylabel('Frequency (Hz)')
    ylim([0,60])
end

figure(5)
title('EEG DFT-CLMS Spectrogram \gamma = 0')
figure(6)
title('EEG DFT-CLMS Spectrogram \gamma = 0.001')
