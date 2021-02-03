%% 1.4 Spectrum of Autoregressive Processes
clear all; close all; clc
%% b)
fs = 1;
nSamples = 1000;
nTrans = 500;
coeff = {2.76, -3.81, 2.65, -0.92};

Mdl = arima('AR', coeff, 'Constant', 0, 'Variance', 1);
Y = simulate(Mdl,nSamples);
Y = Y(nTrans+1:end);
n = length(Y);

[h, w] = freqz(1,[1 -cell2mat(coeff)],n);
psdTrue = abs(h).^2;

figure(1)
hold on
plot(w/pi,mag2db(psdTrue),'k','linewidth',2)
for p = 2:14
    ToEstMdl = arima(p,0,0);
    ToEstMdl.Constant = 0;
    EstMdl = estimate(ToEstMdl,Y);
    coeffEst{p} = EstMdl.AR;
    varEst(p) = EstMdl.Variance;
    h = freqz(sqrt(EstMdl.Variance),[1 -cell2mat(EstMdl.AR)],n);
    psdEst{p} = abs(h).^2;
    psdErr(p) = mean((psdTrue-psdEst{p}).^2);
    plot(w/pi, mag2db(psdEst{p}),'linewidth',1.1)
end
legend('Ground','p=2','p=3','p=4','p=5','p=6','p=7','p=8','p=9','p=10','p=11','p=12','p=13','p=14')

%plot of best psd estimates
figure(2)
subplot(2,1,1)
hold on
plot(w/pi,mag2db(psdTrue),'k','linewidth',2)
plot(w/pi,mag2db(psdEst{4}), 'linewidth', 1.1) 
xlim([0.2 0.32])
ylim([50 100])
title('PSD Estimates Using AR Models')
xlabel('Normalised Frequency (\pi rad/sample)')
ylabel('PSD (dB)')
legend('Ground','p=4')
subplot(2,1,2)
hold on
plot(mag2db(varEst))
plot(mag2db(psdErr))
xlim([2 14])


%% c)
nSamples = 1e4;

YHRes = simulate(Mdl,nSamples);
YHRes = YHRes(nTrans+1:end);
nH = length(YHRes);

[h, w] = freqz(1,[1 -cell2mat(coeff)],nH);
psdTrueHRes = abs(h).^2;

figure(3)
hold on
plot(w/pi,mag2db(psdTrueHRes),'k','linewidth',2)

for p = 2:14
    ToEstMdlHRes = arima(p,0,0);
    ToEstMdlHRes.Constant = 0;
    EstMdlHRes = estimate(ToEstMdlHRes,YHRes);
    coeffEstHRes{p} = EstMdlHRes.AR;
    varEstHRes(p) = EstMdlHRes.Variance;
    h = freqz(sqrt(EstMdlHRes.Variance),[1 -cell2mat(EstMdlHRes.AR)],nH);
    psdEstHRes{p} = abs(h).^2;
    psdErrHRes(p) = mean((psdTrueHRes-psdEstHRes{p}).^2);
    plot(w/pi, mag2db(psdEstHRes{p}),'linewidth',1.1)
end

figure(4)
subplot(2,1,1)
hold on
plot(w/pi,mag2db(psdTrueHRes),'k','linewidth',2)
title('High Resolution PSD Estimates Using AR Models')
xlabel('Normalised Frequency (\pi rad/sample)')
ylabel('PSD (dB)')
legend('Ground','p=3','p=4','p=8')

subplot(2,1,2)
hold on
plot(mag2db(varEstHRes))
plot(mag2db(psdErrHRes))

