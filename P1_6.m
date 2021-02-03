%% 1.6 Robust Regression
clear all; close all; clc;

load PCR/PCAPCR.mat
addpath '/Users/Kevin/Documents/MATLAB/Year 4/Adaptive Signal Processing and Machine Intelligence/PCR'

%% a)
s_X = svd(X);
s_Xnoise = svd(Xnoise);

r = rank(X);

figure(1)
subplot(2,1,1)
hold on
stem(s_Xnoise)
stem(s_X)
title('Singular Values of Clean and Noisy Input Training Data')
xlabel('Input Variable')
ylabel('Singular Value')
legend('Noisy', 'Clean')

sqrErr = (s_X - s_Xnoise).^2;

subplot(2,1,2)
stem(sqrErr,'r')
title('Square Error')
xlabel('Input Variable')
ylabel('Error')

%% b)
[U_Xnoise, S_Xnoise, V_Xnoise] = svd(Xnoise);
Xdenoise = U_Xnoise(:,1:r) * S_Xnoise(1:r,1:r) * V_Xnoise(:,1:r).';

errNoise = abs(vecnorm(X - Xnoise));
errDenoise = abs(vecnorm(X - Xdenoise));

figure(2)
hold on
stem(errNoise, 'r')
stem(errDenoise, 'Color', '#A2142F')
title('Difference Error') %Difference Between Clean Input Training Data to its Noisy and Denoised Versions
xlabel('Input Variable')
ylabel('Error')
legend('Noisy', 'Denoised')

%% c)
B_OLS = inv(Xnoise.' * Xnoise) * X.' * Y;
B_PCR = V_Xnoise(:,1:r) * inv(S_Xnoise(1:r,1:r)) * U_Xnoise(:,1:r).' * Y;

Y_OLS = Xnoise * B_OLS;
Y_PCR = Xdenoise * B_PCR;
Ytest_OLS = Xtest * B_OLS;

[U_Xtest, S_Xtest, V_Xtest] = svd(Xtest);
XtestPCA = U_Xtest(:,1:r) * S_Xtest(1:r,1:r) * V_Xtest(:,1:r).';

Ytest_PCR = XtestPCA * B_PCR;

errOLS = abs(vecnorm(Y-Y_OLS));
errPCR = abs(vecnorm(Y-Y_PCR));
errOLStest = abs(vecnorm(Y-Ytest_OLS));
errPCRtest = abs(vecnorm(Y-Ytest_PCR));


figure(3)
subplot(2,1,1)
hold on
stem(errOLS)
stem(errPCR,'x','--')
title('Solution Difference Error for Training Data')
xlabel('Output Variable')
ylabel('Error')
legend('OLS', 'PCR')

subplot(2,1,2)
hold on
stem(errOLStest)
stem(errPCRtest,'x','--')
title('Solution Difference Error for Testing Data')
xlabel('Output Variable')
ylabel('Error')

[YrOLS, YeOLS] = regval(B_OLS);
[YrPCR, YePCR] = regval(B_PCR);

errMSE_OLS = mean(abs(YrOLS - YeOLS)).^2;
errMSE_PCR = mean(abs(YrPCR - YePCR)).^2;

figure(4)
hold on
stem(errMSE_OLS)
stem(errMSE_PCR, 'x','--')
title('Solution MSE for an Ensemble Test Data')
xlabel('Output Variable')
ylabel('Error')
legend('OLS', 'PCR')