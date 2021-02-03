%% 2.1 The Least Mean Square (LMS) Algorithm
clear all; close all; clc;

%% b)
fs = 1;
nSamples = 1000;
nRls = 100;
nTrans = 500; %c)

coeff = [0.1, 0.8];
var = 0.25;
u = [0.05 0.01];

ARmdl = arima('AR', coeff, 'Constant', 0, 'Variance', var);
ARsig = simulate(ARmdl,nSamples,'NumPaths',nRls);


e = zeros(nRls,nSamples);
avg_err = zeros(length(u),nSamples);
sing_err = zeros(length(u),nSamples);
w1 = zeros(nRls,nSamples);
w2 = zeros(nRls,nSamples);
w1_avg = zeros(length(u),nSamples);
w2_avg = zeros(length(u),nSamples);

for k = 1:length(u)
    for i = 1:nRls
        delay1 = [0 ARsig(:,i)'];
        input1 = delay1(1:end-1)';
        delay2 = [0 0 ARsig(:,i)'];
        input2 = delay2(1:end-2)';

        x = [input1 input2]';

        [~,e(i,:),w] = lms(x,ARsig(:,i)',u(k),0); %need to create own LMS filter which allows for multiple inputs
        mse(i) = mean(e(i,(nTrans+1):end).^2); %c)
        w1(i,:) = w(1,:); %d)
        w2(i,:) = w(2,:);
    end
    avg_err(k,:) = mean(e.^ 2);
    sing_err(k,:) = e(1,:).^2;
    ms_err(k,:) = mse; %c)
    w1_avg(k,:) = mean(w1); %d)
    w2_avg(k,:) = mean(w2);
end

figure(1)
hold on
plot(10*log10(sing_err(1,:)))
plot(10*log10(sing_err(2,:)))
title('Error of LMS Adaptive Predictor')
xlabel('Samples')
ylabel('Squared Error (dB)')
legend('\mu = 0.05','\mu = 0.01' )

figure(2)
hold on
plot(10*log10(avg_err(1,:)))
plot(10*log10(avg_err(2,:)))
title('Average Error of 100 LMS Adaptive Predictions')
xlabel('Samples')
ylabel('Squared Error (dB)')
legend('\mu = 0.05','\mu = 0.01' )

%% c)
mse_ensemble = mean(ms_err,2); %average of ensemble(/two) learning curves
emse = mse_ensemble - var;
M = emse/var;

Rxx = [25/27, 25/54; 25/54, 25/27];
M_approx = u/ 2 * trace(Rxx);

%% d)
coeff1005_est = mean(w1_avg(1,(nTrans+1):end),2);
coeff2005_est = mean(w2_avg(1,(nTrans+1):end),2);

coeff1001_est = mean(w1_avg(2,(nTrans+1):end),2);
coeff2001_est = mean(w2_avg(2,(nTrans+1):end),2);

figure(3)
hold on
plot([0 999],[coeff(1) coeff(1)],'Color','#7E2F8E','LineStyle','--')
plot([0 999],[coeff(2) coeff(2)],'r--')
plot(w1_avg(1,:),'Color','#0072BD')
plot(w2_avg(1,:),'Color','#D95319')
title('Coefficient Estimates for \mu = 0.05')
xlabel('Samples')
ylabel('Averaged Weights')
ylim([0 1])

figure(4)
hold on
plot([0 999],[coeff(1) coeff(1)],'Color','#7E2F8E','LineStyle','--')
plot([0 999],[coeff(2) coeff(2)],'r--')
plot(w1_avg(2,:),'Color','#0072BD')
plot(w2_avg(2,:),'Color','#D95319')
title('Coefficient Estimates for \mu = 0.01')
xlabel('Samples')
ylabel('Averaged Weights')
lgd = legend('$$\alpha_1$$','$$\alpha_2$$','$$\hat{\alpha}_1$$','$$\hat{\alpha}_2$$','Interpreter','Latex');
lgd.NumColumns = 4;
ylim([0 1])
%% f)

leak = [0.25 0.5 0.75];

a1_est = cell([1 length(leak)]);
a2_est = cell([1 length(leak)]);

for j = 1:length(leak)
    for k = 1:length(u)
        for i = 1:nRls
            delay1 = [0 ARsig(:,i)'];
            input1 = delay1(1:end-1)';
            delay2 = [0 0 ARsig(:,i)'];
            input2 = delay2(1:end-2)';

            x = [input1 input2]';

            [~,~,w] = lms(x,ARsig(:,i)',u(k),leak(j));
            w1(i,:) = w(1,:); % for 2nd order
            w2(i,:) = w(2,:);
        end
        w1_avg(k,:) = mean(w1); 
        w2_avg(k,:) = mean(w2);
    end
    a1_est{j} = w1_avg;
    a2_est{j} = w2_avg;
end

figure(5)
subplot(3,2,1)
hold on
plot(a1_est{1}(1,:))
plot(a2_est{1}(1,:))
plot(coeff(1) * ones(1,nSamples),'Color','#7E2F8E','LineStyle','--')
plot(coeff(2)* ones(1,nSamples),'r--')
ylim([0 1])
title('Coefficient Estimates for \mu = 0.05, \gamma = 0.25')
subplot(3,2,3)
hold on
plot(a1_est{2}(1,:))
plot(a2_est{2}(1,:))
plot(coeff(1) * ones(1,nSamples),'Color','#7E2F8E','LineStyle','--')
plot(coeff(2)* ones(1,nSamples),'r--')
ylim([0 1])
title('                                         \mu = 0.01, \gamma = 0.5')
ylabel('Averaged Weights')
subplot(3,2,5)
hold on
plot(a1_est{3}(1,:))
plot(a2_est{3}(1,:))
plot(coeff(1) * ones(1,nSamples),'Color','#7E2F8E','LineStyle','--')
plot(coeff(2)* ones(1,nSamples),'r--')
ylim([0 1])
title('                                        \mu = 0.01, \gamma = 0.75')
xlabel('Samples')
subplot(3,2,2)
hold on
plot(a1_est{1}(2,:))
plot(a2_est{1}(2,:))
plot(coeff(1) * ones(1,nSamples),'Color','#7E2F8E','LineStyle','--')
plot(coeff(2)* ones(1,nSamples),'r--')
ylim([0 1])
title('Coefficient Estimates for \mu = 0.01, \gamma = 0.25')
subplot(3,2,4)
hold on
plot(a1_est{2}(2,:))
plot(a2_est{2}(2,:))
plot(coeff(1) * ones(1,nSamples),'Color','#7E2F8E','LineStyle','--')
plot(coeff(2)* ones(1,nSamples),'r--')
ylim([0 1])
title('                                         \mu = 0.01, \gamma = 0.5')
ylabel('Averaged Weights')
subplot(3,2,6)
hold on
plot(a1_est{3}(2,:))
plot(a2_est{3}(2,:))
plot(coeff(1) * ones(1,nSamples),'Color','#7E2F8E','LineStyle','--')
plot(coeff(2)* ones(1,nSamples),'r--')
ylim([0 1])
title('                                        \mu = 0.01, \gamma = 0.75')
xlabel('Samples')
