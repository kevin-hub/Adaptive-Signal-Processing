%% 4 From LMS to Deep Learning
clear; close all; clc;
%% 1)
load 'time-series.mat'

y = y - mean(y);
N = length(y);
u = 1e-5;

ypp = pp(y',4,1);

[p,e,~] = lms(ypp,y,1e-5,0);

mse = mean(abs(e).^2);
Rp = pow2db(var(p)/var(e));

figure(1)
hold on
plot(y,'k')
plot(p)
legend('Zero-mean','Prediction')
title('Zero-mean and One-step Prediction of Time-series Using LMS')
xlabel('Time (samples)')
ylabel('Amplitude')

%% 2)
[p_tanh,e_tanh,~] = tanhptron(ypp,y,1e-5,0,1);

mse_tanh = mean(abs(e_tanh).^2);
Rp_tanh = pow2db(var(p_tanh)/var(e_tanh));

figure(2)
hold on
plot(y,'k')
plot(p_tanh)
legend('Zero-mean','tanh Prediction')
title('Zero-mean and One-step Prediction of Time-series Using tanh-LMS')
xlabel('Time (samples)')
ylabel('Amplitude')

%% 3)
a = 1:100;
mse_tanhamp = zeros(1,length(a));
Rp_tanhamp = zeros(1,length(a));
p_tanhamp = zeros(length(a),N);
for i = a
    [p_tanhamp(i,:),e_tanhamp,~] = tanhptron(ypp,y,1e-5,0,i);

    mse_tanhamp(i) = mean(abs(e_tanhamp).^2);
    Rp_tanhamp(i) = pow2db(var(p_tanhamp(i,:))/var(e_tanhamp));
    
end
figure(3)
title('Error and Gain of Time-series Using Scaled tanh-LMS')
yyaxis left
plot(mse_tanhamp)
ylabel('Squared Error (dB)')
yyaxis right
plot(Rp_tanhamp)
ylabel('Prediction Gain (dB)')

[~,a_opt] = min(mse_tanhamp);

figure(4)
hold on
plot(y,'k')
plot(p_tanhamp(a_opt,:))
legend('Zero-mean','Prediction')
title('Zero-mean and One-step Prediction of Time-series Using Scaled tanh-LMS')
xlabel('Time (samples)')
ylabel('Amplitude')

%% 4)
load 'time-series.mat'

ypp = pp(y',4,1);
augy = [ones(1,length(ypp)); ypp];

for i = a
    [p_bias(i,:),e_bias,~] = tanhptron(augy,y,1e-5,0,i);

    mse_bias(i) = mean(abs(e_bias).^2);
    Rp_bias(i) = pow2db(var(p_bias(i,:))/var(e_bias));
end

figure(7)
hold on
plot(y,'k')
plot(p_bias(a_opt,:))
legend('Zero-mean','Prediction')
title('Zero-mean and One-step Prediction of Time-series Using Biased tanh-LMS')
xlabel('Time (samples)')
ylabel('Amplitude')

figure(5)
title('Error and Gain of Time-series Using Biased tanh-LMS')
yyaxis left
plot(mse_bias)
ylabel('Squared Error (dB)')
yyaxis right
plot(Rp_bias)
ylabel('Prediction Gain (dB)')

%% 5)

batch = repmat(y(1:20),1,100);
batchpp = pp(batch,4,1);

for i = a
    [~,~,w_pre] = tanhptron(augy,y,1e-5,0,i);
    
    w_init(:,i) = w_pre(:,end); 
    
    [p_init(i,:),e_init,~] = tanhptron(augy,y,1e-5,0,a(i),w_init(:,i));

    mse_init(i) = mean(abs(e_init).^2);
    Rp_init(i) = pow2db(var(p_init(i,:))/var(e_init));
end

figure(6)
title('Error and Gain of Time-series Using Batched tanh-LMS')
yyaxis left
plot(mse_init)
ylabel('Squared Error (dB)')
yyaxis right
plot(Rp_init)
ylabel('Prediction Gain (dB)')

figure(8)
hold on
plot(y,'k')
plot(p_init(a_opt,:))
legend('Zero-mean','Prediction')
title('Zero-mean and One-step Prediction of Time-series Using Batched tanh-LMS')
xlabel('Time (samples)')
ylabel('Amplitude')