%% 2.2 Adaptive Step Sizes
clear all; close all; clc;

%% a)

fs = 1;
nSamples = 1000;
nRps = 100;
nTrans = 500; %c)

coeff = 0.9;
w0 = coeff;
var = 0.5;

u = [0.01 0.1];
u0 = 0;

MAmdl = arima('MA', coeff, 'Constant', 0, 'Variance', var);
[MAsig,innovation] = simulate(MAmdl,nSamples,'NumPaths',nRps);

inputs = innovation'; % white noise i.e. input into the system

e = zeros(nRps,nSamples);
e_avg = zeros(length(u),nSamples);
w2 = zeros(nRps, nSamples);
w_avg = zeros(length(u),nSamples);

benv.name = 'B';
angfar.name = 'AF';
angfar.param = 0.8;
mattxie.name = 'MX';

e_benv = zeros(nRps,nSamples);
e_angfar = zeros(nRps,nSamples);
e_mattxie = zeros(nRps,nSamples);

w2_benv = zeros(nRps, nSamples);
w2_angfar = zeros(nRps, nSamples);
w2_mattxie = zeros(nRps, nSamples);

for k = 1:length(u)
    for i = 1:nRps
        input1 = inputs(i,:);
        input2 = [0 inputs(i,1:end-1)];
        
        x = [input1' input2']';
        
        [~,e(i,:),w] = lms(x,MAsig(:,i)',u(k),0);
    
        w2(i,:) = w(2,:); % only store second coefficient i.e. the delay input with a coeff = 0.9
        %the first order term n(n) has a constant coeff of 1 and thus the
        %adaptive process is not needed/or of interest
        if k == 1
            [~,e_benv(i,:),w_benv] = gass(x,MAsig(:,i)',u0,0,5e-3,benv);
            [~,e_angfar(i,:),w_angfar] = gass(x,MAsig(:,i)',u0,0,5e-3,angfar);
            [~,e_mattxie(i,:),w_mattxie] = gass(x,MAsig(:,i)',u0,0,5e-3,mattxie);
            
            w2_benv(i,:) = w_benv(2,:);
            w2_angfar(i,:) = w_angfar(2,:);
            w2_mattxie(i,:) = w_mattxie(2,:);

        end
    end
    e_avg(k,:) = mean(e.^2);
    w_avg(k,:) = mean(w2);
    
    if k == 1
        w_benv_avg = mean(w2_benv);
        w_angfar_avg = mean(w2_angfar);
        w_mattxie_avg = mean(w2_mattxie);
        e_benv_avg = mean(e_benv.^2);
        e_angfar_avg = mean(e_angfar.^2);
        e_mattxie_avg = mean(e_mattxie.^2);
    end
end

figure(1)
subplot(2,1,1)
hold on
plot(w0 - w_avg(1,:))
plot(w0 - w_avg(2,:))
plot(w0 - w_benv_avg)
plot(w0 - w_angfar_avg)
plot(w0 - w_mattxie_avg)
title('Weight Error of Fixed and Adaptive Learning Rates')
xlabel('Iteration (Samples)')
ylabel('Weight Error')
legend('\mu = 0.01','\mu = 0.1','Benveniste','Ang & Farhang','Matthews & Xie')

subplot(2,1,2)
hold on
plot(10*log10(e_avg(1,:)))
plot(10*log10(e_avg(2,:)))
plot(10*log10(e_benv_avg))
plot(10*log10(e_angfar_avg))
plot(10*log10(e_mattxie_avg))
title('Error of LMS Predictor with Fixed and Adaptive Learning Rates')
xlabel('Iteration (Samples)')
ylabel('Squared Error (dB)')

%% c)

for i = 1:nRps
    input1 = inputs(i,:);
    input2 = [0 inputs(i,1:end-1)];

    x = [input1' input2']';

    [~,e_gnsd(i,:),w_gnsd] = gnsd(x,MAsig(:,i)',2,0,5e-3);
    
    w2_gnsd(i,:) = w_gnsd(2,:);

end

e_gnsd_avg = mean(e_gnsd.^2);
w_gnsd_avg = mean(w2_gnsd);

figure(3)
subplot(2,1,1)
hold on
plot(w0-w_gnsd_avg)
plot(w0-w_benv_avg)
xlim([0 100])
title('Weight Error by GNSD and Benveniste GASS Algorithm')
xlabel('Iteration (Samples)')
ylabel('Weight Error')
legend('GNSD','Benveniste')

subplot(2,1,2)
hold on
plot(10*log10(e_gnsd_avg))
plot(10*log10(e_benv_avg))
xlim([0 300])
title('Error of LMS Predictor Using GNSD and Benveniste GASS Algorithm')
xlabel('Iteration (Samples)')
ylabel('Squared Error (dB)')
