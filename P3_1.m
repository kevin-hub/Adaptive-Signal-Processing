%% 3.1 Complex LMS and Widely Linear Modelling
clear; close all; clc;
%% a)
N = 1000;
nRls = 100;
b = [1.5+1j, 2.5-0.5j];
d = zeros(nRls,N);
w = zeros(nRls,N);
h = zeros(nRls,N);
g = zeros(nRls,N);
e_clms = zeros(nRls,N);
e_aclms = zeros(nRls,N);

for i = 1:nRls
    noise = wgn(N,1,0,'complex');
    x = [0 noise'];
    xpp = pp(noise',2,0);

    d = x(2:end) + b(1)*x(1:end-1)+b(2)*conj(x(1:end-1));
    
    [~,e_clms(i,:),~] = clms(xpp, d, 0.1,0);
    [~,e_aclms(i,:),~,~] = aclms(xpp,d,0.1,0);
    
end
e_clms_avg = mean(abs(e_clms).^2);
e_aclms_avg = mean(abs(e_aclms).^2);

figure(1)
hold on
plot(pow2db(e_clms_avg))
plot(pow2db(e_aclms_avg))
title('Learning Curves of CLMS and ACLMS')
xlabel('Iterations (samples)')
ylabel('Squared Error (dB)')
legend('CLMS','ACLMS')
% scatter(real(e_clms_avg),imag(e_clms_avg),'filled')
% scatter(real(e_aclms_avg),imag(e_aclms_avg),'filled')

%% b)

low = load('wind-dataset/low-wind.mat');
medium = load('wind-dataset/medium-wind.mat');
high = load('wind-dataset/high-wind.mat');

wind(1,:) = (low.v_east + low.v_north .* 1j)';
wind(2,:) = (medium.v_east + medium.v_north .* 1j)';
wind(3,:) = (high.v_east + high.v_north .* 1j)';


figure(2)
subplot(3,1,1)
scatter(low.v_east, low.v_north,'fill')
title('Low Wind Regime Scatter Diagram \rho = 0.16')
subplot(3,1,2)
scatter(medium.v_east, medium.v_north,'fill')
title('Medium Wind Regime Scatter Diagram \rho = 0.45')
ylabel('\Im')
subplot(3,1,3)
scatter(high.v_east, high.v_north,'fill')
title('High Wind Regime Scatter Diagram \rho = 0.62')
xlabel('\Re')


N = length(wind);
u = [0.1, 0.01, 0.001];
e_clms = zeros(1,N);
e_aclms = zeros(1,N);
se_clms = zeros(3,30);
se_aclms = zeros(3,30);
for i = 1:3
    for M = 1:30
        x = pp(wind(i,:),M,1);
        
        [~,e_clms(M,:),~] = clms(x, wind(i,:), u(i),0);
        [~,e_aclms(M,:),~,~] = aclms(x,wind(i,:),u(i),0);
        
    end
    se_clms(i,:) = mean(abs(e_clms).^2,2);
    se_aclms(i,:) = mean(abs(e_aclms).^2,2);
end

figure(3)
subplot(3,1,1)
hold on
plot(pow2db(se_clms(1,:)))
plot(pow2db(se_aclms(1,:)))
title('Low Wind Regime Learning Curves')
legend('CLMS','ACLMS')
subplot(3,1,2)
hold on
plot(pow2db(se_clms(2,:)))
plot(pow2db(se_aclms(2,:)))
title('Medium Wind Regime Learning Curves')
ylabel('Squared Error (dB)')
subplot(3,1,3)
hold on
plot(pow2db(se_clms(3,:)))
plot(pow2db(se_aclms(3,:)))
title('High Wind Regime Learning Curves')
xlabel('Order/M')
%% c)
L = 3; % number of phases
N = 1000;
f0 = 50;
fs = 10000;
t = (0:N-1)/fs;

% balanced
amp = ones(L,1); 
phi = [0; -2/3*pi; 2/3*pi]; %phase shift
vabc = amp .* cos(2*pi*f0*t+phi);
v0ab = ct(vabc);
vclarke_bal = v0ab(2,:) + 1j * v0ab(3,:);

figure(4)
hold on
scatter(real(vclarke_bal),imag(vclarke_bal),'fill','k')

%unbalanced
da = 0.1:0.1:0.5;
dp = pi*(-0.2: 0.1: 0.2);

for i = 1:length(da)
    vabc = (amp + [-da(i); 0; da(i)]) .* cos(2*pi*f0*t+phi);
    v0ab = ct(vabc);
    vclarke = v0ab(2,:) + 1j * v0ab(3,:);
    
    scatter(real(vclarke),imag(vclarke))
    rho_da(i) = circularity(vclarke); 
end
xlim([-2 2])
ylim([-2 2])
title('Circularity Diagram of Balanced Power System')
xlabel('\Re')
ylabel('\Im')
legend({'Balanced','\Delta V= 0.1 \rho = 0.12','\Delta V= 0.2 \rho = 0.23','\Delta V= 0.3 \rho = 0.34','\Delta V= 0.4 \rho = 0.44','\Delta V= 0.5 \rho = 0.53'},'Location','southeast');

figure(5)
hold on
scatter(real(vclarke_bal),imag(vclarke_bal),'fill','k')

for i = 1:length(dp)
    vabc = amp .* cos(2*pi*f0*t+phi+[0; dp(i); dp(i)]);
    v0ab = ct(vabc);
    vclarke = v0ab(2,:) + 1j * v0ab(3,:);
    if i ~= 3
        scatter(real(vclarke),imag(vclarke))
        rho_dp(i) = circularity(vclarke); 
    end
end
xlim([-2 2])
ylim([-2 2])
title('Circularity Diagram of Unbalanced Power System')
xlabel('\Re')
ylabel('\Im')
legend({'Balanced','\Delta\phi = -0.2\pi \rho = 0.41','\Delta\phi = -0.1\pi \rho = 0.21','\Delta\phi = 0.1\pi \rho = 0.21','\Delta\phi = 0.2\pi \rho = 0.41'},'Location','southeast');

%% e)
fs = 1000; %lower sampling frequency
t = (0:N-1)/fs;
amp = ones(L,1); 
phi = [0; -2/3*pi; 2/3*pi]; %phase shift
vabc = amp .* cos(2*pi*f0*t+phi);
v0ab = ct(vabc);
vclarke_bal = v0ab(2,:) + 1j * v0ab(3,:);
[ppvclarke_bal] = pp(vclarke_bal,1,1);
[~,eab_bal_clms,hab_bal_clms] = clms(ppvclarke_bal,vclarke_bal,0.05,0);
[~,eab_bal_aclms,hab_bal_aclms,gab_bal_aclms] = aclms(ppvclarke_bal,vclarke_bal,0.05,0);
f0_bal_clms = abs(fs/(2*pi) * atan(imag(hab_bal_clms)./real(hab_bal_clms)));
f0_bal_aclms = abs(fs/(2*pi) * atan(sqrt((imag(hab_bal_aclms)).^2-abs(gab_bal_aclms).^2)./real(hab_bal_aclms)));

figure(6)
subplot(2,2,1)
hold on
plot([1 1000],[50 50],'k')
plot(f0_bal_clms)
plot(f0_bal_aclms)
title('Balanced Nominal Frequency Estimation') 
ylim([0 100])
ylabel('Frequency (Hz)');
legend('True','CLMS','ACLMS')

subplot(2,2,2)
hold on
plot(pow2db(abs(eab_bal_clms.^2)))
plot(pow2db(abs(eab_bal_aclms.^2)))
title('Balanced Learning Curves')
legend('CLMS','ACLMS')

da = 0.1; dp = 0.2;
vabc = (amp + [-da; 0; da]) .* cos(2*pi*f0*t+phi+[0; dp; dp]);
v0ab = ct(vabc);
vclarke_unbal = v0ab(2,:) + 1j * v0ab(3,:);
[ppvclarke_unbal] = pp(vclarke_unbal,1,1);
[~,eab_unbal_clms,hab_unbal_clms] = clms(ppvclarke_unbal,vclarke_unbal,0.05,0);
[~,eab_unbal_aclms,hab_unbal_aclms,gab_unbal_aclms] = aclms(ppvclarke_unbal,vclarke_unbal,0.05,0);
f0_unbal_clms = abs(fs/(2*pi) * atan(imag(hab_unbal_clms)./real(hab_unbal_clms)));
f0_unbal_aclms = abs(fs/(2*pi) * atan(sqrt((imag(hab_unbal_aclms)).^2-abs(gab_unbal_aclms).^2)./real(hab_unbal_aclms)));

subplot(2,2,3)
hold on
plot([1 1000],[50 50],'k')
plot(f0_unbal_clms)
plot(f0_unbal_aclms)
title('Unbalanced Nominal Frequency Estimation') 
ylim([0 100])
xlabel('Time (sample)');
ylabel('Frequency (Hz)');

subplot(2,2,4)
hold on
plot(pow2db(abs(eab_unbal_clms.^2)))
plot(pow2db(abs(eab_unbal_aclms.^2)))
title('Unbalanced Learning Curves') 
xlabel('Time (sample)');