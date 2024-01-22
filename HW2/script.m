%% 1.1
close all; clear; clc;
beta = 0.5;
T = 1;
Fs = 10;
dt = T/Fs;



% Scen1: Ideal
t = -6*T:dt:6*T;
x_t1 = sinc(t / T) .* cos(pi * beta * t / T) ./ (1 - (2 * beta * t / T).^2);
x_t1(+T/(2*beta)/dt + ceil(length(t)/2)) = pi/4 * sinc(1/(2*beta));
x_t1(-T/(2*beta)/dt + ceil(length(t)/2)) = pi/4 * sinc(1/(2*beta));

% Scen2: e = 0.1T
e = 0.1;
t = (-6-e)*T:dt:(6-e)*T;
x_t2 = sinc(t / T) .* cos(pi * beta * t / T) ./ (1 - (2 * beta * t / T).^2);
x_t2(T/(2*beta)/dt + ceil(length(t)/2) + Fs*e) = pi/4 * sinc(1/(2*beta));

% Scen2: e = 0.2T
e = 0.2;
t = (-6-e)*T:dt:(6-e)*T;
x_t3 = sinc(t / T) .* cos(pi * beta * t / T) ./ (1 - (2 * beta * t / T).^2);
x_t3(-T/(2*beta)/dt + ceil(length(t)/2) + Fs*e) = pi/4 * sinc(1/(2*beta));


figure;
plot(x_t1);
hold on;
plot(x_t2);
hold on;
plot(x_t3);
title('Raised Cosine Pulse');
xlabel('Time');
ylabel('Amplitude');
legend('\beta=1 T=1 \epsilon=0','\beta=1 T=1 \epsilon=0.1T', '\beta=1 T=1 \epsilon=0.2T')
grid on;


%% 1.2
N = 10e5;
bits = randi([0,1],[1,N]);
modulated_symbols = 2*bits - 1;
L = T*Fs;
Zero_Pad_MS = zeros(L*length(bits), 1);
Zero_Pad_MS(1:L:end) = modulated_symbols;

transmitted_signal1 = conv(Zero_Pad_MS, x_t1);
transmitted_signal2 = conv(Zero_Pad_MS, x_t2);
transmitted_signal3 = conv(Zero_Pad_MS, x_t3);

figure;
plot(transmitted_signal1(1:1000));
hold on;
plot(transmitted_signal2(1:1000));
hold on;
plot(transmitted_signal3(1:1000));
title('Modulated Signal (First 1000 Symbols)');
xlabel('N_{Symbol}');
ylabel('Amplitude');
legend('\epsilon=0','\epsilon=0.1T', '\epsilon=0.2T')
grid on;
%% 1.3
SNR_dB = 0:1:10;

transmitted_signal_power1 = sum(abs(transmitted_signal1).^2) / length(transmitted_signal1);
transmitted_signal_power2 = sum(abs(transmitted_signal2).^2) / length(transmitted_signal2);
transmitted_signal_power3 = sum(abs(transmitted_signal3).^2) / length(transmitted_signal3);

noise_power1 = transmitted_signal_power1 ./ (10.^(SNR_dB./10));
noise_power2 = transmitted_signal_power2 ./ (10.^(SNR_dB./10));
noise_power3 = transmitted_signal_power3 ./ (10.^(SNR_dB./10));

noise1 = sqrt(noise_power1/2) .* randn(size(transmitted_signal1));
noise2 = sqrt(noise_power2/2) .* randn(size(transmitted_signal2));
noise3 = sqrt(noise_power3/2) .* randn(size(transmitted_signal3));

received_signal1 = transmitted_signal1 + noise1;
received_signal2 = transmitted_signal2 + noise2;
received_signal3 = transmitted_signal3 + noise3;


figure;
subplot(3,1,1);
plot(received_signal1(1:1000, 1));
title('Received Signal (First 1000 Symbols) For SNR=10 and \epsilon = 0');
xlabel('N_{Symbol}');
ylabel('Amplitude');
grid on;

subplot(3,1,2);
plot(received_signal2(1:1000, 1));
title('Received Signal (First 1000 Symbols) For SNR=10 and \epsilon = 0.1T');
xlabel('N_{Symbol}');
ylabel('Amplitude');
grid on;

subplot(3,1,3);
plot(received_signal3(1:1000, 1));
title('Received Signal (First 1000 Symbols) For SNR=10 and \epsilon = 0.2T');
xlabel('N_{Symbol}');
ylabel('Amplitude');
grid on;
%% 1.4
T_sampling = 6*L+1 : L : (N+6-1)*L+1;
Threshold = 0;
detected_symbols1 = received_signal1(T_sampling,:) > Threshold;
detected_symbols2 = received_signal2(T_sampling,:) > Threshold;
detected_symbols3 = received_signal3(T_sampling,:) > Threshold;


%% 1.5
errors1 = sum(abs((detected_symbols1 - bits')))./N;
errors2 = sum(abs((detected_symbols2 - bits')))./N;
errors3 = sum(abs((detected_symbols3 - bits')))./N;

figure;
semilogy(errors1, '-o', 'LineWidth',3);
hold on;
semilogy(errors2, '-o', 'LineWidth',3);
hold on;
semilogy(errors3, '-o', 'LineWidth',3);
grid on;
title('Binary PAM Detection Error: \beta = 1');
xlabel('SNR(dB)');
ylabel('log(Error)');
legend('\epsilon=0','\epsilon=0.1T', '\epsilon=0.2T')
grid on;




%% 1.6
clear; clc;
beta = 0;
N = 10e5;
T = 1;
Fs = 10;
dt = T/Fs;
t = -6*T:dt:6*T;
x_t1 = sinc(t / T) .* cos(pi * beta * t / T) ./ (1 - (2 * beta * t / T).^2);
% x_t1(+T/(2*beta)/dt + ceil(length(t)/2)) = pi/4 * sinc(1/(2*beta));
% x_t1(-T/(2*beta)/dt + ceil(length(t)/2)) = pi/4 * sinc(1/(2*beta));

probs = rand(1, N);

modulated_symbols = zeros(1, N);
modulated_symbols(probs >= 0    & probs < 0.1) = -3;
modulated_symbols(probs >= 0.1  & probs < 0.5) = -1;
modulated_symbols(probs >= 0.5  & probs < 0.9) = 1;
modulated_symbols(probs >= 0.9  & probs <=  1) = 3;
L = T*Fs;
Zero_Pad_MS = zeros(L*length(probs), 1);
Zero_Pad_MS(1:L:end) = modulated_symbols;

transmitted_signal = conv(Zero_Pad_MS', x_t1);

figure;
plot(transmitted_signal(1:1000));
title('4-PAM Modulated Signal (First 1000 Symbols)');
xlabel('N_{Symbol}');
ylabel('Amplitude');
grid on;

SNR_dB = 0:1:10;
transmitted_signal_power = sum(abs(transmitted_signal).^2) / length(transmitted_signal);
noise_power = transmitted_signal_power ./ (10.^(SNR_dB./10));
noise = sqrt(noise_power'/2) .* randn(size(transmitted_signal));

received_signal = transmitted_signal + noise;

figure;
plot(received_signal(6, 1:1000));
title('Received Signal (First 1000 Symbols) For SNR = 6');
xlabel('N_{Symbol}');
ylabel('Amplitude');
grid on;




%% 1.7
T_sampling = 6*L+1 : L : (N+6-1)*L+1;
samples = received_signal(:,T_sampling);

detected_symbols_ML = zeros(11, N);
detected_symbols_ML(samples <= -2) = -3;
detected_symbols_ML(samples > -2   &   samples <= 0) = -1;
detected_symbols_ML(samples > 0    &   samples <= 2) = +1;
detected_symbols_ML(samples >  +2) = +3;

detected_symbols_MAP = zeros(11, N);
MAP_diff_thr = noise_power.^2*log(3)/2;
detected_symbols_MAP(samples <= -2-MAP_diff_thr') = -3;
detected_symbols_MAP(samples > -2-MAP_diff_thr'   &   samples <= 0) = -1;
detected_symbols_MAP(samples > 0    &   samples <= 2+MAP_diff_thr') = +1;
detected_symbols_MAP(samples >  +2+MAP_diff_thr') = +3;


errors_ML  = sum(abs((detected_symbols_ML'  - modulated_symbols')))./N;
errors_MAP = sum(abs((detected_symbols_MAP' - modulated_symbols')))./N;

figure;
semilogy(errors_ML, '-o', 'LineWidth',3);
hold on;
semilogy(errors_MAP, '-o', 'LineWidth',3);
title('4-PAM Detection Error: ML and MAP Detector');
xlabel('SNR(dB)');
ylabel('log(Error)');
legend('ML','MAP')
grid on;

