%% 1.1
close all; clear; clc;
rng default
N_bits = 10e4;
N_vecs = 20;
N_msgs = N_bits/N_vecs;
vectors = randi([0 1],N_msgs, N_vecs);



%% 2
intIndices  = randperm(N_msgs);
turboenc = comm.TurboEncoder('InterleaverIndices', intIndices);
N_codded = 3*N_msgs + 12;
codedVectors = [];
for i=1:N_vecs   
    codedVectors = [codedVectors, turboenc(vectors(:,i))];
end

%% 1.2
qpskQI = [+1+1i, +1-1i, -1+1i, -1-1i] / sqrt(2);
bitPairs = reshape(codedVectors, 2, []);
constellationMap = bi2de(bitPairs', 'left-msb');
qpskSymbols = zeros(1, N_codded);
for i=1:N_codded * N_vecs/2
   qpskSymbols(i) = qpskQI(constellationMap(i) + 1);
end
transmittedSignals = reshape(qpskSymbols, [N_codded/2, N_vecs]);

%% 1.3
nfft  = N_codded/2; 
cplen = 16;  
OFDMSignals = ofdmmod(transmittedSignals,nfft,cplen);

%% 1.4
SNR_dB = 10:0.5:17;
TxSigPower = sum(abs(OFDMSignals).^2) / length(OFDMSignals);
noise_power = TxSigPower ./ (10.^(SNR_dB./10));
noise = sqrt(noise_power*2) .* randn(size(OFDMSignals));
RxSignal = [];
for i=1:length(SNR_dB)
    RxSignal = [RxSignal, OFDMSignals + noise(:,i)];
end


figure;
sgtitle('Received Signal (First 1000 Symbols) For SNR=1 and SNR=10');
subplot(2,1,1);
plot(abs(RxSignal(1:1000,1)));
xlabel('N_{Symbol}');
ylabel('Amplitude');
grid on;
hold on;
subplot(2,1,2);
plot(abs(RxSignal(1:1000,10)));
hold off;
xlabel('N_{Symbol}');
ylabel('Amplitude');
grid on;

%% 1.5
y_1 = ofdmdemod(RxSignal,nfft,cplen);

quadrature = real(y_1)>0;
inphase    = imag(y_1)>0;
positions  = 3 - 2*quadrature - inphase;
y_2 = zeros(N_codded, N_vecs, length(SNR_dB));
turboDec = comm.TurboDecoder('InterleaverIndices', intIndices);
y_3 = zeros(N_msgs, N_vecs, length(SNR_dB));
for i=1:N_vecs
    for j=1:length(SNR_dB)
        temp = de2bi(positions(:,i,j), 'left-msb');
        y_2(:,i,j) = reshape(temp',[N_codded,1]);
        y_3(:,i,j) = turboDec(y_2(:,i,j));
    end
end

%% 1.6
error = sum(abs((y_3 - vectors)))./N_msgs;
RxError = zeros(length(SNR_dB), 1);
for i=1:length(SNR_dB)
    RxError(i) = mean(error(1,:,i));
end

figure;
semilogy(SNR_dB, RxError, '-o', 'LineWidth',3);
grid on;
title('Received Signal Errors');
xlabel('SNR(dB)');
ylabel('P(Error)');
grid on;




%% 1.7
bitPairs = reshape(vectors, 2, []);
constellationMap = bi2de(bitPairs', 'left-msb');
qpskSymbols = zeros(1, N_msgs);
for i=1:N_msgs * N_vecs/2
   qpskSymbols(i) = qpskQI(constellationMap(i) + 1);
end
transmittedSignals = reshape(qpskSymbols, [N_msgs/2, N_vecs]);

nfft = N_msgs/2;
OFDMSignals = ofdmmod(transmittedSignals,nfft,cplen);

SNR_dB2 = N_msgs/N_codded * SNR_dB;
TxSigPower = sum(abs(OFDMSignals).^2) / length(OFDMSignals);
noise_power = TxSigPower ./ (10.^(SNR_dB2./10));
noise = sqrt(noise_power*2) .* randn(size(OFDMSignals));
RxSignal = [];
for i=1:length(SNR_dB2)
    RxSignal = [RxSignal, OFDMSignals + noise(:,i)];
end

y_1 = ofdmdemod(RxSignal,nfft,cplen);
quadrature = real(y_1)>0;
inphase    = imag(y_1)>0;
positions  = 3 - 2*quadrature - inphase;
y_2 = zeros(N_msgs, N_vecs, length(SNR_dB2));
for i=1:N_vecs
    for j=1:length(SNR_dB2)
        temp = de2bi(positions(:,i,j), 'left-msb');
        y_2(:,i,j) = reshape(temp',[N_msgs,1]);
    end
end

error = sum(abs((y_2 - vectors)))./N_msgs;
RxError2 = zeros(length(SNR_dB2), 1);
for i=1:length(SNR_dB2)
    RxError2(i) = mean(error(1,:,i));
end
figure;
semilogy(SNR_dB2, RxError2, '-o', 'LineWidth',3);
grid on;
title('Received Signal Errors');
xlabel('SNR(dB)');
ylabel('P(Error)');
grid on;

%% Comparrision
figure;
semilogy(RxError, '-o', 'LineWidth',3);
hold on;
semilogy(RxError2, '-o', 'LineWidth',3);
hold off;
grid on;
title('Codding Effect Comparision');
xlabel('E_b / N_0');
ylabel('P(Error)');
legend('Pe Coding', 'Pe No-Coding');
grid on;