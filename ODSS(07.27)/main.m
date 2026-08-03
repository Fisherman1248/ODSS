clc;
clear;
close all;

%% ODSS parameter settings
q_list = [1.5 2.0];

numCases = length(q_list);

BER_results = cell(1, numCases);
SNR_results = cell(1, numCases);
folder_results = cell(1, numCases);

%% Run ODSS transmitter and receiver
for ii = 1:numCases

    q = q_list(ii);

    fprintf('\n');
    fprintf('============================================\n');
    fprintf('Start ODSS AWGN simulation: q = %.3g\n', q);
    fprintf('============================================\n');

    %% Generate transmitter data
    fprintf('\n=== Generate transmitter data: q = %.3g ===\n', q);
    folderName = TX(q);
    folder_results{ii} = folderName;

    %% Run receiver
    fprintf('\n=== Run receiver: q = %.3g ===\n', q);
    [BER_LS, SNR_dB] = RX(folderName);
    BER_results{ii} = BER_LS;
    SNR_results{ii} = SNR_dB;

    fprintf('Completed simulation for q = %.3g\n', q);
end

%% Compare BER curves for different q values
figure;
hold on;

for ii = 1:numCases

    semilogy( ...
        SNR_results{ii}, ...
        BER_results{ii}, ...
        'o-', ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('q = %.2f', q_list(ii)));
end

grid on;
ylim([1e-6 1]);

xlabel('Time-domain waveform SNR (dB)');
ylabel('Bit Error Rate');

title('ODSS BER over pure AWGN');
legend('Location', 'southwest');

%% Save summary results
save('ODSS_AWGN_All_Results.mat', ...
    'q_list', ...
    'BER_results', ...
    'SNR_results', ...
    'folder_results');

fprintf('\n');
disp('ALL ODSS AWGN DATA GENERATED!');