% Universidade de Brasília (UnB)
% Projeto 2 - Comunicações móveis
% Discente: Aureliano Magalhães de Sousa Neto
% Data: 17/11/2025

% -------------------------------------------------------------------------
% Parte 1: Simulação da Probabilidade de Interrupção - Envoltória Rayleigh
% -------------------------------------------------------------------------

% Parâmetros de simulação
amostras = 1e8;
gamma_th_dB = -30:1:30;
gamma_th = 10.^(gamma_th_dB/10);

gamma_bar_dB = [-20, 0, 20];
gamma_bar = 10.^(gamma_bar_dB/10);

Pout_simulado = zeros(length(gamma_bar), length(gamma_th));
Pout_analitico = zeros(length(gamma_bar), length(gamma_th));

for i = 1:length(gamma_bar)

    sigma = 1/sqrt(2);

    z1 = gaussian_boxmuller(amostras); % Geração de gaussianas via Box-Muller
    z2 = gaussian_boxmuller(amostras); % Geração de gaussianas via Box-Muller
    x = sigma * z1;                    % Geração de VA 1
    y = sigma * z2;                    % Geração de VA 2
    beta = sqrt(x.^2 + y.^2);          % Envoltória Rayleigh

    gamma_s = gamma_bar(i) .* (beta.^2);
    
    for j = 1:length(gamma_th)
        Pout_simulado(i,j) = mean(gamma_s < gamma_th(j));
        Pout_analitico(i,j) = 1 - exp(-gamma_th(j)/gamma_bar(i));
    end
end

% -------------------------------
% Plot
% -------------------------------
figure;
set(gcf, 'Color', 'w');
colors = ['b', 'r', 'g'];

for i = 1:length(gamma_bar)
    semilogy(gamma_th_dB, Pout_analitico(i,:), 'Color', colors(i), ...
        'LineWidth', 2, 'DisplayName', sprintf('Analítico - \\gamma_s = %.0f dB', gamma_bar_dB(i)));
    hold on;

    semilogy(gamma_th_dB, Pout_simulado(i,:), 'o', 'Color', colors(i), ...
        'MarkerSize', 4, 'LineWidth', 1, ...
        'DisplayName', sprintf('Simulado - \\gamma_s = %.0f dB', gamma_bar_dB(i)));
end

grid on;
xlabel('Limiar de SNR (\gamma_{th}) [dB]', 'FontSize', 12);
ylabel('Probabilidade de Outage (P_{out})', 'FontSize', 12);
title('Probabilidade de Outage - Canal Rayleigh', 'FontSize', 14);
legend('Location', 'southeast', 'FontSize', 10);
ylim([1e-5 1]);
xlim([-30 30]);
set(gca, 'FontSize', 11);





% -----------------------------------------------------------
% Parte 2: Simulação da Probabilidade de Interrupção - Envoltória Rice
% -----------------------------------------------------------


clearvars -except gaussian_boxmuller marcum_q_custom lcg_random;
clc;

% -------------------------------
% Parâmetros
% -------------------------------
amostras = 1e8;
gamma_th_dB = -30:1:30;
gamma_th = 10.^(gamma_th_dB/10);

gamma_bar_dB = [-20, 0, 20];
gamma_bar = 10.^(gamma_bar_dB/10);

K_R = [0.1, 1, 10];  % NÃO ALTERAR (a pedido do usuário)

% -------------------------------
% Loop sobre KR
% -------------------------------
for indice_Rice = 1:length(K_R)
    KR_value = K_R(indice_Rice);

    Pout_simulado = zeros(length(gamma_bar), length(gamma_th));
    Pout_analitico = zeros(length(gamma_bar), length(gamma_th));

    % Parâmetros da distribuição Rice (potência unitária)
    A = sqrt(KR_value / (KR_value + 1));  
    sigma = 1/sqrt(2*(KR_value + 1));

    % -------------------------------
    % Geração Gaussianas via Box-Muller
    % -------------------------------
    z1 = gaussian_boxmuller(amostras);
    z2 = gaussian_boxmuller(amostras);

    % Componentes do canal
    x = A + sigma * z1;
    y = sigma * z2;

    % Envoltória Rice
    envoltoria = sqrt(x.^2 + y.^2);

    % -------------------------------
    % Calcular para cada SNR média
    % -------------------------------
    for i = 1:length(gamma_bar)

        SNR_inst = gamma_bar(i) .* (envoltoria.^2);

        for j = 1:length(gamma_th)

            % Simulação
            Pout_simulado(i,j) = mean(SNR_inst < gamma_th(j));

            % Parâmetros da Marcum-Q
            a_param = sqrt(2*KR_value);
            b_param = sqrt(2*(KR_value+1)*gamma_th(j)/gamma_bar(i));

            Pout_analitico(i,j) = 1 - marcum_q_custom(a_param, b_param);
        end
    end

    % -------------------------------
    % Plot para este KR
    % -------------------------------
    figure;
    set(gcf, 'Color', 'w');
    cores = ['b', 'r', 'g'];
    estilos = {'-', '--', ':'};

    for i = 1:length(gamma_bar)
        semilogy(gamma_th_dB, Pout_analitico(i,:), ...
            'Color', cores(i), 'LineStyle', estilos{i}, ...
            'LineWidth', 2, ...
            'DisplayName', sprintf('Analítico - \\gamma_s = %.0f dB', gamma_bar_dB(i)));
        hold on;

        semilogy(gamma_th_dB, Pout_simulado(i,:), ...
            'o', 'Color', cores(i), 'MarkerSize', 5, ...
            'LineWidth', 1.5, ...
            'DisplayName', sprintf('Simulado - \\gamma_s = %.0f dB', gamma_bar_dB(i)));
    end

    grid on;
    xlabel('Limiar de SNR (\gamma_{th}) [dB]');
    ylabel('Probabilidade de Interrupção (P_{out})');
    title(sprintf('Probabilidade de Interrupção - Rice (K_R = %.1f)', KR_value));
    legend('Location', 'southeast');
    ylim([1e-5 1]);
    xlim([-30 30]);
    set(gca, 'FontSize', 11);

end
