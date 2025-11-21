% Universidade de Brasília (UnB)
% Projeto 2 - Comunicações móveis
% Discente: Aureliano Magalhães de Sousa Neto
% Data: 17/11/2025

% -------------------------------------------------------------------------
% Parte 1: Simulação da Probabilidade de Interrupção - Envoltória Rayleigh
% -------------------------------------------------------------------------

% Parâmetros de simulação
amostras_comparacao = 1e8;             % número de amostras (Monte Carlo)
gamma_th_dB = -30:1:30;                % limiar de SNR em dB (eixo X)
gamma_th = 10.^(gamma_th_dB/10);       % conversão p/ escala linear
gamma_bar_dB = [-20, 0, 20];           % valores médios de SNR (em dB)
gamma_bar = 10.^(gamma_bar_dB/10);     % valores médios de SNR em escala linear
Pout_simulado = zeros(length(gamma_bar), length(gamma_th));
Pout_analitico = zeros(length(gamma_bar), length(gamma_th));

for i = 1:length(gamma_bar)
    x1 = randn(1, amostras_comparacao) * sqrt(0.5); % Primeira VA Gaussiana ~N(0,1/√2)
    x2 = randn(1, amostras_comparacao) * sqrt(0.5); % Segunda VA Gaussiana ~N(0,1/√2)
    beta = sqrt(x1.^2 + x2.^2);          % Envoltória Rayleigh: β = √(x1² + x2²)
    gamma_s = gamma_bar(i) * (beta.^2);  % SNR instantânea: γ = γ̄ * β²

    for j = 1:length(gamma_th)
        Pout_simulado(i,j) = mean(gamma_s < gamma_th(j));         % Estimativa por simulação
        Pout_analitico(i,j) = 1 - exp(-gamma_th(j)/gamma_bar(i)); % Expressão analítica
    end
end

% -----------------------------------------------------------
% Plot dos resultados solicitados na etapa 1 do roteiro
% -----------------------------------------------------------

figure;
set(gcf, 'Color', 'w');
colors = ['b', 'r', 'g'];

for i = 1:length(gamma_bar)

% Curva analítica (linha contínua)
    semilogy(gamma_th_dB, Pout_analitico(i,:), 'Color', colors(i), 'LineWidth', 2, 'DisplayName', sprintf('Analítico - \\gamma_s = %.0f dB', gamma_bar_dB(i)));
    hold on;
    
% Curva simulada (marcadores)
    semilogy(gamma_th_dB, Pout_simulado(i,:), 'o', 'Color', colors(i), 'MarkerSize', 4, 'LineWidth', 1, 'DisplayName', sprintf('Simulado - \\gamma_s = %.0f dB', gamma_bar_dB(i)));
end

grid on;
xlabel('Limiar de SNR (\gamma_{th}) [dB]', 'FontSize', 12);
ylabel('Probabilidade de Outage (P_{out})', 'FontSize', 12);
title('Probabilidade de Outage - Canal Rayleigh', 'FontSize', 14);
legend('Location', 'southeast', 'FontSize', 10);
ylim([1e-5 1]);
xlim([-30 30]);
set(gca, 'FontSize', 11);

% ----------------------------------------------------------------------
% Plot comparativo para diferentes quantidades de amostras (Monte Carlo)
% ----------------------------------------------------------------------

amostras_comparacao = [1e6, 1e7, 1e8];   % Diferentes tamanhos de amostra
gamma_bar_green_index = 3;     % Índice da curva verde (γ̄ = 20 dB)
gamma_bar_green_dB = gamma_bar_dB(gamma_bar_green_index);
gamma_bar_green = gamma_bar(gamma_bar_green_index);
Pout_simulado_comparacao = zeros(length(amostras_comparacao), length(gamma_th));

for n = 1:length(amostras_comparacao)
    x1 = randn(1, amostras_comparacao(n)) * sqrt(0.5);
    x2 = randn(1, amostras_comparacao(n)) * sqrt(0.5);
    beta = sqrt(x1.^2 + x2.^2);  
    gamma_s = gamma_bar_green * (beta.^2);  % SNR instantânea para γ̄ = 20 dB

% Cálculo da probabilidade de outage para cada limiar
    for j = 1:length(gamma_th)
        Pout_simulado_comparacao(n,j) = mean(gamma_s < gamma_th(j));
    end
end

Pout_analitico_green = 1 - exp(-gamma_th/gamma_bar_green);
figure;
set(gcf, 'Color', 'w');

semilogy(gamma_th_dB, Pout_analitico_green, 'k-', 'LineWidth', 3, 'DisplayName', sprintf('Analítico - \\gamma_s = %.0f dB', gamma_bar_green_dB));
hold on;

colors_nt = ['r', 'b', 'g'];  % Cores diferentes para cada Nt
markers = ['o', 's', '^'];    % Marcadores diferentes para cada Nt

for n = 1:length(amostras_comparacao)
    semilogy(gamma_th_dB, Pout_simulado_comparacao(n,:), 'Color', colors_nt(n), 'Marker', markers(n), 'MarkerSize', 5, 'LineWidth', 1.5, 'MarkerFaceColor', colors_nt(n), 'DisplayName', sprintf('Simulado - Nt = %.0e', amostras_comparacao(n)));
end

grid on;
xlabel('Limiar de SNR (\gamma_{th}) [dB]', 'FontSize', 12);
ylabel('Probabilidade de Outage (P_{out})', 'FontSize', 12);
title(sprintf('Comparação de P_{out} para Diferentes Nt (\\gamma_s = %.0f dB)', gamma_bar_green_dB), 'FontSize', 14);
legend('Location', 'southeast', 'FontSize', 10);
ylim([1e-5 1]);
xlim([-30 30]);
set(gca, 'FontSize', 11);


