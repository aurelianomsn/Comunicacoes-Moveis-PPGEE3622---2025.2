% Universidade de Brasília (UnB)
% Projeto 2 - Comunicações móveis
% Discente: Aureliano Magalhães de Sousa Neto
% Data: 17/11/2025

% -------------------------------------------------------------------------
% Parte 1: Simulação da Probabilidade de Interrupção - Envoltória Rayleigh
% -------------------------------------------------------------------------

amostras = 1e6;                     % número de amostras (Monte Carlo)
gamma_th_dB = -30:1:30;             % limiar de SNR em dB (eixo X)
gamma_th = 10.^(gamma_th_dB/10);    % conversão p/ escala linear
gamma_bar_dB = [-20, 0, 20];        % valores médios de SNR (em dB)
gamma_bar = 10.^(gamma_bar_dB/10);  % valores médios de SNR em escala linear

Pout_simulado = zeros(length(gamma_bar), length(gamma_th));
Pout_analitico = zeros(length(gamma_bar), length(gamma_th));

for i = 1:length(gamma_bar)

    x1 = randn(1, amostras) * sqrt(0.5); % Primeira VA Gaussiana ~N(0,1/2)
    x2 = randn(1, amostras) * sqrt(0.5); % Segunda VA Gaussiana ~N(0,1/2)
    beta = sqrt(x1.^2 + x2.^2);          % Envoltória Rayleigh
    gamma_s = gamma_bar(i) * (beta.^2);  % SNR instantânea: γ = γ̄ * β²

    for j = 1:length(gamma_th)
        Pout_simulado(i,j) = mean(gamma_s < gamma_th(j));         % Simulação
        Pout_analitico(i,j) = 1 - exp(-gamma_th(j)/gamma_bar(i)); % Analítico
    end
end

% Plot exigido pelo roteiro

figure;
set(gcf, 'Color', 'w');
colors = ['b', 'r', 'g'];

for i = 1:length(gamma_bar)
% Curva analítica
    semilogy(gamma_th_dB, Pout_analitico(i,:), 'Color', colors(i), 'LineWidth', 2, 'DisplayName', sprintf('Analítico - \\gamma_s = %.0f dB', gamma_bar_dB(i)));
    hold on;

% Curva simulada
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

% Plot comparativo para diferentes quantidades de amostras (Monte Carlo)

amostras_comparativas = [1e6, 1e7, 1e8];          % Diferentes tamanhos de amostra
gamma_bar_green_index = 3;           % Índice da curva verde (gamma_bar = 20 dB)
gamma_bar_green_dB = gamma_bar_dB(gamma_bar_green_index);
gamma_bar_green = gamma_bar(gamma_bar_green_index);

Pout_simulado_comparacao = zeros(length(amostras_comparativas), length(gamma_th));

for n = 1:length(amostras_comparativas)

    x1 = randn(1, amostras_comparativas(n)) * sqrt(0.5);
    x2 = randn(1, amostras_comparativas(n)) * sqrt(0.5);

    beta = sqrt(x1.^2 + x2.^2);
    gamma_s = gamma_bar_green * (beta.^2);  % SNR instantânea para gamma_bar=20 dB

    for j = 1:length(gamma_th)
        Pout_simulado_comparacao(n,j) = mean(gamma_s < gamma_th(j));
    end
end

Pout_analitico_green = 1 - exp(-gamma_th/gamma_bar_green);

figure;
set(gcf, 'Color', 'w');

semilogy(gamma_th_dB, Pout_analitico_green, 'k-', 'LineWidth', 3, 'DisplayName', sprintf('Analítico - \\gamma_s = %.0f dB', gamma_bar_green_dB));
hold on;

colors_nt = ['r', 'b', 'g'];    % Cores diferentes
markers = ['o', 's', '^'];      % Marcadores diferentes

for n = 1:length(amostras_comparativas)
    semilogy(gamma_th_dB, Pout_simulado_comparacao(n,:), 'Color', colors_nt(n), 'Marker', markers(n), 'MarkerSize', 5, 'LineWidth', 1.5, 'MarkerFaceColor', colors_nt(n), 'DisplayName', sprintf('Simulado - Nt = %.0e', amostras_comparativas(n)));
end

grid on;
xlabel('Limiar de SNR (\gamma_{th}) [dB]', 'FontSize', 12);
ylabel('Probabilidade de Outage (P_{out})', 'FontSize', 12);
title(sprintf('Comparação de P_{out} para Diferentes Nt (\\gamma_s = %.0f dB)', gamma_bar_green_dB), 'FontSize', 14);

legend('Location', 'southeast', 'FontSize', 10);
ylim([1e-5 1]);
xlim([-30 30]);
set(gca, 'FontSize', 11);


% --------------------------------------------------------------------
% Parte 2: Simulação da Probabilidade de Interrupção - Envoltória Rice
% --------------------------------------------------------------------

K_R = [0.1, 1, 10]; % fatores de Rice

for indice_fator_K = 1:length(K_R)

    K = K_R(indice_fator_K);

    % Parâmetros da distribuição Rice (potência unitária)
    A = sqrt(K / (K + 1));
    sigma = 1 / sqrt(2 * (K + 1));

    Pout_sim = zeros(length(gamma_bar), length(gamma_th)); % Probabilidade de interrupção (simulado)
    Pout_ana = zeros(length(gamma_bar), length(gamma_th)); % Probabilidade de interrupção (analítico)

    for i = 1:length(gamma_bar)

        x = A + sigma * randn(1, amostras);
        y = sigma * randn(1, amostras);
        rice_env = sqrt(x.^2 + y.^2); % Geração da envoltória Rice

        % SNR instantânea
        gamma_inst = gamma_bar(i) * (rice_env.^2);

          for j = 1:length(gamma_th)
            Pout_sim(i,j) = mean(gamma_inst < gamma_th(j)); 

            a = sqrt(2*K);
            b = sqrt(2*(K+1)*gamma_th(j)/gamma_bar(i));
            Pout_ana(i,j) = 1 - marcum_q_custom(a, b);      
        end
    end

% Plot exigido pelo roteiro
    figure;
    set(gcf, 'Color', 'w');
    cores = ['b', 'r', 'g'];

    for i = 1:length(gamma_bar)
        semilogy(gamma_th_dB, Pout_ana(i,:), 'Color', cores(i), 'LineWidth', 2, 'DisplayName', sprintf('Analítico  -  \\gamma_s = %.0f dB', gamma_bar_dB(i)));
        hold on;

        semilogy(gamma_th_dB, Pout_sim(i,:), 'o', 'Color', cores(i), 'MarkerSize', 5, 'LineWidth', 1.2, 'DisplayName', sprintf('Simulado  -  \\gamma_s = %.0f dB', gamma_bar_dB(i)));
    end

    grid on;
    xlabel('Limiar de SNR (\gamma_{th}) [dB]', 'FontSize', 12);
    ylabel('P_{out}', 'FontSize', 12);
    title(sprintf('Outage - Envoltória Rice (K_R = %.2f)', K), 'FontSize', 14);

    legend('Location', 'southeast', 'FontSize', 10);
    ylim([1e-5 1]);
    xlim([-30 30]);
end



% --------------------------------------------------------------------
% Parte 3: Probabilidade de erro (constelação M-QAM)
% --------------------------------------------------------------------

Valores_M = [4, 16, 64]; % Ordens de modulação
Cores = ['b', 'r', 'g'];
Marcadores = ['o', 's', '^'];

figure('Position', [100, 100, 1200, 800]);
set(gcf, 'Color', 'w');

for indice_M = 1:length(Valores_M)
    M = Valores_M(indice_M);
    bits_por_simbolo = log2(M);  % bits por símbolo
    
% Geração da constelação M-QAM
    raiz_M = sqrt(M);
    constelacao = modulacao_MQAM_manual(0:M-1, M, 'UnitAveragePower', true);

% Arrays para resultados
    Pe_simulado_rayleigh = zeros(size(gamma_th_dB));
    Pe_teorico_rayleigh = zeros(size(gamma_th_dB));
    Pe_teorico_awgn = zeros(size(gamma_th_dB));
    
    for indice_SNR = 1:length(gamma_th_dB)
        SNR_media = gamma_th(indice_SNR);
        
% Simulação do canal Rayleigh
        dados = randi([0 M-1], 1, amostras); % Geração de símbolos aleatórios
        simbolos_transmitidos = modulacao_MQAM_manual(dados, M, 'UnitAveragePower', true);
        canal = (randn(1, amostras) + 1j*randn(1, amostras)) / sqrt(2); % Geração do canal Rayleigh (potência unitária)        
        potencia_ruido = 1 / SNR_media;  % Potência do sinal é 1 (UnitAveragePower)
        ruido = sqrt(potencia_ruido/2) * (randn(1, amostras) + 1j*randn(1, amostras)); % Ruído AWGN        
        simbolos_recebidos_rayleigh = canal .* simbolos_transmitidos + ruido; % Sinal recebido (canal Rayleigh)        
        simbolos_equalizados = simbolos_recebidos_rayleigh ./ canal; % Equalização ZF
        
% Demodulação
        dados_recebidos_rayleigh = demodulacao_MQAM_manual(simbolos_equalizados, M, 'UnitAveragePower', true);      
        erros_rayleigh = sum(dados ~= dados_recebidos_rayleigh);
        Pe_simulado_rayleigh(indice_SNR) = erros_rayleigh / amostras; % Cálculo da probabilidade de erro
        
% Probabilidade de erro teórica para Rayleigh (fórmula fornecida no roteiro)
        Termo_C = sqrt(1.5 * SNR_media / (M - 1 + 1.5 * SNR_media));
        termo1 = 2 * (raiz_M-1)/raiz_M * (1 - Termo_C);
        termo2 = ((raiz_M-1)/raiz_M)^2 * (1 - (4/pi) * Termo_C * atan(1/Termo_C));
        Pe_teorico_rayleigh(indice_SNR) = termo1 - termo2;
        Pe_teorico_awgn(indice_SNR) = 4 * (1 - 1/raiz_M) * qfunc(sqrt(3 * SNR_media / (M - 1))); % Probabilidade de erro teórica para AWGN
    end
    
    subplot(2,1,1);
    semilogy(gamma_th_dB, Pe_teorico_rayleigh, '--', 'Color', Cores(indice_M), 'LineWidth', 2, 'DisplayName', sprintf('Teórico Rayleigh M=%d', M));
    hold on;
    semilogy(gamma_th_dB, Pe_simulado_rayleigh, Marcadores(indice_M), 'Color', Cores(indice_M), 'MarkerSize', 6, 'LineWidth', 1, 'DisplayName', sprintf('Simulado Rayleigh M=%d', M));
    
    subplot(2,1,2);
    semilogy(gamma_th_dB, Pe_teorico_awgn, ':', 'Color', Cores(indice_M), 'LineWidth', 2, 'DisplayName', sprintf('AWGN M=%d', M));
    hold on;
end

subplot(2,1,1);
grid on;
xlabel('SNR Média [dB]', 'FontSize', 12);
ylabel('Probabilidade de Erro de Símbolo (P_e)', 'FontSize', 12);
title('Desempenho M-QAM em Canal Rayleigh', 'FontSize', 14);
legend('Location', 'southwest', 'FontSize', 10);
ylim([1e-4, 1]);
xlim([-30, 30]);

subplot(2,1,2);
grid on;
xlabel('SNR Média [dB]', 'FontSize', 12);
ylabel('Probabilidade de Erro de Símbolo (P_e)', 'FontSize', 12);
title('Comparação: Desempenho M-QAM em Canal AWGN', 'FontSize', 14);
legend('Location', 'southwest', 'FontSize', 10);
ylim([1e-4, 1]);
xlim([-30, 30]);

% Gráfico combinado para comparação direta
figure('Position', [100, 100, 1000, 600]);
set(gcf, 'Color', 'w');

for indice_M = 1:length(Valores_M)
    M = Valores_M(indice_M);
    
    Pe_teorico_rayleigh_combinado = zeros(size(gamma_th_dB));
    Pe_teorico_awgn_combinado = zeros(size(gamma_th_dB));
    
    for indice_SNR = 1:length(gamma_th_dB)
        SNR_media = gamma_th(indice_SNR);
        
% Cálculo teórico Rayleigh
        Termo_C = sqrt(1.5 * SNR_media / (M - 1 + 1.5 * SNR_media));
        termo1 = 2 * (sqrt(M)-1)/sqrt(M) * (1 - Termo_C);
        termo2 = ((sqrt(M)-1)/sqrt(M))^2 * (1 - (4/pi) * Termo_C * atan(1/Termo_C));
        Pe_teorico_rayleigh_combinado(indice_SNR) = termo1 - termo2;
        
% Cálculo teórico AWGN
        Pe_teorico_awgn_combinado(indice_SNR) = 4 * (1 - 1/sqrt(M)) * qfunc(sqrt(3 * SNR_media / (M - 1)));
    end
    
% Plot Rayleigh
    semilogy(gamma_th_dB, Pe_teorico_rayleigh_combinado, '-', 'Color', Cores(indice_M), 'LineWidth', 2, 'DisplayName', sprintf('Rayleigh M=%d', M));
    hold on;
    
% Plot AWGN
    semilogy(gamma_th_dB, Pe_teorico_awgn_combinado, '--', 'Color', Cores(indice_M), 'LineWidth', 2, 'DisplayName', sprintf('AWGN M=%d', M));
end

grid on;
xlabel('SNR Média [dB]', 'FontSize', 12);
ylabel('Probabilidade de Erro de Símbolo (P_e)', 'FontSize', 12);
title('Comparação Rayleigh vs AWGN para M-QAM', 'FontSize', 14);
legend('Location', 'southwest', 'FontSize', 10);
ylim([1e-4, 1]);
xlim([-30, 30]);