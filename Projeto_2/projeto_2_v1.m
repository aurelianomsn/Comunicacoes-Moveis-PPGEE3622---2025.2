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
    x1 = randn(1, amostras_comparacao) * sqrt(0.5); % Primeira VA Gaussiana ~N(0,1/sqrt(2))
    x2 = randn(1, amostras_comparacao) * sqrt(0.5); % Segunda VA Gaussiana ~N(0,1/sqrt(2))
    beta = sqrt(x1.^2 + x2.^2);          % Envoltória Rayleigh: Beta = sqrt(x1² + x2²)
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

% ----------------------------------------------------------------------
% Plot comparativo para diferentes quantidades de amostras (Monte Carlo)
% ----------------------------------------------------------------------

amostras_comparacao = [1e6, 1e7, 1e8];   % Diferentes tamanhos de amostra
gamma_bar_green_index = 3;     % Índice da curva verde (gamma_bar = 20 dB)
gamma_bar_green_dB = gamma_bar_dB(gamma_bar_green_index);
gamma_bar_green = gamma_bar(gamma_bar_green_index);
Pout_simulado_comparacao = zeros(length(amostras_comparacao), length(gamma_th));

for n = 1:length(amostras_comparacao)
    x1 = randn(1, amostras_comparacao(n)) * sqrt(0.5);
    x2 = randn(1, amostras_comparacao(n)) * sqrt(0.5);
    beta = sqrt(x1.^2 + x2.^2);  
    gamma_s = gamma_bar_green * (beta.^2);  % SNR instantânea para gamma_bar= 20 dB

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

colors_nt = ['r', 'b', 'g'];  % Cores diferentes para cada quantidade de amostras diferente
markers = ['o', 's', '^'];    % Marcadores diferentes para cada quantidade de amostras diferente

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

% -----------------------------------------------------------
% Parte 2: Simulação da Probabilidade de Interrupção - Envoltória Rice
% -----------------------------------------------------------

% Parâmetros de simulação
amostras = 1e6;                 % número de amostras (Monte Carlo)
Limiar_SNR_dB = -30:1:30;              % limiar de SNR em dB (eixo X)
Limiar_SNR = 10.^(Limiar_SNR_dB/10);   % conversão p/ escala linear
SNR_Media_dB = [-20, 0, 20];           % Valores médios de SNR (em dB)
SNR_Media = 10.^(SNR_Media_dB/10);     % Valores médios de SNR em escala linear
K_R = [0.1, 1, 10];                    % Fatores de Rice

% Loop sobre cada fator de Rice
for indice_Rice = 1:length(K_R)
    Fator_Rice = K_R(indice_Rice);
    
    % Inicialização das variáveis de saída
    Prob_Interrupcao_Simulada = zeros(length(SNR_Media), length(Limiar_SNR));
    Prob_Interrupcao_Analitica = zeros(length(SNR_Media), length(Limiar_SNR));
    
    % Loop sobre cada valor médio de SNR
    for i = 1:length(SNR_Media)
        % Geração MANUAL da envoltória Rice conforme teoria
        % Para distribuição Rice: h = x + jy, onde:
        % x ~ N(μ_I, σ²), y ~ N(μ_Q, σ²)
        % com μ_I² + μ_Q² = A² (componente LOS)
        
        % Parâmetros da distribuição Rice (potência unitária)
        % E[β²] = A² + 2σ² = 1
        A = sqrt(Fator_Rice / (Fator_Rice + 1));  % Componente LOS
        sigma = 1/sqrt(2*(Fator_Rice + 1));       % Componente disperso
        
        % Geração de variáveis Gaussianas usando Box-Muller
        u1 = rand(1, amostras);
        u2 = rand(1, amostras);
        u3 = rand(1, amostras);
        u4 = rand(1, amostras);
        
        % Transformação Box-Muller
        z1 = sqrt(-2 * log(u1)) .* cos(2 * pi * u2);
        z2 = sqrt(-2 * log(u3)) .* cos(2 * pi * u4);
        
        % Componentes em fase e quadratura
        % Assumindo componente LOS apenas na fase (μ_I = A, μ_Q = 0)
        x = A + sigma * z1;  % Componente em fase com LOS
        y = sigma * z2;      % Componente em quadratura
        
        % Envoltória Rice
        envoltoria = sqrt(x.^2 + y.^2);
        
        % Verificação da potência média
        potencia_media = mean(envoltoria.^2);
        
        % SNR instantânea
        SNR_Instantanea = SNR_Media(i) * (envoltoria.^2);
        
        % Cálculo da probabilidade de interrupção
        for j = 1:length(Limiar_SNR)
            % Estimativa por simulação (Monte Carlo)
            Prob_Interrupcao_Simulada(i,j) = mean(SNR_Instantanea < Limiar_SNR(j));
            
            % Expressão analítica CORRETA da probabilidade de interrupção Rice
            % P_out(γ_th) = 1 - Q1(√(2K_R), √(2(K_R+1)γ_th/γ̄))
            a_param = sqrt(2 * Fator_Rice);
            b_param = sqrt(2 * (Fator_Rice + 1) * Limiar_SNR(j) / SNR_Media(i));
            
            % Usar função Marcum Q
            Prob_Interrupcao_Analitica(i,j) = 1 - marcum_q_custom(a_param, b_param);
        end
    end
    
    % -----------------------------------------------------------
    % Plot dos resultados para cada Fator_Rice
    % -----------------------------------------------------------
    figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [100, 100, 800, 600]);
    
    cores = ['b', 'r', 'g'];
    estilos = {'-', '--', ':'};
    
    for i = 1:length(SNR_Media)
        % Curva analítica (linha contínua)
        semilogy(Limiar_SNR_dB, Prob_Interrupcao_Analitica(i,:), ...
                 'Color', cores(i), 'LineStyle', estilos{i}, 'LineWidth', 2, ...
                 'DisplayName', sprintf('Analítico - \\gamma_s = %.0f dB', SNR_Media_dB(i)));
        hold on;
        
        % Curva simulada (marcadores)
        semilogy(Limiar_SNR_dB, Prob_Interrupcao_Simulada(i,:), ...
                 'o', 'Color', cores(i), 'MarkerSize', 5, 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('Simulado - \\gamma_s = %.0f dB', SNR_Media_dB(i)));
    end
    
    grid on;
    xlabel('Limiar de SNR (\gamma_{th}) [dB]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Probabilidade de Interrupção (P_{out})', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Probabilidade de Interrupção - Envoltória Rice (K_R = %.1f)', Fator_Rice), ...
          'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'southeast', 'FontSize', 10);
    ylim([1e-5 1]);
    xlim([-30 30]);
    
    % Melhorar aparência
    set(gca, 'FontSize', 11, 'GridAlpha', 0.3);
    set(gca, 'YMinorGrid', 'on', 'XMinorGrid', 'on');
end

% -----------------------------------------------------------
% Plot adicional: Comparação entre diferentes K_R para γ̄ = 0 dB
% -----------------------------------------------------------
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [100, 100, 800, 600]);

SNR_Media_index = 2; % γ̄ = 0 dB
cores_kr = ['r', 'g', 'b'];
nomes_kr = {'K_R = 0.1', 'K_R = 1', 'K_R = 10'};

for indice_Rice = 1:length(K_R)
    Fator_Rice = K_R(indice_Rice);
    
    % Recalcular curvas analíticas para plot comparativo
    Pout_kr = zeros(1, length(Limiar_SNR));
    for j = 1:length(Limiar_SNR)
        a_param = sqrt(2 * Fator_Rice);
        b_param = sqrt(2 * (Fator_Rice + 1) * Limiar_SNR(j) / SNR_Media(SNR_Media_index));
        Pout_kr(j) = 1 - marcum_q_custom(a_param, b_param);
    end
    
    semilogy(Limiar_SNR_dB, Pout_kr, ...
             'Color', cores_kr(indice_Rice), 'LineWidth', 2, ...
             'DisplayName', nomes_kr{indice_Rice});
    hold on;
end

grid on;
xlabel('Limiar de SNR (\gamma_{th}) [dB]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Probabilidade de Interrupção (P_{out})', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Comparação: P_{{out}} para Diferentes K_R (\\gamma_s = %.0f dB)', SNR_Media_dB(SNR_Media_index)), ...
      'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 10);
ylim([1e-5 1]);
xlim([-30 30]);

set(gca, 'FontSize', 11, 'GridAlpha', 0.3);
set(gca, 'YMinorGrid', 'on', 'XMinorGrid', 'on');

