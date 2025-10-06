% Universidade de Brasília (UnB)
% Projeto 1 - Comunicações móveis
% Discente: Aureliano Magalhães de Sousa Neto
% Data: 06/10/2025

%% Atraso multipercurso - espalhamento de atraso eficaz
M  = 1;                           % Quantidade de amostras
fc = 3;                           % [GHz]

% Cálculo das estatísticas - (passo 2 - 10/55)
media_uma_nlos  = -0.204*log10(1+fc) - 6.28; % media do espalhamento [log]
desvio_uma_nlos = 0.39;                      % desvio padrão do espalhamento [log]

% Geração de uma amostra aleatória - (passo 3 - 10/55)
espalhamento_multip_log = normrnd (media_uma_nlos, desvio_uma_nlos, [M,1]); % espalhamento de atraso eficaz

% Conversão para linear - (passo 4 - 10/55)
espalhamento_multip_linear = 10.^espalhamento_multip_log;                         % espalhamento de atraso eficaz [seg]

% Dados para o gráfico 0.5 GHz <= fc <= 100 GHz
fc_q1c         = 0.5:1:100;                      % [GHz] para traçar gráfico 0.5 GHz <= fc <= 100 GHz
tamanho_desvio_umx_q1c = length(fc_q1c);         % tamanho do vetor  
desvio_umx_q1c = desvio_uma_nlos * ones(1,tamanho_desvio_umx_q1c);   
media_umx_q1c  = -0.204*log10(1+fc_q1c) - 6.28;  % media do espalhamento [log]
media_umx_q1c  = 10.^media_umx_q1c;              % segundos
media_umx_q1c_us = media_umx_q1c * 10^6;         % [µs]

% Média e desvio padrão do espalhamento de atraso vs Frequência (escala linear): Gráfico parte 1, letra C

figure;
subplot(1,2,1)
plot(fc_q1c, media_umx_q1c_us, '-', 'LineWidth', 1.5, 'DisplayName', 'Média');
hold on;
plot(fc_q1c, desvio_umx_q1c, '--', 'LineWidth', 1.5, 'DisplayName', 'Desvio Padrão');
xlabel('Frequência (GHz)');
ylabel('Espalhamento de atraso eficaz (\mus)');
title('Média e desvio padrão do espalhamento de atraso vs Frequência (escala linear)');
grid on;
ytickformat('%.2f');   % mostra com duas casas decimais
legend('Location','best');

% Média do Espalhamento de atraso vs Frequência (escala logarítmica): Gráfico parte 1, letra C
subplot(1,2,2)
semilogx(fc_q1c, media_umx_q1c_us, '-', 'LineWidth', 1.5,  'DisplayName', 'Média');  % eixo x logarítmico
hold on;
semilogx(fc_q1c, desvio_umx_q1c, '--', 'LineWidth', 1.5, 'DisplayName', 'Desvio Padrão');  % eixo x logarítmico
xlabel('Frequência (GHz)');
ylabel('Média do Espalhamento de atraso eficaz (\mus)');
title('Média do Espalhamento de atraso vs Frequência (escala logarítmica)');
grid on;
ytickformat('%.2f');   % mostra com duas casas decimais
legend('Location','best');


%% Atraso multipercurso - componentes de atraso multipercurso
N = 100;                          % Número de multipercursos
fator_proporcionalidade_rt = 2.3; % Fator de proporcionalidade (slide 12/55)

% Cálculo média ut (mu) - (passo 2 - 15/55)
ut = fator_proporcionalidade_rt * espalhamento_multip_linear;

% Geração de atrasos exponenciais -  (passo 3 - 15/55)
componente_multipercurso_exponencial = exprnd (ut, [N,1]);

% Normalização do atraso - (passo 4 - 15/55)
componente_multipercurso_normalizado = componente_multipercurso_exponencial - min(componente_multipercurso_exponencial);

% Ordenação - (passo 5 - 15/55)
componente_multipercurso_ordenado = sort(componente_multipercurso_normalizado);


%% Parte: Potência multipercurso - perfil de potência multipercurso 

% Geração dos termos de sombreamento - passo 1 (19/55)
desvio_padrao_sombreamento = 6; % dB
sombreamento = normrnd(0, desvio_padrao_sombreamento, [N,1]);

% Cálculo da potência preliminar - passo 2 (19/55)
pot_preliminar  = exp(-componente_multipercurso_ordenado*(fator_proporcionalidade_rt - 1)/ut).*10.^(-sombreamento/10);

% Fator de Rice - passo 3 (19/55)
Kr = 0; % NLoS 

% Normalização das potências dispersas
ganho_canal = 0;

for x = 2:N
    ganho_canal = ganho_canal + pot_preliminar(x);
end

pot_dispersa_normalizada = (1/(Kr+1))*(pot_preliminar/ganho_canal);

figure;
h = stem(componente_multipercurso_ordenado*10^6, pot_dispersa_normalizada, 'k', 'LineWidth', 1.5);
set(h, 'Marker', '^');

set(gca, 'YScale', 'log'); % eixo y em semilog
xlabel('Domínio de Atraso [\mus]');
ylabel('PDP');
title('Perfil de atraso de potência considerando N = 100 componentes multipercurso');
grid on;

% Recálculo do espalhamento de atraso do canal de acordo com a sua definição
somatorio_1 = 0;
somatorio_2 = 0;

for x = 1:N
    somatorio_1 = somatorio_1 + pot_preliminar(x)*componente_multipercurso_ordenado(x);
end

atraso_medio = (1/ganho_canal)*somatorio_1;

for x = 1:N
    somatorio_2 = somatorio_2 + pot_preliminar(x)*(componente_multipercurso_ordenado(x)-atraso_medio)^2;
end

espalhamento_atraso_recalculado = sqrt ((1/ganho_canal)*somatorio_2);


%% Parte: Direções de chegada

% Ângulo em azimute
% Gerar uma amostra do espalhamento angular azimutal - (passo 1 27/55)
media_espalhamento_ang_azimutal           = -0.27*log10(fc)+2.08; %[unid logaritmica]
desvio_padrao_espalhamento_ang_azimutal   = 0.11;                 %[unid logaritmica]
amostra_espalhamento_angular_azimutal_log = normrnd (media_espalhamento_ang_azimutal, desvio_padrao_espalhamento_ang_azimutal, [M,1]);

% Converter a amostra para escala linear e radianos - (passo 2 27/55)
amostra_espalhamento_angular_azimutal_dg = 10.^amostra_espalhamento_angular_azimutal_log;
amostra_espalhamento_angular_azimutal_rad = amostra_espalhamento_angular_azimutal_dg*(pi/180);

% Gerar ângulos iniciais - (passo 3 27/55)
angulos_azimute_iniciais = 1.42.*amostra_espalhamento_angular_azimutal_rad.*(-log((pot_dispersa_normalizada)./(max(pot_dispersa_normalizada)))).^(1/2);

% Gere sinais aleatórios - (passo 4 27/55)
u_n_azimute = randsample([-1,1], N, true);
u_n_azimute = u_n_azimute';

% Geração de flutuações aleatórias - (passo 5 27/55) 
y_n_azimute = normrnd (0,amostra_espalhamento_angular_azimutal_rad/7, [N,1]);

% Cálculo dos ângulos finais - (passo 6 27/55) teta_n equação 7
angulos_azimute_finais = u_n_azimute .* angulos_azimute_iniciais + y_n_azimute; %[radianos]

% Converter radianos para graus
angulos_azimute_final_deg = angulos_azimute_finais * 180/pi;

angulos_azimute_final_deg = mod(rad2deg(angulos_azimute_finais), 360);
figure;

 subplot(1,2,1)
 h=stem(angulos_azimute_final_deg, pot_dispersa_normalizada, 'k', 'filled', 'MarkerSize', 3);
 set(h, 'Marker', '^');
 xlim([0, 360]);
 xticks(0:15:360);               
 xlabel('Ângulo de chegada em azimute (graus)', 'FontSize', 12);
 ylabel('Potência Dispersa Normalizada', 'FontSize', 12);
 title('Espectro de potência angular (ângulo azimutal)', 'FontSize', 14);
 grid on
 set(gca, 'YScale', 'log'); % Ajustar eixo y para log
 
% Criar gráfico polar

subplot(1,2,2)
polarplot(deg2rad(angulos_azimute_final_deg), pot_dispersa_normalizada, 'b.', 'MarkerSize', 15)

title('Ângulo de chegada em azimute - Visualização Polar')
rlim([0 max(pot_dispersa_normalizada)])

hold on
 for i = 1:length(angulos_azimute_final_deg)
     polarplot([0 deg2rad(angulos_azimute_final_deg(i))], [0 pot_dispersa_normalizada(i)], 'b-', 'LineWidth', 1)
 end


% Ângulo em elevação

% Geração do espalhamento angular em elevação - (passo 1 32/55)
media_espalhamento_ang_elevacao           = -0.3236*log10(fc)+1.512;
desvio_padrao_espalhamento_ang_elevacao   = 0.16;
amostra_espalhamento_angular_elevacao_log = normrnd (media_espalhamento_ang_elevacao, desvio_padrao_espalhamento_ang_elevacao, [M,1]);

% Conversão para unidaades lineares - (passo 2 32/55)
amostra_espalhamento_angular_elevacao_dg = 10.^amostra_espalhamento_angular_elevacao_log;
amostra_espalhamento_angular_elevacao_rad = amostra_espalhamento_angular_elevacao_dg*(pi/180);

% Geração dos ângulos de chegada inicial - (passo 3 32/55)
angulos_elevacao_iniciais = -amostra_espalhamento_angular_elevacao_rad.*log((pot_dispersa_normalizada)./(max(pot_dispersa_normalizada)));

% Geração de amostras aleatórias no conjunto discreto [-1, 1] - (passo 4 32/55)
u_n_elevacao = randsample([-1,1], N, true);
u_n_elevacao = u_n_elevacao'; % Comentado às 15h19 de 28/09/2025

% Geração de flutuações aleatórias - (passo 5 33/55) 
y_n_elevacao = normrnd (0,amostra_espalhamento_angular_elevacao_rad/7, [N,1]);

% Escolha arbitrariamente um angulo entre 0 e pi/2 para representar angulo de elevação média - (passo 6 33/55)
ang_elevacao_medio = 45; %escolhido 45º

% Cálculo dos ângulos finais - (passo 7 33/55) teta_n equação 7
angulos_elevacao_finais = u_n_elevacao .* angulos_elevacao_iniciais + y_n_elevacao;

% Converter radianos para graus
angulos_elevacao_final_deg = angulos_elevacao_finais * 180/pi;

angulos_elevacao_final_deg = mod(rad2deg(angulos_elevacao_final_deg), 360);

subplot (1,2,1)
h = stem(angulos_elevacao_final_deg, pot_dispersa_normalizada, 'k', 'filled', 'MarkerSize', 3);
set(h, 'Marker', '^');
xlim([0, 360]);
xticks(0:15:360);                
xlabel('Ângulo de chegada em elevação (graus)', 'FontSize', 12)
ylabel('Potência Dispersa Normalizada', 'FontSize', 12)
title('Espectro de potência angular (ângulo de elevação)', 'FontSize', 14)
grid on
set(gca, 'YScale', 'log');

% Criar gráfico polar
subplot(1,2,2)
polarplot(deg2rad(angulos_elevacao_final_deg), pot_dispersa_normalizada, 'b.', 'MarkerSize', 15)
title('Ângulo de chegada em elevação - Visualização Polar')
rlim([0 max(pot_dispersa_normalizada)])

% Adicionar linhas stem (da origem até cada ponto)
hold on
for i = 1:length(angulos_elevacao_final_deg)
    polarplot([0 deg2rad(angulos_elevacao_final_deg(i))], [0 pot_dispersa_normalizada(i)], 'b-', 'LineWidth', 1)
end


% Geração dos vetores de direção de chegada (slide 35/55)

r_n = [cosd(angulos_azimute_final_deg).*sind(angulos_elevacao_final_deg); sind(angulos_azimute_final_deg).*sind(angulos_elevacao_final_deg);cosd(angulos_elevacao_final_deg)];

% Número de vetores (cada 3 elementos formam um vetor 3D)
N = length(r_n)/3;

% Reorganizar em matriz 3xN (cada coluna é um vetor)
r_mat = reshape(r_n, [3, N]);

% Plotar vetores 3D saindo da origem
figure;
quiver3(zeros(1,N), zeros(1,N), zeros(1,N), ...
        r_mat(1,:), r_mat(2,:), r_mat(3,:), 0.5, 'b');

xlabel('x'); ylabel('y'); zlabel('z');
title('Vetores direção de chegada das componentes multipercurso');
axis equal;
grid on;


%%  Parte: Desvio Doppler

% Fixe uma velocidade escalar para móvel - (passo 1 38/55)
vrx_1 = 5;             % [m/s]
vrx_2 = 50;            % [m/s]

% Calcule o comprimento de onda  - (passo 2 38/55)
lambda_fc = 3*10^8/(fc*10^9); % [m]

% Fixe um angulo azimutal (entre 0 e 2pi) e angulo de elevacao (entre 0 e pi)  - (passo 3 38/55)
ang_azim_doppler = pi/2;
ang_elev_doppler = pi/2;

% Cálculo do vetor velocidade do receptor  - (passo 4 38/55)
vetor_vel_rx_1 = vrx_1*[cos(ang_azim_doppler)*sin(ang_elev_doppler); sin(ang_azim_doppler)*sin(ang_elev_doppler);cos(ang_elev_doppler)];   %38/55 (passo 4) - 3x1
vetor_vel_rx_2 = vrx_2*[cos(ang_azim_doppler)*sin(ang_elev_doppler); sin(ang_azim_doppler)*sin(ang_elev_doppler);cos(ang_elev_doppler)];   %38/55 (passo 4) - 3x1

qtd_blocos = length(r_n)/3;
r_blocos = reshape(r_n, 3, qtd_blocos);   % 3 x 100

vel = vetor_vel_rx_1(:);  % garante 3x1
desvio_doppler_1 = (1/lambda_fc) * (vel' * r_blocos)';  % 100x1

figure;
subplot(1,2,1)
h = stem(desvio_doppler_1, pot_dispersa_normalizada, 'k', 'filled', 'MarkerSize', 3);
set(h, 'Marker', '^');
xlabel('Desvio Doppler (Hz)', 'FontSize', 12);
ylabel('Potência Dispersa Normalizada', 'FontSize', 12);
title('Dispersão da potência no domínio de desvio Doppler (Vrx = 5m/s)', 'FontSize', 14);
grid on;
set(gca, 'YScale', 'log'); % eixo y em semilog

vel = vetor_vel_rx_2(:);  % garante 3x1
desvio_doppler_2 = (1/lambda_fc) * (vel' * r_blocos)';  % 100x1

subplot(1,2,2)
h = stem(desvio_doppler_2, pot_dispersa_normalizada, 'k', 'filled', 'MarkerSize', 3);
set(h, 'Marker', '^');
xlabel('Desvio Doppler (Hz)', 'FontSize', 12);
ylabel('Potência Dispersa Normalizada', 'FontSize', 12);
title('Dispersão da potência no domínio de desvio Doppler (Vrx = 50m/s)', 'FontSize', 14);
grid on;
set(gca, 'YScale', 'log'); % eixo y em semilog


%%  Parte: Determinação do sinal recebido considerando a caracterização multipercurso gerada nos passos anteriores

% Gere três pulsos diferentescom δt ∈ {10^-7, 10^-5, 10^−3} s

deltas = [1e-7, 1e-5, 1e-3]; % três larguras de pulso

fc = fc*10^9; %[GHz]

% Gera e armazena os sinais transmitidos (letra A da parte 5 do projeto)
for k = 1:length(deltas)
    delta_t = deltas(k);
    [t, s_tx] = gerar_pulso(0, delta_t, N);  % pulso transmitido
    t = t(:).';                               % garante vetor linha
    sinal_recebido = zeros(1, N);

    dt = t(2) - t(1);

    for n = 1:length(componente_multipercurso_ordenado)
        % elementos escalares
        tau_n = componente_multipercurso_ordenado(n);  % escalar
        alpha_n = sqrt(pot_dispersa_normalizada(n));   % escalar
        nu_n = desvio_doppler_1(n);                    % escalar

        % termo de fase
        fase = exp(-1j*2*pi*((fc + nu_n)*tau_n - nu_n*t));  % 1×Nt

        %  pulso atrasado
        pulso_atrasado = min(round(tau_n / dt), N);  % nunca maior que Nt
        sinal_atrasado = circshift(s_tx, [0 pulso_atrasado]);
        if pulso_atrasado > 0
            sinal_atrasado(1:pulso_atrasado) = 0;
        end

    % Determinação do sinal recebido (letra B da parte 5 do projeto)
        sinal_recebido = sinal_recebido + alpha_n * fase .* sinal_atrasado;  
    end

    % Plot comparativo: sinal tx x sinal rx (letra C da parte 5 do projeto)
    figure;
    plot(t/delta_t, abs(s_tx), 'b', 'LineWidth', 1.5); hold on;
    plot(t/delta_t, abs(sinal_recebido), 'r--', 'LineWidth', 1.5);
    xlabel('Tempo / \delta t');
    ylabel('|r (t)|');
    legend('Sinal Transmitido', 'Sinal Recebido');
    title(['NLoS - \delta t = ', num2str(delta_t), ' s, \sigma_{\tau} = ', num2str(espalhamento_multip_linear*1e9), ' ns']);
    grid on;
end

%% Banda de coerência e Tempo de coerência

%% === AUTOCORRELAÇÃO DO CANAL NORMALIZADA (Cenário NLoS) ===
ganho_canal=sum(pot_preliminar);

%% --- Parâmetros de varredura (somente positivos) ---
kappa_MHz = linspace(0, 2, 2001);   % varre até 2 MHz
kappa = kappa_MHz * 1e6;            % converter para Hz

sigma_ms = linspace(0, 10, 2001);   % até 10 ms
sigma = sigma_ms * 1e-3;            % converter para segundos

%% --- Função anônima para ρ_TT ---
rho_TT = @(k,s,pot,tau,nu,Omega) ...
    (1/Omega) * sum(pot(:) .* exp(-1j*2*pi*tau(:)*k + 1j*2*pi*nu(:)*s));

%% === (b) Autocorrelação em frequência ρ_TT(κ; 0) ===
rho_k = arrayfun(@(k) rho_TT(k, 0, pot_preliminar, ...
    componente_multipercurso_ordenado, desvio_doppler_1, ganho_canal), kappa);
rho_k_mag = abs(rho_k);

figure;
plot(kappa_MHz, rho_k_mag, 'b', 'LineWidth', 1.5); hold on;
xlabel('\kappa [MHz]');
ylabel('| \rho_{TT}(\kappa, 0) |');
title('Autocorrelação em frequência ρ_{TT}(\kappa,0)');
grid on;

% --- Linhas tracejadas azuis com anotação ---
thresholds = [0.95, 0.90];
banda_coerencia_MHz = zeros(size(thresholds));

for ii = 1:length(thresholds)
    th = thresholds(ii);
    cross_idx = find(rho_k_mag < th, 1, 'first');
    if isempty(cross_idx)
        Tc_kHz = kappa_MHz(end);
    else
        Tc_kHz = kappa_MHz(cross_idx);
    end
    banda_coerencia_MHz(ii) = Tc_kHz;

    % Linha vertical tracejada azul
    xline(Tc_kHz, '--', 'Color', [0 0.45 0.9], ...
        'Label', sprintf('\\rho_B=%.2f → %.3f MHz', th, Tc_kHz), ...
        'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'aligned');
end

% 🔹 Zoom visual para 0–0.5 MHz
xlim([0 0.5]);
ylim([0 1.05]);
set(gca, 'XMinorGrid', 'on');


%% === (c) Autocorrelação temporal ρ_TT(0; σ) ===
v_list = [5, 50];
nu_list = {desvio_doppler_1, desvio_doppler_2};
tempo_coerencia_ms = zeros(length(v_list), length(thresholds));

figure;
for v_idx = 1:length(v_list)
    vrx = v_list(v_idx);
    nu = nu_list{v_idx};
    
    rho_s = arrayfun(@(s) rho_TT(0, s, pot_preliminar, ...
        componente_multipercurso_ordenado, nu, ganho_canal), sigma);
    rho_s_mag = abs(rho_s);
    
    subplot(length(v_list), 1, v_idx);
    plot(sigma_ms, rho_s_mag, 'LineWidth', 1.5);
    xlabel('\sigma [ms]');
    ylabel('| \rho_{TT}(0,\sigma) |');
    title(sprintf('Autocorrelação temporal ρ_{TT}(0,σ) (v_{rx} = %d m/s)', vrx));
    grid on; hold on;

    for tt = 1:length(thresholds)
        th = thresholds(tt);
        cross_idx = find(rho_s_mag < th, 1, 'first');
        if isempty(cross_idx)
            Tc_ms = sigma_ms(end);
        else
            Tc_ms = sigma_ms(cross_idx);
        end
        tempo_coerencia_ms(v_idx, tt) = Tc_ms;

        xline(Tc_ms, '--', ...
            sprintf('\\rho_T=%.2f → %.3f ms', th, Tc_ms), ...
            'Color', [0 0.6 0.9], 'LabelVerticalAlignment', 'bottom', ...
            'LabelOrientation', 'aligned');
    end

    % 🔹 Zoom específico por velocidade
    if vrx == 50
        xlim([0 0.5]);   % foco no intervalo 0–0.5 ms
    else
        xlim([0 4]);     % foco mais amplo (0–4 ms)
    end
    ylim([0 1.05]);
    set(gca, 'XMinorGrid', 'on');
end



