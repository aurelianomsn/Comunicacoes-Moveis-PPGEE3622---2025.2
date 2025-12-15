% Universidade de Brasília (UnB)
% Projeto 3 - Comunicações móveis
% Discente: Aureliano Magalhães de Sousa Neto
% Data: 03/12/2025

clear; clc;

Nbc       = 100;               % número de blocos de coerência por rede
Ncf       = 300;               % número total de redes avaliadas
fc        = 3 * 10^9;          % freq da portadora [Hz]
Bw        = 20 * 10^6;         % largura de banda [Hz]
Fn_dB     = 9;                 % figura de ruído [dB]
Fn_linear = 10^(Fn_dB/10);     % figura de ruído [linear]
h_ap      = 15;                % altura dos APs [m]
h_ue      = 1.65;              % altura dos UEs [m]
T0        = 296.15;            % temperatura local [K]
L_x       = 1000;              % dimensão X da área [m]
L_y       = 1000;              % dimensão Y da área [m]
Pp        = 0.2;               % potência das sequências piloto
Pdl       = 0.2;               % potência de transmissão em downlink
tau_p     = 50;                % comprimento das sequências piloto
k0        = 1.381 * 10^-23;    % constante de Boltzmann [J/K]
M_values  = [100, 150, 200];   % valores de APs possíveis
K_values  = [10, 20, 30];      % valores de UEs possíveis
pl_fs     = 32.44 + 20*log10(fc/10^6) + 20*log10(0.001); % atenuação no espaço livre para 1m
sigma_2_w = k0 * T0 * Bw * Fn_linear;                    % potência de ruído

% Estruturas para armazenar resultados
sinr_estatistica_M = cell(length(M_values), 1); % SINR com conhecimento estatístico para diferentes M. K=20 fixo.
sinr_instantanea_M = cell(length(M_values), 1); % SINR com conhecimento perfeito para diferentes M.
taxa_estatistica_M = cell(length(M_values), 1); % Taxas com conhecimento estatístico para diferentes M (em bits/s).
taxa_instantanea_M = cell(length(M_values), 1); % Taxas com conhecimento perfeito para diferentes M (em bits/s).

sinr_estatistica_K = cell(length(K_values), 1); % SINR com conhecimento estatístico para diferentes K. M=100 fixo.
sinr_instantanea_K = cell(length(K_values), 1); % SINR com conhecimento perfeito para diferentes K.
taxa_estatistica_K = cell(length(K_values), 1); % Taxas com conhecimento estatístico para diferentes K.
taxa_instantanea_K = cell(length(K_values), 1); % Taxas com conhecimento perfeito para diferentes K.

%% SIMULAÇÃO 1: Variação de APs (M) sendo mantida a mesma quantidade de UEs (k= 20)

for m_idx = 1:length(M_values)
    m = M_values(m_idx);
    k = 20;  % Número de usuários
    
    sinr_estatistica_completa = zeros(Ncf, k);
    sinr_instantanea_media    = zeros(Ncf, k);
    taxa_estatistica_completa = zeros(Ncf, k);
    taxa_instantanea_completa = zeros(Ncf, k);
    
   
    pos_ap                  = zeros(3, m);  % geração da matriz para posições dos APs
    pos_ue                  = zeros(3, k);  % geração da matriz para posições dos UEs
    dist_ap_ue              = zeros(m, k);  % geração da matriz para distâncias entre APs e UEs
    omega_dB                = zeros(m, k);  % geração da matriz para desvanecimento em larga escala   
    sinal_de_projecao       = zeros(m, k);
    coeficientes_de_canal_g = zeros(m, k);
    canal_estimado_gmk      = zeros(m, k);
    sinr_instantanea_rede   = zeros(Nbc, k);
    v_mk = zeros(m, k); % geração da matriz para ruído equivalente
    x_sf = zeros(m, k); % geração da matriz para V.A de média zero e desvio padrão 8 dB
    h_i  = zeros(m, k); % geração da matriz para parte real do desvanecimento em pequena escala   
    h_q  = zeros(m, k); % geração da matriz para parte complexa do desvanecimento em pequena escala   
    h    = zeros(m, k); % geração da matriz para desvanecimento em pequena escala   
    
 for rede = 1:Ncf 
        for i = 1:m
            x = (rand * L_x) - L_x/2;
            y = (rand * L_y) - L_y/2;
            pos_ap(:, i) = [x; y; h_ap]; % posições dos APs (já transposta)
        end
        
        for i = 1:k
            x = (rand * L_x) - L_x/2;
            y = (rand * L_y) - L_y/2;
            pos_ue(:, i) = [x; y; h_ue]; % posições dos UEs (já transposta)
        end
        
        for i = 1:m
            for j = 1:k
                dist_ap_ue(i,j) = norm(pos_ap(:, i) - pos_ue(:, j)); % cômputo das distâncias entre os APs e os UEs
                x_sf(i,j) = 8 * randn(1,1); % V.A de média zero e desvio padrão 8 dB
            end
        end

dist_ap_ue = max(dist_ap_ue, 1);
        
        for i = 1:m
            for j = 1:k
                omega_dB(i,j) = pl_fs + 28 * log10(dist_ap_ue(i,j)) + x_sf(i,j); % desvanecimento em larga escala
            end
        end

omega_linear = 10.^(-omega_dB/10);  % conversão do desvanecimento em larga escala para escala linear
        
        for bloco = 1:Nbc % Loop interno: blocos de coerência
            for i = 1:m
                for j = 1:k
                    h_i(i,j) = randn / sqrt(2);        % parte real do desvanecimento em pequena escala
                    h_q(i,j) = randn / sqrt(2);        % parte complexa do desvanecimento em pequena escala
                    h(i,j) = h_i(i,j) + 1j * h_q(i,j); % desvanecimento em pequena escala  
                    coeficientes_de_canal_g(i,j) = sqrt(omega_linear(i,j)) .* h(i,j);  % coeficientes de canal em todos os enlaces
                end
            end
            
v_mk = sqrt(sigma_2_w/2) * (randn(m, k) + 1j * randn(m, k)); % ruído equivalente: 

            for i = 1:m
                for j = 1:k
                    sinal_de_projecao(i,j) = sqrt(tau_p * Pp) * coeficientes_de_canal_g(i,j) + v_mk(i,j); 
                end
            end
            
c_mk = (sqrt(tau_p * Pp) .* omega_linear) ./ (tau_p * Pp .* omega_linear + sigma_2_w); 
canal_estimado_gmk = c_mk .* sinal_de_projecao;      
gamma_mk = sqrt(Pp * tau_p) .* omega_linear .* c_mk; 
eta_m = 1 ./ sum(gamma_mk, 2);
eta_mk = repmat(eta_m, 1, k);
            
            % Calcular SINR para conhecimento instantâneo
            for user = 1:k
                sinal_desejado = 0;
                for ap = 1:m
                    sinal_desejado = sinal_desejado + sqrt(eta_mk(ap, user)) * coeficientes_de_canal_g(ap, user) * conj(canal_estimado_gmk(ap, user));
                end
                sinal_desejado = Pdl * abs(sinal_desejado)^2;
                
                interferencia = 0;
                for outro_user = 1:k
                    if outro_user ~= user
                        soma_interf = 0;
                        for ap = 1:m
                            soma_interf = soma_interf + sqrt(eta_mk(ap, outro_user)) * coeficientes_de_canal_g(ap, user) * conj(canal_estimado_gmk(ap, outro_user));
                        end
                        interferencia = interferencia + Pdl * abs(soma_interf)^2;
                    end
                end
                
                sinr_instantanea_rede(bloco, user) = sinal_desejado / (interferencia + sigma_2_w);
            end
        end
        
        % Calcular SINR para conhecimento estatístico
        for user = 1:k
            numerador = 0;
            for ap = 1:m
                numerador = numerador + sqrt(eta_mk(ap, user)) * gamma_mk(ap, user);
            end
            numerador = Pdl * numerador^2;
            
            denominador = 0;
            for outro_user = 1:k
                for ap = 1:m
                    denominador = denominador + eta_mk(ap, outro_user) * gamma_mk(ap, outro_user) * omega_linear(ap, user);
                end
            end
            denominador = Pdl * denominador + sigma_2_w;
            
            sinr_estatistica_completa(rede, user) = numerador / denominador;
        end
        
        % Calcular médias
        sinr_instantanea_media(rede, :) = mean(sinr_instantanea_rede, 1);
        
        % Calcular taxas
        for user = 1:k
            taxa_estatistica_completa(rede, user) = log2(1 + sinr_estatistica_completa(rede, user));
            
            taxas_instantaneas = zeros(Nbc, 1);
            for bloco = 1:Nbc
                taxas_instantaneas(bloco) = log2(1 + sinr_instantanea_rede(bloco, user));
            end
            taxa_instantanea_completa(rede, user) = mean(taxas_instantaneas);
        end   
 end % Fim do Loop Monte Carlo
    
    % Converter e armazenar resultados
    sinr_estatistica_M{m_idx} = 10*log10(sinr_estatistica_completa(:));
    sinr_instantanea_M{m_idx} = 10*log10(sinr_instantanea_media(:));
    taxa_estatistica_M{m_idx} = taxa_estatistica_completa(:) * Bw;
    taxa_instantanea_M{m_idx} = taxa_instantanea_completa(:) * Bw;
    
end % Fim da simulação 1

%% SIMULAÇÃO 2: Variação de K (UEs) com M=100 fixo

for k_idx = 1:length(K_values)
    m = 100;  
    k = K_values(k_idx);
    
    sinr_estatistica_completa = zeros(Ncf, k);
    sinr_instantanea_media = zeros(Ncf, k);
    taxa_estatistica_completa = zeros(Ncf, k);
    taxa_instantanea_completa = zeros(Ncf, k);
    
    % Arrays temporários
    pos_ap = zeros(3, m);  
    pos_ue = zeros(3, k);  
    dist_ap_ue = zeros(m, k);
    omega_dB = zeros(m, k);
    sinal_de_projecao = zeros(m, k);
    v_mk = zeros(m, k);
    x_sf = zeros(m, k);
    h_i = zeros(m, k);
    h_q = zeros(m, k);
    h = zeros(m, k);
    coeficientes_de_canal_g = zeros(m, k);
    canal_estimado_gmk = zeros(m, k);
    
    % Gerar ruído para projeção dos pilotos
    v_mk = sqrt(sigma_2_w/2) * (randn(m, k) + 1j*randn(m, k));
    
    % Loop Monte Carlo
    for rede = 1:Ncf
        % Distribuição dos APs - já na forma 3 x M
        for i = 1:m
            x = (rand * L_x) - L_x/2;
            y = (rand * L_y) - L_y/2;
            pos_ap(:, i) = [x; y; h_ap];
        end
        
        % Distribuição dos UEs - já na forma 3 x K
        for i = 1:k
            x = (rand * L_x) - L_x/2;
            y = (rand * L_y) - L_y/2;
            pos_ue(:, i) = [x; y; h_ue];
        end
        
        % Distâncias entre os APs e os UEs
        for i = 1:m
            for j = 1:k
                dist_ap_ue(i,j) = norm(pos_ap(:, i) - pos_ue(:, j));
                x_sf(i,j) = 8 * randn(1,1);
            end
        end
        dist_ap_ue = max(dist_ap_ue, 1);
        
        % Desvanecimento em larga escala
        for i = 1:m
            for j = 1:k
                omega_dB(i,j) = pl_fs + 28*log10(dist_ap_ue(i,j)) + x_sf(i,j);
            end
        end
        omega_linear = 10.^(-omega_dB/10);
        
        sinr_instantanea_rede = zeros(Nbc, k);
        
        % Loop interno: blocos de coerência
        for bloco = 1:Nbc
            % Gerar coeficientes de pequena escala
            for i = 1:m
                for j = 1:k
                    h_i(i,j) = randn / sqrt(2);
                    h_q(i,j) = randn / sqrt(2);
                    h(i,j) = h_i(i,j) + 1j*h_q(i,j);
                    coeficientes_de_canal_g(i,j) = sqrt(omega_linear(i,j)) .* h(i,j);
                end
            end
            
            % Estimar canal
            for i = 1:m
                for j = 1:k
                    sinal_de_projecao(i,j) = sqrt(tau_p * Pp) * coeficientes_de_canal_g(i,j) + v_mk(i,j);
                end
            end
            
            c_mk = (sqrt(tau_p * Pp) .* omega_linear) ./ (tau_p * Pp .* omega_linear + sigma_2_w);
            canal_estimado_gmk = c_mk .* sinal_de_projecao;
            gamma_mk = sqrt(Pp * tau_p) .* omega_linear .* c_mk;
            eta_m = 1 ./ sum(gamma_mk, 2);
            eta_mk = repmat(eta_m, 1, k);
            
            % Calcular SINR para conhecimento instantâneo
            for user = 1:k
                sinal_desejado = 0;
                for ap = 1:m
                    sinal_desejado = sinal_desejado + sqrt(eta_mk(ap, user)) * ...
                        coeficientes_de_canal_g(ap, user) * conj(canal_estimado_gmk(ap, user));
                end
                sinal_desejado = Pdl * abs(sinal_desejado)^2;
                
                interferencia = 0;
                for outro_user = 1:k
                    if outro_user ~= user
                        soma_interf = 0;
                        for ap = 1:m
                            soma_interf = soma_interf + sqrt(eta_mk(ap, outro_user)) * ...
                                coeficientes_de_canal_g(ap, user) * conj(canal_estimado_gmk(ap, outro_user));
                        end
                        interferencia = interferencia + Pdl * abs(soma_interf)^2;
                    end
                end
                
                sinr_instantanea_rede(bloco, user) = sinal_desejado / (interferencia + sigma_2_w);
            end
        end
        
        % Calcular SINR para conhecimento estatístico
        for user = 1:k
            numerador = 0;
            for ap = 1:m
                numerador = numerador + sqrt(eta_mk(ap, user)) * gamma_mk(ap, user);
            end
            numerador = Pdl * numerador^2;
            
            denominador = 0;
            for outro_user = 1:k
                for ap = 1:m
                    denominador = denominador + eta_mk(ap, outro_user) * ...
                        gamma_mk(ap, outro_user) * omega_linear(ap, user);
                end
            end
            denominador = Pdl * denominador + sigma_2_w;
            
            sinr_estatistica_completa(rede, user) = numerador / denominador;
        end
        
        % Calcular médias
        sinr_instantanea_media(rede, :) = mean(sinr_instantanea_rede, 1);
        
        % Calcular taxas
        for user = 1:k
            taxa_estatistica_completa(rede, user) = log2(1 + sinr_estatistica_completa(rede, user));
            
            taxas_instantaneas = zeros(Nbc, 1);
            for bloco = 1:Nbc
                taxas_instantaneas(bloco) = log2(1 + sinr_instantanea_rede(bloco, user));
            end
            taxa_instantanea_completa(rede, user) = mean(taxas_instantaneas);
        end
    end
    
    % Converter e armazenar resultados
    sinr_estatistica_K{k_idx} = 10*log10(sinr_estatistica_completa(:));
    sinr_instantanea_K{k_idx} = 10*log10(sinr_instantanea_media(:));
    taxa_estatistica_K{k_idx} = taxa_estatistica_completa(:) * Bw;
    taxa_instantanea_K{k_idx} = taxa_instantanea_completa(:) * Bw;
end

%% Resultados das simulações

% Gráfico 1: ECDF da SINR para diferentes valores de M (K=20 fixo)
figure(1);
set(gcf, 'Position', [100, 100, 800, 600]);
hold on;

cores_est = {'b-', 'g-', 'r-'};
cores_inst = {'b--', 'g--', 'r--'};
legendas = {};

for m_idx = 1:length(M_values)
    m = M_values(m_idx);
    
    % ECDF para conhecimento estatístico (ECSI)
    [F_est, x_est] = ecdf(sinr_estatistica_M{m_idx});
    plot(x_est, F_est, cores_est{m_idx}, 'LineWidth', 1.5);
    legendas{end+1} = sprintf('ECSI - M = %d', m);
    
    % ECDF para conhecimento perfeito (PCSI)
    [F_inst, x_inst] = ecdf(sinr_instantanea_M{m_idx});
    plot(x_inst, F_inst, cores_inst{m_idx}, 'LineWidth', 1.5);
    legendas{end+1} = sprintf('PCSI - M = %d', m);
end

hold off;
grid on;
xlabel('SINR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('ECDF', 'FontSize', 12, 'FontWeight', 'bold');
title('ECDF da SINR para Diferentes Números de APs (K=20 fixo)', 'FontSize', 14, 'FontWeight', 'bold');
legend(legendas, 'Location', 'southeast', 'FontSize', 10);
xlim([-10, 35]);

% Gráfico 2: ECDF da Taxa Alcançável para diferentes valores de M (K=20 fixo)
figure(2);
set(gcf, 'Position', [100, 100, 800, 600]);
hold on;

legendas = {};

for m_idx = 1:length(M_values)
    m = M_values(m_idx);
    
    % ECDF para conhecimento estatístico (ECSI) em Mbit/s
    [F_est, x_est] = ecdf(taxa_estatistica_M{m_idx} / 1e6);
    plot(x_est, F_est, cores_est{m_idx}, 'LineWidth', 1.5);
    legendas{end+1} = sprintf('ECSI - M = %d', m);
    
    % ECDF para conhecimento perfeito (PCSI) em Mbit/s
    [F_inst, x_inst] = ecdf(taxa_instantanea_M{m_idx} / 1e6);
    plot(x_inst, F_inst, cores_inst{m_idx}, 'LineWidth', 1.5);
    legendas{end+1} = sprintf('PCSI - M = %d', m);
end

hold off;
grid on;
xlabel('Taxa Alcançável (Mbit/s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('ECDF', 'FontSize', 12, 'FontWeight', 'bold');
title('ECDF da Taxa Alcançável para Diferentes Números de APs (K=20 fixo)',  'FontSize', 14, 'FontWeight', 'bold');
legend(legendas, 'Location', 'southeast', 'FontSize', 10);

% Gráfico 3: ECDF da SINR para diferentes valores de K (M=100 fixo)
figure(3);
set(gcf, 'Position', [100, 100, 800, 600]);
hold on;

cores_est_K = {'b-', 'g-', 'r-'};
cores_inst_K = {'b--', 'g--', 'r--'};
legendas = {};

for k_idx = 1:length(K_values)
    k = K_values(k_idx);
    
    % ECDF para conhecimento estatístico (ECSI)
    [F_est, x_est] = ecdf(sinr_estatistica_K{k_idx});
    plot(x_est, F_est, cores_est_K{k_idx}, 'LineWidth', 1.5);
    legendas{end+1} = sprintf('ECSI - K = %d', k);
    
    % ECDF para conhecimento perfeito (PCSI)
    [F_inst, x_inst] = ecdf(sinr_instantanea_K{k_idx});
    plot(x_inst, F_inst, cores_inst_K{k_idx}, 'LineWidth', 1.5);
    legendas{end+1} = sprintf('PCSI - K = %d', k);
end

hold off;
grid on;
xlabel('SINR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('ECDF', 'FontSize', 12, 'FontWeight', 'bold');
title('ECDF da SINR para Diferentes Números de UEs (M=100 fixo)', 'FontSize', 14, 'FontWeight', 'bold');
legend(legendas, 'Location', 'southeast', 'FontSize', 10);
xlim([-10, 35]);

% Gráfico 4: ECDF da Taxa Alcançável para diferentes valores de K (M=100 fixo)
figure(4);
set(gcf, 'Position', [100, 100, 800, 600]);
hold on;

legendas = {};

for k_idx = 1:length(K_values)
    k = K_values(k_idx);
    
    % ECDF para conhecimento estatístico (ECSI) em Mbit/s
    [F_est, x_est] = ecdf(taxa_estatistica_K{k_idx} / 1e6);
    plot(x_est, F_est, cores_est_K{k_idx}, 'LineWidth', 1.5);
    legendas{end+1} = sprintf('ECSI - K = %d', k);
    
    % ECDF para conhecimento perfeito (PCSI) em Mbit/s
    [F_inst, x_inst] = ecdf(taxa_instantanea_K{k_idx} / 1e6);
    plot(x_inst, F_inst, cores_inst_K{k_idx}, 'LineWidth', 1.5);
    legendas{end+1} = sprintf('PCSI - K = %d', k);
end

hold off;
grid on;
xlabel('Taxa Alcançável (Mbit/s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('ECDF', 'FontSize', 12, 'FontWeight', 'bold');
title('ECDF da Taxa Alcançável para Diferentes Números de UEs (M=100 fixo)', 'FontSize', 14, 'FontWeight', 'bold');
legend(legendas, 'Location', 'southeast', 'FontSize', 10);

%% Salvar todos os resultados
save('resultados_completos_projeto.mat', ...
    'M_values', 'K_values', ...
    'sinr_estatistica_M', 'sinr_instantanea_M', ...
    'taxa_estatistica_M', 'taxa_instantanea_M', ...
    'sinr_estatistica_K', 'sinr_instantanea_K', ...
    'taxa_estatistica_K', 'taxa_instantanea_K');