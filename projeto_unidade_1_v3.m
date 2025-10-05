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
h = stem(angulos_elevacao_final_deg, pot_dispersa_normalizada, 'k', 'filled', 'MarkerSize', 3)
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
title('Vetores 3D a partir dos ângulos');
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

% Gere três pulsos diferentes com δt ∈ {10^-7, 10^-5, 10^-3} s
deltas = [1e-7, 1e-5, 1e-3]; % três larguras de pulso
Nt = 10^5; % amostras

for k = 1:length(deltas)
    delta_t = deltas(k);

    % Inicializa sinal recebido
    sinal_recebido = zeros(1, Nt);

    % Para simplicidade, defino o sinal transmitido como o pulso gerado sem Doppler
    [t, s_tx] = gerar_pulso(0, delta_t, Nt);   % atraso = 0 → pulso transmitido
    
    % Loop de multipercurso
    for n = 1:length (componente_multipercurso_ordenado)
        % gera pulso com atraso específico
        [t, s] = gerar_pulso(componente_multipercurso_ordenado(n), delta_t, Nt);

        % contribuição desse caminho
        componente = sqrt(pot_dispersa_normalizada(n)) .* ...
             exp(-1j*2*pi*((fc + desvio_doppler_1(n))*componente_multipercurso_ordenado(n) ...
                            - desvio_doppler_1(n)*t)) .* s;


        % soma no sinal recebido total
        sinal_recebido = sinal_recebido + componente;
    end

    % Plota sinal transmitido vs recebido
    figure;
    plot(t/delta_t, real(s_tx), 'b', 'LineWidth', 1.5); hold on;
    plot(t/delta_t, real(sinal_recebido), 'r--', 'LineWidth', 1.5);
    xlabel('Tempo / \delta t');
    ylabel('Amplitude');
    legend('Sinal Transmitido', 'Sinal Recebido');
    title(['Comparação - \delta t = ', num2str(delta_t)]);
    grid on;
end





%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

deltas = [1e-12, 1e-12, 1e-12]; % três larguras de pulso
deltas = [1e-5, 1e-5, 1e-5]; % três larguras de pulso
deltas = [1e-3, 1e-5, 1e-7]; % três larguras de pulso
Nt = 10^5; % amostras

for k = 1:length(deltas)
    delta_t = deltas(k);
    
    % SINAL TRANSMITIDO: usa a função gerar_pulso com atraso = 0
    [t, s_tx] = gerar_pulso(0, delta_t, Nt);
    
    % Inicializa sinal recebido
    sinal_recebido = zeros(1, Nt);

    % Loop de multipercurso (TODOS os caminhos são dispersos em NLoS)
    for n = 1:length(componente_multipercurso_ordenado)
        % Calcula fase estática UMA VEZ
        fase_estatica = 2*pi*(fc + desvio_doppler_1(n)) * componente_multipercurso_ordenado(n);
        
        % Gera pulso com atraso para este caminho usando a função
        [~, s_rx_n] = gerar_pulso(componente_multipercurso_ordenado(n), delta_t, Nt);
        
        % Todos os caminhos usam pot_dispersa_normalizada (sem componente LoS especial)
        amplitude = sqrt(pot_dispersa_normalizada(n));
        
        % Contribuição deste caminho
        componente = amplitude .* exp(-1j * (fase_estatica - 2*pi*desvio_doppler_1(n)*t)) .* s_rx_n;
        
        % Soma no sinal recebido total
        sinal_recebido = sinal_recebido + componente;
    end

    % Plota sinal transmitido vs recebido
    figure;
    plot(t/delta_t, real(s_tx), 'b', 'LineWidth', 1.5); hold on;
    plot(t/delta_t, real(sinal_recebido), 'r--', 'LineWidth', 1.5);
    xlabel('Tempo / \delta t');
    ylabel('Amplitude');
    legend('Sinal Transmitido', 'Sinal Recebido');
    title(['NLoS - \delta t = ', num2str(delta_t), ' s, \sigma_{\tau} = ', num2str(espalhamento_multip_linear*1e9), ' ns']);
    grid on;
end


















%***************************************************





% 
% deltas = [1e-7, 1e-5, 1e-3]; % três larguras de pulso
%  for k = 1:length(deltas)
%     delta_t = deltas(k);
%  end
% 
% % Inicializa sinal recebido
% sinal_recebido = zeros(1, N);
% 
% % Para cada caminho multipercurso
% for n = 1:N
%     % gera pulso com atraso específico
%     [t, s] = gerar_pulso(componente_multipercurso_ordenado(n), delta_t, N);
% 
%     % contribuição desse caminho
%     componente = sqrt(pot_dispersa_normalizada(n)) .* ...
%                  exp(-1j*2*pi*(fc + desvio_doppler_1(n)*componente_multipercurso_ordenado(n) - desvio_doppler_1(n)*t)) .* s;
% 
%     % soma no sinal recebido total
%     sinal_recebido = sinal_recebido + componente;
% end
% 
% % plotar parte real do sinal
% plot(t/delta_t, real(sinal_recebido));
% xlabel('Tempo / \delta t');
% ylabel('Re\{r(t)\}');
% grid on;
% 
% 
% %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ - comentado em 03/10
% % figure;
% % 
% % for k = 1:length(deltas) %garante a execução da função para cada delta_t
% %     delta_t = deltas(k);
% % 
% %     % Gera o pulso transmitido (atraso = 0)
% %     [t, s] = gerar_pulso(espalhamento_multip_linear, delta_t, N);
% % 
% %     % Subplot para cada pulso
% %     subplot(3,1,k);
% %     plot(t/delta_t, s, 'LineWidth', 1.5);
% %     grid on;
% % 
% %     % Eixos e título
% %     xlabel('Tempo / \delta t');
% %     ylabel('s(t)');
% %     title(['Pulso retangular com \delta t = ' num2str(delta_t) ' s']);
% %     xticks(0:1:5);
% %     xticklabels({'0','\delta_t','2\delta_t','3\delta_t','4\delta_t','5\delta_t'});
% % end
% 
% 
% %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% 
% 
% 
% %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ - fazendo por conta
% % figure;
% % 
% % for k = 1:length(deltas) %garante a execução da função para cada delta_t
% %     delta_t = deltas(k);
% % 
% %     % Gera o pulso transmitido (atraso = 0)
% %     [t, s] = gerar_pulso(espalhamento_multip_linear, delta_t, N);
% % 
% %     % Subplot para cada pulso
% %     subplot(3,1,k);
% %     plot(t/delta_t, s, 'LineWidth', 1.5);
% %     grid on;
% % 
% %     % Eixos e título
% %     xlabel('Tempo / \delta t');
% %     ylabel('s(t)');
% %     title(['Pulso retangular com \delta t = ' num2str(delta_t) ' s']);
% %     xticks(0:1:5);
% %     xticklabels({'0','\delta_t','2\delta_t','3\delta_t','4\delta_t','5\delta_t'});
% % end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 









% --- Parâmetros --- Nt = 1e5; % número de amostras deltas = [1e-7, 1e-5, 1e-3]; % três larguras de pulso Tmax = 5e-3; % tempo máximo (cobre o maior delta_t) % Vetor de tempo comum t = linspace(0, Tmax, Nt); figure; hold on; grid on; for k = 1:length(deltas) delta_t = deltas(k); % Pulso começa em t=0 e termina em delta_t s = double(t >= 0 & t <= delta_t); % Plotar no mesmo eixo (tempo absoluto) plot(t*1e3, s, 'LineWidth', 1.5, 'DisplayName', ... ['\delta t = ' num2str(delta_t) ' s']); end xlabel('Tempo (ms)'); ylabel('s(t)'); title('Pulsos retangulares transmitidos'); legend show;

%  atraso = [0, espalhamento_multip_linear];   % representar sinal transmitido e sinal recebido
% delta_t = 10^-7;    % largura do pulso em segundos
% Nt = 10^5;          % número de amostras
% 
% figure; hold on; grid on;
% 
% for k = 1:length(atraso)
%     [t, s] = gerar_pulso(atraso(k), delta_t, Nt);
%     plot(t/delta_t, s, 'LineWidth', 2);
% end
% 
% xlabel('Tempo absoluto - t (seg)');
% ylabel('s(t)');
% title('Pulso Transmitido e Recebido');
% grid on;
% xticks(0:1:5);                       % coloca ticks em 0,1,2,3,4,5
% xticklabels({'0','\delta_t','2\delta_t','3\delta_t','4\delta_t','5\delta_t'});
% 
% Implementação do sinal recebido
% sinal_recebido = 0;
% for n = 1:N
%     sinal_recebido = sinal_recebido + sqrt(pot_dispersa_normalizada)*exp(-j*2*pi*((fc+desvio_doppler_1*atraso_multipercurso_ordenado-desvio_doppler_1*t)))*s(t-tn);
% 
% 
%     --- Parâmetros ---
% delta_t = 1e-7;     % largura do pulso
% Nt = 1e5;           % número de amostras
% N = length(componente_multipercurso_ordenado); % nº de caminhos
% 
% --- Inicializar tempo e sinal recebido ---
% t = linspace(0, 5*delta_t, Nt);  
% sinal_recebido = zeros(1, Nt);  
% 
% --- Construir sinal recebido ---
% for n = 1:N
%     Gerar pulso atrasado do caminho n
%     [~, s_n] = gerar_pulso(componente_multipercurso_ordenado(n), delta_t, Nt);
% 
%     Peso complexo (amplitude, fase Doppler, etc.)
%     ganho_n = sqrt(pot_dispersa_normalizada(n)) .* ...
%               exp(-1j*2*pi*((fc+desvio_doppler_1)*componente_multipercurso_ordenado(n) ...
%               - desvio_doppler_1*t));
% 
%     Acumular no sinal recebido
%     sinal_recebido = sinal_recebido + ganho_n .* s_n;
% end
% 
% --- Plotar sinais ---
% figure; hold on; grid on;
% 
% Pulsos individuais
% for k = 1:length(atraso)
%     [t_p, s_p] = gerar_pulso(atraso(k), delta_t, Nt);
%     plot(t_p/delta_t, s_p, 'LineWidth', 1.5);
% end
% 
% Sinal resultante (real)
% plot(t/delta_t, real(sinal_recebido), 'k', 'LineWidth', 2);
% 
% xlabel('Tempo / \delta_t');
% ylabel('s(t)');
% legend('Pulsos individuais', 'Sinal recebido');
% title('Pulso Transmitido e Sinal Recebido');
