function [t, s] = gerar_pulso(atraso, delta_t, Nt)
% gerar_pulso - Gera um pulso transmitido s(t)
%
% Entradas:
%   atraso   - atraso inicial do pulso (s)
%   delta_t  - largura do pulso (s)
%   Nt       - número de amostras do domínio temporal
%
% Saídas:
%   t        - vetor tempo absoluto t ∈ [0, 5*delta_t]
%   s        - vetor do pulso transmitido s(t)

    % --- 1) Definir o vetor tempo ---
    t = linspace(0, 5*delta_t, Nt);   % tempo absoluto

    % --- 2) Gerar o pulso ---
    % Aqui usamos um pulso retangular centrado no atraso inicial
    % O pulso vale 1 durante [atraso, atraso + delta_t] e 0 fora
    s = double(t >= atraso & t <= atraso + delta_t);

end


