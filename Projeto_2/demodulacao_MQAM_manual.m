function dados_demod = demodulacao_MQAM_manual(simbolos, M, varargin)
% Entrada: 
%   simbolos = vetor complexo com símbolos recebidos
%   M = ordem da modulação
%   varargin = parâmetros opcionais: 'UnitAveragePower', true/false
% Saída: dados_demod = vetor de inteiros entre 0 e M-1

    UnitAveragePower = true;
    
    for i = 1:2:length(varargin)
        if strcmpi(varargin{i}, 'UnitAveragePower')
            UnitAveragePower = varargin{i+1};
        end
    end
    
    raiz_M = sqrt(M);
    
    constelacao = zeros(1, M); % Cria a mesma constelação usada na função de modulação
    indice = 1;
    
    for i = 1:raiz_M
        for j = 1:raiz_M
            I = (2*i - 1 - raiz_M);
            Q = (2*j - 1 - raiz_M);
            constelacao(indice) = I + 1j*Q;
            indice = indice + 1;
        end
    end
    
    % Aplica mesma normalização usada na função de modulação
    if UnitAveragePower
        potencia_media = mean(abs(constelacao).^2);
        constelacao = constelacao / sqrt(potencia_media);
    end
    
    % Decisão por mínima distância manual
    dados_demod = zeros(1, length(simbolos));
    
    for k = 1:length(simbolos)
        distancias = abs(simbolos(k) - constelacao).^2; % Calcula distâncias para todos os pontos da constelação        
        [~, indice_min] = min(distancias); % Encontra ponto mais próximo
        dados_demod(k) = indice_min - 1;
    end
end
