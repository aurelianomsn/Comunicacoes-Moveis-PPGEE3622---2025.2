function simbolos = modulacao_MQAM_manual(dados, M, varargin)
% Entrada: 
%   dados = vetor de inteiros entre 0 e M-1
%   M = ordem da modulação (4, 16, 64, ...)
%   varargin = parâmetros opcionais: 'UnitAveragePower', true/false
% Saída: simbolos = vetor complexo com símbolos M-QAM

    UnitAveragePower = true;
    
    for i = 1:2:length(varargin)
        if strcmpi(varargin{i}, 'UnitAveragePower')
            UnitAveragePower = varargin{i+1};
        end
    end
    
    raiz_M = sqrt(M);
    
    constelacao = zeros(1, M); % Cria constelação M-QAM
    indice = 1;
    
    for i = 1:raiz_M
        for j = 1:raiz_M
            I = (2*i - 1 - raiz_M); % Coordenadas da constelação
            Q = (2*j - 1 - raiz_M); % Coordenadas da constelação
            
            constelacao(indice) = I + 1j*Q;
            indice = indice + 1;
        end
    end
    
    % Normalização para potência unitária
    if UnitAveragePower
        potencia_media = mean(abs(constelacao).^2);
        constelacao = constelacao / sqrt(potencia_media);
    end
    
    simbolos = constelacao(dados + 1);  
end
