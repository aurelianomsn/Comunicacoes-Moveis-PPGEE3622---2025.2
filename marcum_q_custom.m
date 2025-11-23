function Q = marcum_q_custom(a, b)
    % Q1(a,b) = ∫_b^∞ x * exp(-(x^2 + a^2)/2) * I0(a*x) dx

    if b > 50
        Q = 0; % Para "b" muito grande, a integral é praticamente zero.
        return;
    end

    if a < 1e-10
        Q = exp(-b^2/2); % Para "a" muito pequeno, há aproximação analítica.
        return;
    end

    N = 2000;  % número de pontos na integração
    x_max = b + 10 * sqrt(a^2 + b^2);
    if x_max <= b
        Q = 0;
        return;
    end

    x = linspace(b, x_max, N);
    dx = x(2) - x(1);

    integrando = x .* exp(-(x.^2 + a^2)/2) .* besseli(0, a.*x);

    % Regra de Simpson
    Q = (dx/3) * (integrando(1) + integrando(end) +  4*sum(integrando(2:2:end-1)) + 2*sum(integrando(3:2:end-2)));
    Q = max(0, Q);
end

