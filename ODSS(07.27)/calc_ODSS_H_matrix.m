function H = calc_ODSS_H_matrix(tau_path, A_path, a_path, q, W, T, pulseType, t, M_scale)
% Calculate H_{n,m}[n,m] according to Eq. (65).
% The coefficients are placed on the diagonal of H.

    tau_path = tau_path(:);
    A_path = A_path(:);
    a_path = a_path(:);
    t = t(:);
    M_scale = M_scale(:).';

    Path_num = numel(A_path);
    Nscale = numel(M_scale);
    M_tot = sum(M_scale);

    if numel(tau_path) ~= Path_num || numel(a_path) ~= Path_num
        error('tau_path, A_path, and a_path must have the same length.');
    end

    alpha_path = 1 + a_path;

    if any(alpha_path <= 0)
        error('All scale factors alpha_i = 1+a_i must be positive.');
    end

    g_tx_base = odssTransmitPulse(t, q, W, T, pulseType);
    g_tx_base = g_tx_base(:);

    g_rx_fun = @(u) odssTransmitPulse(u, q, W, T, pulseType);

    H = complex(zeros(M_tot,M_tot));

    idx = 1;

    for n = 0:Nscale-1
        Mn = M_scale(n+1);

        for m = 0:Mn-1
            H_nm = 0;

            for p = 1:Path_num
                tau_argument = q^n*((m/W)*(alpha_path(p)-1) - alpha_path(p)*tau_path(p));
                scale_argument = 1/alpha_path(p);

                A_value = delay_scale_ambiguity(g_tx_base, t, g_rx_fun, tau_argument, scale_argument);

                H_nm = H_nm + A_path(p)*A_value;
            end

            H(idx,idx) = H_nm;
            idx = idx + 1;
        end
    end
end