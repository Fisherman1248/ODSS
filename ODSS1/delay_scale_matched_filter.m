function A = delay_scale_matched_filter(r, t, g_rx, tau, alpha)
%DELAY_SCALE_MATCHED_FILTER
% Compute:
%
% A(tau,alpha) = integral conj(g_rx(alpha*(t-tau))) * sqrt(alpha) * r(t) dt
%
% Inputs:
%   r     : received signal samples
%   t     : time axis
%   g_rx  : receiver pulse function handle
%   tau   : delay hypothesis
%   alpha : scale hypothesis
%
% Output:
%   A     : matched-filter output

    r = r(:);
    t = t(:);

    if numel(r) ~= numel(t)
        error('r and t must have the same length.');
    end

    if alpha <= 0
        error('alpha must be positive.');
    end

    % Delay-scale matched-filter template
    template = sqrt(alpha) .* g_rx(alpha .* (t - tau));

    template = template(:);

    % MATLAB ' is conjugate transpose
    A = trapz(t, conj(template) .* r);

    % For uniform sampling, this is equivalent to:
    % Ts = t(2) - t(1);
    % A = Ts * (template' * r);
end