function A = delay_scale_ambiguity(r, t, g_rx_fun, tau, alpha)
%DELAY_SCALE_AMBIGUITY
%
% Computes:
%
% A_{g_rx,r}(tau,alpha)
%   = integral conj(g_rx(alpha*(t-tau))) ...
%              * sqrt(alpha) * r(t) dt
%
% g_rx_fun must accept arbitrary time arguments.

    r = r(:);
    t = t(:);

    if numel(r) ~= numel(t)
        error('r and t must have the same number of samples.');
    end

    if numel(t) < 2
        error('t must contain at least two samples.');
    end

    if alpha <= 0
        error('alpha must be positive.');
    end

    % Construct the delay-scale receive template directly.
    % No interp1 is required because g_rx_fun can be evaluated
    % directly at arbitrary time arguments.
    template = sqrt(alpha) .* g_rx_fun(alpha .* (t - tau));

    template = template(:);

    % Reuse the existing matched-filter function
    A = matchedFilterSample(r, t, template);
end



function A = matchedFilterSample(r, t, template)
%MATCHEDFILTERSAMPLE
% Compute one sampled matched-filter output:
%
%   A = integral conj(template(t)) * r(t) dt
%
% For uniformly sampled signals:
%
%   A approximately equals Ts * template^H * r

    r = r(:);
    t = t(:);
    template = template(:);

    if numel(r) ~= numel(t)
        error('r and t must have the same number of samples.');
    end

    if numel(template) ~= numel(r)
        error('template and r must have the same number of samples.');
    end

    if numel(t) < 2
        error('The time vector must contain at least two samples.');
    end

    TsLocal = t(2)-t(1);

    % Matched-filter inner product
    % MATLAB '' performs conjugate transpose.
    A = TsLocal*(template'*r);

end