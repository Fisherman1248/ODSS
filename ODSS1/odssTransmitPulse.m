function g_tx = odssTransmitPulse(tLocal, q, W, T, pulseType)
% odssTransmitPulse generates the transmit pulse g_tx(tLocal)
%
% Inputs:
%   tLocal    - local pulse-time argument, scalar or vector
%   q         - geometric scale ratio
%   W         - basic chirplet bandwidth in Hz
%   T         - basic chirplet duration in seconds
%   pulseType - 'rect' or 'phydyas'
%
% Output:
%   g_tx      - transmit pulse g_tx(tLocal)

    % Pulse is zero outside [0,T)
    g_tx = complex(zeros(size(tLocal)));

    valid = (tLocal >= 0) & (tLocal < T);

    if ~any(valid)
        return;
    end

    tValid = tLocal(valid);

    %% Basic chirplet g_0(t), Eq. (71)

    % Convert normalized frequencies 1/sqrt(q), sqrt(q)
    % into physical frequencies in Hz.
    %% Basic chirplet, Eq. (71)

% Normalized frequencies exactly as defined in the paper
f1Normalized = 1/sqrt(q);
f2Normalized = sqrt(q);

% Convert the normalized chirplet bandwidth to W Hz
frequencyScale = W / (f2Normalized - f1Normalized);

f1 = frequencyScale * f1Normalized;
f2 = frequencyScale * f2Normalized;

chirpRate = (f2-f1)/T;

g0 = exp(1j*2*pi .* ...
    (f1*tValid + 0.5*chirpRate*tValid.^2));

    %% Pulse-shaping window g_w(t), Eq. (72)

    switch lower(pulseType)

        case 'rect'

            gw = ones(size(tValid));

        case 'phydyas'

            K = 3;
            A = [0.91143783, 0.41143783];

            gw = ones(size(tValid));

            for r = 1:K-1
                gw = gw + 2*(-1)^r*A(r) .* cos(2*pi*r*tValid/(K*T));
            end

        otherwise

            error('Unknown pulse type: %s', pulseType);
    end

    %% g_tx(t) = g_w(t)g_0(t)

    g_tx(valid) = gw .* g0;

end