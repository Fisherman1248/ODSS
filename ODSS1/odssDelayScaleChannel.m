function [r, rs, w, tRx, Ps, Pn, pathSignals] = ...
    odssDelayScaleChannel(sTx, Fs, taup, alphap, hp, SNRdB, fc)
% odssDelayScaleChannel
%
% Implements the wideband delay-scale channel
%
%   rs(t) = sum_p hp(p) * sqrt(alphap(p)) ...
%           * sTx(alphap(p) * (t - taup(p)))
%
% Inputs:
%   sTx       - transmitted waveform samples
%   Fs        - sampling frequency, Hz
%   taup      - path delays, seconds
%   alphap    - path time-scale factors, alpha_p > 0
%   hp        - complex path gains
%   SNRdB     - received waveform SNR in dB; use Inf for no noise
%   fc        - optional carrier frequency for a LOWPASS baseband envelope
%
% Outputs:
%   r           - noisy received waveform
%   rs          - noise-free received waveform
%   w           - complex AWGN
%   tRx         - receiver sampling times
%   Ps          - measured received signal power
%   Pn          - noise power
%   pathSignals - contribution from every channel path
%
% Important:
%   If sTx already contains the physical carrier/chirp frequencies,
%   set fc = 0.
%
%   If sTx is a complex lowpass envelope centered at 0 Hz, set fc to
%   the physical carrier frequency.

    if nargin < 7 || isempty(fc)
        fc = 0;
    end

    if nargin < 6 || isempty(SNRdB)
        SNRdB = Inf;
    end

    % Use column vectors consistently
    sTx    = sTx(:);
    taup   = taup(:);
    alphap = alphap(:);
    hp     = hp(:);

    Np = numel(hp);

    if numel(taup) ~= Np || numel(alphap) ~= Np
        error('taup, alphap, and hp must have the same length.');
    end

    if any(alphap <= 0)
        error('Every scale factor alphap must be positive.');
    end

    if any(taup < 0)
        error(['This implementation assumes nonnegative relative delays. ' ...
               'Set the first arrival to tau = 0.']);
    end

    if Fs <= 0
        error('Fs must be positive.');
    end

    %% Transmit-side time grid

    Ntx = numel(sTx);
    Ttx = Ntx / Fs;

    tTx = (0:Ntx-1).' / Fs;

    %% Receiver observation duration
    %
    % Path p is active when
    %
    %   0 <= alpha_p(t - tau_p) < Ttx
    %
    % therefore
    %
    %   tau_p <= t < tau_p + Ttx/alpha_p.

    tEnd = max(taup + Ttx ./ alphap);

    Nrx = max(1, ceil(tEnd * Fs));
    tRx = (0:Nrx-1).' / Fs;

    %% Apply the delay-scale channel

    pathSignals = complex(zeros(Nrx, Np));
    rs = complex(zeros(Nrx, 1));

    for p = 1:Np

        % Argument of transmitted waveform:
        %
        % u_p(t) = alpha_p (t - tau_p)
        u = alphap(p) * (tRx - taup(p));

        % Fractional delay and arbitrary time scaling are both handled
        % through interpolation.
        %
        % Samples outside the support of sTx are set to zero.
        sScaledDelayed = interp1( ...
            tTx, sTx, u, 'linear', 0);

        % When sTx is a LOWPASS complex envelope, time scaling of the
        % physical carrier produces the extra baseband phase:
        %
        % exp(j*2*pi*fc*(alpha_p(t-tau_p) - t)).
        %
        % When sTx already contains its physical frequency, fc should be 0.
        if fc ~= 0
            carrierPhase = exp( ...
                1i * 2*pi*fc .* (u - tRx));

            sScaledDelayed = sScaledDelayed .* carrierPhase;
        end

        % Paper's sqrt(alpha_p) energy-normalization factor
        pathSignals(:,p) = ...
            hp(p) * sqrt(alphap(p)) * sScaledDelayed;

        rs = rs + pathSignals(:,p);
    end

    %% Add complex AWGN

    Ps = mean(abs(rs).^2);

    if isinf(SNRdB)
        Pn = 0;
        w = complex(zeros(size(rs)));
    else
        Pn = Ps * 10^(-SNRdB/10);

        w = sqrt(Pn/2) .* ...
            (randn(size(rs)) + 1i*randn(size(rs)));
    end

    r = rs + w;
end