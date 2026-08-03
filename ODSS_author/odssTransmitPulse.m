%% Local function: basic chirplet and optional PHYDYAS window
function g_tx = odssTransmitPulse(tLocal, q, W, T, pulseType)

g_tx = complex(zeros(size(tLocal)));

valid = (tLocal >= 0) & (tLocal < T);

if ~any(valid)
    return;
end

tv = tLocal(valid);

freqScale = 1;% W/(sqrt(q)-1/sqrt(q));
f1 = freqScale/sqrt(q);
f2 = freqScale*sqrt(q);
kappa = (f2-f1)/T;

g0 = exp(1j*2*pi*(f1*tv + 0.5*kappa*tv.^2));

switch lower(pulseType)

    case 'rect'
        window = ones(size(tv));

    case 'phydyas'
        K = 3;
        A = [0.91143783, 0.41143783];

        window = ones(size(tv));

        for kk = 1:K-1
            window = window + 2*(-1)^kk*A(kk).*cos(2*pi*kk*tv/(K*T));
        end

    otherwise
        error('Unknown pulseType: %s', pulseType);
end

g_tx(valid) = window.*g0;

end