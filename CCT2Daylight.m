function [spd, lambda] = CCT2Daylight(cct, lambda)
arguments
    cct
    lambda (:,1) = nan
end
    % CCT2Daylight Returns Daylight SPD of various Color Tempuratures according
    % to CIE015:2004
    %
    % [spd, lambda] = CCT2Daylight(cct)
    % *cct*, A vector of CCT values to calculate Daylight SPDs for. 
    % *spd*, an M x N matrix of SPDs. M is the number of wl samples in lambda and
    % N is the number of elements in CCT.
    % *lambda*, the wavelength index for spd. 
    
    if ~isvector(cct)
        error('CCT must be a vector!');
    end
    persistent wl S0 S1 S2;
    if isempty(wl)
        load('daylightVectors.mat', 'wl', 'S0', 'S1', 'S2');
    end

    
    cctL = cct >= 4000 & cct < 7000;
    cctH = cct >= 7000 & cct <= 25000;
    
    Xd = nan(numel(cct), 1);
    Xd(cctL) =...
        - 4.6070e9  ./  cct(cctL).^3 ...
        + 2.9678e6  ./  cct(cctL).^2 ...
        + 0.09911e3 ./  cct(cctL)... 
        + 0.244063;
    Xd(cctH) =...
        - 2.0064e9  ./  cct(cctH).^3 ...
        + 1.9018e6  ./  cct(cctH).^2 ...
        + 0.24748e3 ./  cct(cctH)... 
        + 0.237040;
    
    Yd = -3.*(Xd.^2) + (2.87.*Xd) - 0.275;
    
    M1 = (-1.3515 - 1.7703 .* Xd + 5.9114 .* Yd) ./...
        (0.0241 + 0.2562 .* Xd - 0.7341 .* Yd);
    
    M2 = (0.03 - 31.4424 .* Xd + 30.0717 .* Yd) ./...
        (0.0241 + 0.2562 .* Xd - 0.7341 .* Yd);
    
    spd = S0 + S1 * M1' + S2 * M2';
    
    if isnan(lambda)
        lambda = wl;
    else
        spd = interp1(wl, spd, lambda, 'spline', 0);
    end
    
    if(any(cct > 25000))
        warning('Daylight Series Illuminant Model is not valid for CCT > 25000k, Reverting to Planckian');
        spd(:, cct>25000) = CCT2Planckian(cct(cct>=25000));
    end
end