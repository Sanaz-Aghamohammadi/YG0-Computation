function [spd, wl] = CCT2RefIlluminant(cct,wl)
    %CCT2Refilluminant returms TM-30 Reference Illuminants. cct can be a vector
    %of values. 
    cct = cct(:)';
    
    if nargin<2
        wl = (380:780)';
    else
        wl = wl(:); % make column vector
    end
    spd = zeros(numel(wl), numel(cct));
    
    % cct<4000
    if(cct(cct<4000))
        [spec, lambda] = CCT2Planckian(cct(cct<4000));
        spd(:, cct<4000) = interp1(lambda, spec, wl, 'spline', 0);
    end
    
    % cct >= 4000 & cct < 5000
    if(~isempty(cct(cct >= 4000 & cct < 5000)))   
        specP = CCT2Planckian(cct(cct >= 4000 & cct < 5000), wl);
        specD = CCT2Daylight(cct(cct >= 4000 & cct < 5000),  wl);
    
        XYZ = spec2XYZ([specP specD], wl);
            
        specP = 100 * specP ./ XYZ(1:size(specP,2),2)';
        specD = 100 * specD ./ XYZ(size(specP,2)+1:end,2)';
            
        spd(:, cct >= 4000 & cct < 5000) = ((5000 - cct(cct >= 4000 & cct < 5000)) ./ 1000) .* specP ...
            + (1-((5000 - cct(cct >= 4000 & cct < 5000)) ./ 1000)) .* specD;
    end
    
    if(~isempty(cct(cct>=5000)))
        % cct > 5000 
        spd(:, cct >= 5000) = CCT2Daylight(cct(cct >= 5000), wl);
    end

    XYZ = spec2XYZ(spd, wl, '10deg');
    spd = (100 * spd) ./ XYZ(:, 2)';
end


