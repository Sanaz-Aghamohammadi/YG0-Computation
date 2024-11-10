function [cct, duv] = XYZ2CCT(XYZ)
    % XYZ2CCT, cacluate CCT of some color. Uses Ohno 2014 algorithm, specified
    % by TM30. Multiple XYZ input not supported.
    %
    % Triangular method from﻿Ohno, Y. (2014) ‘Practical use and calculation of
    % CCT and Duv’, LEUKOS - Journal of Illuminating Engineering Society of
    % North America, 10(1), pp. 47–55. doi: 10.1080/15502724.2014.839020.
    % Very dense data table, worst case error < approx. 0.5K
    if any(size(XYZ, 2) ~= 3)
        error('Input must be nx3 mat XYZ');
    end
    
    if (all(XYZ == [0 0 0]))
       cct = 0;
       duv = inf;
       return
    end
    
    persistent T uv1960Table;
    if isempty(T)
        load('PlanckLocus.mat', 'cct', 'uvY1960');
        T = cct;
        uv1960Table = uvY1960(:, 1:2);
    end
    
    cct = zeros(size(XYZ, 1), 1);
    duv = zeros(size(XYZ, 1), 1);
    
    for j = 1:size(XYZ, 1)
        
        if ~any(isnan(XYZ(j,:))) % skip NaN rows
        
        uv = XYZ2uvY1960(XYZ(j, :));
        uv = uv(:, 1:2);
        dist = sum((uv1960Table - uv).^2, 2).^.5;
        m = find(dist == min(dist));
        
        %% CCT is at the table limit
        if m == size(uv1960Table,1) || m == 1
            cct(j) = T(m);
            duv(j) = nan;
            continue;
        end
        
        %% Triangular Approach
        l = sum((uv1960Table(m-1,:) - uv1960Table(m+1,:)).^2, 2).^.5;
        x = (dist(m-1)^2 - dist(m+1)^2 + l^2) / (2*l);
        cct(j) = T(m-1) + (T(m+1) - T(m-1)) * (x / l) ;
        
        vTx = uv1960Table(m-1,2) + (uv1960Table(m+1,2) - uv1960Table(m-1,2))*(x/l);
        duv(j) = (dist(m-1)^2 - x^2).^.5 * sign(uv(2) - vTx);
        
        %% Parabolic Approach
        
        if abs(duv(j)) > 0.002
            
            X = (T(m+1) - T(m)) * (T(m-1) - T(m+1)) * (T(m) - T(m-1));
            a = ( (T(m-1)  * (dist(m+1) -  dist(m))) ...
                + (T(m)    * (dist(m-1) -  dist(m+1))) ...
                + (T(m+1)  * (dist(m)   -  dist(m-1)))) ...
                * X^-1;
            b = - ((T(m-1)^2  * (dist(m+1) - dist(m))    ) ...
                +  (T(m)  ^2  * (dist(m-1) - dist(m+1))  ) ...
                +  (T(m+1)^2  * (dist(m)   -   dist(m-1))  )) ...
                * X^-1;
            
            c = - ((dist(m-1)  * (T(m+1) - T(m))    * T(m)   * T(m+1) ) ...
                +  (dist(m)    * (T(m-1) - T(m+1))  * T(m-1) * T(m+1) )  ...
                +  (dist(m+1)  * (T(m)   - T(m-1))  * T(m-1) * T(m)  )) ...
                * X^-1;
            
            cct(j) = -(b/(2*a));
            duv(j) = sign(uv(2) - vTx) * (a * cct(j)^2 + b * cct(j) + c);
        end
        
        if duv(j) >= 0.05
            duv(j) = nan;
            continue;
        end
        
        else % if any NaN
            cct(j) = nan;
            duv(j) = nan;
        end
    end
end