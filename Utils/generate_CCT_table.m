clear 
cct(1) = 1000;

while (cct(end) < 25000)
    cct(end+1) = cct(end) * 1.0025;
end

cct(cct > 24950) = [];
cct = [cct 24990:25000];

[spec, wl] = CCT2Planckian(cct);
XYZ = spec2XYZ(spec, wl);
uvY = XYZ2uvY(XYZ);
uvY1960 = XYZ2uvY1960(XYZ);
xyY = XYZ2xyY(XYZ);

if 0
    save('../PlanckLocus.mat', 'uvY1960', 'cct')
end

clear XYZ2CCT

%% Test
temps = rand(1000,1) * 24000 + 1000;

err = zeros(100, numel(temps));
for T = temps(:)'
    samples(:, 2) = linspace(-0.05, 0.05, 100);
    samples(:, 1) = T;
    
    xy = CCTduv2xy(samples(:, 1), samples(:, 2));
    
    CalculatedTemps = xy2CCT(xy);
    err(:, T==temps) = abs(CalculatedTemps - T);
end

maxerr = quantile(err(:), .99)
if maxerr > .9;
    warning('Maximum error too high!'); 
end

[~, highErrTemps] = ind2sub(size(err), find(err(:) > 1));

if ~isempty(highErrTemps)
    fprintf('High Error Temp: %fK\n', temps(unique(highErrTemps)));
end