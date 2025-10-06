function validateLambda(wl)
if ~isvector(wl) && any(~isreal(wl))
    error("Wavelength must be a vector of real numbers");
end
if any(diff(wl, 2) ~= 0)
    error("Wavelength must be equal spacing!");
end
end