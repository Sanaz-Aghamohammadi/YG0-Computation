function XYZ = reflectance2XYZ(reflectance, varargin)
% reflectance2XYZ converts spectral reflectance data to XYZ values.
% If no other param is given, CIE 1931 2 deg observer under D65
% illuminant is set to be default.
%
% All data along wavelength is recommended to align in the same
% direction. For example, if SPD is wavelength-by-samples, then
% illuminant SPD (if provided) should be wavelength-by-1.
% 
%
% Usage:
% XYZ = reflectance2XYZ(reflectance)
% XYZ = reflectance2XYZ(reflectance, Name, Value)
%
% @param: reflectance, n*m or m*n matrix
%                      SPD n wl interval and m samples
% 
% @return: XYZ, 3*m or m*3 matrix
%               tristimulus values, default with whiteY = 100
%
% NameValue pairs:
% 
% 'wavelength' ---- vector of wavelength values
%  | guessWL(reflectance) (default), wavelength guessed based on spec
%  | [wl], n*1 or 1*n wl vector;
%
% 'illuminant' ---- illuminant of source spd
%  | 'D65'(default)
%  | 'D50'
%  | 'A'
%  | [illuminant], n*1 or 1*n vector; source SPD
%  | CCT, scalar value; TM-30 reference illuminant will be used
%  | [Xn, Yn, Zn], XYZ values of white; XYZ2CCT() will be used
%
% 'observer' ---- color matching functions
%  | '1931' (dafault) | '2deg' | 'xyz' | 'xyz2'; CIE 1931 2 deg observer.
%  | '1964' | '10deg' | 'xyz10'; CIE 1964 10 deg observer.
%  | [CMF], 3*n or n*3 vector; costume color matching functions
%
% 'whiteY' ---- Y value of PRD
%  | 100 (default)
%  | Yn, non-neg scalar value
%
%
% See also guessWL, CCT2RefIlluminant, XYZ2CCT
%

in = parseInputs(varargin{:}); 

% WAVELENGTH
persistent wvl lastWvl
if isempty(wvl) || ~isequal(lastWvl,in.wavelength) || ~any(length(wvl)==size(reflectance)) % last conditional is fix for situation where you change the wavelength sampling of the reflectances between function calls, without specify wavelength range
    if isempty(in.wavelength) % if none given, guess
        warning('Assuming wavelength sample based on size of reflectance.')
        wvl = guessWL(reflectance(:, 1));
    else
        wvl = in.wavelength;
    end
    
    lastWvl = in.wavelength;
end

% ILLUMINANT
persistent ill lastIll
if isempty(ill) || ~isequal(lastIll,in.illuminant) || ~any(size(ill)==length(wvl))
    if ischar(in.illuminant) % illuminant name given
        ill = getIlluminant(in.illuminant,wvl);
    elseif isscalar(in.illuminant) % CCT given
        ill = CCT2RefIlluminant(in.illuminant,wvl);
    elseif length(in.illuminant)==3 % XYZ given
        ill = CCT2RefIlluminant(XYZ2CCT(in.illuminant),wvl);
    else % vector given
        ill = in.illuminant;
    end
    lastIll = in.illuminant;
end
% orientation
ill = ill(:);

% OBSERVER
persistent obs lastObs
if isempty(obs) || ~isequal(lastObs,in.observer) || ~any(size(obs)==length(wvl))
    if ischar(in.observer)
        obs = getObserverFunction(in.observer,wvl);
    else
        obs = in.observer;
    end
    
    lastObs = in.observer;
end
% orientation
if size(obs,1)~=length(wvl)
    obs=obs';
end
if size(obs,1)~=length(wvl)
    error('Observer must match wavelength in at least one dimension')
end

% Reflectance shape
if size(reflectance,1) ~= length(wvl)
    reflectance = reflectance';
    if size(reflectance,1) ~= length(wvl)
        error('Reflectance must match wvl in at least one dimension.')
    end
    switched = 1;
else
    switched = 0;
end
if size(reflectance,2) == size(reflectance,1)
    warning('Square reflectance matrix assumed to be wvlxN orientation.')
end

% calculate XYZ
XYZ = obs'*diag(ill)*reflectance * (in.whiteY / (obs(:,2)'*ill));

% output shape
if switched
    XYZ = XYZ';
end

end
    
    

function options = parseInputs(varargin)
    persistent parser; 
    if isempty(parser)
        % Set up parser
        parser = inputParser();
        parser.FunctionName = mfilename;
        
        parser.addParameter('wavelength', [], ...
            @(x) isvector(x) && all(isreal(x)));

        parser.addParameter('illuminant', 'D65');
        
        parser.addParameter('observer', '1931');
        
        parser.addParameter('whiteY', 100, ...
            @(x) isreal(x) && x>0);
    end
    parser.parse(varargin{:})
    options = parser.Results;
end
