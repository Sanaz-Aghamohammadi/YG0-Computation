function [illuminant, lambda] = getIlluminant(name, varargin)
    % GETOBSERVERFUNCTION Get a tristimulus observer function
    % 
    % [observer, lambda] = getObserevrFunction(NAME) gets the observer with the
    % name NAME at it's default wavelength spacing.
    % 
    % [observer, lambda] = getObserevrFunction(NAME, WL) gets the observer with the
    % name NAME at each wavelength in the vector WL. Interpolates from the
    % observers default wavelength vector to the requested vector using
    % pchip. Does not return NAN. 
    % 
    % WL must be a vector of real numbers. No
    % other checks are done, but WL SHOULD also be monotonic for most use
    % cases. 
    
    in = parseInputs(varargin{:});    
    
    % Persistent + isEmpty to only load files once
    persistent CIE_Illuminant_A CIE_Illuminant_D50 CIE_Illuminant_D65
    switch lower(name)
        case 'a'
            if isempty(CIE_Illuminant_A)
                load('CIE_Illuminants.mat', 'CIE_Illuminant_A');
            end
            illuminant = CIE_Illuminant_A.data;
            lambda = CIE_Illuminant_A.lambda;
            
        case 'd50'
            if isempty(CIE_Illuminant_D50)
                load('CIE_Illuminants.mat', 'CIE_Illuminant_D50');
            end
            illuminant = CIE_Illuminant_D50.data;
            lambda = CIE_Illuminant_D50.lambda;
        case 'd65'
            if isempty(CIE_Illuminant_D65)
                load('CIE_Illuminants.mat', 'CIE_Illuminant_D65');
            end
            illuminant = CIE_Illuminant_D65.data;
            lambda = CIE_Illuminant_D65.lambda;
    end
    
    if ~isempty(in.lambda)
        illuminant = interp1(lambda, illuminant, in.lambda, 'linear', 0);
        lambda = in.lambda;
    end
    
    function options = parseInputs(varargin)
        narginchk(0,1); % doc nargincheck
        persistent parser; 
        if isempty(parser)
            % Set up parser
            parser = inputParser();
            parser.FunctionName = mfilename;
            
            parser.addOptional('lambda', [], ...
                @(x) isvector(x) && all(isreal(x)));
        end
        parser.parse(varargin{:})
        options = parser.Results;
    end
end