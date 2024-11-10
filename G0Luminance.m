function YG0 = G0Luminance(xy_input, varargin)

% G0Luminance computes the G0 luminance (YG0) for the given chromaticity coordinates.
% If no additional parameters are specified, the function uses the CIE 1931 2-degree observer 
% under the D65 illuminant as the default.
%
% Usage:
% YG0 = G0Luminance(xy_input)
% YG0 = G0Luminance(xy_input, Name, Value)
%
% @param: 
%   xy_input: Chromaticity coordinates, provided as an m*2 or 2*m matrix,
%             where each row (or column) represents the x and y chromaticity values.
% 
% @return:
%   YG0: A column vector (m*1 or 1*m) of the computed YG0 luminance values.
%
% Name-Value pairs:
%
% 'illuminant' ---- Illuminant for the source SPD.
%   | 'D65' (default)
%   | 'D50'
%   | 'A'
%   | [illuminant], an n*1 or 1*n vector representing the SPD of the illuminant.
%   | CCT, a scalar value; TM-30 reference illuminant will be used.
%   | [Xn, Yn, Zn], the XYZ values of white; the function XYZ2CCT() will be used to compute the illuminant.
%
% 'illuminantWavelength' ---- Wavelengths corresponding to the custom illuminant SPD.
%   | A numeric vector specifying the wavelengths (in nm) for the custom illuminant.
%   | Required if 'illuminant' is provided as a vector.
%
% 'observer' ---- Color matching functions (CMFs) to use.
%   | '1931' (default) or '2deg' or 'xyz' or 'xyz2'; for the CIE 1931 2-degree observer.
%   | '1964' or '10deg' or 'xyz10'; for the CIE 1964 10-degree observer.
%   | [CMF], a custom 3*n or n*3 matrix for user-defined CMFs.
%
% 'observerWavelength' ---- Wavelengths corresponding to the custom observer CMFs.
%   | A numeric vector specifying the wavelengths (in nm) for the custom observer.
%   | Required if 'observer' is provided as a vector.
%
% 'whiteY' ---- The Y value of the reference white.
%   | 100 (default)
%   | Yn, a non-negative scalar value.
%
%
% Note: Since the optimal color dataset has a 1 nm wavelength step, it is recommended that 
% all wavelength-dependent inputs (illuminants and CMFs) also have a 1 nm interval. Otherwise, interpolation will be used.
%
% See also: CCT2RefIlluminant, XYZ2CCT, getIlluminant,
% getObserverFunction,reflectance2XYZ, XYZ2xyY.

    in = parseInputs(varargin{:});
    wvl = (400:1:700)';

    % ILLUMINANT
    persistent ill lastIll
    if isempty(ill) || ~isequal(lastIll, in.illuminant) || ~any(size(ill) == length(wvl))
        if ischar(in.illuminant) % illuminant name given
            ill = getIlluminant(in.illuminant, wvl);
        elseif isscalar(in.illuminant) % CCT given
            ill = CCT2RefIlluminant(in.illuminant, wvl);
        elseif length(in.illuminant) == 3 % XYZ given
            ill = CCT2RefIlluminant(XYZ2CCT(in.illuminant), wvl);
        else % vector given
            if isempty(in.illuminantWavelength)
                error('When providing a custom illuminant vector, you must also provide ''illuminantWavelength''.');
            end
            if length(in.illuminantWavelength) ~= length(in.illuminant)
                error('The length of ''illuminantWavelength'' (%d) must match the length of the illuminant vector (%d).', length(in.illuminantWavelength), length(in.illuminant));
            end
            if any(in.illuminantWavelength < 400) || any(in.illuminantWavelength > 700)
                warning('Custom illuminant wavelengths are outside the 400-700 nm range. Extrapolation may be inaccurate.');
            end
            if length(unique(in.illuminantWavelength)) < 2
                error('''illuminantWavelength'' must contain at least two unique values for interpolation.');
            end
            % Interpolate the illuminant to fit the default wavelength range
            warning('Interpolating the illuminant to match the 400-700 nm range with 1 nm interval.');
            ill = interp1(in.illuminantWavelength, in.illuminant, wvl, 'spline', 'extrap');
        end
        lastIll = in.illuminant;
    end
    % Ensure column vector
    ill = ill(:);

    % OBSERVER
    persistent obs lastObs
    if isempty(obs) || ~isequal(lastObs, in.observer) || ~any(size(obs) == length(wvl))
        if ischar(in.observer)
            obs = getObserverFunction(in.observer, wvl);
        else
            if isempty(in.observerWavelength)
                error('When providing a custom observer vector, you must also provide ''observerWavelength''.');
            end
            if size(in.observer, 2) == 3 || size(in.observer, 1) == 3
                % Assuming CMFs are provided as a 3xn or nx3 matrix
                if length(in.observerWavelength) ~= length(in.observer)
                    error('The length of ''observerWavelength'' (%d) must match the length of the observer vector (%d).', ...
                          length(in.observerWavelength), length(in.observer));
                end
                if any(in.observerWavelength < 400) || any(in.observerWavelength > 700)
                    warning('Custom observer wavelengths are outside the 400-700 nm range. Extrapolation may be inaccurate.');
                end
                if length(unique(in.observerWavelength)) < 2
                    error('''observerWavelength'' must contain at least two unique values for interpolation.');
                end
                % Interpolate the observer to fit the default wavelength range
                warning('Interpolating the observer to match the 400-700 nm range with 1 nm interval.');
                if size(in.observer, 1) == 3
                    obs = interp1(in.observerWavelength, in.observer', wvl, 'spline', 'extrap')';
                else
                    obs = interp1(in.observerWavelength, in.observer, wvl, 'spline', 'extrap');
                end
            lastObs = in.observer;
            end
        end
    end


    % Ensure correct orientation
    if size(obs, 1) ~= length(wvl)
        obs = obs';
    end

    % xy input
    if size(xy_input, 2) ~= 2
        xy_input = xy_input';
        if size(xy_input, 2) ~= 2
            error('Matrix must be 2*m or m*2.');
        end
        switched = true;
    else
        switched = false;
    end

    if size(xy_input, 2) == size(xy_input, 1)
        warning('Square xy matrix assumed to be m*2 orientation.');
    end

    % Load lookup table
    persistent cachedResultsTable
    if isempty(cachedResultsTable)
        if ~exist('lookupTable.mat', 'file')
            error('lookupTable.mat not found. Please ensure it is in the MATLAB path.');
        end
        tmp = load('lookupTable.mat', 'resultsTable');
        if ~isfield(tmp, 'resultsTable')
            error('lookupTable.mat does not contain ''resultsTable''.');
        end
        cachedResultsTable = tmp.resultsTable;
    end
    resultsTable = cachedResultsTable;


    XYZ = reflectance2XYZ(resultsTable.Reflectance, 'wavelength', wvl, 'illuminant', ill, 'observer', obs, 'whiteY', in.whiteY);
    xyY = XYZ2xyY(XYZ);

    % NEAREST NEIGHBOR SEARCH
    diff_x = xyY(:,1) - xy_input(:,1)'; 
    diff_y = xyY(:,2) - xy_input(:,2)'; 
    distances_sq = diff_x.^2 + diff_y.^2; 
    
    [~, nearest_indices] = min(distances_sq, [], 1); 
    
    YG0 = XYZ(nearest_indices, 2);
    
    % Adjust output orientation if necessary
    if switched
        YG0 = YG0';
    end

end


function options = parseInputs(varargin)
    persistent parser; 
    if isempty(parser)
        % Set up parser
        parser = inputParser();
        parser.FunctionName = mfilename;

        parser.addParameter('illuminant', 'D65');

        parser.addParameter('illuminantWavelength', [], @(x) isnumeric(x) && isvector(x));

        parser.addParameter('observer', '1931');

        parser.addParameter('observerWavelength', [], @(x) isnumeric(x) && isvector(x));

        parser.addParameter('whiteY', 100, ...
            @(x) isreal(x) && x > 0);
    end
    parser.parse(varargin{:})
    options = parser.Results;
end
