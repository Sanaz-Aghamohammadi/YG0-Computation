function xyY = XYZ2xyY( XYZ )
    % XYZ2xyY: compute xyY (chromaticity & luminance) from XYZ
    %
    % usage:  xyY = mjmXYZ2xyY( XYZ )
    %
    % input:  XYZ  (Nx3 or 3xN) XYZ tristimulus values (normalized or not; doesn't matter)
    %
    % output: xyY  (Nx3 or 3xN) xy chromaticity and Y luminance
    %
    % MJMurdoch 20160802
    
    % error check
    if nargin < 1
        help mfilename
    end
    
    % check orientation
    switched = size(XYZ,2)~=3;
    if switched
        XYZ = XYZ';
    end
    if size(XYZ,2) == size(XYZ,1)
        warning('Square matrix assumed to be Nx3 orientation.')
    end
    
    % compute xy from XY, preserve luminance Y
    xyY = [ XYZ(:,1:2)./sum(XYZ,2) XYZ(:,2) ];
    
    if switched
        xyY=xyY';
    end
end