function uvY = XYZ2uvY1960(XYZ)
    u = (4 .* XYZ(:,1)) ./ (XYZ(:,1) + 15*XYZ(:,2) + 3*(XYZ(:,3)));
    v = (6 .* XYZ(:,2)) ./ (XYZ(:,1) + 15*XYZ(:,2) + 3*(XYZ(:,3)));
    uvY = [u v XYZ(:,2)];
end