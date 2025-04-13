function doubleArow(ax, point1, point2, infostr, position)
quiver(ax, point1(1), point1(2), point2(1) - point1(1), point2(2) - point1(2), 1, '-k');
quiver(ax, point2(1), point2(2), point1(1) - point2(1), point1(2) - point2(2), 1, '-k');

% text position 
delta = 0.001;
x = 0.5 * (point1(1) + point2(1));
y = 0.5 * (point1(2) + point2(2));
switch position
    case "top"
        y = y + 0.005;
        allign = "center";
    case "down"
        y = y - 0.005;
        allign = "center";
    case "left"
        x = x - delta * ceil(length(infostr)/2);
        allign = "right";
    case "right"
        x = x + delta * ceil(length(infostr)/2);
        allign = "left";
end

% place text info
text(ax, x, y, infostr, "HorizontalAlignment", allign, 'FontSize',8);
end