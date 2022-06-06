function [value,isterminal,direction] = borderFunc(t,x, xLow, xUp, yLow, yUp)
isterminal = 1;
direction = 0;
if x(1) > xUp || x(1) < xLow
    value = 0;
elseif x(2) > yUp || x(2) < yLow
    value = 0;
else
    value = 1;
end
end