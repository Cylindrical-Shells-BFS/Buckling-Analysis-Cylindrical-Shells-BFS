function [ oV ] = linearScalarMap( iV, iMin, iMax, oMin, oMax )
% iV:   Input value in between the input range (iMin =< iV =< iMax) 
% oV:   Output value in between the output range (oMin =< oV =< oMax) 
% iMin: Minimum value of the input range
% iMax: Maximum value of the input range
% oMin: Minimum value of the input range
% oMax: Minimum value of the input range

li = iMax-iMin; % Length of the input range
lo = oMax-oMin; % Length of the output range

li_current = iV-iMin; % Length of iV along the input range

oV = li_current/li*lo + oMin; % Calculating the current output value
end

