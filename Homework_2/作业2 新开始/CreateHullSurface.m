function y = CreateHullSurface(L, B, D, x, z)
    y = 0.5 * B * (1 - (x / (L * 0.5))^2) * (1 - (z / D)^2)
%myFun - Description
%
% Syntax: output = myFun(input)
%
% Long description
    
end