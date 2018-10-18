function [ L ] = eulers2L( eulers )
%EULERS2L Summary of this function goes here
%   Detailed explanation goes here

    a = eulers;
    L =            [cos(a(1)) * cos(a(3)) - sin(a(1)) * sin(a(2)) * sin(a(3)), sin(a(1)) * cos(a(3)) + cos(a(1)) * sin(a(2)) * sin(a(3)), -cos(a(2)) * sin(a(3)) ;
                    - sin(a(1)) * cos(a(2))                                  , cos(a(1)) * cos(a(2))                                    , sin(a(2))              ;
                    cos(a(1)) * sin(a(3)) + sin(a(1)) * sin(a(2)) * cos(a(3)), sin(a(1)) * sin(a(3)) - cos(a(1)) * sin(a(2)) * cos(a(3)), cos(a(2)) * cos(a(3)) ];
end

