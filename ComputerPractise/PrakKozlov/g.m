function [ g ] = g( phi, h) 
c2=0.005302;
c3=0.000007;
c4=0.00014;
ge=9.78030;
a=6378137;
g=ge*(1+c2*(sin(phi))^2-c3*(sin(2*phi))^2-2*h/a)-c4;
end