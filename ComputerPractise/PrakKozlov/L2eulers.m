function [ eulers ] = L2eulers(L)
    gam = atan2(-L(1,3), L(3,3));
    th = atan2(L(2,3), sqrt(L(2,1)^2 + L(2,2)^2));
    csi = atan2(-L(2,1), L(2,2));
    eulers = [-csi, gam, th];
%L2EULERS Summary of this function goes here
%   Detailed explanation goes here


end

% def orient_matrix_to_eulers(L):
%     gam = np.arctan(-L[0,2]/L[2,2])
%     th = np.arctan(L[1,2]/np.sqrt(L[1,0]**2 + L[1,1]**2))
%     csi = np.arctan(-L[1,0]/L[1,1])
%     eulers = np.array([csi, th, gam])
%     return eulers