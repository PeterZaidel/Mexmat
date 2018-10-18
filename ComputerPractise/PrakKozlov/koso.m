function [ k ] = koso( q )
k=zeros(3);
k=[0 q(3) -q(2);-q(3) 0 q(1);q(2) -q(1) 0];
end