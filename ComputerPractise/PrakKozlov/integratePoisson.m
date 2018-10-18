function A=integratePoisson(A_,w,t)
% (1-cos(norm(w)*t))/norm(w)/norm(w)*skewsym(w)*skewsym(w)
A=(eye(3)+sin(norm(w)*t)/norm(w)*koso(w)+(1-cos(norm(w)*t))/norm(w)/norm(w)*koso(w)*koso(w))*A_;
end