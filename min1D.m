function [f,p]= min1D (lambda,v0,dir,t,Q,n)
ones_vector=ones(n+1,1);
[R,p] = chol(Q+diag(v0+lambda*dir));      
determ=(det(R))^2;
if p==0
    f = t*(ones_vector'*(v0+lambda*dir))-log(determ);
else
f=inf;
end

end