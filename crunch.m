function [delx,dely,delz,dels,bnd] = crunch(G,Ae,Ai,x,y,z,s,rQ,rA,rC,rsz)

n=size(G,1);
ne=size(Ae,1);
ni=size(Ai,1);

S = sparse(0*speye(ni));
invZ = sparse(0*speye(ni));
for i=1:ni, S(i,i)=s(i); invZ(i,i)=1.0/z(i); end

mat = sparse(0*speye(n+ne+ni));
mat(1:n,1:n) = G;
mat((n+1):(n+ne),1:n)=Ae;
mat(1:n,(n+1):(n+ne))=Ae';
mat((n+ne+1):(n+ne+ni),1:n)=Ai;
mat(1:n,(n+ne+1):(n+ne+ni))=Ai';

%mat((n+ne+1):(n+ne+ni),(n+ne+1):(n+ne+ni)) = -invZ*S;
for i=1:ni,
    mat(n+ne+i,n+ne+i)=-1.0/z(i)*s(i);
end

rhs = [-rQ;-rA;-rC-invZ*rsz];

del=mat\rhs;
delx = del(1:n,1);
%dely = -del((n+1:n+ne),1);
dely = del((n+1:n+ne),1);
delz = -del((n+ne+1:n+ne+ni),1);
dels = invZ*(-rsz-S*delz);

bnd = bound(s,dels,1);
bnd = bound(z,delz,bnd);
