function [delx,dely,delz,dels,bnd] = crunch2(G,Ae,Ai,x,y,z,s,rQ,rA,rC,rsz)

n=size(G,1);
ne=size(Ae,1);
ni=size(Ai,1);

S = sparse(0*speye(ni));
invZ = sparse(0*speye(ni));
for i=1:ni, S(i,i)=s(i); invZ(i,i)=1.0/z(i); end

invS = sparse(0*speye(ni));
Z = sparse(0*speye(ni));
for i=1:ni, invS(i,i)=1.0/s(i); Z(i,i)=z(i); end
CtZSinvC = Ai'*Z*invS*Ai;

mat2 = sparse(0*speye(n+ne));
mat2(1:n,1:n) = G+CtZSinvC;
mat2((n+1):(n+ne),1:n)=Ae;
mat2(1:n,(n+1):(n+ne))=Ae';

rhs2 = [-rQ-Ai'*invS*(Z*rC+rsz);-rA];
del2=mat2\rhs2;
delx2 = del2(1:n,1);
%dely2 = -del2((n+1:n+ne),1);
dely2 = del2((n+1:n+ne),1);

ZinvS = sparse(0*speye(ni));
for i=1:ni, ZinvS(i,i)=z(i)/s(i); end

delz2 = ZinvS*(-Ai*delx2-rC-invZ*rsz);

dels2 = invZ*(-rsz-S*delz2);

delx=delx2;
dely=dely2;
delz=delz2;
dels=dels2;

bnd = bound(s,dels,1);
bnd = bound(z,delz,bnd);


