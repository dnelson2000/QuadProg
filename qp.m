
clear
Gijtest70

scale = 10.0;

G=scale*G+1e-6*speye(size(G,1),size(G,1));


Ae=sparse(0,size(Ai,2));
be=[];
ne=0;
ni=size(Ai,1);
ni2=0;

%bi=-bi/10000*0.1;
thicknesses = -0*bi;

bi = -0*bi;
bi2 = zeros(0,1); %bi2 = zeros(size(Ai,1),1);
Ai2=sparse(size(Ai,1),size(Ai,2));
for ci=1:size(Ai,1),
    nodes=find(Ai(ci,:)~=0);
    if ( size(nodes,2) == 2 )
        Ai2(1+ni2,:)=Ai(ci,:);
        bi2(1+ni2,1)=thicknesses(ci,1);
        ni2=ni2+1;
    end    
end

for i=1:size(lo),
    if ( abs(lo(i)-up(i)) < 1e-6 ),
        Ae(1+ne,i)=1;
        be(1+ne,1)=lo(i);
        ne=ne+1;
      else
          if ( lo(i) > -10000000 ),
              Ai(1+ni,i)=-1;
              bi(1+ni,1)=-lo(i);
              ni=ni+1;
              Ai2(1+ni2,i)=-1;
              bi2(1+ni2,1)=-lo(i);
              ni2=ni2+1;
          end
          if ( up(i) < 10000000 ),
              Ai(1+ni,i)=1;
              bi(1+ni,1)=up(i);
              ni=ni+1;     
              Ai2(1+ni2,i)=1;
              bi2(1+ni2,1)=up(i);
              ni2=ni2+1;     
          end
    end
end

Ai=Ai2;
bi=bi2;


% these are flipped so instead of Ai*x<bi we have -Ai*x>-bi
Ai=-Ai;
bi=-bi;
xorig=lux(:,3);

n = size(G,1);
ne = size(Ae,1);
ni = size(Ai,1);
c=zeros(n,1);
tsig=3;

nrm = max(max(G));
componentnrm = max(max(Ae));
nrm = max(componentnrm,nrm);
componentnrm = max(max(Ai));
nrm = max(componentnrm,nrm);

componentnrm = max(abs(bi));
nrm = max(componentnrm,nrm);
componentnrm = max(abs(be));
nrm = max(componentnrm,nrm);
componentnrm = max(abs(c));
nrm = max(componentnrm,nrm);

sdatanrm = sqrt(nrm);
alpha = sdatanrm;
beta = sdatanrm;

x=zeros(n,1);
y=zeros(ne,1);
s=alpha*ones(ni,1);
z=beta*ones(ni,1);
%
if (0),
x=xorig;
s=-(Ai*xorig-bi);
coeff=[Ae' Ai'];
yz=inv(coeff'*coeff+0.001*eye(size(coeff,2)))*coeff'*(G*xorig);
y=-yz(1:size(Ae,1));
z=yz((size(Ae,1)+1):size(yz));
end

rQ = c + G*x - Ae'*y - Ai'*z;
rA = Ae * x - be;
rC = Ai*x - bi - s;
rsz = s.*z;

[delx,dely,delz,dels,bnd] = crunch2(G,Ae,Ai,x,y,z,s,rQ,rA,rC,rsz);

s=s+dels;
z=z+delz;
x=x+delx;
y=y+dely;

viol = -min(min(s),min(z));
shift = 1.e3 + 2*viol;
s=s+shift;
z=z+shift;

mu=z'*s/ni;

tic

% loop
for iter=1:30,
    
rQ = c + G*x - Ae'*y - Ai'*z;
rA = Ae * x - be;
rC = Ai*x - bi - s;

% predictor

rsz = s.*z;
[delx,dely,delz,dels,bnd] = crunch2(G,Ae,Ai,x,y,z,s,rQ,rA,rC,rsz);

alpha = bnd;
muaff = (z + alpha*delz)'*(s + alpha*dels) / ni;
sigma =(muaff/mu)^tsig;

% corrector

rsz = rsz+(dels.*delz)+(-sigma*mu);
[delx,dely,delz,dels,bnd] = crunch2(G,Ae,Ai,x,y,z,s,rQ,rA,rC,rsz);
alpha = 0.995*bnd;

s=s+alpha*dels;
x=x+alpha*delx;
y=y+alpha*dely;
z=z+alpha*delz;

mu=z'*s/ni;
iter
x(1:10)

end

toc

figure(3)
xxdraw=x;
hold off
for i=1:size(quads,1),
    tris(i*2-1,:)=quads(i,1:3);
    tris(i*2,:)=quads(i,[1 3 4]);
end
for i=1:size(xxdraw,1),
    px(i)=nodeToXY(i,1);
    py(i)=nodeToXY(i,2);
    pz(i)=xxdraw(i);
end
trimesh(tris,px,py,pz);

find(Ai*x<bi-1e-11~=0)
find(Ae*x-be>1e-5)
find(Ae*x-be<-1e-5)
find(x>up+1e-11~=0)
find(x<lo-1e-11~=0)
