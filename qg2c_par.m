clear

f0=8.64; % f0
%beta=1.728e-3;
global L W nx ny a11 a12 a21 a22 b11 b12 b21 b22 massfac Hfree del F1 F2 rd
W=2; % width
L=2*W; % length
rd= 0.00001; % linear damping pv drag


beta = 8;
%Rd=60; % deformation radius (km, day units) gives dif
del=1; % H1/H2

nx=128;ny=64; % grid points/ fourier modes

k0x=2*pi/L; % lowest wavenumber in x
k0y=pi/W; % lowest wavenumber in y
[k,l]=meshgrid([0:nx/2,-nx/2+1:-1]*k0x,[1:ny-1]*k0y); % k,l wavenumber
l0=[0:ny-1]*k0y; % wavenumber set for mean flow

F = 10;
F1 = F*(1+del); % usual
F2 = del*F1;

Rd = sqrt(1/F);

wv2=(k.*k+l.*l); % K^2
det=wv2.*(wv2+F1+F2); % determinant of inversion matrix
a11=-(wv2+F2)./det; % components of inversion matrix
a12=-F1./det;
a21=-F2./det;
a22=-(wv2+F1)./det;

a11(:,1)=zeros(ny-1,1); % take care of singularities
a12(:,1)=zeros(ny-1,1);
a21(:,1)=zeros(ny-1,1);
a22(:,1)=zeros(ny-1,1);

wv20=(l0.*l0)'; % mean flow K^2
det=wv20.*(wv20+F1+F2); % same as above for mean flow
det(1)=1;
b11=-(wv20+F2)./det;
b12=-F1./det;
b21=-F2./det;
b22=-(wv20+F1)./det;

b11(1)=0;
b12(1)=0;
b21(1)=0;
b22(1)=0;

[x,y]=meshgrid([1/2:1:nx]/nx*L,[1/2:1:ny]/ny*W); % set up grid;
%dU = U0;
delta = 0.1;
U1=0*cos(2*pi*y/W);

%indi = find(U1(:,1) == min(U1(:,1)));

%U1(ny/2+1:end,:) = U1(ny/2:-1:1,:);
U2=U1*0.4;%(0)*exp(-(y-1).^2/(delta^2));

topo = f0*exp(-(x-1).^2/0.04);%f0*cos(2*pi*y/W);
%topo(1:ny/2,:) = topo(1:ny/2,:)/2;

U1(:,129) = U1(:,1);
U2(:,129) = U2(:,1);
U1y = diffxy(y,U1(:,1:end-1),1,1);
U2y = diffxy(y,U2(:,1:end-1),1,1);
U1yy = diffxy(y,U1(:,1:end-1),1,2);
U2yy = diffxy(y,U2(:,1:end-1),1,2);
% U1yy = -real(ifft2(l.^2.*fft2(U2)));
% U2yy = -real(ifft2(l.^2.*fft2(U2)));
% U1yy(isnan(U1yy)) = 0;
% U2yy(isnan(U2yy)) = 0;

beta1=beta+F1*(U1(:,1:end-1)-U2(:,1:end-1)) - U1yy;

%ind = interp1(beta1(1:ny/4,1),y(1:ny/4,1),0);

beta2=beta-F2*(U1(:,1:end-1)-U2(:,1:end-1)) - U2yy + topo;
%topo1 = interp1(y(:,1),beta2(:,1),ind);
%beta2 = beta2-topo;

plot(beta1(:,1),y(:,1),beta2(:,1),y(:,1),zeros*y(:,1),y(:,1))

% beta1=beta+F1.*(U1-U2);
% beta2=beta-F2.*(U1-U2);
by=beta*(y-W/2); % set up beta field
 %beta1y = beta1.*(y-W/2);
 %beta2y = beta2.*(y-W/2);

%amp=4e4;
amp = 10;
qforce=amp*cos(pi*y/W); % wind

y0=[0:ny]'*W/ny; % grid for mean flow
Hfree=cosh((y0-W/2)/Rd); % calculating zero pv mode with unit value on boundary
massfac=[1/2,ones(1,ny-1),1/2];
Hfree=Hfree/sum(massfac*Hfree)/(1+del)-topo/f0; % free solution correction, impose mass conservation

% remember to add topography back into the above solver!!!
dx=L/nx;
k0=2*pi/nx/dx;
if ~exist('q1')
    amp=0.0001;
    q1=0;
    q2=0;
    for ik=[1,2,5,8]
        for il=[-3,-1,0,1,3]
            q1 = q1 + amp*rand(1,1)*cos(ik*k0*x+il*k0*y+2*pi*rand(1,1));
            q2 = q2 + amp*rand(1,1)*cos(ik*k0*x+il*k0*y+2*pi*rand(1,1));
        end
    end
end;

% q1(30:45,30:45) = 1;
% q2(30:45,30:45) = 1;


qg2c_step(q1, q2, by, U1, U2, beta1, beta2,qforce)

return
