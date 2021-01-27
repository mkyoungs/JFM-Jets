dt=1/128; % time step
tmax=200; % max time
tpl=1/dt; % time plot in number of time steps


dx=L/nx; % grid spacing
dy=W/ny; % grid spacing y

[xh,yh]=meshgrid([0:nx]*dx,[0:ny]*dy); % extended grid with end points

t=0; % time
tc=0; % count

q1_p=q1;q2_p=q2;


zy=zeros(1,nx);

psimax=[];

ts=[];
ubs=[];
qps=[];
while t<=tmax+dt/2
      
    % second order Adams-Bashforth
  tmp=q1;q1=1.5*q1-0.5*q1_p;q1_p=tmp; % extrapolation
  tmp=q2;q2=1.5*q2-0.5*q2_p;q2_p=tmp;
% p1 = psi1, p2 = psi2 total psi
  [p1,p2]=invert(q1,q2,a11,a12,a21,a22,b11,b12,b21,b22,massfac,Hfree,del);
  

  u1=-p(p1,dy); % differentiate
  

  v1=p(p1',dx)';

  
  u2=-p(p2,dy);
  v2=p(p2',dx)';
 

  if(rem(tc,tpl)==0) % plot t_0 time point
    ff=[rs(q2_p);zy;zy;rs(q1_p)];
    imagesc(ff);axis('xy');title(num2str(t));drawnow();
    ts=[ts,t];
    ubs=[ubs,mean(u1')'];
    qp=q1_p-mean(q1_p')'*ones(1,nx);
    qps=[qps,max(qp(:))];
    
    p1t(:,:,t+1) = 0.25*(p1(1:end-1,1:end-1)+ p1(2:end,1:end-1) + p1(2:end,2:end) + p1(1:end-1,2:end));
    u1t(:,:,t+1) = (u1(:,1:end-1) + u1(:,2:end))/2;
    v1t(:,:,t+1) = (v1(1:end-1,:) + v1(2:end,:))/2;
    v2t(:,:,t+1) = (v2(1:end-1,:) + v2(2:end,:))/2;
  end

  qt1=q1+by; % add background PV
  qt2=q2+by;

%   Fy=flux(q1,v1,dy,dt,0); % compute northward flux of total pv vq
%   Fy(1,:)=0; Fy(ny+1,:)=0;
  Fx=flux(q1',U1',dx,dt,1)'; % compute eastward flux of total pv uq
  dq1dt=-p(Fx',dx)'-beta1.*v1(2:end,:)-r*q1-r*F1*qforce; % PV equation including linear damping and forcing

%   Fy=flux(q2,v2,dy,dt,0);
%   Fy(1,:)=0; Fy(ny+1,:)=0;
  Fx=flux(q2',U2',dx,dt,1)';
  dq2dt=-p(Fx',dx)'-beta2.*v2(2:end,:)-r*q2+r*F2*qforce; % other layer
  
  dq2dt = dq2dt- mean(dq2dt,2);
  dq1dt = dq1dt- mean(dq1dt,2);

%  keyboard;
  q1=q1_p+dt*dq1dt; % give new PV
  q2=q2_p+dt*dq2dt;

  tc=tc+1; % time step
  t=tc*dt;

end

KEflux = -mean(u1t(:,:,end-100:end).*v1t(:,:,end-100:end),3).*U1y;
ketotal = sum(KEflux(:));
PEflux = F1.*(U1(:,1:end-1)-U2(:,1:end-1)).*mean(v2t(:,:,end-100:end).*p1t(:,:,end-100:end),3);
petotal = sum(PEflux(:));

figure;
plot(ts,max(ubs));
figure;
contourf(KEflux)
titK = sprintf('KE Flux %d',ketotal);
title(titK)
figure;
contourf(PEflux)
titP = sprintf('PE Flux %d',petotal);
title(titP)

figure;
plot(1:64,beta1(:,1),1:64,beta2(:,1))



function df=p(f,dy) % compute derivative in y
  df=(1/dy)*diff(f);
  return
end

function qp=rs(q) % rescaling
  del=max(max(q))-min(min(q));
  if(del==0)
    qp=q;
    return;
  end
  qp=(q-min(min(q)))/del;
end

% compute psi, q sep. and compute correct advection, remove zonal mean from
% dq/dt
function fa=flux(f,v,dy,dt,bcf) % compute fluxes with flux correction scheme (limit magnitude)
  if(bcf==0) 
    f=[-f(1,:);f;-f(end,:)];
  else
    f=[f(end,:);f;f(1,:)];
  end
  n=size(f,1);nm1=n-1;
  fbar=0.5*(f(1:nm1,:)+f(2:n,:));
  delf=f(2:n,:)-f(1:nm1,:);
  absv=abs(v);
  fup=v.*fbar-0.5*absv.*delf;
  flw=v.*fbar-0.5*dt/dy*absv.^2.*delf;
%  fbar=shift(f,1);
  delfp=circshift(delf,1);
  delfm=circshift(delf,-1);
  r=((v>=0).*delfp+(v<0).*delfm)./delf;
  if(sum(sum(isnan(r)))>0)
    r(isnan(r))=0;
  end
  if(sum(sum(isinf(r)))>0)
    r(isinf(r))=1e20*sign(r(isinf(r)));
  end
% van leer
  psi=(r+abs(r))./(1+abs(r));
% superbee
%  psi=max(0,max(min(1,2*r),min(2,r)));
%  if(sum(sum(isnan(psi)))>0)
%    psi(isnan(psi))=0;
%  endif
  fa=fup+psi.*(flw-fup);
%  fa=flw;
end


% compute PV inversion
function [p1,p2]=invert(q1,q2,a11,a12,a21,a22,b11,b12,b21,b22,massfac,Hfree,del) 
  nx=size(q1,2);
  ny=size(q1,1);
  z1b=mean(q1')'; % compute means
  z2b=mean(q2')'; 
  q1=q1-z1b*ones(1,nx);
  q2=q2-z2b*ones(1,nx); % interpolates to corners of box
  z1=(q1(1:ny-1,:)+q1(1:ny-1,[nx,1:nx-1])...
    +q1(2:ny,:)+q1(2:ny,[nx,1:nx-1]))/4;
  z2=(q2(1:ny-1,:)+q2(1:ny-1,[nx,1:nx-1])...
    +q2(2:ny,:)+q2(2:ny,[nx,1:nx-1]))/4;

  z1b=dct(z1b); % inversions for mean for cosine transform
  z2b=dct(z2b);
  p1b=idct(b11.*z1b+b12.*z2b);
  p2b=idct(b21.*z1b+b22.*z2b); % average back because of cosine
  p1b=[0.5*p1b(1)+0.5*p1b(2);0.5*(p1b(1:ny-1)+p1b(2:ny));0.5*p1b(ny-1)+0.5*p1b(ny)];
  p2b=[0.5*p2b(1)+0.5*p2b(2);0.5*(p2b(1:ny-1)+p2b(2:ny));0.5*p2b(ny-1)+0.5*p2b(ny)];
  delmass=massfac*(p1b-p2b); % conserve mass
  p1b=p1b-delmass*Hfree;
  p2b=p2b+del*delmass*Hfree;

  z1=fft(dst(z1)')'; % Fourier transform for perturbations
  z2=fft(dst(z2)')';
  p1=a11.*z1+a12.*z2;
  p2=a21.*z1+a22.*z2;
  p1=real([zeros(1,nx);...
       idst(ifft(p1')');...
       zeros(1,nx)]);
  p2=real([zeros(1,nx);...
       idst(ifft(p2')');...
       zeros(1,nx)]);
  p1=[p1,p1(:,1)]+p1b*ones(1,nx+1);
  p2=[p2,p2(:,1)]+p2b*ones(1,nx+1);
  return
end
