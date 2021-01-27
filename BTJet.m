function [qymin, u1max, psiratio, kideal,petotout, ketotout,petotout2,ketotout2,ketotal,petotal,kvec,cvec] = BTJet(beta,alpha,F)
% NDJet function to compute various linear stability parameters for Gaussian jets
% beta0 = beta delta^2/U_1
% alpha = ratio top and bottom velocities
% F = delta^2/Rd^2
% del = ratio of top layer depth to bottom layer depth
% gamma = magnitude of topography
del = 1;

% set up grid
% domain is -5delta to 5delta in width and 20 in length with jet width 1
W = 2;
L = 24*2;
dy=W/512*4;
ny=256/4;
ny=2*ny;dy=dy/2;
y=[1/2:1:ny]'/ny*W-W/2;
[xgrid, ygrid] = meshgrid(dy:dy:L,y);



% set F in each layer
F1 = F*(del + 1)/del;
F2 = F1*del;

% set up velocity Gaussian jet with half-width delta
U1=1/2 + 1/2*cos(pi*y);
U2=alpha*U1;

topo = 8*cos(pi*y);





% U1mat =  1/2 + 1/2*exp(-(ygrid).^2);
% U2mat = 0*exp(-(ygrid).^2);

% compute derivatives of velocity in matrix
U1y = diffxy(y,U1,1,1);
U2y = diffxy(y,U2,1,1);


% set up differentiation matrix for solver
ddy2=[-1,zeros(1,ny-1);eye(ny);zeros(1,ny-1),-1];
ddy2=diff(diff(ddy2));
ddy2 = ddy2/dy^2;

% set up differentiation matrix for pv derivative, different boundary
% conditions
ddy2n=[-1,zeros(1,ny-1);eye(ny);zeros(1,ny-1),-1];
ddy2n=diff(diff(ddy2n));
ddy2n(1,1) = 2; ddy2n(1,2) = -5; ddy2n(1,3) = 4; ddy2n(1,4) = -1;
ddy2n(end,end) = 2; ddy2n(end,end-1) = -5; ddy2n(end,end-2) = 4; ddy2n(end,end-3) = -1;
ddy2n = ddy2n/dy^2;

% compute background PV gradient
Q1y=beta-ddy2n*U1+F1*(U1-U2);
Q2y=beta-ddy2n*U2+F2*(U2-U1)+topo;

ind = interp1(Q1y(1:ny/2,1),y(1:ny/2,1),0);

topo1 = interp1(y(:,1),Q2y(:,1),ind);

Q2y = Q2y-topo1;

u1max = max(U1y);

qymin = Q2y(ny/2);

% plot background PV gradient
plot(Q1y,y,Q2y,y,y*0,y);
drawnow;

h12 = figure;
h112 = axes; hold on;
plot(U1,y,U2,y,'linewidth',2)
text(-0.07,0.9,'(a)','fontsize',16)
xlim([-0.1 1])
xlabel('Velocity')
ylabel('X')
box on
grid on
legend('Upper Layer','Lower Layer')
%title('Velocity')

h122 = axes; hold on
plot(Q1y,y,Q2y,y,'linewidth',2);
text(-18,0.9,'(b)','fontsize',16,'backgroundcolor','w')
grid on
xlabel('Q_y')
box on
%legend('Upper Layer','Lower Layer')
%title('Potential Vorticity Gradient')



set(h12,'position',[100 100 1000 500],'paperpositionmode','auto')
set(h112,'position',[0.06 0.11 0.44 0.85],'fontsize',16)
set(h122,'position',[0.54 0.11 0.44 0.85],'fontsize',16)

print -depsc2 ~/'Dropbox (MIT)'/Work/GFD/SetUpBT


% set up matrices for solver
M0=[ddy2-F1*eye(ny),F1*eye(ny);F2*eye(ny),ddy2-F2*eye(ny)];
DDY2 = [ddy2,zeros(ny);zeros(ny),ddy2];
U=[diag(U1),zeros(ny);zeros(ny),diag(U2)];
Qy=[diag(Q1y),zeros(ny);zeros(ny),diag(Q2y)];

% initiate ks
ks=[2:2:64*6]*pi/(L);

% initiate vectors for parameters

eigenvec = [];
kvec = [];
cvec = [];
petotal = [];
ketotal = [];
petotal2 = [];
ketotal2 = [];
for k=ks
    k; %,fflush(1);
    M=U+Qy/(M0-k*k*eye(2*ny));
    if sum(isnan(M(:))) == 0 && sum(isinf(M(:))) == 0
        [V,D] = eig(M);
        %vec = zeros(2*ny,1);
        c = diag(D);
        if (sum(imag(c) ~= 0)) ~= 0
            [~, ind] = sort(imag(c),'descend');
            for i = 1:min(sum(imag(c)>0),1)
                vec = V(:,ind(i));
                
                
                psis = (M0-k*k*eye(2*ny))\vec;
                psi1 = real(psis(1:ny)*exp(1i*k*(dy:dy:L)));
                psi2 = real(psis(ny+1:end)*exp(1i*k*(dy:dy:L)));
                
                v1 = diffxy(xgrid,psi1,2);
                v2 = diffxy(xgrid,psi2,2);
                u1 = -diffxy(ygrid,psi1,1);
                u2 = -diffxy(ygrid,psi2,1);
                
                totalenergy = 1/4*sum(v1.^2 + u1.^2,2)*dy/L + ...
                    1/4*sum(v2.^2 + u2.^2,2)*dy/L + F*sum((psi1-psi2).^2,2)*dy/L/2;
                totalenergy = sum(totalenergy)*dy;
                
                vec = vec/sqrt(totalenergy);
                
                psis = (M0-k*k*eye(2*ny))\vec;
                psi1 = real(psis(1:ny)*exp(1i*k*(dy:dy:L)));
                psi2 = real(psis(ny+1:end)*exp(1i*k*(dy:dy:L)));
                
                
                
                
                % buoyancy/ mass fluxes and KE fluxes
                
                KEflux = 1/2*sum(diffxy(xgrid,psi1,2).*diffxy(ygrid,psi1,1),2)*dy/L.*U1y+...
                    1/2*sum(diffxy(xgrid,psi2,2).*diffxy(ygrid,psi2,1),2)*dy/L.*U2y;
                
                PEflux = F.*(U1-U2).*sum(diffxy(xgrid,psi2,2).*psi1,2)*dy/L;
                
                
                ketotal = [ketotal, sum(KEflux)/(sum(KEflux)+sum(PEflux))];
                petotal = [petotal,sum(PEflux)/(sum(KEflux)+sum(PEflux))];
                ketotal2 = [ketotal2, sum(KEflux)];
                petotal2 = [petotal2, sum(PEflux)];
                
                eigenvec = [eigenvec,vec];
                kvec = [kvec, k];
                cvec = [cvec, c(ind(i))];
            end
        end
    end
    %  c = eig(M);
    
    
end
%

if ~isempty(cvec) && length(cvec)> 3
    % plot things
    titlebase = sprintf('BT F %g beta %g alpha %g',F,beta,alpha);
    
    kmax = max(ks);
    
%     h1 = figure;
%     scatter(kvec,real(cvec),20,kvec,'filled')
%     xlabel('Zonal Wavenumber k')
%     grid on
%     xlim([0 kmax])
%     ylabel('Real Phase Speed (c_r)')
%     tit = ['Real Phase Speed ',titlebase];
%     title(tit)
%     SaveFigureGFD(h1,titlebase,tit);
%     
%     h2 = figure;
%     scatter(kvec,imag(cvec),20,kvec,'filled')
%     xlabel('Zonal Wavenumber k')
%     xlim([0 kmax])
%     grid on
%     ylabel('Imaginary Phase Speed (c_i)')
%     tit = ['Imaginary Phase Speed ',titlebase];
%     title(tit)
%     SaveFigureGFD(h2,titlebase,tit);
%     
%     h3 = figure;
%     scatter(kvec,kvec.*imag(cvec),20,kvec,'filled')
%     xlabel('Zonal Wavenumber k')
%     xlim([0 kmax])
%     grid on
%     ylabel('Growth Rate (kc_i)')
%     tit = ['Growth Rate ',titlebase];
%     title(tit)
%     SaveFigureGFD(h3,titlebase,tit);
%     
%     h4 = figure;
%     scatter(kvec,ketotal,20,kvec,'filled')
%     xlabel('Zonal Wavenumber k')
%     xlim([0 kmax])
%     grid on
%     totalgrowthKE = sum(ketotal);
%     ylabel('KE Conversion')
%     tit = [sprintf('KE Conversion =%g    ',totalgrowthKE),titlebase];
%     title(tit)
%     SaveFigureGFD(h4,titlebase,tit);
%     
%     h5 = figure;
%     scatter(kvec,petotal,20,kvec,'filled')
%     xlabel('Zonal Wavenumber k')
%     totalgrowth = sum(petotal);
%     xlim([0 kmax])
%     grid on
%     ylabel('PE Conversion')
%     tit = [sprintf('PE Conversion =%g   ',totalgrowth),titlebase];
%     title(tit)
%     SaveFigureGFD(h5,titlebase,tit);
%     
%     h6 = figure;
%     plot(Q1y,y,Q2y,y,y*0,y);
%     xlabel('Potential Vorticity Gradient')
%     ylabel('Y')
%     legend('Q_y^1','Q_y^2')
%     tit = ['Potential Vorticity Gradient  ', titlebase];
%     title(tit)
%     SaveFigureGFD(h6,titlebase,tit)
    
    [~,locs] = findpeaks(kvec.*imag(cvec),'SortStr','descend');
    %sigs=[sigs;sig];
    
    if ~isempty(locs)
        loc = locs(1);
        psis = (M0-kvec(loc)*kvec(loc)*eye(2*ny))\eigenvec(:,loc);
        psi1 = real(psis(1:ny)*exp(1i*kvec(loc)*(dy:dy:L)));
        psi2 = real(psis(ny+1:end)*exp(1i*kvec(loc)*(dy:dy:L)));
        %         psistot = real(psis*exp(1i*kvec(loc)*(dy:dy:L)));
        %
        %         qs = real(eigenvec(:,loc)*exp(1i*kvec(loc)*(dy:dy:L)));
        psiratio = sum(sum(abs(psi2),2)*dy^2)/sum(sum(abs(psi1),2)*dy^2);
        %     q1 = real(eigenvec(1:ny,locs(1))*exp(1i*kvec(locs(1))*(dy:dy:L)));
        %     q2 = real(eigenvec(ny+1:end,locs(1))*exp(1i*kvec(locs(1))*(dy:dy:L)));
        v1 = diffxy(xgrid,psi1,2);
        v2 = diffxy(xgrid,psi2,2);
        u1 = -diffxy(ygrid,psi1,1);
        u2 = -diffxy(ygrid,psi2,1);
        
        % compute ageostrophic component of velocities
        
        
        kideal = kvec(locs(1));
        
        if length(locs)>1
            kideal = min(kvec(locs(1)),kvec(locs(2)));
        end
        %         alpha = 0.03;
        %         kideal = 1/delta*(1+alpha*sqrt(beta));
        %        bcgrowth = F1*kvec(locs(1))/(kvec(locs(1))^2 + pi^2/4);
        
        petotout = petotal(locs(1));
        ketotout = ketotal(locs(1));
        petotout2 = petotal2(locs(1));
        ketotout2 = ketotal2(locs(1));
        
        v1ave = F1*sum(v1.*(psi2-psi1),2)*dy/L;
        reynolds1 = -diffxy(y,sum(u1.*v1,2)*dy/L,1);
        du1dt = v1ave + reynolds1;
        
        v2ave = -F2*sum(v2.*(psi2-psi1),2)*dy/L;
        reynolds2 = -diffxy(y,sum(u2.*v2,2)*dy/L,1);
        du2dt = v2ave + reynolds2;
        
        % where does the factor of 4 come from?
        perturbation = 4*[diffxy(y,du1dt);diffxy(y,du2dt)];
        dpsidt = M0\perturbation;
        
        dudt1 = diffxy(y,dpsidt(1:ny));
        dudt2 = diffxy(y,dpsidt(ny+1:end));
        
        %    dudt = sum(du1dt)*dy + sum(du2dt)*dy;
        %    dudt = du1dt(ny/2) + du2dt(ny/2);
%         
%         h71 = figure;
%         plot(dudt1-reynolds1,y,'.g',reynolds1,y,'--g',dudt1,y,'-g',v1ave,'og',...
%             dudt2-reynolds2,y,'.m',reynolds2,y,'--m',dudt2,y,'-m',v2ave,'om',dudt1+dudt2,y,'k')
%         grid on
%         legend('Layer 1 - v^a','Layer 1 - uv','Layer 1 - total','Layer 1 - vq',...
%             'Layer 2 - v^a','Layer 2 - uv', 'Layer 2 - total','Layer 2 - vq','Sum')
%         xlabel('Eddy PV flux');
%         ylabel('y')
%         titl = ['Eddy PV flux  ',titlebase];
%         title(titl)
%         SaveFigureGFD(h71,titlebase,titl)
%         
        %         h72 = figure;
        %         plot(imag(psis(1:ny)),y,'g--',real(psis(1:ny)),y,'g.',...
        %             imag(psis(ny+1:end)),y,'m--', real(psis(ny+1:end)),y,'m.')
        %         legend('Imag - Layer 1','Real - Layer 1','Imag - Layer 2','Real - Layer 2')
        %         grid on
        %         ylabel('y')
        %         xlabel('psi')
        %         titl = ['Streamfunction  ',titlebase];
        %         title(titl)
        %         SaveFigureGFD(h72,titlebase,titl)
        % % %
        
        
        
        
        
        %                 h7 = figure;
        %                 plot(dudt,y,y*0,y,dudt1,y,dudt2,y);
        %                 xlabel('sum dudt')
        %                 ylabel('Y')
        %                 tit = [sprintf('Mean Velocity Change %g ',usum), titlebase];
        %                 title(tit)
        %                 SaveFigureGFD(h7,titlebase,tit)
        %
        %                 figure;
        %                 plot(-e1,y,-e2,y)
       % close all
    else
        
        petotout = NaN;
        ketotout = NaN;
        petotout2 = NaN;
        ketotout2 = NaN;
        bcgrowth = NaN;
        kideal = NaN;
        psiratio = NaN;
   
        
    end
else
    petotout = NaN;
    ketotout = NaN;
    petotout2 = NaN;
    ketotout2 = NaN;
    bcgrowth = NaN;
    kideal = NaN;
    psiratio = NaN;
end

end




