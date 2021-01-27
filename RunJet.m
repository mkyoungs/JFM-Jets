%% run NDJet
clear

newd = 0.05;
%Fs = 2:2:60;
Fs = [3:6:40]*(0.2^2/newd^2);
betas = [0:5:30]*(0.2^2/newd^2);
deltas = newd;%*[0.4 0.5 0.75 1 1.2];
reynolds = NaN(length(betas),length(Fs),length(deltas));
reynolds2 = reynolds;
pe = reynolds;
pe2 = reynolds;
%dudt = reynolds;
qymin = reynolds;
BTmin = reynolds;
psiratio = reynolds;
l = length(deltas);
l2 = length(Fs);

for i = 1:length(betas)
    disp(i)
    for k = 1:l2
        for j = 1:l
            [qymin(i,k,j),BTmin(i,k,j),~,~,pe(i,k,j),reynolds(i,k,j),pe2(i,k,j),reynolds2(i,k,j)] = NDJet(betas(i),deltas(j),Fs(k));
        end
 
    end
end

%cd ~/'Dropbox (MIT)'/Work/GFD/Output/
%save outputdelta2 qymin pe reynolds pe2 reynolds2 Fs betas deltas
%%

qymin = squeeze(qymin);
pe = squeeze(pe);
reynolds = squeeze(reynolds);
pe2 = squeeze(pe2);
reynolds2 = squeeze(reynolds2);
%copt = squeeze(kopt);
%bcgrowth = squeeze(bcgrowth);
BTmin = squeeze(BTmin);
psiratio = squeeze(psiratio);

[betamat, Fmat] = meshgrid(betas,Fs);


%     kideal = 0.55/deltas + deltas*(0.26*betamat + 0.7*Fmat);
%     bcgrowth = (Fmat.*kideal./(kideal.^2 + pi^2/4))';
% cutoff = 0;%btgrowth(1,1);

% h1 = figure; hold on
titlebase = sprintf('Background U0 delta %g',deltas);
% contourf(betas,deltas,dudt','linestyle','none')
% contour(betas,deltas,dudt',[0 0],'k','linewidth',2)
% contour(betas,deltas,qymin',[0 0],'w','linewidth',2)
% %plot(betas,1./sqrt(sqrt(2)*F).*(1-betas./sqrt(2)./F))
% %plot(1./sqrt(F - 1./(deltas.^2)),deltas)
% xlabel('Beta')
% ylabel('Delta')
% colorbar
% box on
% set(gca,'tickdir','out')
% titl = ['dudt 1 ',titlebase];
% title(titl)
% SaveFigureGFD(h1,titlebase,titl)

h3 = figure; 
h2a = axes; 
hold on
cmap = brewermap([],'RdBu');
colormap(cmap);
contourf(betas*deltas^2,Fs*deltas^2,pe',[-0.1:0.05:1.1],'linestyle','none'); shading flat
contour(betas*deltas^2,Fs*deltas^2,pe',[0 0],':k','linewidth',2)
contour(betas*deltas^2,Fs*deltas^2,BTmin',[0 0],'--g','linewidth',2)
contour(betas*deltas^2,Fs*deltas^2,reynolds',[0 0],'--k','linewidth',2)
contour(betas*deltas^2,Fs*deltas^2,qymin',[0 0],'k','linewidth',2)
%contour(betas,Fs,bcgrowth'-btgrowth',[cutoff cutoff],'m','linewidth',2)
caxis([-1.5 1.5])
text(0.05,1.45,'(a)','backgroundcolor','w','fontsize',14)
%plot(betas,betas/2)
%plot(betas,1./sqrt(sqrt(2)*F).*(1-betas./sqrt(2)./F))
%plot(1./sqrt(F - 1./(deltas.^2)),deltas)
%xlabel('\beta \delta^2','fontsize',24)
ylabel('F \delta^2','fontsize',24)
% hy1 = colorbar;
% set(hy1,'ylim',[-0.1 1.1])
box on
set(gca,'tickdir','out')
hy1 = colorbar;
set(hy1,'ylim',[-0.1 1.1])
ylabel(hy1,'Normalized Energy Conversion')
titl = ['PE Conversion Normalized ',titlebase];
%title('Normalized APE Conversion')
%SaveFigureGFD(h2,titlebase,titl)

h3a = axes; hold on
colormap(cmap)
contourf(betas*deltas^2,Fs*deltas^2,reynolds',[-0.1:0.05:1.1],'linestyle','none')
caxis([-1.5 1.5])
text(0.05,1.45,'(b)','backgroundcolor','w','fontsize',14)
% text(0.3,1.3,'1','fontsize',20,'backgroundcolor','w')
% text(0.2,0.9,'2','fontsize',20,'backgroundcolor','w')
% text(0.7,1,'3','fontsize',20,'backgroundcolor','w')
% text(0.5,0.7,'4','fontsize',20,'backgroundcolor','w')
% text(0.6,0.23,'5','fontsize',20,'backgroundcolor','w')
%text(0.9,0.2,'6','fontsize',14,'backgroundcolor','w')
contour(betas*deltas^2,Fs*deltas^2,reynolds',[0 0],'--k','linewidth',2)
contour(betas*deltas^2,Fs*deltas^2,pe',[0 0],':k','linewidth',2)
contour(betas*deltas^2,Fs*deltas^2,qymin',[0 0],'k','linewidth',2)
%contour(betas,Fs,bcgrowth'-btgrowth',[cutoff cutoff],'m','linewidth',2)
%plot(betas,1./sqrt(sqrt(2)*F).*(1-betas./sqrt(2)./F))
%plot(1./sqrt(F - 1./(deltas.^2)),deltas)
xlabel('\beta \delta^2','fontsize',24)
ylabel('F \delta^2','fontsize',24)
%ylabel('F')
hy1 = colorbar;
set(hy1,'ylim',[-0.1 1.1])
ylabel(hy1,'Normalized Energy Conversion')
box on
set(gca,'tickdir','out')
titl = ['KE Conversion Normalized ',titlebase];
%title('Normalized KE Conversion')

% h3b = axes; hold on
% plot(Fs*deltas^2,reynolds(2,:),'r','linewidth',2)
% text(0.05,0.95,'(c)','fontsize',14)
% Fint = interp1(qymin(2,:),Fs,0);
% plot([Fint Fint]*deltas^2,[0 1],'k','linewidth',2)
% set(gca,'tickdir','out')
% xlabel('F \delta^2','fontsize',24)
% ylabel('Normalized Energy Conversion','fontsize',24)
% xlim([0 1.6])

%plot(h3a,[betas(2) betas(2)]*deltas^2,[min(Fs) max(Fs)]*deltas^2,'r','linewidth',2)

set(h2a,'position',[0.11 0.72 0.72 0.28],'fontsize',16)
set(h3a,'position',[0.11 0.4 0.72 0.28],'fontsize',16)
%set(h3b,'position',[0.11 0.06 0.72 0.28],'fontsize',16)
set(h3,'position',[100 100 550 1000])
SaveFigureGFD(h3,titlebase,titl)

h4 = figure; hold on

contourf(betas,Fs,real(copt)','linestyle','none')
contour(betas,Fs,pe2',[0 0],'k','linewidth',2)
contour(betas,Fs,qymin',[0 0],'k','linewidth',2)
%contour(betas,Fs,bcgrowth'-btgrowth',[cutoff cutoff],'m','linewidth',2)
%plot(betas,1./sqrt(sqrt(2)*F).*(1-betas./sqrt(2)./F))
%plot(1./sqrt(F - 1./(deltas.^2)),deltas)
xlabel('\beta')
%ylabel('F')
colorbar
box on
set(gca,'tickdir','out')
titl = ['PE Conversion ',titlebase];
title(titl)
SaveFigureGFD(h4,titlebase,titl)

h5 = figure; hold on

contourf(betas,Fs,abs(copt)',40,'linestyle','none')
%caxis([-0.1 1.1])
contour(betas,Fs,reynolds2',[0 0],'k','linewidth',2)
contour(betas,Fs,qymin',[0 0],'w','linewidth',2)
%contour(betas,Fs,bcgrowth'-btgrowth',[cutoff cutoff],'m','linewidth',2)
%plot(betas,1./sqrt(sqrt(2)*F).*(1-betas./sqrt(2)./F))
%plot(1./sqrt(F - 1./(deltas.^2)),deltas)
xlabel('Beta')
ylabel('F')
colorbar
box on
set(gca,'tickdir','out')
titl = ['KE Conversion ',titlebase];
title(titl)
SaveFigureGFD(h5,titlebase,titl)