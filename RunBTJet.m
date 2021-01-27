%% run BTJet
clear


%Fs = 2:2:60;
Fs = [2:2:18];
betas = [0:3:15];
alphas = 0;
%deltas = newd;%*[0.4 0.5 0.75 1 1.2];
reynolds = NaN(length(betas),length(alphas),length(Fs));
reynolds2 = reynolds;
pe = reynolds;
pe2 = reynolds;
%dudt = reynolds;
qymin = reynolds;
btgrowth = reynolds;
psiratio = reynolds;
l = length(alphas);
l2 = length(Fs);

parfor i = 1:length(betas)
    tic
    for k = 1:l2
 
        for j = 1:l
            [qymin(i,j,k),btgrowth(i,j,k),psiratio(i,j,k),kopt(i,j,k),pe(i,j,k),reynolds(i,j,k),pe2(i,j,k),reynolds2(i,j,k)] = BTJet(betas(i),alphas(j),Fs(k));
        end
 
    end
    toc
end

% cd ~/'Dropbox (MIT)'/Work/GFD/Output/
% save outputalpha0.1 qymin pe reynolds pe2 reynolds2 Fs betas deltas
%%

qymin = squeeze(qymin);
pe = squeeze(pe);
reynolds = squeeze(reynolds);
pe2 = squeeze(pe2);
reynolds2 = squeeze(reynolds2);
%bcgrowth = squeeze(bcgrowth);
btgrowth = squeeze(btgrowth);
psiratio = squeeze(psiratio);

[betamat, Fmat] = meshgrid(betas,Fs);


%     kideal = 0.55/deltas + deltas*(0.26*betamat + 0.7*Fmat);
%     bcgrowth = (Fmat.*kideal./(kideal.^2 + pi^2/4))';
% cutoff = 0;%btgrowth(1,1);

% h1 = figure; hold on
titlebase = sprintf('BT alpha %g',alphas);
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
cmap = brewermap(40,'RdBu');
colormap(cmap);
contourf(betas,Fs,pe',[-0.1:0.05:1.1],'linestyle','none'); shading flat
contour(betas,Fs,pe',[0 0],':k','linewidth',2)
%contour(betas,Fs,reynolds',[0 0],'--k','linewidth',2)
contour(betas,Fs,qymin',[0 0],'k','linewidth',2)
%contour(betas,Fs,bcgrowth'-btgrowth',[cutoff cutoff],'m','linewidth',2)
caxis([-1.5 1.5])
text(0.5,17,'(a)','backgroundcolor','w','fontsize',14)
text(2.5,4,'1','fontsize',14,'backgroundcolor','w')
text(1,7,'2','fontsize',14,'backgroundcolor','w')
text(5,12,'3','fontsize',14,'backgroundcolor','w')
%plot(betas,betas/2)
%plot(betas,1./sqrt(sqrt(2)*F).*(1-betas./sqrt(2)./F))
%plot(1./sqrt(F - 1./(deltas.^2)),deltas)
xlabel('\beta')
ylabel('F')
% hy1 = colorbar;
% set(hy1,'ylim',[-0.1 1.1])
box on
set(gca,'tickdir','out')
titl = ['PE Conversion Normalized ',titlebase];
%title('Normalized APE Conversion')
%SaveFigureGFD(h2,titlebase,titl)

h3a = axes; hold on
colormap(cmap)
contourf(betas,Fs,reynolds',[-0.1:0.05:1.1],'linestyle','none')
caxis([-1.5 1.5])
text(0.5,17,'(b)','backgroundcolor','w','fontsize',14)


% text(0.5,0.7,'4','fontsize',14,'backgroundcolor','w')
% text(0.6,0.23,'5','fontsize',14,'backgroundcolor','w')
%text(0.9,0.2,'6','fontsize',14,'backgroundcolor','w')
%contour(betas,Fs,reynolds',[0 0],'--k','linewidth',2)
contour(betas,Fs,pe',[0 0],':k','linewidth',2)
contour(betas,Fs,qymin',[0 0],'k','linewidth',2)
%contour(betas,Fs,bcgrowth'-btgrowth',[cutoff cutoff],'m','linewidth',2)
%plot(betas,1./sqrt(sqrt(2)*F).*(1-betas./sqrt(2)./F))
%plot(1./sqrt(F - 1./(deltas.^2)),deltas)
xlabel('\beta')
%ylabel('F')
hy1 = colorbar;
set(hy1,'ylim',[-0.1 1.1])
ylabel(hy1,'Normalized Energy Conversion')
box on
set(gca,'tickdir','out')
titl = ['KE Conversion Normalized ',titlebase];
%title('Normalized KE Conversion')

set(h2a,'position',[0.06 0.11 0.4 0.85],'fontsize',14)
set(h3a,'position',[0.5 0.11 0.4 0.85],'fontsize',14)
set(h3,'position',[100 100 1000 450])
SaveFigureGFD(h3,titlebase,titl)
h4 = figure; hold on

contourf(betas,Fs,pe2',40,'linestyle','none')
contour(betas,Fs,pe2',[0 0],'k','linewidth',2)
contour(betas,Fs,qymin',[0 0],'w','linewidth',2)
%contour(betas,Fs,bcgrowth'-btgrowth',[cutoff cutoff],'m','linewidth',2)
%plot(betas,1./sqrt(sqrt(2)*F).*(1-betas./sqrt(2)./F))
%plot(1./sqrt(F - 1./(deltas.^2)),deltas)
xlabel('Beta')
ylabel('F')
colorbar
box on
set(gca,'tickdir','out')
titl = ['PE Conversion ',titlebase];
title(titl)
SaveFigureGFD(h4,titlebase,titl)

h5 = figure; hold on

contourf(betas,Fs,reynolds2',40,'linestyle','none')
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