%% create growth rate figures for paper

% region 1 beta 120 F 520
[~,~,~,~,~,~,~,~,ketotal1,petotal1,kvec1,cvec1] = NDJet(120,0.05,520);

% region 2 beta 80 F 360
[~,~,~,~,~,~,~,~,ketotal2,petotal2,kvec2,cvec2,psi2,xgrid,ygrid] = NDJet(80,0.05,360);

% region 3 beta 280 F 400
[~,~,~,~,~,~,~,~,ketotal3,petotal3,kvec3,cvec3,psi3] = NDJet(280,0.05,400);

% region 4 beta 200 F 280
[~,~,~,~,~,~,~,~,ketotal4,petotal4,kvec4,cvec4] = NDJet(200,0.05,280);

% region 5 beta 240 F 200
[~,~,~,~,peoput1,~,~,~,ketotal5,petotal5,kvec5,cvec5] = NDJet(240,0.05,200);

% region 6 beta 360 F 80
%[~,~,~,~,peoput,~,~,~,ketotal6,petotal6,kvec6,cvec6] = NDJet(20,0.05,80);

% F 24 beta 23 delta 50/200
%%
h1 = figure(34);
cmap = [0.6 0 0.9; 0.4 0 0.9; 0 0 0.9;0 0.4 0.6; 0 0.5 0];
%cmap = brewermap(6,'PRGn');
h3 = axes; hold on;
scatter(kvec1(1:2:end),kvec1(1:2:end).*imag(cvec1(1:2:end)),50,cmap(1,:),'filled','o')
scatter(kvec2(1:2:end),kvec2(1:2:end).*imag(cvec2(1:2:end)),50,cmap(2,:),'filled','>')
scatter(kvec3(1:2:end),kvec3(1:2:end).*imag(cvec3(1:2:end)),50,cmap(3,:),'filled','d')
 scatter(kvec4(1:2:end),kvec4(1:2:end).*imag(cvec4(1:2:end)),50,cmap(4,:),'filled','^')
 scatter(kvec5(1:2:end),kvec5(1:2:end).*imag(cvec5(1:2:end)),50,cmap(5,:),'filled','s')
%scatter(kvec6,kvec6.*imag(cvec6),20,cmap(6,:),'filled')
text(2,5.8,'(a)','backgroundcolor','w','fontsize',14)
xlabel('Zonal Wavenumber k')
grid on
box on
ylabel('Growth Rate (kc_i)')
legend('Region 1','Region 2','Region 3','Region 4','Region 5')

h4 = axes; hold on;
text(2,1.05,'(b)','backgroundcolor','w','fontsize',14)
scatter(kvec1(1:2:end),ketotal1(1:2:end),50,cmap(1,:),'filled','o')
scatter(kvec2(1:2:end),ketotal2(1:2:end),50,cmap(2,:),'filled','>')
scatter(kvec3(1:2:end),ketotal3(1:2:end),50,cmap(3,:),'filled','d')
 scatter(kvec4(1:2:end),ketotal4(1:2:end),50,cmap(4,:),'filled','^')
 scatter(kvec5(1:2:end),ketotal5(1:2:end),50,cmap(5,:),'filled','s')
%scatter(kvec6,ketotal6,20,cmap(6,:),'filled')
xlabel('Zonal Wavenumber k')
grid on
box on
ylabel('KE Conversion')
ylim([-0.1 1.1])

h6 = axes; hold on
text(2,0.04,'(d)','backgroundcolor','w','fontsize',14)
scatter(kvec1(1:2:end),ketotal1(1:2:end),50,cmap(1,:),'filled','o')
scatter(kvec2(1:2:end),ketotal2(1:2:end),50,cmap(2,:),'filled','>')
scatter(kvec3(1:2:end),ketotal3(1:2:end),50,cmap(3,:),'filled','d')
scatter(kvec4(1:2:end),ketotal4(1:2:end),50,cmap(4,:),'filled','^')
scatter(kvec5(1:2:end),ketotal5(1:2:end),50,cmap(5,:),'filled','s')
%scatter(kvec6,ketotal6,20,cmap(6,:),'filled')
xlabel('Zonal Wavenumber k')
grid on
box on
ylabel('KE Conversion')
ylim([-0.05 0.05])
xlim([0 70])

h5 = axes; hold on;
text(2,1.05,'(c)','backgroundcolor','w','fontsize',14)
scatter(kvec1(1:2:end),petotal1(1:2:end),50,cmap(1,:),'filled','o')
scatter(kvec2(1:2:end),petotal2(1:2:end),50,cmap(2,:),'filled','>')
scatter(kvec3(1:2:end),petotal3(1:2:end),50,cmap(3,:),'filled','d')
 scatter(kvec4(1:2:end),petotal4(1:2:end),50,cmap(4,:),'filled','^')
 scatter(kvec5(1:2:end),petotal5(1:2:end),50,cmap(5,:),'filled','s')
%scatter(kvec6,petotal6,20,cmap(6,:),'filled')
xlabel('Zonal Wavenumber k')
grid on
box on
ylim([-0.1 1.1])
ylabel('PE Conversion')

h7 = axes; hold on;
text(2,0.04,'(e)','backgroundcolor','w','fontsize',14)
scatter(kvec1(1:2:end),petotal1(1:2:end),50,cmap(1,:),'filled','o')
scatter(kvec2(1:2:end),petotal2(1:2:end),50,cmap(2,:),'filled','>')
scatter(kvec3(1:2:end),petotal3(1:2:end),50,cmap(3,:),'filled','d')
scatter(kvec4(1:2:end),petotal4(1:2:end),50,cmap(4,:),'filled','^')
scatter(kvec5(1:2:end),petotal5(1:2:end),50,cmap(5,:),'filled','s')
%scatter(kvec6,petotal6,20,cmap(6,:),'filled')
xlabel('Zonal Wavenumber k')
grid on
box on
ylabel('PE Conversion')
ylim([-0.05 0.05])
xlim([0 70])

set(h1,'position',[100 100 1500 500],'paperpositionmode','auto')
set(h3,'position',[0.04 0.1 0.29 0.85],'fontsize',14)
set(h4,'position',[0.375 0.4 0.29 0.55],'fontsize',14)
set(h6,'position',[0.37 0.1 0.29 0.25],'fontsize',14)
set(h5,'position',[0.7 0.4 0.29 0.55],'fontsize',14)
set(h7,'position',[0.7 0.1 0.29 0.25],'fontsize',14)

print -depsc2 ~/'Dropbox (MIT)'/Work/GFD/Growth

%%

[~,~,~,~,~,~,~,~,ketotal2,petotal2,kvec2,cvec2,psi2,xgrid,ygrid] = NDJet(10,0.05,360);

% region 3 beta 280 F 400
[~,~,~,~,~,~,~,~,ketotal3,petotal3,kvec3,cvec3,psi3] = NDJet(320,0.05,400);

h9 = figure;

h90 = axes;

plot(1/2 + 1/2*exp(-ygrid(:,1).^2/0.05^2),ygrid(:,1))
xlim([0.4 1])
ylim([-0.4 0.4])
xlabel('U')
ylabel('Y')
text(0.45, 0.35,'(a)','backgroundcolor','w','fontsize',15)


h91 = axes; hold on
cmap2 = flipud(brewermap([],'RdBu'));
colormap(cmap2)
contourf(xgrid',ygrid',psi2',[-0.3:0.025:0.3],'linestyle','none')
contour(xgrid',ygrid',psi2',[0 0],'k','linewidth',1)
xlim([0 0.5])
ylim([-0.4 0.4])
xlabel('X')
%ylabel('Y')
caxis([-0.4 0.4])
%title('Down-Gradient Momentum Fluxes')
text(0.05,0.36,'(b)','backgroundcolor','w','fontsize',15)

h92 = axes; hold on
contourf(xgrid',ygrid',psi3')
contourf(xgrid',ygrid',psi3',[-0.3:0.025:0.3],'linestyle','none')
contour(xgrid',ygrid',psi3',[0 0],'k','linewidth',1)
xlim([0 0.5])
ylim([-0.4 0.4])
xlabel('X')
h1 = colorbar;
%title('Up-Gradient Momentum Fluxes')
ylabel(h1,'Stream Function')
caxis([-0.4 0.4])
text(0.05,0.36,'(c)','backgroundcolor','w','fontsize',15)


set(h9,'position',[100 100 1000 500],'paperpositionmode','auto')
set(h90,'position',[0.06 0.11 0.155 0.85],'fontsize',16)
set(h91,'position',[0.26 0.11 0.3 0.85],'fontsize',16)
set(h92,'position',[0.6 0.11 0.3 0.85],'fontsize',16)

print -depsc2 ~/'Dropbox (MIT)'/Work/GFD/Phi

