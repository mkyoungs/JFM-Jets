%% BT Growth Rate

% region 1 beta 1 F 4
[~,~,~,~,~,~,~,~,ketotal1,petotal1,kvec1,cvec1] = BTJet(1,0,4);

% region 2 beta 1 F 7
[~,~,~,~,~,~,~,~,ketotal2,petotal2,kvec2,cvec2] = BTJet(1,0,7);

% region 3 beta 1 F 10
[~,~,~,~,~,~,~,~,ketotal3,petotal3,kvec3,cvec3] = BTJet(1,0,10);


%%
h1 = figure(34);
cmap = [0.6 0 0.9; 0.4 0 0.9; 0 0 0.9;0 0.4 0.6; 0 0.5 0];
h3 = axes; hold on;
scatter(kvec1,kvec1.*imag(cvec1),50,cmap(1,:),'filled','o')
scatter(kvec2,kvec2.*imag(cvec2),50,cmap(2,:),'filled','s')
scatter(kvec3,kvec3.*imag(cvec3),50,cmap(3,:),'filled','d')
%  scatter(kvec4,kvec4.*imag(cvec4),20,cmap(4,:),'filled')
%  scatter(kvec5,kvec5.*imag(cvec5),20,cmap(5,:),'filled')
%scatter(kvec6,kvec6.*imag(cvec6),20,cmap(6,:),'filled')
text(0.2,0.37,'(a)','backgroundcolor','w','fontsize',14)
xlabel('Zonal Wavenumber k')
grid on
box on
ylabel('Growth Rate (kc_i)')
legend('Region 1','Region 2','Region 3','Region 4','Region 5')

h4 = axes; hold on;
text(0.2,1.1,'(b)','backgroundcolor','w','fontsize',14)
scatter(kvec1,ketotal1,50,cmap(1,:),'filled','o')
scatter(kvec2,ketotal2,50,cmap(2,:),'filled','s')
scatter(kvec3,ketotal3,50,cmap(3,:),'filled','d')
%  scatter(kvec4,ketotal4,20,cmap(4,:),'filled')
%  scatter(kvec5,ketotal5,20,cmap(5,:),'filled')
%scatter(kvec6,ketotal6,20,cmap(6,:),'filled')
xlabel('Zonal Wavenumber k')
grid on
box on
ylabel('KE Conversion')
ylim([-0.1 1.2])

h6 = axes; hold on
scatter(kvec1,ketotal1,50,cmap(1,:),'filled','o')
scatter(kvec2,ketotal2,50,cmap(2,:),'filled','s')
scatter(kvec3,ketotal3,50,cmap(3,:),'filled','d')
% scatter(kvec4,ketotal4,20,cmap(4,:),'filled')
% scatter(kvec5,ketotal5,20,cmap(5,:),'filled')
%scatter(kvec6,ketotal6,20,cmap(6,:),'filled')
xlabel('Zonal Wavenumber k')
grid on
box on
text(0.2,0.085,'(d)','backgroundcolor','w','fontsize',14)
ylabel('KE Conversion')
ylim([-0.1 0.1])
xlim([0 6])

h5 = axes; hold on;
text(0.2,1.1,'(c)','backgroundcolor','w','fontsize',14)
scatter(kvec1,petotal1,50,cmap(1,:),'filled','o')
scatter(kvec2,petotal2,50,cmap(2,:),'filled','s')
scatter(kvec3,petotal3,50,cmap(3,:),'filled','d')
%  scatter(kvec4,petotal4,20,cmap(4,:),'filled')
%  scatter(kvec5,petotal5,20,cmap(5,:),'filled')
%scatter(kvec6,petotal6,20,cmap(6,:),'filled')
xlabel('Zonal Wavenumber k')
grid on
box on
ylim([-0.1 1.2])
ylabel('PE Conversion')

h7 = axes; hold on;
text(0.2,0.085,'(e)','backgroundcolor','w','fontsize',14)
scatter(kvec1,petotal1,50,cmap(1,:),'filled','o')
scatter(kvec2,petotal2,50,cmap(2,:),'filled','s')
scatter(kvec3,petotal3,50,cmap(3,:),'filled','d')
% scatter(kvec4,petotal4,20,cmap(4,:),'filled')
% scatter(kvec5,petotal5,20,cmap(5,:),'filled')
%scatter(kvec6,petotal6,20,cmap(6,:),'filled')
xlabel('Zonal Wavenumber k')
grid on
box on
ylabel('PE Conversion')
ylim([-0.1 0.1])
%xlim([0 70])

set(h1,'position',[100 100 1500 500],'paperpositionmode','auto')
set(h3,'position',[0.04 0.1 0.29 0.85],'fontsize',14)
set(h4,'position',[0.375 0.4 0.29 0.55],'fontsize',14)
set(h6,'position',[0.37 0.1 0.29 0.25],'fontsize',14)
set(h5,'position',[0.7 0.4 0.29 0.55],'fontsize',14)
set(h7,'position',[0.7 0.1 0.29 0.25],'fontsize',14)

print -depsc2 ~/'Dropbox (MIT)'/Work/GFD/BTGrowth

