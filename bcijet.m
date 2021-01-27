clear
close all
sigs=[];

delta = 50; % width of jet
beta0 = 50; % beta*L^2/U

% why is there no stability?
%F = 150;


deltas = [200:-20:120,115:-5:25];
%deltas = [20:20:60];

for Find = 1%:8
% for F = 100
    for beta0 = 10%:10:50
        usumtotal = [];
        kemaxgrowth = [];
        minqy = [];
        for delta = 2000000000000%20:20:200
            disp(delta)
            F = 20*Find;
            % parameters to change
            % F*L^2 = L^2/Rd^2*(Dtotal/Dn), for single layer
            
            % set layer depths equal
            del = 1; % H1/H2
            
            %Rd=40;
            Rd = sqrt(500^2/F*2); % deformation radius for total depth (km, day units)
            f0=8.64;
            beta=1.728e-3;
            L=1000;
            %beta=0;
            
            %del=0.2;
            %del=0;
            
            dy=L/256;
            ny=128;
            ny=2*ny;dy=dy/2;
            W=dy*ny;
            
            r = 0.0;
            %U0 = beta/beta0*W*W;
            
            %dU = U0;
            
            %
            % U0 = 15;
            % u2f = 0.2;
            
            F0 = 1/Rd^2;
            F1=1/Rd^2*(1+del);
            F2=del*F1; % make sure to  double checkchange if I change del
            
            U0 = 2*beta/F2;
            
            y=[1/2:1:ny]'/ny*W-W/2;
            
            [xgrid, ygrid] = meshgrid(1:L,y);
            
            
            U1= U0*ones(size(y)); %U0/2 + U0/2*exp(-(y).^2/(delta^2));
            %U1 = U0*exp(-(y).^2/(delta^2));
            U2=0*exp(-(y).^2/(delta^2));
            
            U1mat =  U0*ones(size(ygrid));%U0/2 + U0/2*exp(-(ygrid).^2/(delta^2));%exp(-(ygrid).^2/(delta^2));
            %U1mat = U0*exp(-(ygrid).^2/(delta^2));
            U2mat = 0*exp(-(ygrid).^2/(delta^2));
            
            U1y = diffxy(ygrid,U1mat,1,1);
            U2y = diffxy(ygrid,U2mat,1,1);
            U1yy = diffxy(ygrid,U1mat,1,2);
            U2yy = diffxy(ygrid,U2mat,1,2);
            %
            % beta1=beta+F1*(U1mat-U2mat) - U1yy;
            %
            %
            % beta2=beta-F2*(U1mat-U2mat) - U2yy;
            
            % beta10=beta+F1*(U1mat-U2mat);
            % beta20=beta-F2*(U1mat-U2mat);
            
            
            ddy2=[-1,zeros(1,ny-1);eye(ny);zeros(1,ny-1),-1];
            ddy2=diff(diff(ddy2));
            %ddy2(1,1) = 2; ddy2(1,2) = -5; ddy2(1,3) = 4; ddy2(1,4) = -1;
            %ddy2(end,end) = 2; ddy2(end,end-1) = -5; ddy2(end,end-2) = 4; ddy2(end,end-3) = -1;
            ddy2 = ddy2/dy^2;
            ddy2n=[-1,zeros(1,ny-1);eye(ny);zeros(1,ny-1),-1];
            ddy2n=diff(diff(ddy2n));
            ddy2n(1,1) = 2; ddy2n(1,2) = -5; ddy2n(1,3) = 4; ddy2n(1,4) = -1;
            ddy2n(end,end) = 2; ddy2n(end,end-1) = -5; ddy2n(end,end-2) = 4; ddy2n(end,end-3) = -1;
            ddy2n = ddy2n/dy^2;
            %     ddy2n=[1,zeros(1,ny-1);eye(ny);zeros(1,ny-1),1];
            %     ddy2n=diff(diff(ddy2n))/dy^2;
            
            Q1y=beta-ddy2n*U1+F1*(U1-U2);
            % Q1y=beta+F1*(U1-U2);
            % Q1y = -ddy2*U1;
            %     Q1y(1) = Q1y(2);
            %     Q1y(end) = Q1y(end-1);
            Q2y=beta-ddy2n*U2+F2*(U2-U1);
            % Q2y=beta+F2*(U2-U1);
            %Q2y=-ddy2*U2;
            %     Q2y(1) = Q2y(2);
            %     Q2y(end) = Q2y(end-1);
            
            plot(Q1y,y,Q2y,y,y*0,y);
            minqy = [minqy, min(Q1y)];
            drawnow;
            %return;
            
            M0=[ddy2-F1*eye(ny),F1*eye(ny);F2*eye(ny),ddy2-F2*eye(ny)];
            %Mz = [-F1*eye(ny),F1*eye(ny);F2*eye(ny),-F2*eye(ny)];
            U=[diag(U1),zeros(ny);zeros(ny),diag(U2)];
            Qy=[diag(Q1y),zeros(ny);zeros(ny),diag(Q2y)];
            
            ks=[1:64]*pi/(4*L);
            %ks=[1:128]/128*0.1;
            om=[];
            eigenvec = [];
            kvec = [];
            cvec = [];
            petotal = [];
            ketotal = [];
            for k=ks
                k; %,fflush(1);
                M=U+Qy*inv(M0-k*k*eye(2*ny));
                [V,D] = eig(M);
                vec = zeros(512,1);
                c = diag(D);
                if (sum(imag(c) ~= 0)) ~= 0
                    %c(abs(real(c)) > 60) = NaN;
                    [~, ind] = sort(imag(c),'descend');
                    for i = 1:min(sum(imag(c)>0),1)
                        vec = V(:,ind(i));
                        eigenvec = [eigenvec,vec];
                        kvec = [kvec, k];
                        cvec = [cvec, c(ind(i))];
                        
                        psis = inv(M0-k*k*eye(2*ny))*vec;
                        psi1 = real(psis(1:256)*exp(1i*k*(1:L)));
                        psi2 = real(psis(257:end)*exp(1i*k*(1:L)));
                        
                        
                        
                        
                        % buoyancy/ mass fluxes
                        
                        KEflux = 1/2*diffxy(xgrid,psi1,2).*diffxy(ygrid,psi1,1).*U1y;
                        ketotal = [ketotal, sum(KEflux(:))];
                        %                 DKM = (-diffxy(ygrid,diffxy(xgrid,psi1,2).*diffxy(ygrid,psi1,1),1)+...
                        %                     diffxy(xgrid,diffxy(ygrid,psi1,1).*diffxy(ygrid,psi1,1),2)).*U1mat;
                        %                 dkmtotal = [dkmtotal, sum(DKM(:))];
                        PEflux = F0.*(U1mat-U2mat).*diffxy(xgrid,psi2,2).*psi1;
                        petotal = [petotal,sum(PEflux(:))];
                        
                    end
                    
                end
                %  c = eig(M);
                om=[om,imag(c)*k];
                
            end
            %
            kvec = kvec.^2/F1;
            cvec = cvec*F2/beta;
            if ~isempty(cvec)
                % plot things
                titlebase = sprintf('Background U0 F %d delta %d beta %d',F,delta,beta0);
                
                kmax = 4;
                
                                h1 = figure;
                                scatter(kvec,real(cvec)-U0/2*F2/beta,20,kvec,'filled')
                                xlabel('Zonal Wavenumber k')
                                grid on
                                xlim([0 kmax])
                                ylabel('Real Phase Speed (c_r)')
                                tit = ['Real Phase Speed ',titlebase];
                                title(tit)
                                SaveFigureGFD(h1,titlebase,tit);
                
                                h2 = figure;
                                scatter(kvec,imag(cvec),20,kvec,'filled')
                                xlabel('Zonal Wavenumber k')
                                xlim([0 kmax])
                                grid on
                                ylabel('Imaginary Phase Speed (c_i)')
                                tit = ['Imaginary Phase Speed ',titlebase];
                                title(tit)
                                SaveFigureGFD(h2,titlebase,tit);
                
                                h3 = figure;
                                scatter(kvec,kvec.*imag(cvec),20,kvec,'filled')
                                xlabel('Zonal Wavenumber k')
                                xlim([0 kmax])
                                grid on
                                ylabel('Growth Rate (kc_i)')
                                tit = ['Growth Rate ',titlebase];
                                title(tit)
                                SaveFigureGFD(h3,titlebase,tit);
                
                                h4 = figure;
                                scatter(kvec,ketotal,20,kvec,'filled')
                                xlabel('Zonal Wavenumber k')
                                xlim([0 kmax])
                                grid on
                                totalgrowthKE = sum(ketotal);
                                ylabel('KE Conversion')
                                tit = [sprintf('KE Conversion =%g    ',totalgrowthKE),titlebase];
                                title(tit)
                                SaveFigureGFD(h4,titlebase,tit);
                
                                h5 = figure;
                                scatter(kvec,petotal,20,kvec,'filled')
                                xlabel('Zonal Wavenumber k')
                                totalgrowth = sum(petotal);
                                xlim([0 kmax])
                                grid on
                                ylabel('PE Conversion')
                                tit = [sprintf('PE Conversion =%g   ',totalgrowth),titlebase];
                                title(tit)
                                SaveFigureGFD(h5,titlebase,tit);
                
                                h6 = figure;
                                plot(Q1y,y,Q2y,y,y*0,y);
                                xlabel('Potential Vorticity Gradient')
                                ylabel('Y')
                                legend('Q_y^1','Q_y^2')
                                tit = ['Potential Vorticity Gradient  ', titlebase];
                                title(tit)
                                SaveFigureGFD(h6,titlebase,tit)
                
                [~,locs] = findpeaks(kvec.*imag(cvec),'SortStr','descend');
                %sigs=[sigs;sig];
                
                psis = inv(M0-kvec(locs(1))*kvec(locs(1))*eye(2*ny))*eigenvec(:,locs(1));
                psi1 = real(psis(1:256)*exp(1i*kvec(locs(1))*(1:L)));
                psi2 = real(psis(257:end)*exp(1i*kvec(locs(1))*(1:L)));
                q1 = real(eigenvec(1:256,locs(1))*exp(1i*kvec(locs(1))*(1:L)));
                q2 = real(eigenvec(257:end,locs(1))*exp(1i*kvec(locs(1))*(1:L)));
                v1 = diffxy(xgrid,psi1,2);
                v2 = diffxy(xgrid,psi2,2);
                
                %                 e1 = sum(q1.*q1,2)./Q1y/L;
                %                 e2 = sum(q2.*q2,2)./Q2y/L;
                %                 etot = -kvec(locs(1))*imag(cvec(locs(1)))/2*(e1+e2);
                
                
                dudt = (sum(v1.*q1,2)+ sum(v2.*q2,2))/2/L;
                dudt1 = (sum(v1.*q1,2))/2/L;
                dudt2 =   (sum(v2.*q2,2))/2/L;
                usum = sum(dudt)*dy/L*2;
                usumtotal = [usumtotal, dudt(128)];
                kemaxgrowth = [kemaxgrowth, ketotal(locs(1))];
                
                
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
                usumtotal = [usumtotal, 0];
                kemaxgrowth = [kemaxgrowth, 0];
            end
            
        end
        titlebase = sprintf('Background U0 F %d beta %d',F,beta0);
        
        h8 = figure(341);
        yyaxis left
        plot(deltas,usumtotal)
        xlabel('Delta')
        ylabel('du/dt at center')
        yyaxis right
        plot(deltas,minqy)
        ylabel('Minimum Q_y')
        grid on
        titl = ['dudt ', titlebase];
        title(titl)
        SaveFigureGFD(h8,titlebase,titl);
        
        h9 = figure(342);
        yyaxis left
        plot(deltas,kemaxgrowth)
        ylabel('KE Conversion')
        grid on
        xlabel('Delta')
        yyaxis right
        plot(deltas,minqy)
        ylabel('Minimum Q_y')
        titl = ['KE Conversion ', titlebase];
        title(titl)
        SaveFigureGFD(h9,titlebase,titl)
    end
end
%%



% [~,locs] = findpeaks(kvec.*imag(cvec),'SortStr','descend');
% %sigs=[sigs;sig];
%
% psis = inv(M0-kvec(locs(1))*kvec(locs(1))*eye(2*ny))*eigenvec(:,locs(1));
% psi1 = real(psis(1:256)*exp(1i*kvec(locs(1))*(1:L)));
% psi2 = real(psis(257:end)*exp(1i*kvec(locs(1))*(1:L)));
% q1 = real(eigenvec(1:256,locs(1))*exp(1i*kvec(locs(1))*(1:L)));
% q2 = real(eigenvec(257:end,locs(1))*exp(1i*kvec(locs(1))*(1:L)));
% v1 = diffxy(xgrid,psi1,2);
% v2 = diffxy(xgrid,psi2,2);
%
% e1 = sum(q1.*q1,2)./Q1y/L;
% e2 = sum(q2.*q2,2)./Q2y/L;
% etot = -kvec(locs(1))*imag(cvec(locs(1)))/2*(e1+e2);
%
%
% dudt = (sum(v1.*q1,2)+ sum(v2.*q2,2))/2/L;
% dudt1 = (sum(v1.*q1,2))/2/L;
% dudt2 =   (sum(v2.*q2,2))/2/L;
% usum = sum(dudt)*dy/L*2;
%
% h7 = figure;
% plot(dudt,y,y*0,y,dudt1,y,dudt2,y,etot,y);
% xlabel('sum dudt')
% ylabel('Y')
% tit = [sprintf('Mean Velocity Change %g ',usum), titlebase];
% title(tit)
% SaveFigureGFD(h7,titlebase,tit)
%
% figure;
% plot(-e1,y,-e2,y)
%
%
% % buoyancy/ mass fluxes

%     KEflux = diffxy(xgrid,psi1,2).*diffxy(ygrid,psi1,1).*U1y;
%     ketotal = sum(KEflux(:));
%     DKM = (-diffxy(ygrid,diffxy(xgrid,psi1,2).*diffxy(ygrid,psi1,1),1)+...
%         diffxy(xgrid,diffxy(ygrid,psi1,1).*diffxy(ygrid,psi1,1),2)).*U1mat;
%     dkmtotal = sum(DKM(:));
%     PEflux = F1.*(U1mat-U2mat).*diffxy(xgrid,psi2,2).*psi1;
%     petotal = sum(PEflux(:));
%
%     figure;
%     contourf(KEflux)
%     titK = sprintf('MKE -> EKE Flux %d Mode1',ketotal);
%     title(titK)
%     figure;
%     contourf(DKM)
%     titKM = sprintf('EKE -> MKE Flux %d Mode1',dkmtotal);
%     title(titKM)
%     figure;
%     contourf(PEflux)
%     titP = sprintf('EAPE -> EKE Flux %d Mode1',petotal);
%     title(titP)
%
%
%     if length(locs)>1
%         psis2 = inv(M0-ks(locs(1))*ks(locs(2))*eye(2*ny))*eigenvec(:,locs(2));
%         psi12 = real(psis(1:256)*exp(1i*ks(locs(2))*(1:L)));
%         psi22 = real(psis(257:end)*exp(1i*ks(locs(2))*(1:L)));
%
%
%
%
%         % buoyancy/ mass fluxes mode 2
%
%         KEflux2 = diffxy(xgrid,psi12,2).*diffxy(ygrid,psi12,1).*U1y;
%         ketotal2 = sum(KEflux2(:));
%         PEflux2 = F1.*(U1mat-U2mat).*diffxy(xgrid,psi22,2).*psi12;
%         petotal2 = sum(PEflux2(:));
%
%         figure;
%         contourf(KEflux2)
%         titK = sprintf('KE Flux %d Mode2',ketotal2);
%         title(titK)
%         figure;
%         contourf(PEflux2)
%         titP = sprintf('PE Flux %d Mode2',petotal2);
%         title(titP)
%
%
%
%         if length(locs)>2
%             psi13 = real(eigenvec(1:256,locs(end-2))*exp(1i*ks(locs(end-2))*(1:L)));
%             psi23 = real(eigenvec(257:end,locs(end-2))*exp(1i*ks(locs(end-2))*(1:L)));
%
%
%
%             % buoyancy/ mass fluxes mode 2
%
%             KEflux3 = -diffxy(xgrid,psi13,2).*diffxy(ygrid,psi13,1).*U1y;
%             ketotal3 = sum(KEflux3(:));
%             PEflux3 = F1.*(U1mat-U2mat).*diffxy(xgrid,psi23,2).*psi13;
%             petotal3 = sum(PEflux3(:));
%
%             figure;
%             contourf(KEflux3)
%             titK = sprintf('KE Flux %d Mode3',ketotal3);
%             title(titK)
%             figure;
%             contourf(PEflux3)
%             titP = sprintf('PE Flux %d Mode3',petotal3);
%             title(titP)
%         end

% figure;
% plot(ks,flipud(sigs));
% legend('d=0.2','d=0')
% title(sprintf('Growth rates (1/d), L = %g, U = %g',L,U0));
% xlabel('k 1/km');

%
