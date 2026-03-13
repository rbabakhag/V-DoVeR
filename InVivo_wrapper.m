clc; clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This is the final code for the post processing of the synthetic data %
%%% Date: Aug / 16 / 2024                                                %
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% main directory
load('D:\Reza\VascularDoVeR\Original\results\images\synthetic.carotid\carotid+demo+redone+20190822.mat');
Dir = 'D:\Reza\VascularDoVeR\Original\results\mouse\';
plotdir  = fullfile(Dir,'comparison');
if ~exist(plotdir,'dir')
    mkdir(plotdir)
end 

%% directories and results
% directories
% cfd_dir   = ['D:\Reza\VascularDoVeR\Original\analysis\synthetic.carotid\rawdata'];
% cfi_dir   = ['D:\Reza\VascularDoVeR\Original\analysis\synthetic.carotid\matfiles'];
% dover_dir = ['D:\Reza\VascularDoVeR\Original\results\images\synthetic.carotid\dover'];
ivfm_dir  = ['D:\Reza\VascularDoVeR\Original\results\mouse\ivfm'];
% vfm_dir   = ['D:\Reza\VascularDoVeR\Original\results\images\synthetic.carotid\vfm'];
% results files: cfd
% cfd_results  = dir(fullfile(cfd_dir, '*.txt'));
% % color doppler images
% cfi_results = dir(fullfile(cfi_dir,'*.mat'));
% % results files: dover
% dov_results  = dir(fullfile(dover_dir, '*.mat'));
% results files: ivfm
ivfm_results = dir(fullfile(ivfm_dir, '*.mat'));
% % results files: vfm
% vfm_results  = dir(fullfile(vfm_dir, '*.mat'));
% Pre-allocate memory for Dover, iVFM, and VFM results
ucfd    = [];
vcfd    = [];
Dov_u   = [];
Dov_v   = [];
ivfm_u  = [];
ivfm_v  = [];
vfm_u   = [];
vfm_v   = [];
%read the results and put them together
%Dover
for i  = 1:size(upsi,3)
        Dov_u = cat(3, Dov_u, upsi(:,:,i));
        Dov_v = cat(3, Dov_v, vpsi(:,:,i));
end
%iVFM
for i = 1:size(ivfm_results)
        IVFM   = load(fullfile(ivfm_dir, ivfm_results(i).name));
        ivfm_u = cat(3, ivfm_u, IVFM.uF);
        ivfm_v = cat(3, ivfm_v, IVFM.vF);
end
%VFM
for i = 1:size(uvfm,3)
        vfm_u = cat(3, vfm_u, uvfm(:,:,i));
        vfm_v = cat(3, vfm_v, vvfm(:,:,i));
end

%% Vorticity Calculation 
%{
% CFD GROUND TRUTH
frRng       =   1:13;
X  = data.x;
Y  = data.y;
dx = data.dx;
dy = data.dy;
mask = data.bcmask;
% DoVeR RESULTS
calDoVeR    =   zeros(size(Dov_u));
omegaDoVeR  =   zeros(size(Dov_v));
for nn = 1:numel(frRng)
    fprintf(['Now Evaluating Frame ',num2str(nn),' of ',num2str(size(Dov_u,3)),' for DoVeR \n']);
    [calDoVeR(:,:,nn),omegaDoVeR(:,:,nn)] = coherent_structure_id(X,...
        Y, Dov_u(:,:,nn), -Dov_v(:,:,nn),dx, dx, mask, 'LambdaCI', 3, 3);
end
% VFM RESULTS
calVFM      =   zeros(size(vfm_u));
omegaVFM    =   zeros(size(vfm_v));
for nn = 1:numel(frRng)
    fprintf(['Now Evaluating Frame ',num2str(nn),' of ',num2str(size(vfm_u)),' for VFM \n']);
    [calVFM(:,:,nn),omegaVFM(:,:,nn)] = coherent_structure_id(X,...
        Y, vfm_u(:,:,nn), -vfm_v(:,:,nn),dx, dx, mask, 'LambdaCI', 3, 3);
end
% iVFM RESULTS
caliVFM     =   zeros(size(ivfm_u));
omegaiVFM   =   zeros(size(ivfm_v));
for nn = 1:numel(frRng)
    fprintf(['Now Evaluating Frame ',num2str(nn),' of ',num2str(size(ivfm_u)),' for iVFM \n']);
    [caliVFM(:,:,nn),omegaiVFM(:,:,nn)] = coherent_structure_id(X,...
        Y, ivfm_u(:,:,nn), -ivfm_v(:,:,nn),dx, dx, mask, 'LambdaCI', 3, 3);
end

% SET COLOR MAP TO RED BLUE
rbmap   = colormap(redblue(40));
% Center plot grids about 0
xp  =   X - mean(X(:));
yp  =   Y - mean(Y(:));
%
%% Plot the reconstruction results for each method
for nn = 1:numel(frRng)
  
    figure(2);
    toDoVeR = omegaDoVeR(:,:,nn); 
    toDoVeR(toDoVeR == 0) = nan;
    h = pcolor(xp(1,:,1)*100,yp(:,1,1)*100,toDoVeR);
    set(h, 'EdgeColor', 'none');
    colormap(rbmap); 
    caxis([-300, 300]);  % Adjusted color axis limits
    hold on;
    imcontour(xp(:,:,1)*100,yp(:,:,1)*100,mask,[-0.05 0.05],...
        'Color',[0 0 0],'LineWidth',1.4);
    quiver(xp(10:10:end,10:30:end,1)*100,yp(10:10:end,10:30:end,1)*100,...
        Dov_u(10:10:end,10:30:end,nn)/3,Dov_v(10:10:end,10:30:end,nn)/3,...
        'Color',[0 0 0],'autoscale','off','LineWidth',1.5,'MaxHeadSize',5);
    caxis([-300, 300]);  % Adjusted color axis limits
    hold off;
    set(gca,'FontSize',28,'FontWeight','bold','YDir','reverse');
    set(gcf,'Position',[100 100 500 300],'Color',[1 1 1]);
    xlabel('Position (cm)');
    ylabel('Position (cm)');
    axis image;
    axis([min(xp(:)) max(xp(:))+0.002 min(yp(:)) max(yp(:))]*100);
    axis off;
    export_fig(2,fullfile(plotdir,sprintf('DoVeR_%03i.tiff',nn)),'-tiff','-a1','-r512');
    
    figure(3);
    toiVFM = omegaVFM(:,:,nn); 
    toiVFM(toiVFM == 0) = nan;
    h = pcolor(xp(1,:,1)*100,yp(:,1,1)*100,toiVFM);
    set(h, 'EdgeColor', 'none');
    colormap(rbmap); 
    caxis([-300, 300]);  % Adjusted color axis limits
    hold on;
    imcontour(xp(:,:,1)*100,yp(:,:,1)*100,mask,[-0.05 0.05],...
        'Color',[0 0 0],'LineWidth',1.4);
    quiver(xp(10:10:end,10:30:end,1)*100,yp(10:10:end,10:30:end,1)*100,...
        vfm_u(10:10:end,10:30:end,nn)/3,-vfm_v(10:10:end,10:30:end,nn)/3,...
        'Color',[0 0 0],'autoscale','off','LineWidth',1.5,'MaxHeadSize',5);
    caxis([-300, 300]);  % Adjusted color axis limits
    hold off;
    set(gca,'FontSize',28,'FontWeight','bold','YDir','reverse');
    set(gcf,'Position',[600 600 500 300],'Color',[1 1 1]);
    xlabel('Position (cm)');
    ylabel('Position (cm)');
    axis image;
    axis([min(xp(:)) max(xp(:))+0.002 min(yp(:)) max(yp(:))]*100);
    axis off;
    export_fig(gcf,fullfile(plotdir,sprintf('VFM_%03i.tiff',nn)),'-tiff','-a1','-r512');
    
    figure(4);
    toiVFM = omegaiVFM(:,:,nn).*mask; 
    toiVFM(toiVFM == 0) = nan;
    h = pcolor(xp(1,:,1)*100,yp(:,1,1)*100,toiVFM);
    set(h, 'EdgeColor', 'none');
    colormap(rbmap); 
    caxis([-300, 300]);  % Adjusted color axis limits
    hold on;
    imcontour(xp(:,:,1)*100,yp(:,:,1)*100,mask.*mask,[-0.05 0.05],...
        'Color',[0 0 0],'LineWidth',1.4);
    caxis([-300, 300]);  % Adjusted color axis limits
    quiver(xp(10:10:end,10:30:end,1)*100,yp(10:10:end,10:30:end,1)*100,...
        ivfm_u(10:10:end,10:30:end,nn)/3.*mask(10:10:end,10:30:end),...
        ivfm_v(10:10:end,10:30:end,nn)/3.*mask(10:10:end,10:30:end),...
        'Color',[0 0 0],'autoscale','off','LineWidth',1.5,'MaxHeadSize',5);
    hold off;
    set(gca,'FontSize',28,'FontWeight','bold','YDir','reverse');
    set(gcf,'Position',[600 100 500 300],'Color',[1 1 1]);
    xlabel('Position (cm)');
    ylabel('Position (cm)');
    axis image;
    axis([min(xp(:)) max(xp(:))+0.002 min(yp(:)) max(yp(:))]*100);
    axis off;
    export_fig(gcf,fullfile(plotdir,sprintf('iVFM_%03i.tiff',nn)),'-tiff','-a1','-r512');
end

% Plot and save the color bar for the vector field reconstruction plots
figure(5);
h = colorbar('EastOutside');
colormap(redblue(40));
caxis([-300, 300]);  % Adjusted color axis limits
ylabel(h,'Vorticity (s^-1)','FontName','Times New Roman','FontSize',20,'FontWeight','bold');
set(gca,'FontName','Times New Roman','FontSize',20,'FontWeight','bold','YTick',-0.5:0.5:0.5);
set(gcf,'Color',[1 1 1]);
export_fig(gcf,fullfile(plotdir,'carotid_quiver_colorbar.tiff'),'-tif','-a1','-r512');

% Pre-allocate memory for computing the magnitude and angle of each vector
mcfd    =   zeros(size(vcfd));
mDoVeR  =   zeros(size(vcfd));
mVFM    =   zeros(size(vcfd));
miVFM   =   zeros(size(vcfd));
acfd    =   zeros(size(vcfd));
aDoVeR  =   zeros(size(vcfd));
aVFM    =   zeros(size(vcfd));
aiVFM   =   zeros(size(vcfd));
% Compute the magnitude and angle of the vectors in the reconstruction
for tt = 1:numel(frRng)
    % Compute the magnitude
    mDoVeR(:,:,tt)  = sqrt(Dov_u(:,:,tt).^2 +   Dov_v(:,:,tt).^2).*...
        imerode(mask(:,:),strel('disk',3));
    mVFM(:,:,tt)    = sqrt(vfm_u(:,:,tt).^2  +   vfm_v(:,:,tt).^2 ).*...
        imerode(mask(:,:),strel('disk',3));
    miVFM(:,:,tt)   = sqrt(ivfm_u(:,:,tt).^2 +   ivfm_v(:,:,tt).^2).*...
        imerode(mask(:,:),strel('disk',3));

    % Compute the angle
    aDoVeR(:,:,tt)  = atan2( Dov_v(:,:,tt), Dov_u(:,:,tt) )*180/pi.*...
        imerode(mask(:,:),strel('disk',3));
    aVFM(:,:,tt)    = atan2(-vfm_v(:,:,tt), vfm_u(:,:,tt) )*180/pi.*...
        imerode(mask(:,:),strel('disk',3));
    aiVFM(:,:,tt)   = atan2( ivfm_v(:,:,tt), ivfm_u(:,:,tt) )*180/pi.*...
        imerode(mask(:,:),strel('disk',3));
end

% Compute the difference of the angle between CFD and each reconstruction
dAngDoVeR   =   acfd - aDoVeR;
dAngDoVeR(dAngDoVeR >    pi) = dAngDoVeR(dAngDoVeR >    pi) - 2*pi;
dAngDoVeR(dAngDoVeR <   -pi) = dAngDoVeR(dAngDoVeR <   -pi) + 2*pi;
dAngVFM     =   acfd - aVFM;
dAngVFM(dAngVFM >    pi)    = dAngVFM(dAngVFM >    pi) - 2*pi;
dAngVFM(dAngVFM <   -pi)    = dAngVFM(dAngVFM <   -pi) + 2*pi;
dAngiVFM    =   acfd - aiVFM;
dAngiVFM(dAngiVFM >    pi)  = dAngiVFM(dAngiVFM >    pi) - 2*pi;
dAngiVFM(dAngiVFM <   -pi)  = dAngiVFM(dAngiVFM <   -pi) + 2*pi;

% Use HISTC to bind the magnitude differences 
binDoVeR    = histc(abs(mDoVeR(mcfd~=0)-mcfd(mcfd~=0))/...
    (max(mcfd(:))-min(mcfd(:)))*100,0:1:100);
binVFM      = histc(abs(mVFM(mcfd~=0)-mcfd(mcfd~=0))/...
    (max(mcfd(:))-min(mcfd(:)))*100,0:1:100);
biniVFM     = histc(abs(miVFM(mcfd~=0)-mcfd(mcfd~=0))/...
    (max(mcfd(:))-min(mcfd(:)))*100,0:1:100);

% Plot the CDF of the magnitude differences
figure(18);
plot(0:1:101,cumsum([0;binVFM])/sum(binVFM)*1,'Color',[0.4660, 0.6740, 0.1880],'LineWidth',4);
hold on;
plot(0:1:101,cumsum([0;biniVFM])/sum(biniVFM)*1,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',4);
plot(0:1:101,cumsum([0;binDoVeR])/sum(binDoVeR)*1,'Color',[0.0000, 0.4470, 0.7410],'LineWidth',4);
hold off;
grid on;
axis([0 100 0 1]);
set(gca,'LineWidth',2,'FontSize',20,'FontWeight','bold');
ylabel('Velocity Magn. nRMSE CDF','FontSize',18,'FontWeight','bold');
xlabel('Relative Error, %','FontSize',18,'FontWeight','bold');
set(gcf,'Position',[10 10 300 300],'Color',[1 1 1]);
export_fig(18,fullfile(plotdir,'carotid_error_CDF_mag.png'),'-a1','-r512','-png');
legend('iVFM_{1D}','iVFM_{2D}','DoVeR','Location','EastOutside');
legend boxoff
export_fig(18,fullfile(plotdir,'carotid_error_CDF_mag+legend.png'),'-a1','-r512','-png');

% Use HISTC to bind the angle differences 
binDoVeR = histc(abs(dAngDoVeR(mcfd~=0)),0:1:100);
binVFM  = histc(abs(dAngVFM(mcfd~=0)),0:1:100);
biniVFM = histc(abs(dAngiVFM(mcfd~=0)),0:1:100);
figure(19);
plot(0:1:101,cumsum([0;binVFM])/sum(binVFM)*1,'Color',[0.4660, 0.6740, 0.1880],'LineWidth',4);

% Plot the CDF of the angle differences
hold on;
plot(0:1:101,cumsum([0;biniVFM])/sum(biniVFM)*1,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',4);
plot(0:1:101,cumsum([0;binDoVeR])/sum(binDoVeR)*1,'Color',[0.0000, 0.4470, 0.7410],'LineWidth',4);
hold off;
grid on;
axis([0 100 0 1]);
set(gca,'LineWidth',2,'FontSize',20,'FontWeight','bold');
ylabel('Velocity Direction RMSE CDF','FontSize',18,'FontWeight','bold');
xlabel('RMS Error, \circ','FontSize',18,'FontWeight','bold');
set(gcf,'Position',[10 10 300 300],'Color',[1 1 1]);
export_fig(19,fullfile(plotdir,'carotid_error_CDF_angle.png'),'-a1','-r512','-png');

% Compute the velocity magnitudes for each reconstruction method
% MASK = min(mask,[],3);
for tt = 1:size(Dov_u,3)
    UDoVeR(:,:,tt)  = sqrt(Dov_u(:,:,tt).^2+Dov_v(:,:,tt).^2).*imerode(mask(:,:),strel('disk',1));
    UVFM(:,:,tt)    = sqrt(vfm_u(:,:,tt).^2+vfm_v(:,:,tt).^2).*imerode(mask(:,:),strel('disk',1));
    UiVFM(:,:,tt)   = sqrt(ivfm_u(:,:,tt).^2+ivfm_v(:,:,tt).^2).*imerode(mask(:,:),strel('disk',1));
end

%% Construct the elements for a Bland-Altman or Tukey plot
label1 = 'DoVeR'; label2 = 'CFD';
threshold = quantile(abs(UCFD(abs(UCFD) > 0)),0.01);
sclv = (max(UCFD(:))-min(UCFD(:)));

meanDiffDoVeR   = (UCFD-UDoVeR)/sclv*100;
meanMeanDoVeR   = (UCFD+UDoVeR)/sclv/2*100;
stdDoVeR    = std((UCFD(abs(UCFD)>threshold)-...
    UDoVeR(abs(UCFD)>threshold))/sclv*100);
meanDoVeR   = mean((UCFD(abs(UCFD)>threshold)-...
    UDoVeR(abs(UCFD)>threshold))/sclv*100);

meanDiffVFM     = (UCFD-UVFM)/sclv*100;
meanMeanVFM     = (UCFD+UVFM)/sclv/2*100;
stdVFM      = std((UCFD(abs(UCFD)>threshold)-...
    UVFM(abs(UCFD)>threshold))/sclv*100);
meanVFM     = mean((UCFD(abs(UCFD)>threshold)-...
    UVFM(abs(UCFD)>threshold))/sclv*100);

meanDiffiVFM    = (UCFD-UiVFM)/sclv*100;
meanMeaniVFM    = (UCFD+UiVFM)/sclv/2*100;
stdiVFM     = std((UCFD(abs(UCFD)>threshold)-...
    UiVFM(abs(UCFD)>threshold))/sclv*100);
meaniVFM    = mean((UCFD(abs(UCFD)>threshold)-...
    UiVFM(abs(UCFD)>threshold))/sclv*100);

%% Display the Bland-Altman plot
figure(21);
scatter(meanMeanVFM(1),meanDiffVFM(1),20,'filled',...
    'MarkerFaceColor',[0.4660, 0.6740, 0.1880],'MarkerEdgeColor','none',...
    'MarkerFaceAlpha',1/1,'LineWidth',1.2);
hold on;
plot(linspace(-250,250,501),...
    linspace(meanVFM,meanVFM,501),'LineStyle','-',...
    'Color',[0.4660, 0.6740, 0.1880],'LineWidth',2);
plot(linspace(-250,250,501),...
    linspace(meanVFM-1.96*stdVFM,meanVFM-1.96*stdVFM,501),'LineStyle','--',...
    'Color',[0.4660, 0.6740, 0.1880],'LineWidth',2);
plot(linspace(-250,250,501),...
    linspace(meanVFM+1.96*stdVFM,meanVFM+1.96*stdVFM,501),'LineStyle','--',...
    'Color',[0.4660, 0.6740, 0.1880],'LineWidth',2);

scatter(meanMeaniVFM(1),meanDiffiVFM(1),20,'filled',...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980],'MarkerEdgeColor','none',...
    'MarkerFaceAlpha',1/1,'LineWidth',1.2);
plot(linspace(-250,250,501),...
    linspace(meaniVFM,meaniVFM,501),'LineStyle','-',...
    'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2);
plot(linspace(-250,250,501),...
    linspace(meaniVFM-1.96*stdiVFM,meaniVFM-1.96*stdiVFM,501),'LineStyle','--',...
    'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2);
plot(linspace(-250,250,501),...
    linspace(meaniVFM+1.96*stdiVFM,meaniVFM+1.96*stdiVFM,501),'LineStyle','--',...
    'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2);

scatter(meanMeanDoVeR(1),meanDiffDoVeR(1),20,'filled',...
    'MarkerFaceColor',[0.0000, 0.4470, 0.7410],'MarkerEdgeColor','none',...
    'MarkerFaceAlpha',1/1,'LineWidth',1.2);
plot(linspace(-250,250,501),...
    linspace(meanDoVeR,meanDoVeR,501),'LineStyle','-',...
    'Color',[0.0000, 0.4470, 0.7410],'LineWidth',2);
plot(linspace(-250,250,501),...
    linspace(meanDoVeR-1.96*stdDoVeR,meanDoVeR-1.96*stdDoVeR,501),'LineStyle','--',...
    'Color',[0.0000, 0.4470, 0.7410],'LineWidth',2);
plot(linspace(-250,250,501),...
    linspace(meanDoVeR+1.96*stdDoVeR,meanDoVeR+1.96*stdDoVeR,501),'LineStyle','--',...
    'Color',[0.0000, 0.4470, 0.7410],'LineWidth',2);

scatter(meanMeanVFM(200:200:end),meanDiffVFM(200:200:end),20,'filled',...
    'MarkerFaceColor',[0.4660, 0.6740, 0.1880],'MarkerEdgeColor','none',...
    'MarkerFaceAlpha',1/2,'LineWidth',1.2);
hold on;
scatter(meanMeaniVFM(200:200:end),meanDiffiVFM(200:200:end),20,'filled',...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980],'MarkerEdgeColor','none',...
    'MarkerFaceAlpha',1/2,'LineWidth',1.2);
scatter(meanMeanDoVeR(200:200:end),meanDiffDoVeR(200:200:end),20,'filled',...
    'MarkerFaceColor',[0.0000, 0.4470, 0.7410],'MarkerEdgeColor','none',...
    'MarkerFaceAlpha',1/2,'LineWidth',1.2);

plot(linspace(-250,250,501),...
    linspace(meanVFM,meanVFM,501),'LineStyle','-',...
    'Color',[0.4660, 0.6740, 0.1880],'LineWidth',2);
plot(linspace(-250,250,501),...
    linspace(meanVFM-1.96*stdVFM,meanVFM-1.96*stdVFM,501),'LineStyle','--',...
    'Color',[0.4660, 0.6740, 0.1880],'LineWidth',2);
plot(linspace(-250,250,501),...
    linspace(meanVFM+1.96*stdVFM,meanVFM+1.96*stdVFM,501),'LineStyle','--',...
    'Color',[0.4660, 0.6740, 0.1880],'LineWidth',2);

plot(linspace(-250,250,501),...
    linspace(meaniVFM,meaniVFM,501),'LineStyle','-',...
    'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2);
plot(linspace(-250,250,501),...
    linspace(meaniVFM-1.96*stdiVFM,meaniVFM-1.96*stdiVFM,501),'LineStyle','--',...
    'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2);
plot(linspace(-250,250,501),...
    linspace(meanDoVeR,meanDoVeR,501),'LineStyle','-',...
    'Color',[0.0000, 0.4470, 0.7410],'LineWidth',2);
plot(linspace(-250,250,501),...
    linspace(meanDoVeR-1.96*stdDoVeR,meanDoVeR-1.96*stdDoVeR,501),'LineStyle','--',...
    'Color',[0.0000, 0.4470, 0.7410],'LineWidth',2);
plot(linspace(-250,250,501),...
    linspace(meanDoVeR+1.96*stdDoVeR,meanDoVeR+1.96*stdDoVeR,501),'LineStyle','--',...
    'Color',[0.0000, 0.4470, 0.7410],'LineWidth',2);
hold off;
grid on;
set(gcf,'Color',[1 1 1]);
axis([0 100 -75 100]);
set(gca,'LineWidth',2,'FontSize',20,'FontWeight','bold')
xlabel('u_{mag} mean, %');
ylabel('u_{mag} difference, %');
% axis square
set(gcf,'Position',[10 10 550 300]);
export_fig(gcf,fullfile(plotdir,'carotid_blandAltman_aggregate.png'),'-a1','-r512','-q400','-png');
legend('iVFM_{1D}',...
    sprintf('Mean: %3.3f',meanVFM),...
    sprintf('LoA_{l}: %5.3f',(meanVFM-1.96*stdVFM)),...
    sprintf('LoA_{u}: %5.3f',(meanVFM+1.96*stdVFM)),...
    'iVFM_{2D}',...
    sprintf('Mean: %3.3f',meaniVFM),...
    sprintf('LoA_{l}: %5.3f',(meaniVFM-1.96*stdiVFM)),...
    sprintf('LoA_{u}: %5.3f',(meaniVFM+1.96*stdiVFM)),...
    'DoVeR',...
    sprintf('Mean: %3.3f',meanDoVeR),...
    sprintf('LoA_{l}: %5.3f',(meanDoVeR-1.96*stdDoVeR)),...
    sprintf('LoA_{u}: %5.3f',(meanDoVeR+1.96*stdDoVeR)),'Location','EastOutside');
legend boxoff;
set(gcf,'Position',[10 10 550 400]);

keyboard
%}
%% 2
% RUN POST-PROCESSING CODES FROM PAPER IMAGES
% RUN VSS CALCULATION
nu  = 3.77e-6;  % Viscosity, m^2/s
rho = 1100;     % Density, kg/m^3
vsscdev = zeros(size(upsi));
vssvfm  = zeros(size(upsi));
cfm = data;
for tt = 1:cfm.s(3)
    fprintf(['Now calculating VSS on frame ',num2str(tt,'%02i'),' \r']);
    vsscdev(:,:,tt) = vss_calc(cfm.x,cfm.y,upsi(:,:,tt),vpsi(:,:,tt),...
        nu,rho,[7 7],'rbf');
    vssvfm(:,:,tt) = vss_calc(cfm.x,cfm.y,uivfm(:,:,tt),vivfm(:,:,tt),...
        nu,rho,[7 7],'rbf');
    vssivfm(:,:,tt) = vss_calc(cfm.x,cfm.y,uivfm(:,:,tt),vivfm(:,:,tt),...
        nu,rho,[7 7],'rbf');
end
%%
% GET SEEDING LOCATIONS
info.voldir = '/Users/brettmeyers/Documents';
offset  = 1E-4;
% [temp.x0,temp.y0] = seeding_locs(cfm.x,cfm.y,upsi(:,:,2),-vpsi(:,:,2),cfm.dx/100,cfm.dy/100,8);
% temp.x0 = linspace(0.001801,0.001801,7)';
% temp.y0 = linspace( 0.00074, 0.00116,7)';

temp.x0 = linspace(0.00100,0.00100,10)';
temp.y0 = cat(2,linspace(0.00038,0.00076,6),linspace(0.00085,0.00105,4))';
% temp.y0 = cat(2,linspace(0.00022,0.00052,6),linspace(0.00086,0.00111,6))';
%{
for tt = 1:cfm.s(end)
%     temp.x0 = linspace(0.0005805,0.0005805,12)';
%     temp.y0 = cat(2,linspace(0.00022,0.00052,6),linspace(0.00086,0.00111,6))';
 figure;   
    plot_streamtrace2D(cfm.x*100,cfm.y*100,upsi(:,:,tt).*cfm.mask*100,-vpsi(:,:,tt).*cfm.mask*100,cfm.dx/100*100,cfm.dy/100*100,...
        1,5E-4,0,0.4,mean(cfm.imgray,4),vsscdev(:,:,tt),temp.x0*100,temp.y0*100);
minr = 1;
maxr = size(upsi,1);
minc = 1;
maxc = size(upsi,2);
%%% REZA
    % plot_streamtrace2D(cfm.x(minr:maxr,minc:maxc)*100,cfm.y(minr:maxr,minc:maxc)*100,...
    %     upsi(minr:maxr,minc:maxc,tt).*cfm.mask(minr:maxr,minc:maxc)*100,...
    %     -vpsi(minr:maxr,minc:maxc,tt).*cfm.mask(minr:maxr,minc:maxc)*100,...
    %     cfm.dx/100*100,cfm.dy/100*100,1,2E-4,0,8E-2,1,...
    %     day(tt).VSScdev(minr:maxr,minc:maxc).*cfm.mask(minr:maxr,minc:maxc),...
    %     temp.x0*100,temp.y0*100)
    caxis([-15 15]);
    h = colorbar;
    set(gca,'FontSize',36,'FontWeight','bold','YDir','Normal');
    set(gcf,'Color',[1 1 1],'Position',[100 100 950 400]);
    ylabel('Position (mm)','FontSize',36,'FontWeight','bold');
    xlabel('Position (mm)','FontSize',36,'FontWeight','bold');
    ylabel(h,'VSS (dynes/cm^2)','FontSize',36,'FontWeight','bold');
    pause(1E-2);
%     export_fig(FIG,['/Volumes/GERI/cDEV/cDEV_paper/new_images/mouse_carotid_new/streamtrace_vss_cdev_',num2str(tt,'%02i'),'.png']);
    export_fig(gcf,fullfile('D:\Reza\VascularDoVeR\Original\results\mouse\comparison\', ...
        ['mouse_carotid_streamtrace_DoVeR_',num2str(tt,'%02i'),'.tiff']),'-tiff','-a1','-r128','-q400');
end

for tt = [4 8]%1:cfm.s(end)
%     close all
%     temp.x0 = cat(1,temp.x0,linspace(0.001801,0.001801,5)');
%     temp.y0 = cat(1,temp.y0,linspace( 0.00074, 0.00116,5)');

%     plot_streamtrace2D(cfm.x*100,cfm.y*100,uvfm(:,:,tt).*cfm.MASK*100,-vvfm(:,:,tt).*cfm.MASK*100,cfm.dx/100*100,cfm.dy/100*100,...
%         1,5E-4,0,0.4,mean(cfm.imgray,4),vssvfm(:,:,tt),temp.x0*100,temp.y0*100);
    plot_streamtrace2D(cfm.x(minr:maxr,minc:maxc)*100,cfm.y(minr:maxr,minc:maxc)*100,...
        uvfm(minr:maxr,minc:maxc,tt).*cfm.mask(minr:maxr,minc:maxc)*100,...
        -vvfm(minr:maxr,minc:maxc,tt).*cfm.mask(minr:maxr,minc:maxc)*100,...
        cfm.dx/100*100,cfm.dy/100*100,1,4E-4,0,8E-2,1,day(tt).VSSvfm(minr:maxr,minc:maxc),...
        temp.x0*100,temp.y0*100)
    caxis([-15 15]);
    h = colorbar;
    set(gca,'FontSize',36,'FontWeight','bold','YDir','Normal');
    set(gcf,'Color',[1 1 1],'Position',[100 100 90 400]);
    ylabel('Position (mm)','FontSize',36,'FontWeight','bold');
    xlabel('Position (mm)','FontSize',36,'FontWeight','bold');
    ylabel(h,'VSS (dynes/cm^2)','FontSize',36,'FontWeight','bold');
    pause(1E-2);
%     export_fig(FIG,['/Volumes/GERI/cDEV/cDEV_paper/new_images/mouse_carotid_new/streamtrace_vss_vfm_',num2str(tt,'%02i'),'.png']);
    export_fig(gcf,fullfile('D:\Reza\VascularDoVeR\Original\results\mouse\comparison\' ...
        ,['mouse_carotid_streamtrace_iVFM_',num2str(tt,'%02i'),'.tiff']),'-tiff','-a1','-r128','-q400');
end
%}
%%
for tt = 3:cfm.s(3)
    figure(41);
    subplot(3,1,1)
    imshow(uint8(flipud(cfm.imgray(:,:,:,tt))),'XData',cfm.x(1,:),'YData',cfm.y(:,1));
    hold on;
    quiver(cfm.x(4:4:end,4:4:end),cfm.y(4:4:end,4:4:end),...
        flipud(upsi(4:4:end,4:4:end,tt))/1000,flipud(vpsi(4:4:end,4:4:end,tt))/1000,'y','autoscale','off');
    hold off;
    axis image
    
    subplot(3,1,2)
    imshow(uint8(flipud(cfm.imgray(:,:,:,tt))),'XData',cfm.x(1,:),'YData',cfm.y(:,1));
    hold on;
    quiver(cfm.x(4:4:end,4:4:end),cfm.y(4:4:end,4:4:end),...
        flipud(uvfm(4:4:end,4:4:end,tt))/1000,flipud(vvfm(4:4:end,4:4:end,tt))/1000,'y','autoscale','off');
    hold off;
    axis image

    subplot(3,1,3)
    imshow(uint8(flipud(cfm.imgray(:,:,:,tt))),'XData',cfm.x(1,:),'YData',cfm.y(:,1));
    hold on;
    quiver(cfm.x(4:4:end,4:4:end),cfm.y(4:4:end,4:4:end),...
        flipud(uivfm(4:4:end,4:4:end,tt))/1000,flipud(vivfm(4:4:end,4:4:end,tt))/1000,'y','autoscale','off');
    hold off;
    axis image
    pause(5E-1);
end
%%
xpos = cfm.x; ypos = cfm.y;
grayimg = cfm.imgray;
mask  = cfm.mask;
ucdev = upsi; vcdev = vpsi;
mvfm = sqrt(uvfm.^2+vvfm.^2); mcdev = sqrt(ucdev.^2+vcdev.^2);
save(fullfile('D:\Reza\VascularDoVeR\Original\results\mouse\comparison'),...
    'xpos','ypos','grayimg','mask','ucdev','vcdev','uvfm','vvfm','mcdev','mvfm','-v7');
%%
%{
for tt = 1:cfm.s(end)
    close all;
    plot_velprofile(cfm.x*1000,flipud(cfm.y*1000),upsi(:,:,tt)*1000,vpsi(:,:,tt)*1000,...
        cfm.mask,vsscdev(:,:,tt),uint8(mean(cfm.imgray,4)))
    pause(1E-2);
    xlabel('Position (mm)','FontSize',20,'FontWeight','bold');
    ylabel('Position (mm)','FontSize',20,'FontWeight','bold');
    set(gcf,'Color',[1 1 1]);
    set(gca,'FontSize',20,'FontWeight','bold');
    set(gcf,'Position',[100 100 500 300]);
    FIG = get(gcf);
%     export_fig(FIG,['/Volumes/GERI/cDEV/cDEV_paper/new_images/mouse_carotid_new/profile_vss_cdev_',...
%         num2str(tt,'%02i')],'-tif','-a1','-r128','-q400');
    
    plot_velprofile(cfm.x*1000,flipud(cfm.y*1000),uvfm(:,:,tt)*1000,vvfm(:,:,tt)*1000,...
        cfm.mask,vssvfm(:,:,tt),uint8(mean(cfm.imgray,4)))
    pause(1E-2);
    xlabel('Position (mm)','FontSize',20,'FontWeight','bold');
    ylabel('Position (mm)','FontSize',20,'FontWeight','bold');
    set(gcf,'Color',[1 1 1]);
    set(gca,'FontSize',20,'FontWeight','bold');
    set(gcf,'Position',[100 100 500 300]);
    FIG = get(gcf);
%     export_fig(FIG,['/Volumes/GERI/cDEV/cDEV_paper/new_images/mouse_carotid_new/profile_vss_vfm_',...
%         num2str(tt,'%02i')],'-tif','-a1','-r128','-q400');
    
    pause(1E-1);
end
%}
%%
%{
close all;
plot_velprofile(cfm.x*1000,flipud(cfm.y*1000),mean(upsi,3)*1000,mean(vpsi,3)*1000,...
    cfm.mask,mean(vsscdev,3),uint8(mean(cfm.imgray,4)))
pause(1E-2);
xlabel('Position (mm)','FontSize',20,'FontWeight','bold');
ylabel('Position (mm)','FontSize',20,'FontWeight','bold');
set(gcf,'Color',[1 1 1]);
set(gca,'FontSize',20,'FontWeight','bold');
set(gcf,'Position',[100 100 500 300]);
FIG = get(gcf);
export_fig(FIG,'D:\Reza\VascularDoVeR\Original\results\mouse\comparison\profile_vss_cdev_mean',...
    '-tif','-a1','-r128','-q400');

plot_velprofile(cfm.x*1000,flipud(cfm.y*1000),mean(uvfm,3)*1000,mean(vvfm,3)*1000,...
    cfm.mask,mean(vssvfm,3),uint8(mean(cfm.imgray,4)))
pause(1E-2);
xlabel('Position (mm)','FontSize',20,'FontWeight','bold');
ylabel('Position (mm)','FontSize',20,'FontWeight','bold');
set(gcf,'Color',[1 1 1]);
set(gca,'FontSize',20,'FontWeight','bold');
set(gcf,'Position',[100 100 500 300]);
FIG = get(gcf);
export_fig(FIG,'D:\Reza\VascularDoVeR\Original\results\mouse\comparison\profile_vss_vfm_mean',...
    '-tif','-a1','-r128','-q400');
%}
%% WALL SHEAR STRESS ON MOUSE DATA
% Input "static" variables first
info.rho    = 1100;     % Density, kg/m^3
info.nu     = 3.77e-6;  % Viscosity, m^2/s
info.mu     = info.rho*info.nu;
info.wres(1) = 3;        % x grid points for RBF wall detection 7 is good
info.wres(2) = 3;        % y grid points for RBF wall detection 7 is good
info.gres   = [1,1];
info.nwalls = 3;         % number of walls for corresponding ROIs
info.methd  = 'Edeform'; % 'timeavg''Edeform''downsample'
info.grad_method = 'rbf'; % either 'lsq' (Least Squares) or 'rbf' (radial basis)

%%
%{
cfm.S = size(cfm.vmap);
cfm.MASK = imerode(cfm.mask,strel('disk',1));
cfm.xpos = cfm.x; cfm.ypos = cfm.y;
% Set CFM BC Mask
cfm.BCMASK = cfm.MASK*7;
% Compute zero-crossings, label all as 1
bnds    = cell2mat(bwboundaries(logical(cfm.MASK)));
inds    = sub2ind([size(cfm.MASK,1) size(cfm.MASK,2)],bnds(:,1),bnds(:,2));
keyboard
cfm.BCMASK(inds) = 4;
% find main vessel inlet points, label as 1
inds = find(cfm.BCMASK(:,423) == 4);
cfm.BCMASK(inds(1:end-0),423) = 1;
% find bottom wall of vessel, label as 2
for nn = 145:423%REZA%53:312
    inds = find(cfm.BCMASK(:,nn) == 4);
    cfm.BCMASK(inds(end),nn) = 2;
end
% label interna outlet as 3
cfm.BCMASK(62:91,91) = 3;
% label externa oulet as 5
cfm.BCMASK(123:148,120) = 5;
% find bottom wall of vessel, label as 2
for nn = 91:423%23:312
    inds = find(cfm.BCMASK(:,nn) == 4);
    cfm.BCMASK(inds(1),nn) = 6;
end
%%

[temp.topy,temp.topx] = find(cfm.BCMASK == 2);
[temp.boty,temp.botx] = find(cfm.BCMASK == 6);
[temp.bncy,temp.bncx] = find(cfm.BCMASK == 4);

cfm.top = [temp.topy(:),temp.topx(:)];
cfm.bot = [temp.boty(:),temp.botx(:)];
cfm.bnc = [temp.bncy(:),temp.bncx(:)];

cfm.topp = [];
cfm.botp = [];
cfm.bncp = [];
for n = 1:numel(temp.topy)
    cfm.topp(n,:) = [cfm.y(temp.topy(n),temp.topx(n)),cfm.x(temp.topy(n),temp.topx(n))];
end
for n = 1:numel(temp.boty)
    cfm.botp(n,:) = [cfm.y(temp.boty(n),temp.botx(n)),cfm.x(temp.boty(n),temp.botx(n))];
end
for n = 1:numel(temp.bncy)
    cfm.bncp(n,:) = [cfm.y(temp.bncy(n),temp.bncx(n)),cfm.x(temp.bncy(n),temp.bncx(n))];
end

temp.bnc = [];
% Find smallest value on y
[~,ind] = min(cfm.bnc(:,2));
I = ind(1);
temp.bnc(1,:) = cfm.bnc(ind,:);
cfm.bnc(ind(1),:) = [];
count = 2;
while numel(cfm.bnc) > 0
    temp.R  = sqrt((cfm.bnc(:,1)-temp.bnc(count-1,1)).^2+(cfm.bnc(:,2)-temp.bnc(count-1,2)).^2);
    [~,ind] = min(abs(temp.R));
    temp.bnc(count,:) = cfm.bnc(ind(1),:);
    cfm.bnc(ind(1),:)  = [];
    count = count + 1;
    I = cat(1,I,ind(1));
end

cfm.bnc = temp.bnc;
cfm.bncp= cfm.bncp(I,:);


%%
for t = 1:cfm.s(end)
    fprintf('Calculating WSS on plane %02i \r',t)
    % Place wall locations into a cell array for plotting wall boundaries
    cfm.wallpos  = {cfm.topp,cfm.botp,cfm.bncp};
    cfm.wallinds = {cfm.top ,cfm.bot ,cfm.bnc };
    
    % Load in shear function
    info.shear_function = @findwallgradients_tps3;
    
    % Spline fit to walls (subpixel interpolation of wall positions)
    for n = 1:numel(cfm.wallinds)
        if n == numel(cfm.wallinds)
            temp.x = cfm.wallinds{n}(:,1);
            temp.xp=  cfm.wallpos{n}(:,1);
            temp.y = cfm.wallinds{n}(:,2);
            temp.yp=  cfm.wallpos{n}(:,2);
        else
            temp.x = cfm.wallinds{n}(:,2);
            temp.xp=  cfm.wallpos{n}(:,2);
            temp.y = cfm.wallinds{n}(:,1);
            temp.yp=  cfm.wallpos{n}(:,1);
        end
        
        if n == numel(cfm.wallinds)
            temp.Xs = temp.x-mean([temp.x(1) temp.x(end)]);
            temp.Ys = temp.y-mean([temp.y(1) temp.y(end)]);
            [temp.tht,temp.rho] = cart2pol(temp.Xs,temp.Ys);
            info.spln       = splinefit(temp.tht,temp.rho,15,'r');
            temp.rhospln    = fnval(info.spln,temp.tht)';                  % spline values still in pixels
            [temp.xtf,temp.ytf] = pol2cart(temp.tht,temp.rhospln');
            temp.x      = movmean(temp.xtf + mean([temp.x(1) temp.x(end)]),5);
            temp.yspln  = movmean(temp.ytf + mean([temp.y(1) temp.y(end)]),5);
        else
            info.spln   = splinefit(temp.x,temp.y,15,'r');
            temp.yspln  = fnval(info.spln,temp.x)';                  % spline values still in pixels
        end

        temp.slope = gradient(temp.yspln(:))./gradient(temp.x(:)); % Find Slope along curve
        
        if n == numel(cfm.wallinds)
            temp.theta = atan2(temp.slope',1);
        else
            temp.theta = atan2(temp.slope',1);
        end
        
        cfm.wallftr{n} = [temp.slope(:),temp.theta(:)];
        
        cfm.wallind{n} = [];
        cfm.wallpts{n} = [];
        if n == numel(cfm.wallinds)
            cfm.wallind{n}(:,1) = temp.x;
            cfm.wallind{n}(:,2) = temp.yspln;
            
            cfm.wallpts{n}(:,1) = temp.xp;
            cfm.wallpts{n}(:,2) = temp.yspln*cfm.dx+min(cfm.x(:))-cfm.dx;
        else
            cfm.wallind{n}(:,1) = temp.yspln;
            cfm.wallind{n}(:,2) = temp.x;
            
            cfm.wallpts{n}(:,1) = temp.yspln*cfm.dy+min(cfm.y(:))-cfm.dy;
            cfm.wallpts{n}(:,2) = temp.xp;
        end
    end
%     keyboard
    % Generate Grid (integer values)
    [cfm.X,cfm.Y] = meshgrid(1:cfm.s(2),1:cfm.s(1));
    
    % Calculate Wall Shear Stress
    [cfm.VSScdev, cfm.WSScdev, cfm.vsscdev, cfm.UXcdev, cfm.UYcdev, cfm.VXcdev, cfm.VYcdev] = shear_calc(cfm.x,cfm.y,cfm.X,cfm.Y,...
        UsDoVeR(:,:,t),VsDoVeR(:,:,t),info.nu,info.rho,info.wres,cfm.wallinds,cfm.wallftr,info.grad_method);
    [cfm.VSSvfm, cfm.WSSvfm, cfm.vssvfm, cfm.UXvfm, cfm.UYvfm, cfm.VXvfm, cfm.VYvfm] = shear_calc(cfm.x,cfm.y,cfm.X,cfm.Y,...
        UsiVFM(:,:,t),VsiVFM(:,:,t),info.nu,info.rho,info.wres,cfm.wallinds,cfm.wallftr,info.grad_method);
    
    day(t).wsscdev  = cfm.VSScdev;
    day(t).wssvfm   = cfm.VSSvfm;
    day(t).VSScdev  = cfm.vsscdev;
    day(t).VSSvfm   = cfm.vssvfm;
    day(t).xpos     = cfm.x;
    day(t).ypos     = cfm.y;
    day(t).ucdev    = UsDoVeR(:,:,t);
    day(t).vcdev    = VsDoVeR(:,:,t);
    day(t).uvfm     = UsiVFM(:,:,t);
    day(t).vvfm     = VsiVFM(:,:,t);
    day(t).mask     = cfm.mask;
    day(t).wallpts    = cfm.wallpts;
    day(t).wallind    = cfm.wallind;

    day(t).uxcdev      = cfm.UXcdev;
    day(t).vxcdev      = cfm.VXcdev;
    day(t).uxvfm       = cfm.UXvfm;
    day(t).vxvfm       = cfm.VXvfm;
    
    day(t).uycdev      = cfm.UYcdev;
    day(t).vycdev      = cfm.VYcdev;
    day(t).uyvfm       = cfm.UYvfm;
    day(t).vyvfm       = cfm.VYvfm;
    
    day(t).WSScdev = cfm.WSScdev;
    day(t).WSSvfm = cfm.WSSvfm;
end

%%
for n = 1:3
    if n == 3
        [arclen,seglen] = arclength(cfm.wallpts{n}(:,1),cfm.wallpts{n}(:,2));
        al_bnc = cat(1,0,cumsum(seglen))/8;
    elseif n == 2
        [arclen,seglen] = arclength(cfm.wallpts{n}(:,1),cfm.wallpts{n}(:,2));
        al_bot = cat(1,0,cumsum(seglen));
    else
        [arclen,seglen] = arclength(cfm.wallpts{n}(:,1),cfm.wallpts{n}(:,2));
        al_top = cat(1,0,cumsum(seglen))/2;
    end
end
%%
for t = 1:13
    for n = 3
        if n == 1
            alen = al_top;
        elseif n == 2
            alen = al_bot;
        elseif n == 3
            alen = al_bnc;
        end
    end
    if t == 1
        wssDoVeRMean = -day(t).WSScdev{n}(1,5:5:end)/13;
        wssiVFMMean = -day(t).WSSvfm{n}(1,5:5:end)/13;
    else
        wssDoVeRMean = wssDoVeRMean-day(t).WSScdev{n}(1,5:5:end)/13;
        wssiVFMMean = wssiVFMMean-day(t).WSSvfm{n}(1,5:5:end)/13;
    end
end
%%
% Compute pressures
% Set physical parameters that will be used to non-dimensionalize
% velocities for reconstruction-based pressure estimations
mu  = 3.5E-3;   % Units in Pa*s, or N*s/m^2
rho = 1060;     % Units of kg/m^3
nu  = mu/rho;   % Units of m^2/s

% Set characteristic velocity, length, and time-scales
U0  = 1;
L0  = 1;
T0  = L0/U0;
Re  = L0*U0/nu;
P0  = rho*U0^2;

% Run POD (use ELF method) on DoVeR processing results
[UsDoVeR,VsDoVeR,~,~] = POD_MCD(upsi,vpsi,4,0.9,0.1);

% Run POD (use ELF method) on DoVeR processing results
[UsiVFM ,VsiVFM ,~,~] = POD_MCD(uvfm,vvfm,4,0.9,0.1);

% Non-dimensionalize velocity, length, and time-scales
UnDoVeR     = UsDoVeR/U0.*cfm.mask;
VnDoVeR     = VsDoVeR/U0.*cfm.mask;
UniVFM      = UsiVFM/U0.*cfm.mask;
VniVFM      = VsiVFM/U0.*cfm.mask;
Xn          = cfm.x/L0;
Yn          = cfm.y/L0;
dtn         = 1E-2/T0;

% Find locations for mask in order to truncate pressure estimate region
[r,c] = find(max(cfm.mask,[],3)~=0);
minr = min(r)-1;
maxr = max(r)+1;
minc = min(c)-1;
maxc = max(c)+1;

% Run pressure estimation on DoVeR velocity information
pDoVeR = pressure_lines_9_3(Xn(minr:maxr,minc:maxc)',Yn(minr:maxr,minc:maxc)',...
    permute(UnDoVeR(minr:maxr,minc:maxc,:),[2 1 3]),...
    permute(VnDoVeR(minr:maxr,minc:maxc,:),[2 1 3]),...
    Re,dtn,'standard',permute(cfm.mask(minr:maxr,minc:maxc,:),[2 1 3]));
PDoVeR = zeros(size(UnDoVeR));
PDoVeR(minr:maxr,minc:maxc,:) = P0/133.3*permute(pDoVeR,[2 1 3]);

% Run pressure estimation on iVFM velocity information
piVFM = pressure_lines_9_3(Xn(minr:maxr,minc:maxc)',Yn(minr:maxr,minc:maxc)',...
    permute(UniVFM(minr:maxr,minc:maxc,:),[2 1 3]),...
    permute(VniVFM(minr:maxr,minc:maxc,:),[2 1 3]),...
    Re,dtn,'standard',permute(cfm.mask(minr:maxr,minc:maxc,:),[2 1 3]));
PiVFM = zeros(size(UniVFM));
PiVFM(minr:maxr,minc:maxc,:) = P0/133.3*permute(piVFM,[2 1 3]);
%%
cfm.MASK(cfm.MASK == 0) = nan;
figure(101);
h1 = pcolor((cfm.x(1,minc:maxc)-mean(cfm.x(1,minc:maxc)))*100,...
    (cfm.y(minr:maxr,1)-mean(cfm.y(minr:maxr,1)))*100,...
    (PDoVeR(minr:maxr,minc:maxc,4)-mean(PDoVeR(26:61,312,4))).*cfm.MASK(minr:maxr,minc:maxc));
set(h1, 'EdgeColor', 'none');
caxis([-1 1])
set(gca,'FontSize',36,'FontWeight','bold','YDir','Normal','XTick',-0.1:0.05:0.1);
set(gcf,'Color',[1 1 1],'Position',[100 100 900 400]);
ylabel('Position (mm)','FontSize',36,'FontWeight','bold');
xlabel('Position (mm)','FontSize',36,'FontWeight','bold');
h = colorbar;
ylabel(h,'Pressure (mmHg)','FontSize',36,'FontWeight','bold');
export_fig(gcf, ['D:\Reza\VascularDoVeR\Original\results\mouse\comparison\mouse_carotid_pressure_peaksystole_DoVeR.tiff'],...
            '-tif','-a1','-r128','-q400');

figure(102);
h1 = pcolor((cfm.x(1,minc:maxc)-mean(cfm.x(1,minc:maxc)))*100,...
    (cfm.y(minr:maxr,1)-mean(cfm.y(minr:maxr,1)))*100,...
    (PiVFM(minr:maxr,minc:maxc,4)-mean(PiVFM(26:61,312,4))).*cfm.MASK(minr:maxr,minc:maxc));
set(h1, 'EdgeColor', 'none');
caxis([-1 1])
set(gca,'FontSize',36,'FontWeight','bold','YDir','Normal','XTick',-0.1:0.05:0.1);
set(gcf,'Color',[1 1 1],'Position',[100 100 900 400]);
ylabel('Position (mm)','FontSize',36,'FontWeight','bold');
xlabel('Position (mm)','FontSize',36,'FontWeight','bold');
h = colorbar;
ylabel(h,'Pressure (mmHg)','FontSize',36,'FontWeight','bold');
export_fig(gcf, ['D:\Reza\VascularDoVeR\Original\results\mouse\comparison\mouse_carotid_pressure_peaksystole_iVFM.tiff'],...
            '-tif','-a1','-r128','-q400');
        
figure(103);
h1 = pcolor((cfm.x(1,minc:maxc)-mean(cfm.x(1,minc:maxc)))*100,...
    (cfm.y(minr:maxr,1)-mean(cfm.y(minr:maxr,1)))*100,...
    (PDoVeR(minr:maxr,minc:maxc,8)-mean(PDoVeR(26:61,312,8))).*cfm.MASK(minr:maxr,minc:maxc));
set(h1, 'EdgeColor', 'none');
caxis([-1 1])
set(gca,'FontSize',36,'FontWeight','bold','YDir','Normal','XTick',-0.1:0.05:0.1);
set(gcf,'Color',[1 1 1],'Position',[100 100 900 400]);
ylabel('Position (mm)','FontSize',36,'FontWeight','bold');
xlabel('Position (mm)','FontSize',36,'FontWeight','bold');
h = colorbar;
ylabel(h,'Pressure (mmHg)','FontSize',36,'FontWeight','bold');
export_fig(gcf, ['D:\Reza\VascularDoVeR\Original\results\mouse\comparison\mouse_carotid_pressure_peakdiastole_DoVeR.tiff'],...
            '-tif','-a1','-r128','-q400');

figure(104);
h1 = pcolor((cfm.x(1,minc:maxc)-mean(cfm.x(1,minc:maxc)))*100,...
    (cfm.y(minr:maxr,1)-mean(cfm.y(minr:maxr,1)))*100,...
    (PiVFM(minr:maxr,minc:maxc,8)-mean(PiVFM(26:61,312,8))).*cfm.MASK(minr:maxr,minc:maxc));
set(h1, 'EdgeColor', 'none');
caxis([-1 1])
set(gca,'FontSize',36,'FontWeight','bold','YDir','Normal','XTick',-0.1:0.05:0.1);
set(gcf,'Color',[1 1 1],'Position',[100 100 900 400]);
ylabel('Position (mm)','FontSize',36,'FontWeight','bold');
xlabel('Position (mm)','FontSize',36,'FontWeight','bold');
h = colorbar;
ylabel(h,'Pressure (mmHg)','FontSize',36,'FontWeight','bold');
export_fig(gcf, ['D:\Reza\VascularDoVeR\Original\results\mouse\comparison\mouse_carotid_pressure_peakdiastole_iVFM.tiff'],...
            '-tif','-a1','-r128','-q400');
%%
figure(105);
plot((1E-2:1E-2:13*1E-2)*1E3,squeeze(mean(PDoVeR(109:137,53,:),1)-mean(PDoVeR(26:61,312,:),1)),...
    '-o','LineWidth',1.6,'MarkerEdgeColor',[0.4 0.4 0.4],...
    'Color',[0.4 0.4 0.4],'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',8);
hold on;
plot((1E-2:1E-2:13*1E-2)*1E3,squeeze(mean(PiVFM(109:137,53,:),1)-mean(PiVFM(26:61,312,:),1)),...
    '->','LineWidth',1.6,'MarkerEdgeColor',[0.0 0.0 0.0],...
    'Color',[0.0 0.0 0.0],'MarkerFaceColor',[0.0 0.0 0.0],'MarkerSize',8);
hold off;
legend('DoVer','iVFM');
grid on;
set(gca,'LineWidth',2,'FontSize',28,'FontWeight','bold');
axis([0 140 -4 4]);
ylabel('\DeltaP_{Internal - Main} (mmHg)');
xlabel('Time (ms)');
set(gcf,'Color',[1 1 1],'Position',[100 100 900 400]);
export_fig(gcf, ['D:\Reza\VascularDoVeR\Original\results\mouse\comparison\mouse_carotid_interna_diff_pressure.tiff'],...
            '-tif','-a1','-r128','-q400');

figure(106);
plot((1E-2:1E-2:13*1E-2)*1E3,squeeze(mean(PDoVeR(51:78,24,:),1)-mean(PDoVeR(26:61,312,:),1)),...
    '-o','LineWidth',1.6,'MarkerEdgeColor',[0.4 0.4 0.4],...
    'Color',[0.4 0.4 0.4],'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',8);
hold on;
plot((1E-2:1E-2:13*1E-2)*1E3,squeeze(mean(PiVFM(51:78,24,:),1)-mean(PiVFM(26:61,312,:),1)),...
    '->','LineWidth',1.6,'MarkerEdgeColor',[0.0 0.0 0.0],...
    'Color',[0.0 0.0 0.0],'MarkerFaceColor',[0.0 0.0 0.0],'MarkerSize',8);
hold off;
legend('DoVer','iVFM');
grid on;
set(gca,'LineWidth',2,'FontSize',28,'FontWeight','bold');
axis([0 140 -4 4]);
ylabel('\DeltaP_{External - Main} (mmHg)');
xlabel('Time (ms)');
set(gcf,'Color',[1 1 1],'Position',[100 100 900 400]);
export_fig(gcf, ['D:\Reza\VascularDoVeR\Original\results\mouse\comparison\mouse_carotid_externa_diff_pressure.tiff'],...
            '-tif','-a1','-r128','-q400');
        %%
        field = day;
        nnn = [2 1 3];
% Plot WSS Time Profiles
for nn = 1:3
    clear WSScdev WSSvfm
    
    position = nnn(nn);
    srslng = zeros(1);
    for t = 1:13
        srslng(t) = numel(field(t).WSScdev{position});
    end
    minlng  = min(srslng);

    WSScfd  = [];
    WSScdev = [];
    WSSvfm  = [];
    for t = 1:13
        Xvect = 1:numel(field(t).WSScdev{position});
        xvect = linspace(1,numel(field(t).WSScdev{position}),minlng);
        WSScdev(:,t)    = interp1(Xvect,field(t).WSScdev{position},xvect);
        WSSvfm(:,t)     = interp1(Xvect,field(t).WSSvfm{position},xvect);
    end
    
    if position == 3
        al = sqrt((field(t).wallind{position}(:,1)-field(t).wallind{position}(1,1)).^2 + ...
            (field(t).wallind{position}(:,2)-field(t).wallind{position}(1,2)).^2);
        diff_al = diff(al);
        arclength = cumsum(al/30);
    else
        arclength = sqrt((field(t).wallind{position}(:,1)-field(t).wallind{position}(1,1)).^2 + ...
            (field(t).wallind{position}(:,2)-field(t).wallind{position}(1,2)).^2);
    end

    figure(11);
    subplot(3,1,nn)
    plot(arclength(1:minlng)/arclength(minlng),smooth(mean(WSScdev,2),0.03,'loess'),'-','LineWidth',2,'Color',[0.8500    0.3250    0.0980])
    hold on;
    plot(arclength(1:minlng)/arclength(minlng),smooth(mean(WSSvfm,2),0.03,'loess'),'-','LineWidth',2,'Color',[0.9290    0.6940    0.1250])
    hold off;
    
    set(gca,'FontSize',20,'FontWeight','bold');
    ylabel('WSS, dynes/cm^2','FontSize',20,'FontWeight','bold');
    xlabel('Normalized wall position','FontSize',20,'FontWeight','bold');
    ylim([-40 40]);
    xlim([0 1]);
end
set(gcf,'Position',[10 10 400 800],'Color',[1 1 1]);

subplot(3,1,2)
lgd = legend({'DoVeR','iVFM'},'FontSize',16,'TextColor','black','Location','northwest');
export_fig(gcf,['D:\Reza\VascularDoVeR\Original\results\mouse\comparison\mouse_time_averaged_wss.tiff'],'-tif','-a1','-r128','-q400');

%}

%% Common Figure Settings
plot_settings = {'LineWidth', 2, 'MarkerSize', 8};
legend_settings = {'FontSize', 18, 'Location', 'best', 'Box', 'off', 'FontName', 'Times New Roman'};
axis_settings = {'FontSize', 16, 'FontWeight', 'Bold', 'FontName', 'Times New Roman', ...
                 'LineWidth', 1.5, 'TickDir', 'out', 'XMinorTick', 'on', 'YMinorTick', 'on'};
xlabel_settings = {'FontSize', 25, 'FontWeight', 'Bold', 'FontName', 'Times New Roman'};
ylabel_settings = {'FontSize', 25, 'FontWeight', 'Bold', 'FontName', 'Times New Roman'};
title_settings = {'FontSize', 22, 'FontWeight', 'Bold', 'FontName', 'Times New Roman'};
grid_settings = {'GridLineStyle', '--', 'Color', [1,1,1]};


%% Additional Velocity Field Comparison Plot
xp = cfm.x;
yp = cfm.y;
selected_time_point = size(Dov_u,3);

for i = 1:selected_time_point
    % Figure and subplot settings
    figure;
    
    % CFD Velocity Field with Colormap
    subplot(1, 3, 1);
    % Plot velocity magnitude as a colormap
    contourf(xp(:,:), yp(:,:), sqrt(vfm_u(:,:,selected_time_point).^2 + vfm_v(:,:,selected_time_point).^2), 'LineColor', 'none');
    hold on;
    % Plot quiver with reduced density and black color
    quiver(xp(1:5:end, 1:5:end), yp(1:5:end, 1:5:end), vfm_u(1:5:end, 1:5:end, selected_time_point), -vfm_v(1:5:end, 1:5:end, selected_time_point), 'k');
    hold off;
    % title('CFD', title_settings{:});
    xlabel('X Position', xlabel_settings{:});
    ylabel('Y Position', ylabel_settings{:});
    set(gca, axis_settings{:});
    axis equal;
    grid on;
    colorbar; % Add colorbar for velocity magnitude
    
    % DoVeR Velocity Field with Colormap
    subplot(1, 3, 2);
    % Plot velocity magnitude as a colormap
    contourf(xp(:,:), yp(:,:), sqrt(Dov_u(:,:,selected_time_point).^2 + Dov_v(:,:,selected_time_point).^2), 'LineColor', 'none');
    hold on;
    % Plot quiver with reduced density and black color
    quiver(xp(1:5:end, 1:5:end), yp(1:5:end, 1:5:end), Dov_u(1:5:end, 1:5:end, selected_time_point), Dov_v(1:5:end, 1:5:end, selected_time_point), 'k');
    hold off;
    % title('DoVeR', title_settings{:});
    xlabel('X Position', xlabel_settings{:});
    ylabel('Y Position', ylabel_settings{:});
    set(gca, axis_settings{:});
    axis equal;
    grid on;
    colorbar; % Add colorbar for velocity magnitude
    
    % iVFM Velocity Field with Colormap
    subplot(1, 3, 3);
    % Plot velocity magnitude as a colormap
    contourf(xp(:,:), yp(:,:), sqrt(ivfm_u(:,:,selected_time_point).^2 + ivfm_v(:,:,selected_time_point).^2), 'LineColor', 'none');
    hold on;
    % Plot quiver with reduced density and black color
    quiver(xp(1:5:end, 1:5:end), yp(1:5:end, 1:5:end), ivfm_u(1:5:end, 1:5:end, selected_time_point), ivfm_v(1:5:end, 1:5:end, selected_time_point), 'k');
    hold off;
    % title('iVFM', title_settings{:});
    xlabel('X Position', xlabel_settings{:});
    ylabel('Y Position', ylabel_settings{:});
    set(gca, axis_settings{:});
    axis equal;
    grid on;
    colorbar; % Add colorbar for velocity magnitude
    
    set(gcf, 'Color', 'w'); % Set figure background to white
    % saveas(gcf, fullfile(plotdir, sprintf('velocity_field_comparison_time_%02d.png', selected_time_point)));

end
%% Demo 
%%
% FIND VORTEX FEATURES USING COHERENT STRUCTURE ID FOR ALL FRAME
% DoVeR RESULTS
% calDoVeR    =   zeros(size(vpsi));
% omegaDoVeR  =   zeros(size(vpsi));
% for nn = 1:data.size(end)
%     fprintf(['Now Evaluating Frame ',num2str(nn),' of ',num2str(data.size(end)),' for DoVeR \n']);
%     [calDoVeR(:,:,nn),omegaDoVeR(:,:,nn)] = coherent_structure_id(data.x,...
%         data.y, -upsi(:,:,nn), vpsi(:,:,nn),...
%         data.dx, data.dy, data.mask, 'LambdaCI', 7, 7);
% end
% % VFM RESULTS
% calVFM      =   zeros(size(vpsi));
% omegaVFM    =   zeros(size(vpsi));
% for nn = 1:data.size(end)
%     fprintf(['Now Evaluating Frame ',num2str(nn),' of ',num2str(data.size(end)),' for VFM \n']);
%     [calVFM(:,:,nn),omegaVFM(:,:,nn)] = coherent_structure_id(data.x,...
%         data.y, -uvfm(:,:,nn), vvfm(:,:,nn),...
%         data.dx, data.dy, data.mask, 'LambdaCI', 7, 7);
% end
% iVFM RESULTS
caliVFM     =   zeros(size(vpsi));
omegaiVFM   =   zeros(size(vpsi));
for nn = 1:data.size(end)
    fprintf(['Now Evaluating Frame ',num2str(nn),' of ',num2str(data.size(end)),' for iVFM \n']);
    [caliVFM(:,:,nn),omegaiVFM(:,:,nn)] = coherent_structure_id(data.x,...
        data.y, -ivfm_u(:,:,nn), ivfm_v(:,:,nn),...
        data.dx, data.dy, data.mask, 'LambdaCI', 7, 7);
end
%%
% SET COLOR MAP TO RED BLUE
rbmap   = colormap(redblue(40));
for nn = 1:13%4%1:data.size(end)
    figure(1);
    imshow(uint8(data.imreduced(:,:,:,nn)));
    axis image;
    export_fig(1,fullfile(plotdir,sprintf('rep_image_%03i.tiff',nn)),'-tiff','-a1','-r512');
    
    figure(5);
    imshow(uint8(data.img(260:610,1:85,:,nn)));
    axis image;
    export_fig(5,fullfile(plotdir,sprintf('rep_cscale_%03i.tiff',nn)),'-tiff','-a1','-r512');
    
    figure(2);
    toDoVeR = omegaDoVeR(:,:,nn); toDoVeR(data.mask == 0) = nan;
    hold on;
    imagesc(data.x(1,:,1)*1000,data.y(:,1,1)*1000,uint8(data.imgray(:,:,:,nn)));
    imagesc(data.x(1,1:65,1)*1000,data.y(1:35,1,1)*1000,255*ones([data.size(1)-35 data.size(2)-65 data.size(3)]));
    h = pcolor(data.x(1,:,1)*1000,data.y(:,1,1)*1000,toDoVeR);
    set(h, 'EdgeColor', 'none');
    colormap(rbmap); caxis([-1000,1000]);
%     imcontour(data.x(:,:,1)*1000,data.y(:,:,1)*1000,calDoVeR(:,:,nn),[-0.3 0.3],...
%         'Color',[0.75, 0, 0.75],'LineStyle','-','LineWidth',2);
%     imcontour(data.x(:,:,1)*1000,data.y(:,:,1)*1000,data.mask,[-0.05 0.05],...
%         'Color',[0 0 0],'LineWidth',1.4);
    quiver(data.x(12:8:end,12:24:end,1)*1000,data.y(12:8:end,12:24:end,1)*1000,...
        upsi(12:8:end,12:24:end,nn),-vpsi(12:8:end,12:24:end,nn).*data.mask(12:8:end,12:24:end),...
        'Color',[0 0 0],'autoscale','off','LineWidth',0.75,'MaxHeadSize',5);
    hold off;
    set(gca,'FontSize',28,'FontWeight','bold','YDir','normal');
    set(gcf,'Position',[400 100 350 500],'Color',[1 1 1]);
    xlabel('Position (cm)');
    ylabel('Position (cm)');
    axis image;
%     axis([0.75 4.5 1 5]);
    axis off;
    hold on;
    quiver(0.05,1.75,0.2,0,'Color',[0 0 0],'autoscale','off','LineWidth',0.75,'MaxHeadSize',5);
    annotation('textbox',[0.12 0.51,0.1 0.1],'String','0.4 m/s',...
    'EdgeColor','none','FontSize',10,'FontWeight','bold');
    hold off;
    export_fig(2,fullfile(plotdir,sprintf('DoVeR_%03i.tiff',nn)),'-tiff','-a1','-r512');
    %     clf

    figure(3);
    toiVFM = omegaVFM(:,:,nn); toiVFM(data.mask == 0) = nan;
    hold on;
    imagesc(data.x(1,:,1)*1000,data.y(:,1,1)*1000,uint8(data.imgray(:,:,:,nn)));
    h = pcolor(data.x(1,:,1)*1000,data.y(:,1,1)*1000,toiVFM);
    set(h, 'EdgeColor', 'none');
    colormap(rbmap); caxis([-1000,1000]);
    hold on;
%     imcontour(data.x(:,:,1)*1000,data.y(:,:,1)*1000,calVFM(:,:,nn),[-0.3 0.3],...
%         'Color',[0.75, 0, 0.75],'LineStyle','-','LineWidth',2);
%     imcontour(data.x(:,:,1)*1000,data.y(:,:,1)*1000,data.mask,[-0.05 0.05],...
%         'Color',[0 0 0],'LineWidth',0.75.4);
    quiver(data.x(12:8:end,12:24:end,1)*1000,data.y(12:8:end,12:24:end,1)*1000,...
        uvfm(12:8:end,12:24:end,nn),-vvfm(12:8:end,12:24:end,nn).*data.mask(12:8:end,12:24:end),...
        'Color',[0 0 0],'autoscale','off','LineWidth',0.75,'MaxHeadSize',4);
    hold off;
    set(gca,'FontSize',28,'FontWeight','bold','YDir','normal');
    set(gcf,'Position',[700 100 350 500],'Color',[1 1 1]);
    xlabel('Position (cm)');
    ylabel('Position (cm)');
        axis image;
%     axis([0.75 4.5 1 5]);
    axis off;
    export_fig(3,fullfile(plotdir,sprintf('VFM_%03i.tiff',nn)),'-tiff','-a1','-r512');
    %     clf

    figure(4);
    toiVFM = omegaiVFM(:,:,nn).*data.mask; toiVFM(data.mask == 0) = nan;
    hold on;
    imagesc(data.x(1,:,1)*1000,data.y(:,1,1)*1000,uint8(data.imgray(:,:,:,nn)));
    h = pcolor(data.x(1,:,1)*1000,data.y(:,1,1)*1000,toiVFM);
%     h = pcolor(data.x(1,:,1)*1000,data.y(:,1,1)*1000,toiVFM);
    set(h, 'EdgeColor', 'none');
    colormap(rbmap); caxis([-1000,1000]);
    hold on;
%     imcontour(data.x(:,:,1)*1000,data.y(:,:,1)*1000,caliVFM(:,:,nn).*data.mask,[-0.3 0.3],...
%         'Color',[0.75, 0, 0.75],'LineStyle','-','LineWidth',2);
%     imcontour(data.x(:,:,1)*1000,data.y(:,:,1)*1000,data.mask.*data.mask,[-0.05 0.05],...
%         'Color',[0 0 0],'LineWidth',0.75.4);
    quiver(data.x(12:8:end,12:24:end,1)*1000,data.y(12:8:end,12:24:end,1)*1000,...
        ivfm_u(12:8:end,12:24:end,nn).*sin(TXY(12:8:end,12:24:end)).*data.mask(12:8:end,12:24:end)+...
        -ivfm_u(12:8:end,12:24:end,nn).*cos(TXY(12:8:end,12:24:end)).*data.mask(12:8:end,12:24:end),...
        -ivfm_u(12:8:end,12:24:end,nn).*cos(TXY(12:8:end,12:24:end)).*data.mask(12:8:end,12:24:end)+...
        -ivfm_v(12:8:end,12:24:end,nn).*sin(TXY(12:8:end,12:24:end)).*data.mask(12:8:end,12:24:end),...
        'Color',[0 0 0],'autoscale','off','LineWidth',0.75,'MaxHeadSize',4);
    hold off;
    set(gca,'FontSize',28,'FontWeight','bold','YDir','normal');
    set(gcf,'Position',[1000 100 350 500],'Color',[1 1 1]);
    xlabel('Position (cm)');
    ylabel('Position (cm)');
        axis image;
%     axis([0.75 4.5 1 5]);
    axis off;
    export_fig(4,fullfile(plotdir,sprintf('iVFM_%03i.tiff',nn)),'-tiff','-a1','-r512');
    %     clf
end
%%
figure(6);
rbmap   = colormap(redblue(40));
h = colorbar('eastoutside'); caxis([-1000 1000]);
set(gca,'FontSize',20,'FontWeight','bold','YDir','normal','FontName','Times New Roman');
set(gcf,'Position',[0 0 490 600],'Color',[1 1 1]);
xlabel(h,'Vorticity, [s^{-1}]','FontSize',20,'FontWeight','bold','FontName','Times New Roman');
export_fig(6,fullfile(plotdir,sprintf('carotid_demo_colorbar_%03i.tiff',nn)),'-tiff','-a1','-r512');