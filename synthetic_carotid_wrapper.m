close all; clear; clc;
% keyboard
% Specify directory where all data is stored
basepath = pwd;
dataDir = fullfile(basepath,'analysis','synthetic.carotid');
addpath(basepath,'src','export_fig');
addpath(basepath,'src','discrete_differencing');
addpath(basepath,'src','proper_orthogonal_decomposition');
addpath(basepath,'src','2D_pressure_lines');
addpath(basepath,'src','cDEV');
addpath(basepath,'src','VFM');
addpath(basepath,'src','iVFM');

% Set directory for rawdata
rawDir  = fullfile(dataDir,'rawdata');
% Find all text files in directory
rawInfo = dir(fullfile(rawDir,'*.txt'));
% % Set directory for matfiles
matDir  = fullfile(dataDir,'matfiles');
% Find all text files in directory
matInfo = dir(fullfile(matDir,'*.mat'));
% Read in initial mat-file to know how to interpolate raw data
temp = load(fullfile(matDir,matInfo(2).name));
% Preserve X-axis resolution, interpolate y-axis resolution to match
dx_CFI  = abs((temp.x_axis_B(2)-temp.x_axis_B(1)));
y_axis_CFI = min(temp.r_axis_B):dx_CFI:max(temp.r_axis_B);
clear temp;

% Set time points for each frame sample
frRng = 3;%3:16:190;
vCFI  = [];
for n = 1:numel(frRng)
    fprintf('Evaluating frame %03i of %03i ... \r',n,numel(frRng));
    % Read in current frame mat-file
    temp = load(fullfile(matDir,matInfo(n).name));
    % keyboard
    x_axis_CFI  = temp.x_axis_CFI;
    r_axis_CFI  = temp.r_axis_CFI;
    % Interpolate color flow data to match final image resolution
    % First step, we interpolate to the resolution of the B-mode image
    [x_CFI, z_CFI]  = meshgrid(x_axis_CFI, r_axis_CFI);
    [x_B  , z_B  ]  = meshgrid(temp.x_axis_B, temp.r_axis_B);
    vel_B   = interp2(x_CFI, z_CFI, temp.vel, x_B, z_B);
    CFI     = vel_B;
    % For frame 2, the data must be de-aliased. I'm just using a simple
    % identification and unwrapping scheme
    if n == 2
        bw = CFI;
        bw(bw == 0)  = inf;
        bw(bw < 12)   = 1; %%% why 12?
        bw(bw >= 12)  = 0;
        bw  = imdilate(bwselect(imerode(bw,strel('disk',1)),130,530),strel('disk',1));
        CFI = CFI + bw*60; %%%Why 130 and 530? %%%why multiplied by 60
    end
    % Second step, we interpolate to an isotropic resolution
    [xCFI,yCFI] = meshgrid(temp.x_axis_B,y_axis_CFI);
    vCFI = cat(3,vCFI,medfilt2(interp2(x_B,z_B,CFI,xCFI,yCFI),[3 3])/100);
end

% Initialize counter
count   =   1;
% Pre-allocate memory for CFD u & v velocities
ucfd    =   [];
vcfd    =   [];
for n = 5%3:16:190
    fprintf('Evaluating frame %03i of 200 ... \r',n);
    % Using fscanf, load in each frame, interpolate the pressure data
    FID         = fopen(fullfile(rawDir,rawInfo(n).name));
    datacell    = textscan(FID, '%f%f%f%f%f%f%f','Delimiter',',', 'HeaderLines', 1, 'CollectOutput', 1);
    fclose(FID);
    rawdata     = datacell{1,1};
    % Rotate points around Z-axis
    theta   = -pi/9;    %%%Why?
    rotMat  = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
    rotdata(:,1)    = rawdata(:,1);
    rotdata(:,2:4)  = (rotMat*rawdata(:,2:4)')'*100;
    rotdata(:,5:7)  = (rotMat*rawdata(:,5:7)')'*100;
    rotdata(:,3)    = rotdata(:,3)+2;   %%%why?
    % Interpolate x-velocity & y-velocity on CFI grid (2D)
    UF      =   scatteredInterpolant(rotdata(:,2),rotdata(:,3),rotdata(:,4),...
        rotdata(:,5),'linear','none');
    UCFD    =   UF(xCFI,yCFI,zeros(size(xCFI)));
    UCFD(isnan(UCFD))= 0;
    ucfd    =   cat(3,ucfd,UCFD/100);
    VF      =   scatteredInterpolant(rotdata(:,2),rotdata(:,3),rotdata(:,4),...
        rotdata(:,6),'linear','none');
    VCFD    =   VF(xCFI,yCFI,zeros(size(xCFI)));
    VCFD(isnan(VCFD))= 0;
    vcfd    =   cat(3,vcfd,-VCFD/100);
    count = count + 1;
end

% Create a line vector, from which we will zero out values
% (This was done by hand to try and find the normal vector to the vessel)
yTrunc = 1:120; %%%what is this ytrunc
xTrunc = round(yTrunc*0.3478 + (104.9565));
% Make sure velocities, axes, and scales are set to cm
X   = xCFI/100;
Y   = yCFI/100;
dx  = abs(X(2,2)-X(1,1));
dt  = 1/numel(frRng);    % Frame seperation time, in seconds
% Create mask from CFD data
mask = ucfd.^2+vcfd.^2; mask(mask ~=0) = 1;
bwMask = mask;
for t = 1:numel(frRng)
    for n = 1:numel(yTrunc)
        mask(yTrunc(n),xTrunc(n):end,t) = 0;
    end
end

bcmask = zeros(size(ucfd));
% build boundary mask points
for n = 1:numel(frRng)
    bcMask = mask(:,:,n)*7;     %why 7?
    % Compute zero-crossings, label all as 1
    bnds    = cell2mat(bwboundaries(logical(mask(:,:,n))));
    inds    = sub2ind([size(mask,1) size(mask,2)],bnds(:,1),bnds(:,2));
    bcMask(inds)    = 5;
    % find main vessel inlet points, label as 1
    inds = find(bcMask(:,1) == 5);
    bcMask(inds(1:end-0),1) = 1;
    % find interna outlet points, label as 3
    inds = find(bcMask(:,end) == 5);
    bcMask(inds(1)-1:inds(end)+1,end) = 3;
    % find bottom wall points, label as 2
    [r,c] = find(bcMask(151:end,:) == 5);   %%%Why 151?
    inds  = sub2ind([size(mask,1) size(mask,2)],r+150,c);
    bcMask(inds)    = 2;
    % find bifurcation wall points, label as 4
    [r,c] = find(bcMask(119:143,108:end) == 5);
    inds  = sub2ind([size(mask,1) size(mask,2)],r+118,c+107);
    bcMask(inds)    = 4;
    % find top wall points, label as 6
    [r,c] = find(bcMask(1:150,1:137) == 5);
    inds  = sub2ind([size(mask,1) size(mask,2)],r,c);
    bcMask(inds)    = 6;
    bcmask(:,:,n)   = bcMask;
end

temp.gradientx = zeros(size(bcmask));
temp.gradienty = zeros(size(bcmask));
for nn = 1:numel(frRng)
    fprintf('Generate normals to mask for frame %03i . . . \r',nn);
    % Find Unit Normals and Angles (for imposed velocity directions)
    mthswt = 5;
    temp.mask = mask(:,:,nn);
    temp.bcmask = bcmask(:,:,nn);
    if mthswt == 1
        [temp.mag,temp.grad] = imgradient(temp.mask,'Sobel');
    elseif mthswt == 2
        [temp.mag,temp.grad] = imgradient(temp.mask,'Prewitt');
    elseif mthswt == 3
        [temp.mag,temp.grad] = imgradient(temp.mask,'CentralDifference');
    elseif mthswt == 4
        [temp.mag,temp.grad] = imgradient(temp.mask,'IntermediateDifference');
    elseif mthswt == 5
        [temp.mag,temp.grad] = imgradient(temp.mask,'Roberts');
    end
    temp.mag(temp.bcmask == max(temp.bcmask(:))) = 0;
    temp.grad(temp.bcmask == max(temp.bcmask(:))) = 0;
    temp.mag(temp.bcmask == 0) = 0;
    temp.grad(temp.bcmask == 0) = 0;
    
    temp.gradientx(:,:,nn) = -temp.mag.*cosd(temp.grad);
    temp.gradienty(:,:,nn) =  temp.mag.*sind(temp.grad);
    
    temp.normx = -temp.gradientx(:,:,nn);
    temp.normy = -temp.gradienty(:,:,nn);
    
    temp.len = sqrt(temp.normx .^ 2 + temp.normy .^ 2);
    temp.xnorm(:,:,nn) = temp.normx./temp.len;
    temp.ynorm(:,:,nn) = temp.normy./temp.len;
    
    temp.anglemat(:,:,nn) = atan2d(-temp.gradienty(:,:,nn),temp.gradientx(:,:,nn));
end

fsigm = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)));
psiBC       =   zeros(size(bcmask));
uPholder    =   zeros(size(bcmask));
temp.Uvfm   =   zeros(size(bcmask));
% Compute boundary conditions                               %%% Last Read
for n = 1:numel(frRng)
    % Store current boundary conditions
    psibc       = bcmask(:,:,n);
    % Compute mean velocity, which will be used to compute psi as plug flow
    inds        = find(psibc == 1);
    uMean       = mean(ucfd(inds,1,n));
    uP          = uPholder(:,:,n);
    uVals       = linspace(max(ucfd(inds,1,n)),max(ucfd(inds,1,n)),numel(inds)).*tukeywin(numel(inds),0.05)';
    uP(inds)    = uVals;
    % Integrate and store as plug flow along intlet
    psivectin   = cumtrapz(uVals)*dx;
    psibc(inds) = psivectin;
    % Set bottom wall to last inlet psi-value
    psibc(psibc == 2) = psivectin(end);
    % Interna will have 65% of mass flux outflowing
    % Start by making it a plug boundary
    [inds,~] = find(psibc == 3);
    maxpsi  = psivectin(end);
    tempVals    = (inds-mean(inds))*dx;
    pFit        = fit([tempVals(1);0;tempVals(end)],[0;1*max(ucfd(inds,end,n));0],'poly2');
    uVals       = feval(pFit,tempVals);
    %     uVals       = linspace(max(ucfd(inds,end,n)),max(ucfd(inds,end,n)),numel(inds)).*tukeywin(numel(inds),0.4)';
    psivectout  = cumtrapz(uVals)*dx;
    psivectout  = psivectout + (maxpsi-psivectout(end));
    minpsi      = psivectout(1);
    psibc(psibc == 3) = psivectout;
    uP(inds,end) = uVals;
    % Set bifurcation wall to last value of interna outlet
    psibc(psibc == 4) = minpsi;
    % Externa will have remaining 35% of mass flux outflowing
    % Start by making it a plug boundary
    inds    = find(psibc == 5);
    % Take psivectout, scale to fit other outlet
    tempVals    = linspace(1,numel(psivectout),numel(inds));
    uVals       = interp1(1:numel(psivectout),psivectout,tempVals);
    uVals       = (uVals - min(uVals))*(min(uVals)/max(uVals));
    psibc(psibc == 5) = uVals;
    uP(inds)    = socdiff(uVals,dx,2);
    % Set top wall to 0
    psibc(psibc == 6) = 0;
    % Set psi inside to 0
    psibc(psibc == 7) = 0;
    psiBC(:,:,n) = psibc;
    uPholder(:,:,n) = uP;
    figure(11);
    imagesc(psibc); colorbar;
    pause(1E-1);
end

skp = 8;
% Generate grid information for the iVFM reconstruction method
x           = X - (max(X(:))+min(X(:)))/2; y = -Y;
[TXY,RXY]   = cart2pol(x,y);
[THETA,RHO] = meshgrid(linspace(deg2rad(-145),deg2rad(-35),size(vCFI,2)),...
    linspace(min(RXY(:)),max(RXY(:)),size(vCFI,1)));
[XRT,YRT]   = pol2cart(THETA,RHO);
dr          = RHO(2,2) - RHO(1,1);
dth         = THETA(2,2) - THETA(1,1);
temp.Uvfm   = zeros(size(bcmask));

% Pre-allocate memory for variables
omegaInit   =   zeros(size(bcmask));
ucdev       =   zeros(size(bcmask));
vcdev       =   zeros(size(bcmask));
psi         =   zeros(size(bcmask));
omegapsi    =   zeros(size(bcmask));
uvfm        =   zeros(size(bcmask));
vvfm        =   zeros(size(bcmask));
wvfm        =   zeros(size(bcmask));
uvfm1       =   zeros(size(bcmask));
vvfm1       =   zeros(size(bcmask));
wvfm1       =   zeros(size(bcmask));
uivfm       =   zeros(size(bcmask));
vivfm       =   zeros(size(bcmask));
wivfm       =   zeros(size(bcmask));

% Do vector field reconstruction
for n = 1:numel(frRng)
    % Run median filter on data to reject any outlier vectors
    bwmask  = bwMask(:,:,n); 
    bcMask  = bcmask(:,:,n);
    vcfi    = vCFI(:,:,n);
    vcfi(bcMask < 7) = 0;
    
    % Construct Gaussian filter to smooth data for DoVeR processing
    sigma = 3;
    h1  = fspecial('gaussian',...
        [2*ceil(2*sigma)+1 2*ceil(2*sigma)+1],sigma);
    vdover = vcfi; vdover(vdover == 0) = nan;
    vdover = nanmedfilt2(vdover, [3 3]);
    vdover = nanconv(vdover,h1);
    vdover(isnan(vdover)) = 0;
    dvdx = socdiff_bc(vdover,dx,2,mask(:,:,n));
    % Compute vorticity on smoothed data
    omegaInit(:,:,n) = -2*dvdx;
    
    % Run DoVeR on carotid data
    %%%Reza
    % [ucdev(:,:,n),vcdev(:,:,n),psi(:,:,n),omegapsi(:,:,n)] =...
    %     cdev_cartesian_reza(uPholder(:,:,n).*mask(:,:,n),-vdover.*mask(:,:,n),dx,dx,...
    %     omegaInit(:,:,n),psiBC(:,:,n),bcmask(:,:,n),0.15*max(abs(-vdover(vdover~=0))),'LU');
     [ucdev(:,:,n),vcdev(:,:,n),psi(:,:,n),omegapsi(:,:,n)] =...
        cdev_cartesian_Org(uPholder(:,:,n).*mask(:,:,n),-vdover.*mask(:,:,n),dx,dx,...
        omegaInit(:,:,n),psiBC(:,:,n),bcmask(:,:,n),0.15*max(abs(-vdover(vdover~=0))),'LU');
   
    figure(101);
    imagesc(X(1,skp:skp:end)/dx,Y(skp:skp:end,1)/dx,ucdev(:,:,n)); colorbar; caxis([-0.3 1]);
    hold on;
    quiver(X(skp:skp:end,skp:skp:end)/dx,Y(skp:skp:end,skp:skp:end)/dx,...
        (ucdev(skp:skp:end,skp:skp:end,n))*30,...
        (vcdev(skp:skp:end,skp:skp:end,n))*30,'y','autoscale','off');
    hold off;
    axis image;
    set(gcf,'Position',[100 100 600 600]);
    
    % Construct Gaussian filter to smooth data for VFM processing
    sigma = 3;
    h1 = fspecial('gaussian',...
        [2*ceil(2*sigma)+1 2*ceil(2*sigma)+1],sigma);
    vcfi = imfilter(medfilt2(vcfi,[3 3]),h1);
    dvdy = socdiff_bc(vcfi,dx,1,mask(:,:,n));
    
    [temp.r,temp.c] = find(bcmask(:,:,n) < max(bcmask(:)) &  bcmask(:,:,n) > 0);
    for nn = 1:numel(temp.r)
        temp.Uvfm(temp.r(nn),temp.c(nn),n) =...
            vcfi(temp.r(nn),temp.c(nn))*cosd(-temp.anglemat(temp.r(nn),temp.c(nn),n));
    end
    
    % Run VFM reconstruction with original formulation BCs
    [uvfm(:,:,n),vvfm(:,:,n),wvfm(:,:,n)] = vector_flow_mapping(...
        (temp.Uvfm(:,:,n)+0*uPholder(:,:,n)).*mask(:,:,n),vcfi,dx,dx,dvdy,mask(:,:,n));
    
    figure(102);
    imagesc(X(1,skp:skp:end)/dx,Y(skp:skp:end,1)/dx,uvfm(:,:,n)); colorbar; caxis([-0.3 1]);
    hold on;
    quiver(X(skp:skp:end,skp:skp:end)/dx,Y(skp:skp:end,skp:skp:end)/dx,...
        (uvfm(skp:skp:end,skp:skp:end,n))*30,...
        -(vvfm(skp:skp:end,skp:skp:end,n))*30,'y','autoscale','off');
    hold off;
    axis image;
    set(gcf,'Position',[700 100 600 600]);
    pause(1E-5);
    
    % Run VFM reconstruction with DoVeR formulation BCs
    [uvfm1(:,:,n),vvfm1(:,:,n),wvfm1(:,:,n)] = vector_flow_mapping(...
        (1*temp.Uvfm(:,:,n)+0*uPholder(:,:,n)).*mask(:,:,n),-vcfi,dx,dx,dvdy,mask(:,:,n));
    
    figure(104);
    imagesc(X(1,skp:skp:end)/dx,Y(skp:skp:end,1)/dx,uvfm1(:,:,n)); colorbar; caxis([-0.3 1]);
    hold on;
    quiver(X(skp:skp:end,skp:skp:end)/dx,Y(skp:skp:end,skp:skp:end)/dx,...
        (uvfm1(skp:skp:end,skp:skp:end,n))*30,...
        (vvfm1(skp:skp:end,skp:skp:end,n))*30,'y','autoscale','off');
    hold off;
    axis image;
    set(gcf,'Position',[700 100 600 600]);
    pause(1E-5);
    %%%ReZA
    % % Run iVFM reconstruction (method that uses least squares)
    % [u,v] = ivfm(flip((temp.Uvfm(:,:,n)+uPholder(:,:,n)).*mask(:,:,n),1),flip(-vcfi.*mask(:,:,n),1),...
    %     x,y,XRT,YRT,RXY,TXY,RHO,THETA,dr,dth);
    %New iVFM code
    ang             = -15;
    [uivfm(:,:,n),vivfm(:,:,n)] = ivfm_cartesian_reza_v05(...
    (temp.Uvfm(:,:,n)+0*uPholder(:,:,n)).*bwmask(:,:,n),-vcfi,X,Y,dx,dx,bwmask(:,:,n),ang);

    % uivfm2(:,:,n) = flip(u,1);
    % vivfm2(:,:,n) = flip(v,1);
    wivfm(:,:,n) = socdiff_bc(vivfm(:,:,n),dx,2,mask(:,:,n)) -...
        socdiff_bc(uivfm(:,:,n),dx,1,mask(:,:,n)) ;
    
 end

% Set image output directory
plotdir  = fullfile(basepath,'results','images','synthetic.carotid');
if ~exist(plotdir,'dir')
    mkdir(plotdir)
end 
% CFD GROUND TRUTH
calCFD      =   zeros(size(vcdev));
omegaCFD    =   zeros(size(vcdev));
for nn = 1:numel(frRng)
    fprintf(['Now Evaluating Frame ',num2str(nn),' of ',num2str(20),' for CFD \n']);
    [calCFD(:,:,nn),omegaCFD(:,:,nn)] = coherent_structure_id(X,...
        Y, ucfd(:,:,nn), -vcfd(:,:,nn),dx, dx, mask(:,:,nn), 'LambdaCI', 3, 3);
end
% DoVeR RESULTS
calDoVeR    =   zeros(size(vcdev));
omegaDoVeR  =   zeros(size(vcdev));
for nn = 1:numel(frRng)
    fprintf(['Now Evaluating Frame ',num2str(nn),' of ',num2str(20),' for DoVeR \n']);
    [calDoVeR(:,:,nn),omegaDoVeR(:,:,nn)] = coherent_structure_id(X,...
        Y, ucdev(:,:,nn), -vcdev(:,:,nn),dx, dx, mask(:,:,nn), 'LambdaCI', 3, 3);
end
% VFM RESULTS
calVFM      =   zeros(size(vcdev));
omegaVFM    =   zeros(size(vcdev));
for nn = 1:numel(frRng)
    fprintf(['Now Evaluating Frame ',num2str(nn),' of ',num2str(20),' for VFM \n']);
    [calVFM(:,:,nn),omegaVFM(:,:,nn)] = coherent_structure_id(X,...
        Y, uvfm(:,:,nn), -vvfm(:,:,nn),dx, dx, mask(:,:,nn), 'LambdaCI', 3, 3);
end
% iVFM RESULTS
caliVFM     =   zeros(size(vcdev));
omegaiVFM   =   zeros(size(vcdev));
for nn = 1:numel(frRng)
    fprintf(['Now Evaluating Frame ',num2str(nn),' of ',num2str(20),' for iVFM \n']);
    [caliVFM(:,:,nn),omegaiVFM(:,:,nn)] = coherent_structure_id(X,...
        Y, -uivfm(:,:,nn), vivfm(:,:,nn),dx, dx, mask(:,:,nn), 'LambdaCI', 3, 3);
end

% SET COLOR MAP TO RED BLUE
rbmap   = colormap(redblue(40));
% Center plot grids about 0
xp  =   X - mean(X(:));
yp  =   Y - mean(Y(:));
keyboard
% Plot the reconstrction results for each method
for nn = 1:numel(frRng)
    figure(1);
    toCFD = omegaCFD(:,:,nn).*mask(:,:,nn); toCFD(toCFD == 0) = nan;
    h = pcolor(xp(1,:,1)*100,yp(:,1,1)*100,toCFD);
    set(h, 'EdgeColor', 'none');
    colormap(rbmap); caxis([-750,750]);
    hold on;
    imcontour(xp(:,:,1)*100,yp(:,:,1)*100,mask(:,:,nn),[-0.05 0.05],...
        'Color',[0 0 0],'LineWidth',1.4);
    quiver(xp(10:10:end,10:30:end,1)*100,yp(10:10:end,10:30:end,1)*100,...
        ucfd(10:10:end,10:30:end,nn)/3.*mask(10:10:end,10:30:end,nn),...
        -vcfd(10:10:end,10:30:end,nn)/3.*mask(10:10:end,10:30:end,nn),...
        'Color',[0 0 0],'autoscale','off','LineWidth',1.5,'MaxHeadSize',5);
    hold off;
    set(gca,'FontSize',28,'FontWeight','bold','YDir','reverse');
    set(gcf,'Position',[100 600 500 300],'Color',[1 1 1]);
    xlabel('Position (cm)');
    ylabel('Position (cm)');
    axis image;
    axis([min(xp(:)) max(xp(:))+0.002 min(yp(:)) max(yp(:))]*100);
    axis off;
    hold on;
    quiver(min(xp(:))*100,min(yp(:))*100+1.1,1/3,0,'Color',[0 0 0],'autoscale','off','LineWidth',1.5,'MaxHeadSize',5);
    annotation('textbox',[0.315 0.6,0.1 0.1],'String','1 m/s',...
        'EdgeColor','none','FontSize',16,'FontWeight','bold');
    hold off;
    export_fig(1,fullfile(plotdir,sprintf('cfd_%03i.tiff',nn)),'-tiff','-a1','-r512');
    
    figure(2);
    toDoVeR = omegaDoVeR(:,:,nn); toDoVeR(toDoVeR == 0) = nan;
    h = pcolor(xp(1,:,1)*100,yp(:,1,1)*100,toDoVeR);
    set(h, 'EdgeColor', 'none');
    colormap(rbmap); caxis([-500,500]);
    hold on;
    imcontour(xp(:,:,1)*100,yp(:,:,1)*100,mask(:,:,nn),[-0.05 0.05],...
        'Color',[0 0 0],'LineWidth',1.4);
    quiver(xp(10:10:end,10:30:end,1)*100,yp(10:10:end,10:30:end,1)*100,...
        ucdev(10:10:end,10:30:end,nn)/3,vcdev(10:10:end,10:30:end,nn)/3,...
        'Color',[0 0 0],'autoscale','off','LineWidth',1.5,'MaxHeadSize',5);
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
    toiVFM = omegaVFM(:,:,nn); toiVFM(toiVFM == 0) = nan;
    h = pcolor(xp(1,:,1)*100,yp(:,1,1)*100,toiVFM);
    set(h, 'EdgeColor', 'none');
    colormap(rbmap); caxis([-500,500]);
    hold on;
    imcontour(xp(:,:,1)*100,yp(:,:,1)*100,mask(:,:,nn),[-0.05 0.05],...
        'Color',[0 0 0],'LineWidth',1.4);
    quiver(xp(10:10:end,10:30:end,1)*100,yp(10:10:end,10:30:end,1)*100,...
        uvfm(10:10:end,10:30:end,nn)/3,-vvfm(10:10:end,10:30:end,nn)/3,...
        'Color',[0 0 0],'autoscale','off','LineWidth',1.5,'MaxHeadSize',5);
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
    toiVFM = omegaiVFM(:,:,nn).*mask(:,:,nn); toiVFM(toiVFM == 0) = nan;
    h = pcolor(xp(1,:,1)*100,yp(:,1,1)*100,toiVFM);
    set(h, 'EdgeColor', 'none');
    colormap(rbmap); caxis([-500,500]);
    hold on;
    imcontour(xp(:,:,1)*100,yp(:,:,1)*100,mask(:,:,nn).*mask(:,:,nn),[-0.05 0.05],...
        'Color',[0 0 0],'LineWidth',1.4);
    quiver(xp(10:10:end,10:30:end,1)*100,yp(10:10:end,10:30:end,1)*100,...
        uivfm(10:10:end,10:30:end,nn)/3.*mask(10:10:end,10:30:end,nn),...
        vivfm(10:10:end,10:30:end,nn)/3.*mask(10:10:end,10:30:end,nn),...
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
caxis([-750 750]);
ylabel(h,'Vorticity (s^-1)','FontSize',18,'FontWeight','bold');
set(gca,'FontSize',18,'FontWeight','bold','YTick',-0.5:0.5:0.5);
set(gcf,'Color',[1 1 1]);
export_fig(gcf,fullfile(plotdir,'carotid_quiver_colorbar.tiff'),...
    '-tif','-a1','-r512');

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
mcfd(:,:,tt)    = sqrt(ucfd(:,:,tt).^2  +   vcfd(:,:,tt).^2 ).*...
    imerode(mask(:,:,tt),strel('disk',3));
mDoVeR(:,:,tt)  = sqrt(ucdev(:,:,tt).^2 +   vcdev(:,:,tt).^2).*...
    imerode(mask(:,:,tt),strel('disk',3));
mVFM(:,:,tt)    = sqrt(uvfm(:,:,tt).^2  +   vvfm(:,:,tt).^2 ).*...
    imerode(mask(:,:,tt),strel('disk',3));
miVFM(:,:,tt)   = sqrt(uivfm(:,:,tt).^2 +   vivfm(:,:,tt).^2).*...
    imerode(mask(:,:,tt),strel('disk',3));

    % Compute the angle
acfd(:,:,tt)    = atan2(	-vcfd(:,:,tt),   ucfd(:,:,tt)    )*180/pi.*...
    imerode(mask(:,:,tt),strel('disk',3));
aDoVeR(:,:,tt)  = atan2(	 vcdev(:,:,tt),  ucdev(:,:,tt)   )*180/pi.*...
    imerode(mask(:,:,tt),strel('disk',3));
aVFM(:,:,tt)    = atan2(	-vvfm(:,:,tt),   uvfm(:,:,tt)    )*180/pi.*...
    imerode(mask(:,:,tt),strel('disk',3));
aiVFM(:,:,tt)   = atan2(	 vivfm(:,:,tt),  uivfm(:,:,tt)   )*180/pi.*...
    imerode(mask(:,:,tt),strel('disk',3));
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
MASK = min(mask,[],3);
for tt = 1:size(ucfd,3)
    UCFD(:,:,tt)    = sqrt(ucfd(:,:,tt).^2+vcfd(:,:,tt).^2).*imerode(mask(:,:,tt),strel('disk',1));
    UDoVeR(:,:,tt)  = sqrt(ucdev(:,:,tt).^2+vcdev(:,:,tt).^2).*imerode(mask(:,:,tt),strel('disk',1));
    UVFM(:,:,tt)    = sqrt(uvfm(:,:,tt).^2+vvfm(:,:,tt).^2).*imerode(mask(:,:,tt),strel('disk',1));
    UiVFM(:,:,tt)   = sqrt(uivfm(:,:,tt).^2+vivfm(:,:,tt).^2).*imerode(mask(:,:,tt),strel('disk',1));
end

% Construct the elements for a Bland-Altman or Tukey plot
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

% Display the Bland-Altman plot
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
    linspace(meaniVFM+1.96*stdiVFM,meaniVFM+1.96*stdiVFM,501),'LineStyle','--',...
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

export_fig(gcf,fullfile(plotdir,'carotid_blandAltman_aggregate+legend.png'),'-a1','-r512','-q400','-png');