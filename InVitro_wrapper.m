close all; clear; clc;
% each case is using a different velocity field for calculating initial
% vorticity

%===========================================================
% Steady: case 118, Reza, 09/20/2022 (91, frames)
% Pulsatile: Case 97, Reza, 09/25/2022 (13, frames)
%Pulsatile: added womersley flow profile, Reza, 09/30/2022 (13, frames)
% v24: updated the womersley velocity profile at the inlet, Reza,
% 10/04/2022
% 10/30/2022, Reza, adjusted the code to work for the -5 angle. Still needs to be generalized for all the angles. 
% 11/01/2022, Reza, Will change it to run for pulsaile Re=450, A10, test
% case 277
%===========================================================

%{
prompt5 = "When the data collected M/D/Y (e.g. 08192022)? ";
date = input (prompt5);
prompt6 = "Is it Steady flow data? [Steady:1, Pulsatile:0]";
St = input(prompt6);
prompt1 = "What is the Re value? ";
Re = input (prompt1);
prompt7 = "What is test number?";
Test_num = input(prompt7);
prompt2 = "What is the Trial # ? ";
Trl = input (prompt2);
prompt3 = "What is the PRF value? ";
PRF = input (prompt3);
prompt4 = "What is the beam angle? ";
ang = input (prompt4);
prompt8 = "What is the case Number? (C1:Final True Vel., C2:75% of the True velocity, C3: Color Doppler images) ";
caseNum = input (prompt8);
prompt9 = "Do you want to do the reconstruction with ZERO Initialization?";
ZeroInit = input (prompt9);
%}
server          = 1;                           % 0: if you are not running on the server  1: if you are running on the server
if_ang          = 0;                          % specify if the data is for +-5 and +-10                  
replicate       = 2;                           % If you want this processing to be saved in another folder
St              = 0;                            % 0: if it is pulsatile 1: if it is steady
date            = 9172022;
%%%
        Re          = 450;
        Test_num    = 30;
        Trl         = 3;
        PRF         = 20;
        ang         = 0;
        caseNum     = 1;
%%%
ZeroInit        = 1;                          % if you are going to do the reconstruction with zero initialization
n_proc          = 26;                           % Number of the frames that you want to process
frRng           = 1:n_proc;            %44:56   %9:21   %3:16:100

%==========================================================
%Updated by Reza 09-20-2022
% Specify directory where all data is stored
% basepath = pwd;
%%%Brett
% dataDir ='/Users/meyers18/Library/CloudStorage/Box-Box/aether.lab/projects/Echo.Cardio/vascular.DoVeR/Reza/Reza_Results/Trial/8192022/CDI/450';   %pulsatile
% dataDir = fullfile(basepath,'analysis','synthetic.carotid');
%%%Reza
% dataDir = 'C:\Users\rbabakha\Box\aether.lab\projects\Echo.Cardio\vascular.DoVeR\Reza\Reza_Results\Trial\8302022\CDI\1425';
% dataDir = '/Users/rezabg/Library/CloudStorage/Box-Box/aether.lab/projects/Echo.Cardio/vascular.DoVeR/Reza/Reza_Results/Trial/8302022/CDI/1425';     %steady
% dataDir ='/Users/rezabg/Library/CloudStorage/Box-Box/aether.lab/projects/Echo.Cardio/vascular.DoVeR/Reza/Reza_Results/Trial/8192022/CDI/450';         %pulsatile

% final images base directory
if server
Dir = '/home/ladylovelace/b/aether/Projects/DoVeR/Data/FinalImages/DopplerImages/FinalTests/';
Analysis_Dir = '/home/ladylovelace/b/aether/Projects/DoVeR/Data/DataAnalysis/FinalImagesAnalysis/CDI/';
else
Dir = 'Z:\Projects\DoVeR\Data\FinalImages\DopplerImages\FinalTests\';
Analysis_Dir = 'Z:\Projects\DoVeR\Data\DataAnalysis\FinalImagesAnalysis\CDI\';   
end

%==================================================================
% Steady Cases
if St

    if server
        basepath = ['/home/ladylovelace/b/aether/Projects/DoVeR/Data/DataAnalysis/' ...
            'FinalImagesAnalysis/CDI/Steady/A',int2str(ang),...
            filesep,int2str(Re),filesep,int2str(Test_num),'_Trial_',...
            int2str(Trl),'_',int2str(Re),'_PRF_',int2str(PRF),'_CDI'];
    else
        basepath = ['Z:\Projects\DoVeR\Data\DataAnalysis\FinalImagesAnalysis\CDI\Steady\A',int2str(ang),...
            filesep,int2str(Re),filesep,int2str(Test_num),'_Trial_',...
            int2str(Trl),'_',int2str(Re),'_PRF_',int2str(PRF),'_CDI'];
    end
    
    CDI_Dir = [Dir,int2str(date),filesep,'CDI',filesep,'Output',filesep,'Steady',filesep,'A',int2str(ang),...
        filesep,int2str(Re),filesep,int2str(Test_num),'_Trial_',...
        int2str(Trl),'_',int2str(Re),'_PRF_',int2str(PRF),'_CDI'];
    
    CFD_Dir=[Dir,int2str(date),filesep,'CFD',filesep,'Steady',filesep,'A',int2str(ang),filesep,int2str(Re),filesep,...
            int2str(Test_num),'_Trial_',int2str(Trl),'_',...
            int2str(Re),'_PRF_',int2str(PRF),'_CFD',filesep];
    
    plotdir  = fullfile(basepath,filesep,[int2str(replicate),'_results_',int2str(Test_num),'_Trial_',...
        int2str(Trl),'_','Re_',int2str(Re),'_PRF_',int2str(PRF),...
        '_Case_',int2str(caseNum),'_ZeroInit_',int2str(ZeroInit)]);
    
        if ~exist(plotdir,'dir')
            mkdir(plotdir)
        end 

end

%==================================================================
% Pulsatile cases
if ~St
    if server
    basepath = ['/home/ladylovelace/b/aether/Projects/DoVeR/Data/DataAnalysis/' ...
        'FinalImagesAnalysis/CDI/Pulsatile/A',int2str(ang),...
        filesep,int2str(Re),filesep,int2str(Test_num),'_Trial_',...
        int2str(Trl),'_',int2str(Re),'_PRF_',int2str(PRF),'_CDI'];
    else
    basepath = ['Z:\Projects\DoVeR\Data\DataAnalysis\FinalImagesAnalysis\' ...
        'CDI\Pulsatile\A',int2str(ang),filesep,int2str(Re),filesep,int2str(Test_num),'_Trial_',...
        int2str(Trl),'_',int2str(Re),'_PRF_',int2str(PRF),'_CDI'];
    end

    CDI_Dir = [Dir,int2str(date),filesep,'CDI',filesep,'Output',filesep,'Pulsatile',filesep,'A',int2str(ang),...
        filesep,int2str(Re),filesep,int2str(Test_num),'_Trial_',...
        int2str(Trl),'_',int2str(Re),'_PRF_',int2str(PRF),'_CDI'];
    
    CFD_Dir=[Dir,int2str(date),filesep,'CFD',filesep,'Pulsatile',filesep,'A',int2str(ang),filesep,int2str(Re),filesep,...
            int2str(Test_num),'_Trial_',int2str(Trl),'_',...
            int2str(Re),'_PRF_',int2str(PRF),'_CFD',filesep];
    
    plotdir  = fullfile(basepath,filesep,[int2str(replicate),'_results_',int2str(Test_num),'_Trial_',...
        int2str(Trl),'_','Re_',int2str(Re),'_PRF_',int2str(PRF),...
        '_Case_',int2str(caseNum),'_ZeroInit_',int2str(ZeroInit)]);
    
        if ~exist(plotdir,'dir')
            mkdir(plotdir)
        end 

end

%=======================================================================
addpath(fullfile(basepath,'src','export_fig'));
addpath(fullfile(basepath,'src','discrete_differencing'));
addpath(fullfile(basepath,'src','proper_orthogonal_decomposition'));
addpath(fullfile(basepath,'src','2D_pressure_lines'));
addpath(fullfile(basepath,'src','dopper_velocity_reconstruction','cDEV'));
addpath(fullfile(basepath,'src','dopper_velocity_reconstruction','VFM'));
addpath(fullfile(basepath,'src','dopper_velocity_reconstruction','iVFM'));

% Set directory for rawdata
% rawDir  = fullfile(dataDir,'rawdata');
rawDir  = fullfile(CFD_Dir);

% Find all text files in directory
% rawInfo = dir(fullfile(rawDir,'*.txt'));
rawInfo = dir(fullfile(rawDir,'*.mat'));

% % Set directory for matfiles
% matDir  = fullfile(dataDir,'matfiles');
matDir  = fullfile(CDI_Dir);

% Find all text files in directory
% matInfo = dir(fullfile(matDir,'*.mat'));
matInfo = dir(fullfile(matDir,'*.mat'));

% Read in initial mat-file to know how to interpolate raw data
temp = load(fullfile(matDir,matInfo(2).name));
% Preserve X-axis resolution, interpolate y-axis resolution to match
dx_CFI  = abs((temp.x_axis_B(2)-temp.x_axis_B(1)));
y_axis_CFI = min(temp.r_axis_B):dx_CFI:max(temp.r_axis_B);
clear temp;
% Set time points for each frame sample
vCFI  = [];
temp_CFD = load(fullfile(rawDir,rawInfo(2).name));
startF = temp_CFD.startingFrame;
% keyboard
% for n = 1:numel(frRng)
for n = startF:(startF+n_proc-1)           %44:56   %9:21
    fprintf('Evaluating frame %03i of %03i ... \r',n,numel(frRng));
    % Read in current frame mat-file
    temp = load(fullfile(matDir,matInfo(n).name));
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
    %     if n == 2
    %         bw = CFI;
    %         bw(bw == 0)  = inf; %black 0 white 1
    %         bw(bw < 12)   = 1; %%% why 12?
    %         bw(bw >= 12)  = 0;
    %         bw  = imdilate(bwselect(imerode(bw,strel('disk',1)),130,530),strel('disk',1));
    %         CFI = CFI + bw*60; %%%Why 130 and 530? %%%why multiplied by 60
    %     end
        % Second step, we interpolate to an isotropic resolution
    [xCFI,yCFI] = meshgrid(temp.x_axis_B,y_axis_CFI);
    %     vCFI = cat(3,vCFI,medfilt2(interp2(x_B,z_B,CFI,xCFI,yCFI),[3 3]));
    vCFI = cat(3,vCFI,interp2(x_B,z_B,CFI,xCFI,yCFI));
    
    %         vCFI = cat(3,vCFI,medfilt2(interp2(x_B,z_B,CFI,xCFI,yCFI),[3 3])/100);
end

% keyboard
% Initialize counter
count   =   1;

% Pre-allocate memory for CFD u & v velocities
ucfd    =   [];
vcfd    =   [];
if St
    if if_ang
        uVals1  =   [];
        vVals1  =   [];
        uVals3  =   [];
        vVals3  =   [];
    else 
        uVals1  =   [];
        uVals3  =   [];
    end

end
if ~St
    w_vel  =   [];
    uVal   =   [];
end
for n =  1:n_proc 

    fprintf('Evaluating frame %03i of 200 ... \r',n);
    % Using fscanf, load in each frame, interpolate the pressure data
    CFD = load(fullfile(rawDir,rawInfo(n).name));
    ucfd    =   cat(3,ucfd,CFD.ucfd);
    vcfd    =   cat(3,vcfd,CFD.vcfd);
    if St
        if if_ang
            uVals1 = cat(2,uVals1,CFD.uVals1);
            vVals1 = cat(2,vVals1,CFD.vVals1);
            uVals3 = cat(2,uVals3,CFD.uVals3);
            vVals3 = cat(2,vVals3,CFD.vVals3);
        else
            uVals1 = cat(2,uVals1,CFD.uVals1);
            uVals3 = cat(2,uVals3,CFD.uVals3);
        end
    end

    if ~St
        w_vel  =   cat(2,w_vel,CFD.w_vel2);
        uVal   =   cat(2,uVal,CFD.uVal);
    end

    count = count + 1;
end

% keyboard

%================================================================
X   = xCFI; %/100;
Y   = yCFI;% /100;
dx  = abs(X(2,2)-X(1,1));
dy  = abs(Y(2,2)-Y(1,1));

dt  = 1/numel(frRng);    % Frame seperation time, in seconds

% Create mask from CFD data
mask = ucfd.^2+vcfd.^2; mask(mask ~=0) = 1;
% mask = ones(size(ucfd));
bcmask = zeros(size(ucfd));

% build boundary mask points
for n = 1:numel(frRng)

    bcMask = mask(:,:,n)*7; 
    % Compute zero-crossings, label all as 1
    bnds    = cell2mat(bwboundaries(logical(mask(:,:,n))));
    inds    = sub2ind([size(mask,1) size(mask,2)],bnds(:,1),bnds(:,2));
    bcMask(inds)    = 5;
    % find main vessel inlet points, label as 1
    inds = find(bcMask(:,1) == 5);
    bcMask(inds(1:end),1) = 1;
    % find outlet points, label as 3
    inds = find(bcMask(:,end) == 5);
    bcMask(inds(1:end),end) = 3;
    % find top wall points, label as 2
    [r,c] = find(bcMask(1:10,:) == 5);
    inds  = sub2ind([size(mask,1) size(mask,2)],r,c);
    bcMask(inds)    = 2;
    % find bttom wall points, label as 4
    [r,c] = find(bcMask(:,:) == 5);
    inds  = sub2ind([size(mask,1) size(mask,2)],r,c);
    bcMask(inds) = 4;
    
    bcmask(:,:,n)   = bcMask;

end

%===================================================================
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

%=================================================================
if caseNum==1
uVort       =   zeros(size(bcmask));
vVort       =   zeros(size(bcmask));
end

if caseNum==2
uVort       =   zeros(size(bcmask));
vVort       =   zeros(size(bcmask));
end

% Compute boundary conditions 
% keyboard  
for n = 1:numel(frRng)

    if St
       % Store current boundary conditions
        psibc       = bcmask(:,:,n);
        % Compute mean velocity, which will be used to compute psi as plug flow
        inds        = find(psibc == 1);
        uMean       = mean(ucfd(inds,1,n));
        uP          = uPholder(:,:,n);
    
        if caseNum==1
            uVort1       = uVort(:,:,n);
            vVort1       = vVort(:,:,n);
        end
        if caseNum==2
            uVort1       = uVort(:,:,n);
            vVort1       = vVort(:,:,n);
        end  

        uVals = uVals1(:,n);

        if if_ang
            vVals = vVals1(:,n);
            uP(inds)    = sqrt(uVals.^2+vVals.^2);
        else
            uP(inds)    = uVals;
        end
    
    
       if caseNum==1
           vVort1(inds)=uVals;
           uVort1(inds)=uVals;
           for ii=2:size(uVort1,2)
               uVort1 (:,ii) = uVort1(:,1);
               vVort1 (:,ii) = vVort1(:,1);
           end
       end

       if caseNum==2
           vVort1(inds)=uVals/1.5;
           uVort1(inds)=uVals/1.5;
           for ii=2:size(uVort1,2)
               uVort1 (:,ii) = uVort1(:,1);
               vVort1 (:,ii) = vVort1(:,1);
           end
       end
        % Integrate and store as plug flow along intlet
        if if_ang
            psivectin   = cumtrapz(uVals)*dy-cumtrapz(vVals)*dx;  
        else
            psivectin   = cumtrapz(uVals)*dy; 
        end
        psibc(inds) = psivectin;
        % Set top and bottom walls to first and last inlet psi-value
        psibc(psibc == 2) = psivectin(1);
        psibc(psibc == 4) = psivectin(end);
        % Start by making it a plug boundary
       [inds,~] = find(psibc == 3);
        maxpsi  = psivectin(end);

        uVals = uVals3(:,n)';
        if if_ang
            vVals = vVals3(:,n)'; 
            psivectout   = flip(-cumtrapz(flip(uVals,2))*dy,2)-flip(-cumtrapz(flip(vVals,2))*dx,2); 
    %         psivectout   = flip(-cumtrapz(flip(uVals,2))*dy+cumtrapz(flip(vVals,2))*dx,2);
            uP(inds)    = sqrt(uVals3(:,n).^2+vVals3(:,n).^2);
        else 
            uP(inds)    = uVals;
            psivectout   = flip(-cumtrapz(flip(uVals,2))*dy,2);
        end
        psivectout  = psivectout + psivectin(end);
        psibc(psibc == 3) = psivectout;
        uP(inds,end) = uVals;
        % Set psi inside to 0
        psibc(psibc == 7) = 0;
        psiBC(:,:,n) = psibc;
        uPholder(:,:,n) = uP;
    
        if caseNum==1
        uVort(:,:,n )=uVort1;
        vVort (:,:,n)=vVort1;
        end
        if caseNum==2
        uVort(:,:,n )=uVort1;
        vVort (:,:,n)=vVort1;
        end
    
        figure(1);
        imagesc(psibc); colorbar;
                                
    end 

    if ~St
        % Store current boundary conditions
        psibc       = bcmask(:,:,n);
        % Compute mean velocity, which will be used to compute psi as plug flow
        inds        = find(psibc == 1);
        uMean       = mean(ucfd(inds,1,n));
        uP          = uPholder(:,:,n);
    
        if caseNum==1
            uVort1       = uVort(:,:,n);
            vVort1       = vVort(:,:,n);
        end

        if caseNum==2
            uVort1       = uVort(:,:,n);
            vVort1       = vVort(:,:,n);
        end

%         uVals = w_vel(:,n);
        uVals       = uVal(:,n);
        uP(inds)    = uVals;

        if caseNum==1
           vVort1(inds)=uVals;
           uVort1(inds)=uVals;
               for ii=2:size(uVort1,2)
                   uVort1 (:,ii) = uVort1(:,1);
                   vVort1 (:,ii) = vVort1(:,1);
               end
       end

       if caseNum==2
           vVort1(inds)=uVals/1.5;
           uVort1(inds)=uVals/1.5;
               for ii=2:size(uVort1,2)
                   uVort1 (:,ii) = uVort1(:,1);
                   vVort1 (:,ii) = vVort1(:,1);
               end
       end
   
  
        % Integrate and store as plug flow along intlet
        psivectin   = cumtrapz(uVals)*dy;          
        psibc(inds) = psivectin;
        % Set top and bottom walls to first and last inlet psi-value
        psibc(psibc == 2) = psivectin(1);
        psibc(psibc == 4) = psivectin(end);
        % Start by making it a plug boundary
       [inds,~] = find(psibc == 3);
        maxpsi  = psivectin(end);
        uVals = w_vel(:,n)';
    
        uP(inds)    = uVals;
        psivectout   = flip(-cumtrapz(flip(uVals,2))*dy,2);  
        
        psivectout  = psivectout + psivectin(end);
        psibc(psibc == 3) = psivectout;
    %     keyboard
        uP(inds,end) = uVals;
        % Set psi inside to 0
        psibc(psibc == 7) = 0;
        psiBC(:,:,n) = psibc;
        uPholder(:,:,n) = uP;
    
        if caseNum==1
        uVort(:,:,n )=uVort1;
        vVort (:,:,n)=vVort1;
        end
        if caseNum==2
        uVort(:,:,n )=uVort1;
        vVort (:,:,n)=vVort1;
        end
    
        figure();
        imagesc(psibc); colorbar;
    %     keyboard
    end

end
% keyboard

skp = 4;

%================================================================
% Generate grid information for the iVFM reconstruction method
x           = X - (max(X(:))+min(X(:)))/2; y = -Y;
[TXY,RXY]   = cart2pol(x,y);
[THETA,RHO] = meshgrid(linspace(deg2rad(165),deg2rad(15),size(vCFI,2)),...
    linspace(min(RXY(:)),max(RXY(:)),size(vCFI,1)));
% [THETA,RHO] = meshgrid(linspace(deg2rad(-145),deg2rad(-35),size(vCFI,2)),...
%     linspace(min(RXY(:)),max(RXY(:)),size(vCFI,1)));

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

%===================================================================
% Do vector field reconstruction
for n =1:numel(frRng)
    % Run median filter on data to reject any outlier vectors
    bcMask  = bcmask(:,:,n);
    vcfi    = vCFI(:,:,n);
%     vcfi(bcMask < 7) = 0;
    % Construct Gaussian filter to smooth data for DoVeR processing
    sigma = 3;
    h1  = fspecial('gaussian',...
        [2*ceil(2*sigma)+1 2*ceil(2*sigma)+1],sigma);
%     vdover = vcfi; vdover(vdover == 0) = nan;
    vdover = vcfi; vdover(vdover == 0) = nan;
%     vdover = nanmedfilt2(vdover, [3 3]);
%     vdover = nanconv(vdover,h1);
    vdover(isnan(vdover)) = 0;
    vdover(1,:)=0;
    vdover(end,:)=0;
    vdover(:,1)=0;
    vdover(:,end)=0;

%     dvdx = socdiff_bc(vdover,dx,2,mask(:,:,n));
    dvdx = socdiff_bc(uVort(:,:,n),dx,2,mask(:,:,n));
    dudy = socdiff_bc(vVort(:,:,n),dy,1,mask(:,:,n));

    % Compute vorticity on smoothed data
    if caseNum==1
    omegaInit (:,:,n) = dvdx-dudy;
    end

    if caseNum==2
    omegaInit (:,:,n) = dvdx-dudy;
    end

    if caseNum==3
    omegaInit(:,:,n) = 2*dvdx;%2*dvdx;
    end
    t_y = 35;
    b_y = 155;
    % Run DoVeR on carotid data
%%%Reza
%save workspace
% keyboard

save([plotdir,filesep, 'Main_WorkS_Case_',int2str(caseNum),'_Frame_',int2str(n)]);

%% DoVeR Reconstruction
  if ~ZeroInit  
    [ucdev(:,:,n),vcdev(:,:,n),psi(:,:,n),omegapsi(:,:,n)] =...
        cdev_cartesian_reza(uPholder(:,:,n),-vdover,dx,dx,...
        omegaInit(:,:,n),psiBC(:,:,n),bcmask(:,:,n),0.15*max(abs(-vdover(vdover~=0))),'LU',ucfd(:,:,n),n,caseNum,plotdir,ZeroInit);
  end

  if ZeroInit
    [ucdev(:,:,n),vcdev(:,:,n),psi(:,:,n),omegapsi(:,:,n)] =...
        cdev_cartesian_reza(uPholder(:,:,n)*0,-vdover,dx,dx,...
        omegaInit(:,:,n)*0,psiBC(:,:,n),bcmask(:,:,n),0.15*max(abs(-vdover(vdover~=0))),'LU',ucfd(:,:,n),n,caseNum,plotdir,ZeroInit);
  end
    %     [ucdev(:,:,n),vcdev(:,:,n),psi(:,:,n),omegapsi(:,:,n)] =...
%         cdev_cartesian(uPholder(:,:,n).*mask(:,:,n),vdover.*mask(:,:,n),dx,dx,...
%         omegaInit(:,:,n),psiBC(:,:,n),bcmask(:,:,n),0.15*max(abs(vdover(vdover~=0))),'LU');
 %%
    figure('Name','Run DoVeR on carotid data');
    imagesc(X(1,skp:skp:end)/dx,Y(skp:skp:end,1)/dx,ucdev(:,:,n)); colorbar; caxis([0 2.5]);
    hold on;
        quiver(X(skp:skp:end,skp:skp:end)/dx,Y(skp:skp:end,skp:skp:end)/dx,...
        (ucdev(skp:skp:end,skp:skp:end,n))/4500,...
        (vcdev(skp:skp:end,skp:skp:end,n))/4500,'k', 'AutoScale','off');
    set(gcf,'Position',get(0,'ScreenSize'));
    annotation('textbox',[0.7 0.15,0.1 0.1],'String',['U_m_a_x = ',num2str(max(ucdev(:,:,n),[],'all')), 'm/s'],...
        'EdgeColor','none','FontSize',16,'FontName','times','FontAngle','oblique','FontWeight','bold');
    annotation('textbox',[0.7 0.17,0.1 0.1],'String',['V_m_a_x = ',num2str(max(vcdev(:,:,n),[],'all')), 'm/s'],...
        'EdgeColor','none','FontSize',16,'FontName','times','FontAngle','oblique','FontWeight','bold');
    hold off;
    title(strcat('U velocity component contour and velocity Vectors for frame number: ',{''},int2str(n)))
    axis image;
    set(gca,'FontSize',20,'FontName','times','FontAngle','oblique','YTick', [],'XTick', [])
 
%     keyboard
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
    
    figure('Name','Run VFM reconstruction with original formulation BCs');
    imagesc(X(1,skp:skp:end)/dx,Y(skp:skp:end,1)/dx,uvfm(:,:,n)); colorbar; %caxis([-0.3 1]);
    hold on;
    quiver(X(skp:skp:end,skp:skp:end)/dx,Y(skp:skp:end,skp:skp:end)/dx,...
        (uvfm(skp:skp:end,skp:skp:end,n)),...
        -(vvfm(skp:skp:end,skp:skp:end,n)),'y');
    hold off;
    axis image;
    set(gcf,'Position',[700 100 600 600]);
    pause(1E-5);
    
    % Run VFM reconstruction with DoVeR formulation BCs
    [uvfm1(:,:,n),vvfm1(:,:,n),wvfm1(:,:,n)] = vector_flow_mapping(...
        (1*temp.Uvfm(:,:,n)+0*uPholder(:,:,n)).*mask(:,:,n),-vcfi,dx,dx,dvdy,mask(:,:,n));
    
    figure('Name','Run VFM reconstruction with DoVeR formulation BCs');
    imagesc(X(1,skp:skp:end)/dx,Y(skp:skp:end,1)/dx,uvfm1(:,:,n)); colorbar; %caxis([-0.3 1]);
    hold on;
    quiver(X(skp:skp:end,skp:skp:end)/dx,Y(skp:skp:end,skp:skp:end)/dx,...
        (uvfm1(skp:skp:end,skp:skp:end,n)),...
        (vvfm1(skp:skp:end,skp:skp:end,n)),'y');
    hold off;
    axis image;
    set(gcf,'Position',[700 100 600 600]);
    pause(1E-5);
    
    
    % Run iVFM reconstruction (method that uses least squares)
%     [u,v] = ivfm(flip((temp.Uvfm(:,:,n)+uPholder(:,:,n)).*mask(:,:,n),1),flip(-vcfi.*mask(:,:,n),1),...
%         x,y,XRT,YRT,RXY,TXY,RHO,THETA,dr,dth);
%     uivfm(:,:,n) = flip(u,1);
%     vivfm(:,:,n) = flip(v,1);
%     wivfm(:,:,n) = socdiff_bc(vivfm(:,:,n),dx,2,mask(:,:,n)) -...
%         socdiff_bc(uivfm(:,:,n),dx,1,mask(:,:,n)) ;
%     
%     figure('Name','Run iVFM reconstruction (method that uses least squares)');
%     imagesc(X(1,skp:skp:end)/dx,Y(skp:skp:end,1)/dx,uivfm(:,:,n)); colorbar; caxis([-0.3 1]);
%     hold on;
%     quiver(X(skp:skp:end,skp:skp:end)/dx,Y(skp:skp:end,skp:skp:end)/dx,...
%         (uivfm(skp:skp:end,skp:skp:end,n))*30,...
%         (vivfm(skp:skp:end,skp:skp:end,n))*30,'y','autoscale','off');
%     hold off;
%     axis image;
%     set(gcf,'Position',[700 100 600 600]);
    pause(1E-5);
end
%%
% Set image output directory
% plotdir  = fullfile(basepath,'results','images','synthetic.carotid');


% CFD GROUND TRUTH
calCFD      =   zeros(size(vcdev));
omegaCFD    =   zeros(size(vcdev));
dvdx = socdiff_bc(vcfd,dy,2,mask);
dudy = socdiff_bc(ucfd,dy,1,mask);
omegaCFD = dvdx-dudy;

% for nn = 1:numel(frRng)
%     fprintf(['Now Evaluating Frame ',num2str(nn),' of ',num2str(20),' for CFD \n']);
%     [calCFD(:,:,nn),omegaCFD(:,:,nn)] = coherent_structure_id(X,...
%         Y, ucfd(:,:,nn), -vcfd(:,:,nn),dx, dx, mask(:,:,nn), 'LambdaCI', 3, 3);
% end

% DoVeR RESULTS
calDoVeR    =   zeros(size(vcdev));
omegaDoVeR  =   zeros(size(vcdev));

dvdx = socdiff_bc(vcdev,dy,2,mask);
dudy = socdiff_bc(ucdev,dy,1,mask);
omegaDoVeR = dvdx-dudy;
omegaDoVeR =omegapsi;

% for nn = 1:numel(frRng)
%     fprintf(['Now Evaluating Frame ',num2str(nn),' of ',num2str(20),' for DoVeR \n']);
%     [calDoVeR(:,:,nn),omegaDoVeR(:,:,nn)] = coherent_structure_id(X,...
%         Y, ucdev(:,:,nn), -vcdev(:,:,nn),dx, dx, mask(:,:,nn), 'LambdaCI', 3, 3);
% end

% VFM RESULTS
calVFM      =   zeros(size(vcdev));
omegaVFM    =   zeros(size(vcdev));
dvdx = socdiff_bc(vvfm,dy,2,mask);
dudy = socdiff_bc(uvfm,dy,1,mask);
omegaVFM = dvdx-dudy;

% for nn = 1:numel(frRng)
%     fprintf(['Now Evaluating Frame ',num2str(nn),' of ',num2str(20),' for VFM \n']);
%     [calVFM(:,:,nn),omegaVFM(:,:,nn)] = coherent_structure_id(X,...
%         Y, uvfm(:,:,nn), -vvfm(:,:,nn),dx, dx, mask(:,:,nn), 'LambdaCI', 3, 3);
% end

% iVFM RESULTS
% caliVFM     =   zeros(size(vcdev));
% omegaiVFM   =   zeros(size(vcdev));
% dvdx = socdiff_bc(vivfm,dy,2,mask);
% dudy = socdiff_bc(uivfm,dy,1,mask);
% omegaiVFM = dvdx-dudy;

% for nn = 1:numel(frRng)
%     fprintf(['Now Evaluating Frame ',num2str(nn),' of ',num2str(20),' for iVFM \n']);
%     [caliVFM(:,:,nn),omegaiVFM(:,:,nn)] = coherent_structure_id(X,...
%         Y, -uivfm(:,:,nn), vivfm(:,:,nn),dx, dx, mask(:,:,nn), 'LambdaCI', 3, 3);
% end

% SET COLOR MAP TO RED BLUE
rbmap   = colormap(redblue(40));
% Center plot grids about 0
xp  =   X - mean(X(:));
yp  =   Y - mean(Y(:));

% Plot the reconstrction results for each method
for nn = 1:numel(frRng)

%{ 
 %CDI
    figure (nn*100000)
    imagesc(vCFI(:,:,nn));colorbar;caxis([min(vCFI(:)), max(vCFI(:))+0.0001])
    hold on
    quiver(1:4:size(x_axis_CFI,2),1:4:size(y_axis_CFI,2),vcfd(1:4:end,1:4:end,nn),vCFI(1:4:end,1:4:end,nn),'Color','k','LineWidth',2)
    hold off
    title(strcat('Color Doppler Images, Frame: ',int2str(nn)))
    set(gcf,'Position',get(0,'ScreenSize'));
    set(gca,'FontSize',20,'FontName','times','FontAngle','oblique')
    export_fig(nn*100000,fullfile(plotdir,sprintf('CDI_%03i.tiff',nn)),'-tiff','-a1','-r512');
    %}
    
%CFD
    sf = 12000;
    figure (nn*100);
    toCFD = omegaCFD(:,:,nn).*mask(:,:,nn); toCFD(toCFD == 0) = nan;
    h = pcolor(xp(1,:,1),yp(:,1,1),toCFD);
    set(h, 'EdgeColor', 'none');
    colormap(rbmap); caxis([-2000 2000]); colorbar;
    hold on;
    imcontour(xp(:,:,1),yp(:,:,1),mask(:,:,nn),[-0.05 0.05],...
        'Color',[0 0 0],'LineWidth',1.4);
    quiver(xp(10:10:end,10:30:end,1),yp(10:10:end,10:30:end,1),...
        ucfd(10:10:end,10:30:end,nn)/sf.*mask(10:10:end,10:30:end,nn),...
        vcfd(10:10:end,10:30:end,nn)/sf.*mask(10:10:end,10:30:end,nn),...
        'Color',[0 0 0],'LineWidth',1.5,'MaxHeadSize',5,'AutoScale','off');
    hold off;
    title (strcat('Initial Velocity Field, Frame: ',int2str(nn)))
    set(gca,'FontSize',28,'FontWeight','bold','YDir','reverse');
    set(gcf,'Position',[100 600 500 300],'Color',[1 1 1]);
    xlabel('Position (m)');
    ylabel('Position (m)');
    axis image;
    axis([min(xp(:)) max(xp(:))+0.0001 min(yp(:)) max(yp(:))]);
    hold on;
    set(gcf,'Position',get(0,'ScreenSize'));

%     quiver(max(xp(:))+0.0005,min(yp(:))+0.001,max(ucfd(:)),0,'Color',[0 0 0],'LineWidth',1.5,'MaxHeadSize',5);
    annotation('textbox',[0.7 0.10,0.1 0.1],'String',['Maximum Velocity = ',num2str(max(ucfd(:,:,nn),[],'all')), 'm/s'],...
        'EdgeColor','none','FontSize',16,'FontName','times','FontAngle','oblique','FontWeight','bold');
    hold off;
    set(gca,'FontSize',20,'FontName','times','FontAngle','oblique')
    export_fig(nn*100,fullfile(plotdir,sprintf('cfd_%03i.tiff',nn)),'-tiff','-a1','-r512');

%Dover    
    figure(nn*10000);
    toDoVeR = omegaDoVeR(:,:,nn); toDoVeR(toDoVeR == 0) = nan;
    h = pcolor(xp(1,:,1),yp(:,1,1),toDoVeR);
    set(h, 'EdgeColor', 'none');
    hold on;
    imcontour(xp(:,:,1),yp(:,:,1),mask(:,:,nn),[-0.05 0.05],...
        'Color',[0 0 0],'LineWidth',1.4);
    quiver(xp(10:10:end,10:30:end,1),yp(10:10:end,10:30:end,1),...
        ucdev(10:10:end,10:30:end,nn)/sf,vcdev(10:10:end,10:30:end,nn)/sf,...
        'Color',[0 0 0],'LineWidth',1.5,'MaxHeadSize',5,'AutoScale','off');
   colormap(rbmap); caxis([-2000 2000]);colorbar;
   hold off;
    set(gca,'FontSize',28,'FontWeight','bold','YDir','reverse');
    set(gcf,'Position',[100 100 500 300],'Color',[1 1 1]);
    title (strcat('DoVeR, Frame: ',int2str(nn)))
    xlabel('Position (m)');
    ylabel('Position (m)');
    set(gcf,'Position',get(0,'ScreenSize'));
    axis image;
    axis([min(xp(:)) max(xp(:))+0.0001 min(yp(:)) max(yp(:))]);
         hold on;
    annotation('textbox',[0.7 0.13,0.1 0.1],'String',['Maximum Velocity = ',num2str(max(ucdev(:,:,nn),[],'all')), 'm/s'],...
        'EdgeColor','none','FontSize',16,'FontName','times','FontAngle','oblique','FontWeight','bold');
    hold off;
%     set(gcf,'Position',get(0,'ScreenSize'));
    set(gca,'FontSize',20,'FontName','times','FontAngle','oblique')
    export_fig(nn*10000,fullfile(plotdir,sprintf('DoVeR_%03i.tiff',nn)),'-tiff','-a1','-r512');
 
%VFM    
    figure(nn*1000);
    toiVFM = omegaVFM(:,:,nn); toiVFM(toiVFM == 0) = nan;
    h = pcolor(xp(1,:,1),yp(:,1,1),toiVFM);
    set(h, 'EdgeColor', 'none');
    colormap(rbmap); caxis([-2000 2000]);
    hold on;
    imcontour(xp(:,:,1),yp(:,:,1),mask(:,:,nn),[-0.05 0.05],...
        'Color',[0 0 0],'LineWidth',1.4);
    quiver(xp(10:10:end,10:30:end,1),yp(10:10:end,10:30:end,1),...
        uvfm(10:10:end,10:30:end,nn)/4500,-vvfm(10:10:end,10:30:end,nn)/4500,...
        'Color',[0 0 0],'LineWidth',1.5,'MaxHeadSize',5,'AutoScale','off');
    hold off;
    set(gca,'FontSize',28,'FontWeight','bold','YDir','reverse');
    set(gcf,'Position',[600 600 500 300],'Color',[1 1 1]);
    xlabel('Position (cm)');
    ylabel('Position (cm)');
    axis image;
    axis([min(xp(:)) max(xp(:))+0.0001 min(yp(:)) max(yp(:))]);
    axis image;
    hold on;
    set(gcf,'Position',get(0,'ScreenSize'));
    annotation('textbox',[0.7 0.13,0.1 0.1],'String',['Maximum Velocity = ',num2str(max(uvfm(:,:,nn),[],'all')), 'm/s'],...
        'EdgeColor','none','FontSize',16,'FontWeight','bold');
    hold off;
    export_fig(nn*1000,fullfile(plotdir,sprintf('VFM_%03i.tiff',nn)),'-tiff','-a1','-r512');
 
%iVFM
%     figure(4);
%     toiVFM = omegaiVFM(:,:,nn).*mask(:,:,nn); toiVFM(toiVFM == 0) = nan;
%     h = pcolor(xp(1,:,1),yp(:,1,1),toiVFM);
%     set(h, 'EdgeColor', 'none');
%     colormap(rbmap); caxis([-2000,2000]);
%     hold on;
%     imcontour(xp(:,:,1),yp(:,:,1),mask(:,:,nn).*mask(:,:,nn),[-0.05 0.05],...
%         'Color',[0 0 0],'LineWidth',1.4);
%     quiver(xp(10:10:end,10:30:end,1),yp(10:10:end,10:30:end,1),...
%         uivfm(10:10:end,10:30:end,nn).*mask(10:10:end,10:30:end,nn),...
%         vivfm(10:10:end,10:30:end,nn).*mask(10:10:end,10:30:end,nn),...
%         'Color',[0 0 0],'LineWidth',1.5,'MaxHeadSize',5);
%     hold off;
%     set(gca,'FontSize',28,'FontWeight','bold','YDir','reverse');
%     set(gcf,'Position',[600 100 500 300],'Color',[1 1 1]);
%     xlabel('Position (cm)');
%     ylabel('Position (cm)');
%     axis image;
%     axis([min(xp(:)) max(xp(:))+0.0001 min(yp(:)) max(yp(:))]);
%     hold on;
%     annotation('textbox',[0.7 0.12,0.1 0.1],'String',['Maximum Velocity = ',num2str(max(uivfm(:))), 'm/s'],...
%         'EdgeColor','none','FontSize',16,'FontWeight','bold');
%     hold off;
%     set(gcf,'Position',get(0,'ScreenSize'));
%     export_fig(gcf,fullfile(plotdir,sprintf('iVFM_%03i.tiff',nn)),'-tiff','-a1','-r512');
end

% keyboard
% Plot and save the color bar for the vector field reconstruction plots
figure(5);
h = colorbar('EastOutside');
colormap(redblue(40));
caxis([-2000 2000]);
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

function xxxfft = fourier_transform (xxx,N)

xxxfft = fft(xxx)/N;

for i = 1:N
    fprintf ('%f %f\N', real (xxxfft(i)),imag(xxxfft(i)));
end
end