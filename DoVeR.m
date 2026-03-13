function [U,V,psi,omega] = cdev_cartesian_Org(U,V,dx,dy,omega,psi,mask,vThresh,linalg_solver,plotdir,caseNum,FrmeNum)
% PSI-OMEGA SOLVER Iterative solver for vector field reconstruction
% Find all unique Boundary Mask values
masknum = unique(mask);
% Throw Error if Boundary Mask is not multi-level
if masknum < 3
    error('Error. Number of mask entries is not sufficient.')
end
% Required Settings for Solver
solver_switch    = 'SOR';       % Set as either "GS" or "SOR"
isotropic_switch = 'isotropic'; % Set as either "isotropic" or "generalized"
smoothvort       = 'On';        % Turn On or Off Vorticity smoothing
trendplt         = 'Off';       % Turn On or Off Convergence Trend Plotter
SSE     = 1;            % Set SSE value to 1 (initialize value)
TREND   = zeros(1,1);   % Initialize SSE storage vector
count   = 1;            % Initialize Counter for while-loop
beta    = 1.87;         % Set SOR smoothing coefficient (if SOR is used)
maskmax = max(mask(:)); % Find maximum level value of Boundary Mask
sigma   = 4;            % Set standard deviation filter size for Vorticity
% Store initial values for U & V when they enter the algorithm
Uold = U;
Vori = V;
% Compute Mask for Second-Order central difference
derivmask = mask;
derivmask(derivmask > 1) = 1;
% Find indices of mask 
[rin,cin] = find(mask == maskmax);
index = find(mask == maskmax);
if strcmp(smoothvort,'On')
    % Compute Gaussian Smoothing Kernel
    h1 = fspecial('gaussian',[2*ceil(2*sigma)+1 2*ceil(2*sigma)+1],sigma);
end
% If solver is a matrix inversion approach, generate matrix operators
if strcmp(linalg_solver,'direct') || strcmp(linalg_solver,'LU')
    [DivGrad,D,lower,upper,qtemp] = build_operators(psi,rin,cin);
end
tic
% keyboard
% Start Solver Iteration
while SSE >= 1E-4
    if strcmp(linalg_solver,'direct')
        % Construct RHS vector for direct method solvers
        q = omega(index)*dx*dy+qtemp;
        % Run Laplace 2D Solver (Direct Solution)
        psi(index) = DivGrad\q;
        % Compute New X-component velocity
        U = PSIdiffY(U,psi,rin,cin,dy);
        % Compute New Y-Component Velocity
        V = PSIdiffX(V,Vori,vThresh,psi,rin,cin,dx);
        % Compute Vorticity
        omega = compute_Vort(U,V,rin,cin,dx,dy);
    elseif strcmp(linalg_solver,'LU')
        % Construct RHS vector for direct method solvers
        q = omega(index)*dx*dy+qtemp;
        % Run Laplace 2D Solver (LU Decomposition Solver)
        psi(index) = upper\(lower\q);
        % Compute New X-component velocity
        U = PSIdiffY(U,psi,rin,cin,dy);
        % Compute New Y-Component Velocity
        V = PSIdiffX(V,Vori,vThresh,psi,rin,cin,dx);
        % Compute Vorticity
        omega = compute_Vort(U,V,rin,cin,dx,dy);
    elseif strcmp(linalg_solver,'iterative')
        % Run Laplace 2D Solver (Seperated for profiling)
        psi = laplace2D(psi,omega,beta,rin,cin,dx,dy,solver_switch,isotropic_switch);
        % Compute New X-component velocity
        U = PSIdiffY(U,psi,rin,cin,dy);
        % Compute New Y-Component Velocity
        V = PSIdiffX(V,Vori,vThresh,psi,rin,cin,dx);
        % Run Laplace 2D Solver (Seperated for profiling)
        psi = laplace2D(psi,omega,beta,flipud(rin),flipud(cin),dx,dy,solver_switch,isotropic_switch);
        % Compute New X-component velocity
        U = PSIdiffY(U,psi,flipud(rin),flipud(cin),dy);
        % Compute New Y-Component Velocity
        V = PSIdiffX(V,Vori,vThresh,psi,flipud(rin),flipud(cin),dx);
        % Compute Vorticity
        omega = compute_Vort(U,V,rin,cin,dx,dy);
    end
    if strcmp(smoothvort,'On')
        % Median filter Vorticity and Gaussian Smooth
        omega = imfilter(medfilt2(omega,[3 3]),h1).*derivmask;
    end
    % Add SSE Value from previous iteration to Error Trend
    TREND(count) = SSE;
    % Increase Counter Value
    count = count+1;
    % Compute New SSE
    SSE = norm(abs(Uold-U),2)/norm(Uold,2);
%     SSE = norm(omega(index)+DivGrad*psi(index),2)/res0;
    if strcmp(trendplt,'On')
        figure(101); plot(TREND,'LineWidth',2);
        set(gca,'YScale','log');
        set(gca,'XScale','log');
        set(gca,'FontSize',16,'FontWeight','bold');
        ylabel('L2 Norm Error, Psi','FontSize',16,'FontWeight','bold');
        xlabel('Iteration #','FontSize',16,'FontWeight','bold');
        pause(5E-2)
    end
    % Store Current reconstruction for Convergence SSE calculation
    Uold = U;
end
toc
% keyboard
% Compute New Y-Component Velocity
V = PSIdiffX2(V,Vori,vThresh,psi,rin,cin,dx);
% Compute Vorticity
omega = compute_Vort(U,V,rin,cin,dx,dy);
% keyboard
if strcmp(smoothvort,'On')
    % Median filter Vorticity and Gaussian Smooth
    omega = imfilter(medfilt2(omega,[3 3]),h1).*derivmask;
end
% save workspace
save([plotdir,filesep,'cdev_WorkS_Case_',int2str(caseNum),'_frame_',int2str(FrmeNum)]);
end

function psi = laplace2D(psi,omega,beta,r,c,dx,dy,solver_switch,isotropic_switch)
switch solver_switch
    % Gauss-Siedel formulation
    case 'GS'
        % Compute Updated Psi from Vorticity & Psi Field
        for n = 1:numel(r)
            switch isotropic_switch
                % If dx == dy
                case 'isotropic'
                    psi(r(n),c(n)) = 0.25*(psi(r(n)+1,c(n))... % Right of the current cell
                        +psi(r(n)-1,c(n))...             % Left of the current cell
                        +psi(r(n),c(n)+1)...             % Above current cell
                        +psi(r(n),c(n)-1)...             % Below current cell
                        +omega(r(n),c(n))*dx^2);         % vorticity in current cell
                    % If dx ~= dy
                case 'generalized'
                    psi(r(n),c(n)) = ((psi(r(n)+1,c(n))+psi(r(n)-1,c(n)))*dy^2+... % Righ/Left
                        +(psi(r(n),c(n)+1)+psi(r(n),c(n)-1))*dx^2+...        % Above/Below
                        +omega(r(n),c(n))*dx^2*dy^2)/(2*dx^2+dy^2);    % vorticity in current cell
            end
        end
        % Successive Over-Relaxation Formulation
    case 'SOR'
        % Compute Updated Psi from Vorticity & Psi Field
        for n = 1:numel(r)
            switch isotropic_switch
                % If dx == dy
                case 'isotropic'
                    psi(r(n),c(n)) = 0.25*beta*(psi(r(n)+1,c(n))... % Right of the current cell
                        +psi(r(n)-1,c(n))...             % Left of the current cell
                        +psi(r(n),c(n)+1)...             % Above current cell
                        +psi(r(n),c(n)-1)...             % Below current cell
                        +omega(r(n),c(n))*dx^2)+...      % vorticity in current cell
                        (1.0-beta)*psi(r(n),c(n));
                    % If dx ~= dy
                case 'generalized'
                    psi(r(n),c(n)) = beta*((psi(r(n)+1,c(n))+psi(r(n)-1,c(n)))*dy^2+... % Right/Left
                        +(psi(r(n),c(n)+1)+psi(r(n),c(n)-1))*dx^2+...             % Above/Below
                        +omega(r(n),c(n))*dx^2*dy^2)/(2*dx^2+2*dy^2)+...    % vorticity
                        (1.0-beta)*psi(r(n),c(n));
            end
        end
end
end

function U = PSIdiffY(U,psi,r,c,dy)
for n = 1:numel(r)
    U(r(n),c(n)) = (psi(r(n)+1,c(n))-psi(r(n)-1,c(n)))/(2*dy);
end
end

function V = PSIdiffX(V,Vori,vThresh,psi,r,c,dx)
for n = 1:numel(r)
    V(r(n),c(n)) = -(psi(r(n),c(n)+1)-psi(r(n),c(n)-1))/(2*dx);
    if Vori(r(n),c(n)) ~= 0
        V(r(n),c(n)) = Vori(r(n),c(n));
    else
        if abs(V(r(n),c(n))) > vThresh
            V(r(n),c(n)) = vThresh* (V(r(n),c(n)));
        end
    end
end
end
function V = PSIdiffX2(V,Vori,vThresh,psi,r,c,dx)
    for n = 1:numel(r)
        V(r(n),c(n)) = -(psi(r(n),c(n)+1)-psi(r(n),c(n)-1))/(2*dx);
    end
end

function omega = compute_Vort(U,V,r,c,dx,dy)
omega = zeros(size(U));
for n = 1:numel(r)
    omega(r(n),c(n)) = (V(r(n),c(n)+1)-V(r(n),c(n)-1))/(2*dx)-...
        (U(r(n)+1,c(n))-U(r(n)-1,c(n)))/(2*dy);
end
end

function [DivGrad,D,lower,upper,qtemp] = build_operators(psi,rin,cin)
DivGrad = sparse(numel(rin),numel(cin));
qtemp   = zeros(numel(rin),1);
for n = 1:numel(rin)
    DivGrad(n,n) = 4;
    ind = find(rin(n) == rin & cin(n)-1 == cin);
    if ~isempty(ind)
        DivGrad(n,ind) = -1;
    else
        qtemp(n) = qtemp(n)+psi(rin(n),cin(n)-1);
    end
    ind = find(rin(n) == rin & cin(n)+1 == cin);
    if ~isempty(ind)
        DivGrad(n,ind) = -1;
    else
        qtemp(n) = qtemp(n)+psi(rin(n),cin(n)+1);
    end
    ind = find(rin(n)-1 == rin & cin(n) == cin);
    if ~isempty(ind)
        DivGrad(ind,n) = -1;
    else
        qtemp(n) = qtemp(n)+psi(rin(n)-1,cin(n));
    end
    ind = find(rin(n)+1 == rin & cin(n) == cin);
    if ~isempty(ind)
        DivGrad(ind,n) = -1;
    else
        qtemp(n) = qtemp(n)+psi(rin(n)+1,cin(n));
    end
end
D = sparse(diag(diag(DivGrad)));
[lower,upper] = lu(DivGrad);
end