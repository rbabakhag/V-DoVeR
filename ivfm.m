function [UN,VN] = ivfm(U,V,X,Y,XRT,YRT,RXY,TXY,RHO,THETA,dr,dth,inSwitch)
Uf = scatteredInterpolant(TXY(:),RXY(:),U(:),'linear','none');
Vf = scatteredInterpolant(TXY(:),RXY(:),V(:),'linear','none');
% Vtmp = interp2(TXY,RXY,V,THETA,RHO,'linear');
Vtmp = Vf(THETA,RHO);
Vtmp(isnan(Vtmp)) = 0;
Utmp = Uf(THETA,RHO);
Utmp(isnan(Utmp)) = 0;
Ur  = -Vtmp.*sin(THETA) + Utmp.*cos(THETA);
Ut  =  Vtmp.*cos(THETA) + Utmp.*sin(THETA);

[M,N] = size(Ur);

mask = Ur;
mask(mask ~= 0) = 1;
% mask    = imdilate(imerode(mask,strel('disk',3)),strel('disk',2));
bcbnd   = cell2mat(bwboundaries(mask));
bcmask  = 2*mask;
for i = 1:size(bcbnd,1)
    bcmask(bcbnd(i,1),bcbnd(i,2)) = 1;
end

Ur = Ur.*mask;

inds = find(Ur == 0);
kind = 1:numel(Ur);
kind(inds) = [];

% Compute "normals" for flow domain mask
mthswt = 5;
if mthswt == 1
    [mag,grad] = imgradient(bcmask,'Sobel');
elseif mthswt == 2
    [mag,grad] = imgradient(bcmask,'Prewitt');
elseif mthswt == 3
    [mag,grad] = imgradient(bcmask,'CentralDifference');
elseif mthswt == 4
    [mag,grad] = imgradient(bcmask,'IntermediateDifference');
elseif mthswt == 5
    [mag,grad] = imgradient(bcmask,'Roberts');
end
mag(bcmask == max(bcmask(:)))   = 0;
grad(bcmask == max(bcmask(:)))  = 0;
mag(bcmask == 0)    = 0;
grad(bcmask == 0)   = 0;
gradX   = -mag.*cosd(grad);
gradY   =  mag.*sind(grad);
normX   = gradX;
normY   = gradY;
len = sqrt(normX .^ 2 + normY .^ 2);
Xnorm = normX./len;
Ynorm = normY./len;
Ynorm(isnan(Ynorm)) = 0;
Xnorm(isnan(Xnorm)) = 0;

% Build differential operators "D"
Ddot_M          = spdiags([-0.5*ones(M,1) zeros(M,1) 0.5*ones(M,1)],[-1 0 1],M,M);
Ddot_M(1,1)     = -0.5;
Ddot_M(end,end) =  0.5;
% Ddot_M(1,1)     = -1.0; Ddot_M(1,2)         =  1.0;
% Ddot_M(end,end) =  1.0; Ddot_M(end,end-1)   = -1.0;
Dddot_M         = spdiags([1*ones(M,1) -2*ones(M,1) 1*ones(M,1)],[-1 0 1],M,M);
Dddot_M(1,1)    = -1;
Dddot_M(end,end) = -1;
% Dddot_M(1,1)        =  1; Dddot_M(1,2)          = -2; Dddot_M(1,3)          = 1;
% Dddot_M(end,end)    =  1; Dddot_M(end,end-1)    = -2; Dddot_M(end,end-2)    = 1;
Ddot_N          = spdiags([-0.5*ones(N,1) zeros(N,1) 0.5*ones(N,1)],[-1 0 1],N,N);
Ddot_N(1,1)     = -0.5;
Ddot_N(end,end) =  0.5;
% Ddot_N(1,1)     = -1.0; Ddot_N(1,2)         =  1.0;
% Ddot_N(end,end) =  1.0; Ddot_N(end,end-1)   = -1.0;
Dddot_N         = spdiags([1*ones(N,1) -2*ones(N,1) 1*ones(N,1)],[-1 0 1],N,N);
Dddot_N(1,1)    = -1;
Dddot_N(end,end) = -1;
% Dddot_N(1,1)        =  1; Dddot_N(1,2)          = -2; Dddot_N(1,3)          = 1;
% Dddot_N(end,end)    =  1; Dddot_N(end,end-1)    = -2; Dddot_N(end,end-2)    = 1;

% Q0 = [I_MN 0_MN]
Q0 = speye(M*N,M*N);
Q0(inds,:) = []; Q0(:,inds) = [];
Q0 = kron([1 0],Q0);

% Q1 = [(R(:)*II_MN')/dr.*(kron(I_N,Ddot_M)+I_MN, (kron(Ddot_N,I_M))/dth]
q11 = (kron(speye(N,N),Ddot_M).*RHO(:))/dr+speye(M*N,M*N);
q11(inds,:) = []; q11(:,inds) = [];
q12 = kron(Ddot_N,speye(M,M))/dth;
q12(inds,:) = []; q12(:,inds) = [];
Q1 = [q11,q12];

% Q2 = [diag(nr),diag(nth)]
q21 = spdiags(Ynorm(:),0,M*N,M*N);
q21(inds,:) = []; q21(:,inds) = [];
q22 = spdiags(Xnorm(:),0,M*N,M*N);
q22(inds,:) = []; q22(:,inds) = [];
Q2 = [q21,q22];

q31 = (kron(speye(N,N),Dddot_M).*RHO(:).^2)/dr^2;
q31(inds,:) = []; q31(:,inds) = [];
Q31 = kron(speye(2,2),q31);
q32 = (sqrt(2)/(dr*dth))*kron(Ddot_N,Ddot_M).*RHO(:);
q32(inds,:) = []; q32(:,inds) = [];
Q32 = kron(speye(2,2),q32);
q33 = kron(Dddot_N,speye(M,M))/dth^2;
q33(inds,:) = []; q33(:,inds) = [];
Q33 = kron(speye(2,2),q33);
Q3  = [Q31; Q32; Q33];

ur = Ur(:);
ur(inds) = [];

% inds1 = sub2ind(size(Ur),bcbnd(:,1),bcbnd(:,2));
% if inSwitch == 1
%     ur = cat(1,ur,Ut(inds1));
% end

nabla   = [1,1,1];
for k = 1:3
    for i = [1 2 3]
        tic
        lambda  = -10:1:10;
        resNorm = zeros(numel(lambda),1);
        regNorm = zeros(numel(lambda),1);
        ind = zeros(numel(lambda),1);
        for j = 1:numel(lambda)
            lastwarn('')
            if i == 1
                A = (Q0'*Q0) + nabla(1)*(Q1'*Q1) + 1*nabla(2)*(Q2'*Q2) +...
                    10^lambda(j)*(Q3'*Q3);
            elseif i == 2
                A = (Q0'*Q0) + 10^lambda(j)*(Q1'*Q1) + 1*nabla(2)*(Q2'*Q2) +...
                    nabla(3)*(Q3'*Q3);
            elseif i == 3
                A = (Q0'*Q0) + nabla(1)*(Q1'*Q1) + 1*10^lambda(j)*(Q2'*Q2) +...
                    nabla(3)*(Q3'*Q3);
            end
            b = Q0'*ur(:);
            v = mldivide(A,b);
            [warnMsg, ~] = lastwarn;
            if ~isempty(warnMsg)
                ind(j) = 1;
            end
            resNorm(j) = log(norm(ur(:) - Q0*v,2));
            if 1+mod((i-1)+2,3) == 1
                regNorm(j) = log(norm(Q1*v,2));
            elseif 1+mod((i-1)+2,3) == 2
                regNorm(j) = log(norm(Q2*v,2));
            elseif 1+mod((i-1)+2,3) == 3
                regNorm(j) = log(norm(Q3*v,2));
            end
        end
        lambda(ind == 1) = [];
        resNorm(ind == 1) = [];
        regNorm(ind == 1) = [];
        [dx, ddx] = diffCenterVar(resNorm',lambda);
        [dy, ddy] = diffCenterVar(regNorm',lambda);
        kappa = 2*((dx.*ddy - ddx.*dy)./((dx.^2+dy.^2).^(3/2)));
        rmNaN = isnan(kappa);
        kappa(rmNaN == 1)   = [];
        lambda(rmNaN == 1)  = [];
        try
            if numel(kappa) < 3
                ind = 1;
            else
                [~,~,~,~,~,~,ind] = changepoint_x2_old(kappa,1);
            end
        catch
            keyboard
        end
        nabla(1+mod((i-1)+2,3)) = 10^lambda(ind);
    end
%     nabla
    A = (Q0'*Q0) + nabla(1)*(Q1'*Q1) + 1*nabla(2)*(Q2'*Q2) + nabla(3)*(Q3'*Q3);
    b = Q0'*ur(:);
    v = mldivide(A,b);
    UF = zeros(size(Ur));
    VF = zeros(size(Ur));
    UF(kind) = v(size(v)/2+1:end);
    VF(kind) = v(1:size(v)/2);
%     figure(30);
%     imagesc(UF); colorbar;
%     axis image;
%     figure(32);
%     imagesc(VF); colorbar;
%     axis image;
%     pause(1E-5);
end
Uf = scatteredInterpolant(THETA(:),RHO(:),UF(:),'linear','none');
Vf = scatteredInterpolant(THETA(:),RHO(:),VF(:),'linear','none');
uF = Uf(TXY,RXY);
uF(isnan(uF)) = 0;
vF = Vf(TXY,RXY);
vF(isnan(vF)) = 0;
UN  = uF.*sin(THETA) + vF.*cos(THETA);
VN  = uF.*cos(THETA) - vF.*sin(THETA);
end