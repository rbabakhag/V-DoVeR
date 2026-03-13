function [U,V,omega] = vector_flow_mapping(U,V,dx,dy,dvdy,mask,plotdir,n)


% Get dimensions of input Velocity
[s1,s2] = size(V);
% Compute Mask for Second-Order central difference
derivmask = mask;
derivmask(derivmask > 1) = 1;
% Find indices of Boundary Mask that are not equal to zero
[rr,cc]  = find(mask' ~=0);
ind = unique(cc);
% Initialize Left and Right reconstructed velocity matrices
utemp1 = zeros([s1 s2]);
utemp2 = zeros([s1 s2]);
% Initialize Left and Right Weighting matrices
w1 = zeros([s1 s2]);
w2 = zeros([s1 s2]);
% Start Solver
keyboard

for l = 1:numel(ind)
    % Find column indices of current row index
    idx = find(cc == ind(l));
    % If the row is continuous (no breaks)
    if max(diff(rr(idx))) == 1
        % Store gradient values of current row in temporary vector
        temp.vect = dvdy(ind(l),rr(idx(1)):rr(idx(end)));
        % Integrate temporary vector from Left to Right (with U0)
        temp.u    = cumtrapz(temp.vect)*dx+U(ind(l),rr(idx(1)));
        % Store integrated vector into Left velocity matrix
        utemp1(ind(l),rr(idx(1)):rr(idx(end))) = temp.u;
        % Update Left Weighting matrix
        w1(ind(l),rr(idx(1)):rr(idx(end))) = linspace(1,0,numel(temp.u));
        % Integrate temporary vector from Right to Left (with U0)
        temp.u = -cumtrapz(fliplr(temp.vect))*dx+U(ind(l),rr(idx(end)));
        % Store integrated vector into Right velocity matrix
        utemp2(ind(l),rr(idx(1)):rr(idx(end))) = fliplr(temp.u);
        % Update Right Weighting matrix
        w2(ind(l),rr(idx(1)):rr(idx(end))) = linspace(0,1,numel(temp.u));
    % If the row is one entry long    
    elseif isempty(max(diff(rr(idx))))
        % Store gradient values of current row in temporary vector
        temp.vect = dvdy(ind(l),rr(idx(1)):rr(idx(end)));
        % Integrate temporary vector from Left to Right (with U0)
        temp.u    = cumtrapz(temp.vect)*dx+U(ind(l),rr(idx(1)));
        % Store integrated vector into Left velocity matrix
        utemp1(ind(l),rr(idx(1)):rr(idx(end))) = temp.u;
        % Update Left Weighting matrix
        w1(ind(l),rr(idx(1)):rr(idx(end))) = linspace(1,0,numel(temp.u));
        % Integrate temporary vector from Right to Left (with U0)
        temp.u = -cumtrapz(fliplr(temp.vect))*dx+U(ind(l),rr(idx(end)));
        % Store integrated vector into Right velocity matrix
        utemp2(ind(l),rr(idx(1)):rr(idx(end))) = fliplr(temp.u);
        % Update Right Weighting matrix
        w2(ind(l),rr(idx(1)):rr(idx(end))) = linspace(0,1,numel(temp.u));
    % If the row is not continuous (breaks exist)    
    else
        % Find indices where breaks exist
        idxt = find(diff(rr(idx)) > 1);
        % Store indices of breaks
        for n = 1:numel(idxt)
            if n == 1
                idt = [idxt(n) idxt(n)+1];
            else
                idt = cat(2,idt,[idxt(n) idxt(n)+1]);
            end
        end
        IDT = [1 idt numel(rr(idx))];
        idt = IDT;
        for n = 1:2:numel(idt)-1
            % Store gradient values of current row in temporary vector
            temp.vect = dvdy(ind(l),rr(idx(idt(n))):rr(idx(idt(n+1))));
            % Integrate temporary vector from Left to Right (with U0)
            temp.u    = cumtrapz(temp.vect)*dx+U(ind(l),rr(idx(idt(n))));
            % Store integrated vector into Left velocity matrix
            utemp1(ind(l),rr(idx(idt(n))):rr(idx(idt(n+1)))) = temp.u;
            % Update Left Weighting matrix
            w1(ind(l),rr(idx(idt(n))):rr(idx(idt(n+1)))) = linspace(1,0,numel(temp.u));
            
            % Integrate temporary vector from Right to Left (with U0)
            temp.u = cumtrapz(-fliplr(temp.vect))*dx+U(ind(l),rr(idx(idt(n+1))));
            % Store integrated vector into Right velocity matrix
            utemp2(ind(l),rr(idx(idt(n))):rr(idx(idt(n+1)))) = fliplr(temp.u);
            % Update Right Weighting matrix
            w2(ind(l),rr(idx(idt(n))):rr(idx(idt(n+1)))) = linspace(0,1,numel(temp.u));
        end
    end
end
% Compute Reconstructed Velocity Component
U = (w1.*utemp1+w2.*utemp2).*derivmask;
% Compute Vorticity
omega = socdiff_bc(V,dx,2,derivmask)-socdiff_bc(U,dy,1,derivmask);
% save workspace
save([plotdir,filesep,'VFM_WorkS_Case_','_frame_',int2str(n)]);

end