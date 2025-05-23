% Function: Infill optimization
% Author: Jun Wu (j.wu-1@tudelft.nl)
% Version: 2017-06-19
% updated 2021-01-19 for the mma version http://www.smoptit.se/
% updated 2025-4-6 for readability 

% This code requires the method of moving asymptotes (MMA).
% Downlad mma at http://www.smoptit.se/ (Prof. Svanberg at KTH)

% Example: Infill(200,100,[2],600) 
% Comparison:  Infill(200,100,[1,2],600), Infill(200,100,[1],600)

function Infill(nelx,nely,mdof,nloop)
% mdof[1,2]
% [1] global volume constraint
% [2] local volume constraint
% [1,2] both total and local volume constraints are enabled.

close all;
mkdir('images');
vol_global = 0.4;  % global volume fraction (if applicable)
vol_local = 0.6;  % local volume fraction
penal = 3;      % stiffness penalty
p = 16;         % pNorm, for aggregating local volume constraints
r_hat = 18;      % radius for measumeasuring local volume
rmin = 1.6;     % density filter radius, to avoid checkerboard 
move = 0.01;    % limited move for the design variables
beta = 1;       % beta continuation
eta = 0.5;      % projection threshold, fixed at 0.5
vol_local_pNorm = (nelx*nely*vol_local^p)^(1/p);

%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS
iLoad = 3;
Fsparse = sparse(2*(nely+1)*(nelx+1),1);
if iLoad == 0
    Fsparse(2*nelx*(nely+1)+2*(nely/2)+2,1) = -1;
    fixeddofs   = union([1:1:2],[2*(nely):1:2*(nely+1)]); %#ok<*NBRAK>
elseif iLoad == 1
    Fsparse(2,1) = -1;
    fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
elseif iLoad == 2
    Fsparse(2*nelx*(nely+1)+2*(nely/2)+2,1) = -1;
    fixeddofs   = union([2*(nely):1:2*(nely+1)],[2*(nelx+1)*(nely+1)-2:1:2*(nelx+1)*(nely+1)]);
elseif iLoad == 3
    Fsparse(2*(nely+1)*(nelx)+nely+2,1) = -1;
    fixeddofs = union([1:1:2*(nely+1)],[1]);
end
U = zeros(2*(nely+1)*(nelx+1),1);
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);

%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);

%% PREPARE PDE FILTER
edofVecF = reshape(nodenrs(1:end-1,1:end-1),nelx*nely,1);
edofMatF = repmat(edofVecF,1,4)+repmat([0 nely+[1:2] 1],nelx*nely,1);
iKF = reshape(kron(edofMatF,ones(4,1))',16*nelx*nely,1);
jKF = reshape(kron(edofMatF,ones(1,4))',16*nelx*nely,1);
iTF = reshape(edofMatF,4*nelx*nely,1);
jTF = reshape(repmat([1:nelx*nely],4,1)',4*nelx*nely,1);
sTF = repmat(1/4,4*nelx*nely,1);
TF = sparse(iTF,jTF,sTF);

Rmin = r_hat/2/sqrt(3);
KEF = Rmin^2*[4 -1 -2 -1; -1  4 -1 -2; -2 -1  4 -1; -1 -2 -1  4]/6 + ...
             [4  2  1  2;  2  4  2  1;  1  2  4  2;  2  1  2  4]/36;
sKF = reshape(KEF(:)*ones(1,nelx*nely),16*nelx*nely,1);
KF = sparse(iKF,jKF,sKF);
LF = chol(KF,'lower');

%% INITIALIZE ITERATION
x = repmat(vol_global,nely,nelx);
xTilde = x;
xPhys = (tanh(beta*eta) + tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
xold1 = reshape(x,[nely*nelx,1]);
xold2 = reshape(x,[nely*nelx,1]);
low = 0;
upp = 0;

loopbeta = 0;
loop = 0;
change = 1;
%% START ITERATION

% store results
c_hist = zeros(nloop,1);        % compliance
vol_hist = zeros(nloop,1);      % volume
change_hist = zeros(nloop,1);   % maximum design change
sharp_hist = zeros(nloop,1);    % sharpness
cons_hist = zeros(nloop,2);     % constraints


while change > 0.0001 && loop < nloop
    loopbeta = loopbeta+1;
    loop = loop+1;
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    
    U(freedofs) = K(freedofs,freedofs)\Fsparse(freedofs);

    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
    dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;

    dv = ones(nely,nelx);

    x_pde_hat = (TF'*(LF'\(LF\(TF*xPhys(:)))));
    dfdx_pde = (sum(x_pde_hat.^p))^(1/p-1) * x_pde_hat.^(p-1);
  
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    dx = beta * (1-tanh(beta*(xTilde-eta)).*tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
    dc(:) = H*(dc(:).*dx(:)./Hs);
    dv(:) = H*(dv(:).*dx(:)./Hs);

    %% UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    % METHOD OF MOVING ASYMPTOTES (MMA)
    m = size(mdof,2);
    n = nelx*nely;  
    
    df0dx = reshape(dc,[nelx*nely,1]);
    dfdx = zeros(3,nelx*nely);
    dfdx(1,1:nelx*nely) = reshape(dv,[1,nelx*nely])/(nelx*nely*vol_global);
    dfdx(2,1:nelx*nely) = TF'*(LF'\(LF\(TF*dfdx_pde(:))));
    
    ic = 2;
    tmp = reshape(dfdx(ic,:),[nely,nelx]);
    dfdx(ic,:) = reshape(H*(tmp(:).*dx(:)./Hs),[1,nelx*nely]);   
    
    iter = loopbeta;
    xval = reshape(x,[nelx*nely,1]);
    xmin=max(0.0,xval-move);
    xmax=min(1,xval+move);

    f0val = c;
    fval = zeros(2,1);
    fval(1,1) = sum(sum(xPhys)) / (nelx*nely*vol_global) - 1;
    fval(2,1) = (sum(x_pde_hat.^p))^(1/p)- vol_local_pNorm;
    
    a0 = 1;
    a = zeros(m,1);     
    c_ = ones(m,1)*1000;
    d = zeros(m,1);
    [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
        mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2,...
        f0val,df0dx,fval(mdof),dfdx(mdof,:),low,upp,a0,a,c_,d);

    xnew = reshape(xmma,[nely,nelx]);
    xold2 = xold1;
    xold1 = xval;
    
    xTilde(:) = (H*xnew(:))./Hs;
    xPhys = (tanh(beta*eta) + tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));

    change = max(abs(xnew(:)-x(:)));
    x = xnew;

    %% PRINT RESULTS
    disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
        ' Vol.: ' sprintf('%6.3f',sum(sum(xPhys))/(nelx*nely)) ...
        ' Ch.: ' sprintf('%6.3f',change) ...
        ' Cons.: ' sprintf('%6.3f',fval)]);

    %% UPDATE HEAVISIDE REGULARIZATION PARAMETER
    if beta < 100 && (loopbeta >= 40 || change <= 0.001)
        beta = 2*beta;
        loopbeta = 0;
        change = 1;
        fprintf('Parameter beta increased to %g.\n',beta);
    end
    
    %% Store current values
    c_hist(loop,1) = c;
    vol_hist(loop,1) = sum(sum(xPhys))/(nelx*nely);
    change_hist(loop,1) = change;
    cons_hist(loop,:) = fval;
    sharp_hist(loop,1) = 4*sum(sum(xPhys.*(ones(nely, nelx)-xPhys))) / (nely*nelx);

    %% PLOT DENSITIES
    figure(1);
    set(1, 'Position', [100, 450, 540, min(100+540*nely/nelx,540)]);
    colormap(gray); imagesc(-xPhys, [-1 0]); axis equal; axis tight; axis off; drawnow;
    %title('\rho');

    filename1 = sprintf('images\\rho-It%i.png',loop);
    saveas(1,filename1,'png');
end


% Publication
% Jun Wu, Niels Aage, Ruediger Westermann, Ole Sigmund, 
% Infill Optimization for Additive Manufacturing -- Approaching Bone-like Porous Structures
% IEEE Trans. on Visualization and Computer Graphics, 2017

% The code was developed based on the 110 line topology optiziation code, by E. Andreassen etc, 2011