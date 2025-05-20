function top88HSmma_compliant(nelx,nely,volfrac,penal,rmin,maxiter)
% function top88HSmma(nelx,nely,volfrac,penal,rmin,maxiter)
% adapted from top88 to use MMA OPTIMIZER and HS projection
% https://doi.org/10.1007/s00158-010-0602-y
% MMA defaults:
asyinitMMA=0.5; 
cMMAvalue=1000; 
convtol = 1e-4;
movelimit = 1; % set to 1 or higher to deactivate

%% PROBLEM DEFINITION
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
%%

%% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
% F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
% %F = zeros(2*(nely+1)*(nelx+1),1); F(2*(nely+1):(2*(nely+1)):end)=-1;
% U = zeros(2*(nely+1)*(nelx+1),1);
% fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
% alldofs = [1:2*(nely+1)*(nelx+1)];
% freedofs = setdiff(alldofs,fixeddofs);
%% DEFINE LOADS AND SUPPORTS
n_in   = (nelx+1)*(nely+1);   % top-right – load corner
n_fix  = 1;                   % bottom-left – anchor corner
n_pin  = nely+2;              % second node on bottom edge (prevents rotation)

F      = sparse(2*(nely+1)*(nelx+1),1);
F(2*n_in-1) =  1;             % Fx
F(2*n_in  ) = -1;             % Fy

U      = zeros(2*(nely+1)*(nelx+1),1);

fixeddofs = [ 2*n_fix-1  2*n_fix  ...   % anchor node (ux, uy)
              2*n_pin-1 ];              % extra pin (ux only)

alldofs  = 1:2*(nely+1)*(nelx+1);
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

%% INITIALIZE ITERATION
xmin = 0*1e-4;
beta = 4; eta=0.5; % Heaviside projection parameters
HSdenominator = tanh(beta*eta)+tanh(beta*(1-eta));
x = repmat(volfrac,nely,nelx);
xTilde=x; xPhys=x;

% apply density filter:
xTilde(:) = (H*x(:))./Hs;
% apply projection:
xPhys(:) = (tanh(beta*eta)+tanh(beta*(xTilde(:) - eta)))/HSdenominator;

objhist = NaN*ones(1,maxiter);
volhist = NaN*ones(1,maxiter);
m = 1; % nr. constraints
n = nelx*nely;
NrEl=n;

%% MMA INIT        
xminvec  = xmin*ones(n,1);
xmaxvec  = ones(n,1);
low   = xminvec;
upp   = xmaxvec;

%cMMA = 10000*ones(m,1);
cMMA = cMMAvalue*ones(m,1); 
dMMA = zeros(m,1);
a0MMA = 1;
aMMA = zeros(m,1);

xold1 = x;
xold2 = x;        
    
%% START ITERATION
iter = 0;
change = 1;
loopbeta = 0;

while ((change > convtol) && (iter<maxiter))
  iter = iter + 1;
  loopbeta = loopbeta+1;
  %% FE-ANALYSIS
   %% FE-ANALYSIS ----------------------------------------------------------
  sK = reshape( KE(:) .* (Emin + xPhys(:)'.^penal * (E0-Emin) ), ...
                64*nelx*nely , 1 );
  K  = sparse(iK , jK , sK);   K = (K+K.')/2;
  U(freedofs) = K(freedofs,freedofs) \ F(freedofs);

  %% OBJECTIVE  – maximise compliance (minimise –c) ----------------------
  ce = reshape( sum( (U(edofMat)*KE) .* U(edofMat) , 2 ) , nely, nelx );
  c  = -sum( sum( (Emin + xPhys.^penal*(E0-Emin)).*ce ) );   % >>> minus
  dc =  penal*(E0-Emin)*xPhys.^(penal-1).*ce;                % >>> sign
if iter == 1
      f0fac = 10/abs(c);        % >>> positive scale only ONCE
  end
  f0val  = c * f0fac;           % negative value → maximise compliance
  df0dx  = dc(:) * f0fac;
  df0dx2 = zeros(size(df0dx));

  %% VOLUME  – lower bound  v ≥ 1  ---------------------------------------
  v     = sum(xPhys(:)) / (NrEl*volfrac);
  dv    = ones(size(xPhys)) / (NrEl*volfrac);

  fval  = 1 - v;                % >>> g(x) ≤ 0 is violated if v<1
  dfdx  = -dv(:)';              % >>> sign matches fval
  dfdx2 = zeros(size(dfdx));
  

  %% FILTER SENSITIVITIES -------------------------------------------------
  dH = beta*(1 - tanh(beta*(xTilde(:)-eta)).^2)/HSdenominator;
  dc(:) = H * ( (dc(:).*dH) ./ Hs );
  dv(:) = H * ( (dv(:).*dH) ./ Hs );

  %% MMA UPDATE -----------------------------------------------------------
  xval = x(:);
  [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
        mmasub(m,n,iter,xval,xminvec,xmaxvec,xold1,xold2, ...
                f0val,df0dx,df0dx2, ...
                fval ,dfdx ,dfdx2 ,low,upp, ...
                a0MMA,aMMA,cMMA,dMMA,asyinitMMA);

  if movelimit < 1          % optional move limit
      dx = xmma - xval;
      xmma = min( max(xval - movelimit , xmma) , xval + movelimit );
  end

  xold2 = xold1;   xold1 = xval;
  x      = reshape(xmma , nely , nelx);
  change = max(abs(x(:)-xPhys(:)));      % monitor change on physical field

  %% UPDATE PHYSICAL DENSITIES -------------------------------------------
  xTilde(:) = (H*x(:))./Hs;
  xPhys(:)  = (tanh(beta*eta) + tanh(beta*(xTilde(:)-eta)))/HSdenominator;
xPhysold = xPhys;     % for the next plot

  objhist(iter) = f0val;
  volhist(iter) = 100*v;

  %% PLOTS & TEXT ---------------------------------------------------------
  fprintf('It:%4d  Obj:%10.3e  Vol:%6.2f%%  Δ:%6.3f\n', ...
           iter , -f0val/f0fac , 100*v , change);

  %% PLOT DENSITIES  
  %colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
  ss=get(0,'screensize');
  fig=1; if (iter>1 && ishandle(fig)), set(0,'CurrentFigure',fig); else figure(fig); clf; end
  set(fig,'position',[0 ss(4)*0.04 ss(3)/3 0.92*ss(4)]);
  subplot(2,1,1); colormap(gray); imagesc(1-xPhysold); caxis([0 1]); axis equal; axis tight; title(['Area fraction: ',num2str(volhist(iter)),'%']);
  subplot(2,1,2); imagesc(1-reshape(xold1,nely,nelx)); caxis([0 1]); axis equal; axis tight; title(['design field x:']);
  drawnow;     
  % obj sensitivity:
  fig=2; if (iter>1 && ishandle(fig)), set(0,'CurrentFigure',fig); else figure(fig); clf; end
  set(fig,'position',[ss(3)/3 2*ss(4)/4 ss(3)/3 ss(4)/3]); clf;
  imagesc(reshape(df0dx,nely,nelx)); axis equal; axis tight; colorbar; title('Objecitve sensitivity'); drawnow;
  % convergence:
  fig=3; if (iter>1 && ishandle(fig)), set(0,'CurrentFigure',fig); else figure(fig); clf; end
  set(fig,'position',[2*ss(3)/3 2*ss(4)/4 ss(3)/3 ss(4)/3]); 
  subplot(2,1,1); plot(objhist,'b-'); title(['Objective: ',num2str(f0val)]);
  subplot(2,1,2); plot(volhist,'r-'); title('Constraint'); drawnow;
  


												
  %% UPDATE HEAVISIDE REGULARIZATION PARAMETER
  if (beta < 128) && (loopbeta >= 40 || change <= 0.01)
    beta = sqrt(2)*beta;
    loopbeta = 0;
    change = 1;
    fprintf('Parameter beta increased to %g.\n',beta);				  
  end   
end


end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Andreassen, A. Clausen, M. Schevenels,%
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,           %
%  Technical University of Denmark,                                        %
%  DK-2800 Lyngby, Denmark.                                                %
% Please sent your comments to: sigmund@fam.dtu.dk                         %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
% This version is based on earlier 99-line code                            %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

