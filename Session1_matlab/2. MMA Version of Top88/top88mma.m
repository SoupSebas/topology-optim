function top88mma(nelx,nely,volfrac,penal,rmin,optsettings)
% adapted to use MMA OPTIMIZER
% MMA defaults:
asyinitMMA=0.5; 
cMMAvalue=1000; 

opt = optsettings(1); % 0 = FD-check, 1 = MMA
maxiter = optsettings(2);
convtol = 1e-4;

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
%% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
%F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
%F = zeros(2*(nely+1)*(nelx+1),1); F(2*(nely+1):(2*(nely+1)):end)=-1e-2; 
U = zeros(2*(nely+1)*(nelx+1),1);
dofx=zeros(nely+1,nelx+1); dofy=dofx;
if mod(nely,2), error('nely must be a multiple of 2'); end
yhalf=nely/2 +1;
% supports
dofx(:,1)=1; dofy(:,1)=1;
alldof=[dofx(:) dofy(:)]'; alldof=alldof(:);
fixeddofs=find(alldof);
% loads
load_magnitude = 0.01;
Fx=zeros(nely+1,nelx+1); Fy=Fx;
Fy(yhalf,end)=-2*load_magnitude;
if 1
    Fx(1,floor(.5*nelx))=-load_magnitude;
    Fx(end,floor(.5*nelx))=load_magnitude;
    Fx(1+floor(nely/10),floor(.75*nelx))=-load_magnitude;
    Fx(end-floor(nely/10),floor(.75*nelx))=load_magnitude;
    Fy(1+floor(nely/10),floor(.75*nelx))=load_magnitude;
    Fy(end-floor(nely/10),floor(.75*nelx))=load_magnitude;
end
% ---
tmp=[Fx(:) Fy(:)]';
F=sparse(tmp(:));
%fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
%fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1) - 2*floor((nely+1)/2)]); if mod(nely,2), warning('nely should be even!'); end % middle
%fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)-2*nely]); % top
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);


if 1 % standard cantilever
    F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
    %F = zeros(2*(nely+1)*(nelx+1),1); F(2*(nely+1):(2*(nely+1)):end)=-1;
    U = zeros(2*(nely+1)*(nelx+1),1);
    fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
    alldofs = [1:2*(nely+1)*(nelx+1)];
    freedofs = setdiff(alldofs,fixeddofs);
end

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
x = repmat(volfrac,nely,nelx);

%x=x+0.01*rand(nely,nelx)-0.01/2; % add some randome fluctuations
if opt(1)==0 % FD testing
    x=repmat(linspace(1e-3,1,nelx),nely,1);
    x=repmat(0.1+0.9*sin(linspace(0,pi,nely)'),1,nelx);
end
x = reshape(x,nely,nelx);
xPhys=x;

% apply density filter:
xPhys(:) = (H*x(:))./Hs;

objhist = NaN*ones(1,maxiter);
volhist = NaN*ones(1,maxiter);
m = 1; % nr. constraints
n = nelx*nely;
NrEl=n;

switch opt
    case 0
        %% FD SENS CHECK INIT
        hstep = 10.^[-10:-2];        
        nh = length(hstep);
        xold1 = x;
        
    case 1
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
end
    
%% START ITERATION
iter = 0;
change = 1;

while ((change > convtol) && (iter<maxiter))
  iter = iter + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
  dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;  
  v = sum(xPhys(:))/(NrEl*volfrac);
  dv = ones(size(xPhys))/(NrEl*volfrac);  
  %dv = ones(nely,nelx);
  if opt(1)>0
      if iter==1, f0fac=10/c; end; % normalize objective
      if (iter>=20) & (f0fac*c<0.1), disp(sprintf('Rescaling objective (it.%d) by %f',iter,1/(f0fac*c))); f0fac=1/c; end; % renormalize objective
      if (iter>=50) & (fval<0), cMMA=1000; end; % set to default value again
  else
      f0fac=1; % no scaling when doing FD checks
  end
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  dc(:) = H*(dc(:)./Hs);
  dv(:) = H*(dv(:)./Hs);  

  f0val = c*f0fac;
  df0dx = dc(:)*f0fac;
  df0dx2 = 0*df0dx;
  fval = v -1;
  dfdx = dv(:)';
  dfdx2 = 0*dfdx;
  xval = x(:);  
  
  switch opt
      case 0
          %% FD SENS CHECK
          if iter==1
              % first store the nominal case:
              f0nom = f0val;
              fnom = fval;
              X0 = xval;
              % find the largest sensivity component:
              tmp = abs(df0dx(1:n)); mx = max(tmp); mx=mx(1);
              FDi = find(tmp==mx); FDi=FDi(1);      
              %FDi = 1; % manual pick
              NOMSENS = [df0dx(FDi); dfdx(:,FDi)];
              FDSENS = zeros(size(NOMSENS,1),nh);
              
              % prepare 1st perturbation
              disp('Perturbation 1');
              xpert = X0(:);
              xpert(FDi) = xpert(FDi) + hstep(1);
              if ~UseSpecialDesignVector
                  xnew = reshape(xpert(1:n),nely,nelx);
              else
                  xnew=xpert(:);
              end
          else
              if iter<=(nh+1)
                  % process past perturbed results
                  FDSENS(1,iter-1)=(f0val-f0nom)/hstep(iter-1);
                  for i=1:length(fval)
                      FDSENS(i+1,iter-1)=(fval(i)-fnom(i))/hstep(iter-1);
                  end
                  % prepare next perturbation
                  disp(['Perturbation ',num2str(iter)]);
                  xpert = X0(:);
                  xpert(FDi) = xpert(FDi) + hstep(min(iter,nh));
                  if ~UseSpecialDesignVector
                      xnew = reshape(xpert(1:n),nely,nelx);
                  else
                      xnew=xpert(:);
                  end                  
              else
                  % all perturbations done: evaluate results
                  format short e
                  disp('Nominal sensitivities (f0, f1, ...)')
                  disp(NOMSENS')
                  disp('FD sensitivities:');
                  disp(FDSENS')
                  figure;
                  plot(hstep,100*(repmat(NOMSENS,1,nh)-FDSENS)./repmat(NOMSENS,1,nh));
                  axis tight;
                  title('Relative error [%]');
                  set(gca,'xscale','log','ylim',[-10 10]);
                  return
              end
          end

      case 1
          %% MMA DESIGN UPDATE:          
          [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
              mmasub(m,n,iter,xval,xminvec,xmaxvec,xold1,xold2, ...
              f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0MMA,aMMA,cMMA,dMMA,asyinitMMA);
          
          xold2 = xold1;
          xold1 = xval;
          xnew = reshape(xmma,nely,nelx);
  end
  change = max(abs(xnew(:)-x(:)));  
  x = xnew;  
  xPhysold=xPhys; % for plotting
  %% update design:
  xPhys(:) = (H*xnew(:))./Hs;  

  if opt==0, change=1; end %% prevent termination in FD case
  %objhist(iter)=f0val; volhist(iter)=100*mean(xPhys(:));
  objhist(iter)=f0val; volhist(iter)=100*mean(xPhysold(:));

  %% PRINT RESULTS
  % ADJUSTED: now showing consistent data from the same design (not updated
  % x values) 
  %fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',iter,c, mean(xPhys(:)),change); 
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',iter,c, mean(xPhysold(:)),change); 
  %% PLOT DENSITIES  
  %colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
  ss=get(0,'screensize');
  fig=1; if (iter>1 && ishandle(fig)), set(0,'CurrentFigure',fig); else figure(fig); clf; end
  %set(1,'position',[0 2*ss(4)/4-ss(4)/3 ss(3)/3 2*ss(4)/3]); 
  set(fig,'position',[0 ss(4)*0.04 ss(3)/3 0.92*ss(4)]);
  %subplot(2,1,1); colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis tight; title(['Area fraction: ',num2str(volhist(iter)),'%']);
  subplot(2,1,1); colormap(gray); imagesc(1-xPhysold); caxis([0 1]); axis equal; axis tight; title(['Area fraction: ',num2str(volhist(iter)),'%']);
  %subplot(2,1,2); imagesc(1-x); caxis([0 1]); axis equal; axis tight; title(['design field x:']);
  subplot(2,1,2); imagesc(1-reshape(xold1,nely,nelx)); caxis([0 1]); axis equal; axis tight; title(['design field x:']);
  drawnow; 
    
  % obj sensitivity:
  fig=2; if (iter>1 && ishandle(fig)), set(0,'CurrentFigure',fig); else figure(fig); clf; end
  set(fig,'position',[ss(3)/3 2*ss(4)/4 ss(3)/3 ss(4)/3]); clf;
  if 1
        imagesc(reshape(df0dx,nely,nelx)); axis equal; axis tight; colorbar; title('Sens'); drawnow;
  end
  fig=3; if (iter>1 && ishandle(fig)), set(0,'CurrentFigure',fig); else figure(fig); clf; end
  set(fig,'position',[2*ss(3)/3 2*ss(4)/4 ss(3)/3 ss(4)/3]); 
  subplot(2,1,1); plot(objhist,'b-'); title(num2str(f0val));
  subplot(2,1,2); plot(volhist,'r-'); drawnow;
    
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

