function varargout=anyfilter(varargin)
% ANYFILTER
% Flexible filter library (2D regular rectangular domains). 
% Possible uses:
% [filterdata] = anyfilter(type, nelx, nely, parameters)
%        initializes a filterdata struct for a given filter
% [filterdata] = anyfilter(filterdata, type, nelx, nely, parameters)
%        initializes another subsequent filter operation for a given
%        filter. in this way, a multi-step filter can be built.
% [filterdata, x_filtered] = anyfilter(filterdata, x)
%        applies the defined filter to field x. may alter filterdata by
%        storing additional information
% [filterdata, sens_filtered] = anyfilter(filterdata, x, sens)
%        applies the defined filter to correct sensitivities
%%
% available filter types:
% density, dilate, dmean, erode, gmean, mean, mean2, none, HS, HS01, HS01flipped, 
% KSmask, mill, mill2, mill3, mill4, mill5, mill*, mill*2, mill*3, mill*3P, mill*4, mill*4P, 
% oneminusx, Pmask, symm, stdev
%
% CAUTION: do not use multi-step filters on different input design fields.
% intermediate field data is stored internally. For multiple input design
% fields, use multiple multi-step filters.
%
% ML 2017

% filterdata structure:
% type:       filter type
% nelx, nely: field dimensions
% data:       filter-specific data (struct)

% HOW TO ADD A NEW FILTER
% 1. Under '1. FILTER INITIALIZATION', add a switch case entry with the filter name,
%    store the filter-specific input data in fd.data.###
%    and perform any required initializations.
% 2. Under '2. FORWARD FILTERING', add a switch case entry with the filter name,
%    and perform the forward filter operation. Assign output field to 'xf'.
% 3. Under '3. SENSITIVITY FILTERING', add a switch case entry with the filter name,
%    and perform the appropriate operation to correct the sensitivities.
%    Assign new sensitivities to 'sf'.

if isstr(varargin{1})
    % 1. FILTER INITIALIZATION
    % [filterdata] = anyfilter(type, nelx, nely, parameters)
    fd=struct('type',varargin{1},'nelx',varargin{2},'nely',varargin{3},'data',struct,'x_in',[]);
    switch fd.type
        case 'HS' % further arguments: beta, eta
            error(nargchk(5,5,nargin));
            fd.data.beta=varargin{4};
            fd.data.eta=varargin{5};
            fd.data.num=tanh(fd.data.beta*fd.data.eta)+tanh(fd.data.beta*(1-fd.data.eta));
            
        case 'HS01' % further arguments: beta
            % centered at 0.5, guaranteed to remain between 0 and 1 for any
            % inputs.
            error(nargchk(4,4,nargin));
            fd.data.beta=varargin{4};
            
        case 'HS01flipped' % further arguments: beta
            % centered at 0.5, guaranteed to remain between 0 and 1 for any
            % inputs.
            error(nargchk(4,4,nargin));
            fd.data.beta=varargin{4};
            
        case 'density' % further arguments: radius
            error(nargchk(4,4,nargin));
            fd.data.rmin=varargin{4};
            rmin=fd.data.rmin;
            iH = ones(fd.nelx*fd.nely*(2*(ceil(rmin)-1)+1)^2,1);
            jH = ones(size(iH));
            sH = zeros(size(iH));
            k = 0;
            for i1 = 1:fd.nelx
                for j1 = 1:fd.nely
                    e1 = (i1-1)*fd.nely+j1;
                    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),fd.nelx)
                        for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),fd.nely)
                            e2 = (i2-1)*fd.nely+j2;
                            k = k+1;
                            iH(k) = e1;
                            jH(k) = e2;
                            sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
                        end
                    end
                end
            end
            fd.data.H = sparse(iH,jH,sH);
            fd.data.Hs = sum(fd.data.H,2);
            
        case 'dilate' % further arguments: radius, beta
            error(nargchk(5,5,nargin));
            fd.data.rmin=varargin{4};
            fd.data.beta=varargin{5};
            fd=prep(fd);
            
        case 'dmean' % further arguments: radius, beta, eta, a, b
            % combination of local mean and knock-down function
            error(nargchk(8,8,nargin));
            fd.data.rmin=varargin{4};
            fd.data.beta=varargin{5};
            fd.data.eta=varargin{6};
            fd.data.a=varargin{7};
            fd.data.b=varargin{8};
            if 0 % plot the knock-down function                
                figure(44);
                beta=fd.data.beta; eta=fd.data.eta; a=fd.data.a; b=fd.data.b;
                num=tanh(beta*eta)+tanh(beta*(1-eta));
                xmean = linspace(0,1,200);                
                S=(tanh(beta*eta)+tanh(beta*(xmean-eta)))/num;
                f=1-a*S-b*xmean.*S;
                plot(xmean,f); set(gca,'ylim',[0 1]); set(findobj(gca,'type','line'),'linewidth',2);
                title(sprintf('beta=%.1f, eta=%.2f, a=%.3f, b=%.3f',beta,eta,a,b));
                error('stopped');
            end
            rmin=fd.data.rmin;
            iH = ones(fd.nelx*fd.nely*(2*(ceil(rmin)-1)+1)^2,1);
            jH = ones(size(iH));
            sH = zeros(size(iH));
            k = 0;
            for i1 = 1:fd.nelx
                for j1 = 1:fd.nely
                    e1 = (i1-1)*fd.nely+j1;
                    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),fd.nelx)
                        for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),fd.nely)
                            e2 = (i2-1)*fd.nely+j2;
                            k = k+1;
                            iH(k) = e1;
                            jH(k) = e2;
                            %sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
                            sH(k) = 1;
                        end
                    end
                end
            end
            fd.data.H = sparse(iH,jH,sH);
            fd.data.Hs = sum(fd.data.H,2);           
            
        case 'dmeanline' % further arguments: radius, beta, eta, a, b, orientation (- | / \)
            % combination of local mean and knock-down function
            error(nargchk(9,9,nargin));
            fd.data.rmin=varargin{4};
            fd.data.beta=varargin{5};
            fd.data.eta=varargin{6};
            fd.data.a=varargin{7};
            fd.data.b=varargin{8};
            fd.data.orientation=varargin{9};
            if 0 % plot the knock-down function                
                figure(44);
                beta=fd.data.beta; eta=fd.data.eta; a=fd.data.a; b=fd.data.b;
                num=tanh(beta*eta)+tanh(beta*(1-eta));
                xmean = linspace(0,1,200);                
                S=(tanh(beta*eta)+tanh(beta*(xmean-eta)))/num;
                f=1-a*S-b*xmean.*S;
                plot(xmean,f); set(gca,'ylim',[0 1]); set(findobj(gca,'type','line'),'linewidth',2);
                title(sprintf('beta=%.1f, eta=%.2f, a=%.3f, b=%.3f',beta,eta,a,b));
                error('stopped');
            end
            rmin=fd.data.rmin;
            %iH = ones(fd.nelx*fd.nely*(2*(ceil(rmin)-1)+1)^2,1);
            iH = ones(fd.nelx*fd.nely*(2*(ceil(rmin)-1)+1),1); % no ^2
            jH = ones(size(iH));
            sH = zeros(size(iH));
            k = 0;
            for i1 = 1:fd.nelx
                for j1 = 1:fd.nely
                    e1 = (i1-1)*fd.nely+j1;
                    switch fd.data.orientation
                        case '-'
                            for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),fd.nelx)
                                for j2 = j1
                                    e2 = (i2-1)*fd.nely+j2;
                                    k = k+1;
                                    iH(k) = e1;
                                    jH(k) = e2;
                                    sH(k) = 1;
                                end
                            end
                        case '|'
                            for i2 = i1
                                for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),fd.nely)
                                    e2 = (i2-1)*fd.nely+j2;
                                    k = k+1;
                                    iH(k) = e1;
                                    jH(k) = e2;
                                    sH(k) = 1;
                                end
                            end
                        case '/'
                            xsteps=[max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),fd.nelx)]-i1;
                            ysteps=min(fd.nely,max(j1+xsteps, 1));
                            for i2 = xsteps+i1
                                %for j2 = max(j1+xdiff,fd.nely):-1:min(j1-xdiff,1)
                                for j2 = ysteps
                                    e2 = (i2-1)*fd.nely+j2;
                                    k = k+1;
                                    iH(k) = e1;
                                    jH(k) = e2;
                                    sH(k) = 1;
                                end
                            end
                        case '\'
                            xsteps=[max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),fd.nelx)]-i1;
                            ysteps=min(fd.nely,max(j1-xsteps, 1));
                            for i2 = xsteps+i1
                                for j2 = ysteps
                                    e2 = (i2-1)*fd.nely+j2;
                                    k = k+1;
                                    iH(k) = e1;
                                    jH(k) = e2;
                                    sH(k) = 1;
                                end
                            end
                        otherwise
                            error('wrong orientation');
                    end
                end
            end
            fd.data.H = sparse(iH,jH,sH);
            fd.data.Hs = sum(fd.data.H,2);           
            
            
        case 'erode' % further arguments: radius, beta
            error(nargchk(5,5,nargin));
            fd.data.rmin=varargin{4};
            fd.data.beta=varargin{5};
            fd=prep(fd);    
            

        case 'gmean' % further arguments: radius, P, beta, eta, a, b
            % global knock-down function based on P-norm maximum of local mean values
            warning('dmean is simpler and works maybe even better!')
            error(nargchk(9,9,nargin));
            fd.data.rmin=varargin{4};
            fd.data.P=varargin{5};
            fd.data.beta=varargin{6};
            fd.data.eta=varargin{7};
            fd.data.a=varargin{8};
            fd.data.b=varargin{9};
            rmin=fd.data.rmin;
            iH = ones(fd.nelx*fd.nely*(2*(ceil(rmin)-1)+1)^2,1);
            jH = ones(size(iH));
            sH = zeros(size(iH));
            k = 0;
            for i1 = 1:fd.nelx
                for j1 = 1:fd.nely
                    e1 = (i1-1)*fd.nely+j1;
                    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),fd.nelx)
                        for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),fd.nely)
                            e2 = (i2-1)*fd.nely+j2;
                            k = k+1;
                            iH(k) = e1;
                            jH(k) = e2;
                            %sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
                            sH(k) = 1;
                        end
                    end
                end
            end
            fd.data.H = sparse(iH,jH,sH);
            fd.data.Hs = sum(fd.data.H,2);           
            
        case 'grad' % further arguments: output_type, extra_arg
            error(nargchk(4,5,nargin));
            fd.data.out=varargin{4}; % 'mag' (magnitude), 'theta' (angle), 'gx', 'gy', 'dot'
            if strcmp(fd.data.out, 'dot')
                % dot product: next argument is vector
                fd.data.vec=varargin{5};
            end
            nely=fd.nely; nelx=fd.nelx; ne=nelx*nely;
            % Sobel stencils:
            gx=[-1 -1; -1 0; -1 1; 1 -1; 1 0; 1 1];
            gy=[-1 -1; 0 -1; 1 -1; -1 1; 0 1; 1 1];
            wgt=-[.25 .5 .25 -.25 -.5 -.25]';            
            tmpx=repmat(1:nelx,6*nely,1); tmpx=tmpx(:); 
            tmpy=repmat(1:nely,6,1); tmpy=repmat(tmpy(:),nelx,1);
            iHx=tmpx+repmat(gx(:,1),ne,1); iHy=tmpx+repmat(gx(:,2),ne,1);
            jHx=tmpy+repmat(gy(:,1),ne,1); jHy=tmpy+repmat(gy(:,2),ne,1);
            sHx=repmat(wgt,ne,1); sHy=repmat(wgt,ne,1);
            % handle boundaries
            for i=1:4 % NESW
                switch i
                    case 1 % N: no y=0
                        kx=find(jHx==0); ky=find(jHy==0);
                        jHx(kx)=1; jHy(ky)=1; % reflect from inside domain
                        sHy([ky,ky+3])=2*sHy([ky,ky+3]); % double weights, also those already inside the domain
                    case 2 % E: no x=nelx+1
                        kx=find(iHx==(nelx+1)); ky=find(iHy==(nelx+1));
                        iHx(kx)=nelx; iHy(ky)=nelx; % reflect from inside domain
                        sHx([kx,kx-3])=2*sHx([kx,kx-3]); % double weights, also those already inside the domain
                    case 3 % S: no y=nely+1
                        kx=find(jHx==(nely+1)); ky=find(jHy==(nely+1));
                        jHx(kx)=nely; jHy(ky)=nely; % reflect from inside domain
                        sHy([ky,ky-3])=2*sHy([ky,ky-3]); % double weights, also those already inside the domain
                    case 4 % W: no x=0
                        kx=find(iHx==0); ky=find(iHy==0);
                        iHx(kx)=1; iHy(ky)=1; % reflect from inside domain
                        sHx([kx,kx+3])=2*sHx([kx,kx+3]); % double weights, also those already inside the domain
                end
            end
            tmp1=repmat(1:ne,6,1);
            fd.data.Hx = sparse(tmp1(:),nely*(iHx-1)+jHx,sHx);
            fd.data.Hy = sparse(tmp1(:),nely*(iHy-1)+jHy,sHy);
            
        case 'KSmask' % further arguments: P, MaskOffsetsX, MaskOffsetsY
            error(nargchk(6,6,nargin));
            fd.data.P=varargin{4};
            mox=varargin{5}(:); fd.data.mox=mox;
            moy=varargin{6}(:); fd.data.moy=moy;
            fd.data.N=length(mox);
            % determine ghost layer dimensions:
            if 0
                gx=[min(0,min(mox)), max(fd.nelx,fd.nelx+max(mox))]; Lx=diff(gx);
                gy=[min(0,min(moy)), max(fd.nely,fd.nely+max(moy))]; Ly=diff(gy);
            else
                % adjust for flipping: 
                gx=[min([0,min(mox),min(-mox)]), max([fd.nelx,fd.nelx+max(mox),fd.nelx+max(-mox)])]; Lx=diff(gx);
                gy=[min([0,min(moy),min(-moy)]), max([fd.nely,fd.nely+max(moy),fd.nely+max(-moy)])]; Ly=diff(gy);
            end
            if (fd.data.P<0) ghostval=1; else ghostval=0; end % adjust these bounds as required! depends on ranges of input and output data.
            fd.data.gx=gx; fd.data.Lx=Lx; fd.data.x0=max(0,-gx(1)); % 0-based, still '1' to be added for indexing
            fd.data.gy=gy; fd.data.Ly=Ly; fd.data.y0=max(0,-gy(1));
            fd.data.ghostval=ghostval;
            % prepare mask indexing on embedded domain:
            fd.data.mi=moy+mox*Ly;
            % for sensitivities:            
            % flip mask in [0,0]:
            fd.data.mi2=-moy-mox*Ly;
            
        case 'mean' % further arguments: radius
            error(nargchk(4,4,nargin));
            fd.data.rmin=varargin{4};
            fd=prep2(fd);
            
        case 'mean2' % further arguments: radius
            error(nargchk(4,4,nargin));
            fd.data.rmin=varargin{4};
            rmin=fd.data.rmin;
            iH = ones(fd.nelx*fd.nely*(2*(ceil(rmin)-1)+1)^2,1);
            jH = ones(size(iH));
            sH = zeros(size(iH));
            k = 0;
            for i1 = 1:fd.nelx
                for j1 = 1:fd.nely
                    e1 = (i1-1)*fd.nely+j1;
                    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),fd.nelx)
                        for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),fd.nely)
                            e2 = (i2-1)*fd.nely+j2;
                            k = k+1;
                            iH(k) = e1;
                            jH(k) = e2;
                            %sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
                            sH(k) = 1;
                        end
                    end
                end
            end
            fd.data.H = sparse(iH,jH,sH);
            fd.data.Hs = sum(fd.data.H,2);           
            
        case {'maxfs1','maxfs2'} % further arguments: radius, a, b
            % maximum feature size
            error(nargchk(6,6,nargin));
            fd.data.rmin=varargin{4}; rmin=fd.data.rmin;
            fd.data.a=varargin{5};
            fd.data.b=varargin{6};
            iH = ones(fd.nelx*fd.nely*(2*(ceil(rmin)-1)+1)^2,1);
            jH = ones(size(iH));
            sH = zeros(size(iH));
            k = 0;
            for i1 = 1:fd.nelx
                for j1 = 1:fd.nely
                    e1 = (i1-1)*fd.nely+j1;
                    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),fd.nelx)
                        for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),fd.nely)
                            e2 = (i2-1)*fd.nely+j2;
                            k = k+1;
                            iH(k) = e1;
                            jH(k) = e2;
                            %sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
                            sH(k) = 1;
                        end
                    end
                end
            end
            fd.data.H = sparse(iH,jH,sH);
            fd.data.Hs = sum(fd.data.H,2);  
            
        case 'maxfs2line' % further arguments: radius (= half line length), a, b, orientation (- | / \)
            % maximum feature size
            error(nargchk(7,7,nargin));
            fd.data.rmin=varargin{4}; rmin=fd.data.rmin;
            fd.data.a=varargin{5};
            fd.data.b=varargin{6};
            fd.data.orientation = varargin{7};
            %iH = ones(fd.nelx*fd.nely*(2*(ceil(rmin)-1)+1)^2,1);
            iH = ones(fd.nelx*fd.nely*(2*(ceil(rmin)-1)+1),1); % no ^2
            jH = ones(size(iH));
            sH = zeros(size(iH));
            k = 0;
            for i1 = 1:fd.nelx
                for j1 = 1:fd.nely
                    e1 = (i1-1)*fd.nely+j1;
                    switch fd.data.orientation
                        case '-'
                            for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),fd.nelx)
                                for j2 = j1
                                    e2 = (i2-1)*fd.nely+j2;
                                    k = k+1;
                                    iH(k) = e1;
                                    jH(k) = e2;
                                    sH(k) = 1;
                                end
                            end
                        case '|'
                            for i2 = i1
                                for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),fd.nely)
                                    e2 = (i2-1)*fd.nely+j2;
                                    k = k+1;
                                    iH(k) = e1;
                                    jH(k) = e2;
                                    sH(k) = 1;
                                end
                            end
                        case '/'
                            xsteps=[max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),fd.nelx)]-i1;
                            ysteps=min(fd.nely,max(j1+xsteps, 1));
                            for i2 = xsteps+i1
                                %for j2 = max(j1+xdiff,fd.nely):-1:min(j1-xdiff,1)
                                for j2 = ysteps
                                    e2 = (i2-1)*fd.nely+j2;
                                    k = k+1;
                                    iH(k) = e1;
                                    jH(k) = e2;
                                    sH(k) = 1;
                                end
                            end
                        case '\'
                            xsteps=[max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),fd.nelx)]-i1;
                            ysteps=min(fd.nely,max(j1-xsteps, 1));
                            for i2 = xsteps+i1
                                for j2 = ysteps
                                    e2 = (i2-1)*fd.nely+j2;
                                    k = k+1;
                                    iH(k) = e1;
                                    jH(k) = e2;
                                    sH(k) = 1;
                                end
                            end
                        otherwise
                            error('wrong orientation');
                    end
                end
            end
            fd.data.H = sparse(iH,jH,sH);
            fd.data.Hs = sum(fd.data.H,2);  
        
        case 'mill' % further arguments: P
            % milling from above (cumsum approach)
            error(nargchk(4,4,nargin));
            fd.data.P=varargin{4};
            fd.data.N=1;
            
        case 'mill2' % further arguments: direction
            % milling from 1 direction (1..4)
            error(nargchk(4,4,nargin));
            fd.data.dir = varargin{4};
            
        case 'mill3' % further arguments: P
            % milling from 2 directions, P-normed
            error(nargchk(4,4,nargin));
            fd.data.P=-varargin{4}; % note: taking the negative!
            fd.data.N=2;
            
        case 'mill4' % further arguments: P, directions
            % milling from multiple directions (45-multiples)
            error(nargchk(5,5,nargin));
            fd.data.P=-varargin{4}; % note: taking the negative!
            fd.data.dir=varargin{5};
            fd.data.N=length(fd.data.dir);
            fd.data.csi=[1 2 1 2]; % index
            fd.data.csd=[1 1 2 2]; % direction
            fd.data.css={'forward','reverse'}; % strings            
            
        case 'mill5' % further arguments: direction
            % difference from 4: parameterization. This has a design var
            % linked to each mill insertion position (row of elements)
            error(nargchk(4,4,nargin));
            fd.data.dir = varargin{4};            
            M=sparse([]);
            for i=1:length(fd.data.dir)
                tmp=newmilltest1(fd.nelx,fd.nely,fd.data.dir(i));
                M=[M, tmp];
            end
            fd.data.M = M;
            
        case 'mill5a' % further arguments: P, direction
            % difference: no summing, but local P-norming.
            error(nargchk(5,5,nargin));
            fd.data.P = varargin{4}; 
            fd.data.dir = varargin{5};
            fd.data.ndir = length(fd.data.dir);
            fd.data.M=cell(fd.data.ndir,1);
            for i=1:fd.data.ndir 
                fd.data.M{i}=newmilltest1(fd.nelx,fd.nely,fd.data.dir(i));
            end
            
        case 'mill*' % further arguments: P, angles
            % as 4, but with arbitary angles (uses mapping)
            error(nargchk(5,5,nargin));
            fd.data.P = -varargin{4}; % note: taking the negative!
            fd.data.angles = -varargin{5} +90;  % adjusted to give logical angles: 0 = from right, 90 is from top
            % form milling matrices:
            fd.data.N=length(fd.data.angles);
            fd.data.M=cell(fd.data.N,1);
            for i=1:fd.data.N
                [F,B,nemy,nemx]=FieldRotationMappings(fd.nely,fd.nelx,fd.data.angles(i));
                C = CumSumMatrix(nemy,nemx);
                if 0 % some toollength/shape testing:
                    TL1 = ToolLengthMatrix(nemy,nemx,50) + speye(nemx*nemy);
                    TL = ToolShapeMatrix(nemy,nemx,[20],[10]); % L, D
                    TL=TL*TL;
                else
                    TL=speye(nemx*nemy);
                end
                %fd.data.M{i}=B*C*F;
                fd.data.M{i}=B*C*TL*F;
            end
            
            %% NOTE: mill*2 and ToolShapeMatrix2 are failed experiments
        case 'mill*2' % further arguments: P, angles, LD (Nx2 array of Li, Di)
            % as *, but with arbitary angles (uses mapping)
            error(nargchk(6,6,nargin));
            fd.data.P = -varargin{4}; % note: taking the negative!
            fd.data.angles = -varargin{5} +90;  % adjusted to give logical angles: 0 = from right, 90 is from top
            LD=varargin{6}; 
            fd.data.L=LD(:,1);
            fd.data.D=LD(:,2);
            fd.data.Power=3;
            % form milling matrices:
            fd.data.N=length(fd.data.angles);
            fd.data.M=cell(fd.data.N,1);
            for i=1:fd.data.N
                [F,B,nemy,nemx]=FieldRotationMappings(fd.nely,fd.nelx,fd.data.angles(i));
                C = CumSumMatrix(nemy,nemx);
                %T = ToolShapeMatrix(nemy,nemx,fd.data.L,fd.data.D);
                %T = T/mean(sum(T,2)); % normalize a bit wrt filter entries
                T = ToolShapeMatrix3(nemy,nemx,fd.data.L,fd.data.D);
                fd.data.M{i}=B*C*T*F;
                % multipass:
                fd.data.M{i}=B*C*T*C*T*F;
                %fd.data.M{i}=B*C*T*C*T*C*T*C*T*F;
                
                %fd.data.M{i}=B*C*F;
                %fd.data.M{i}=B*C*T'*F;
                %fd.data.M{i}=B*T*C*F;
                %T2 = ToolShapeMatrix2(nemy,nemx,fd.data.L,fd.data.D);
                %fd.data.M{i}=B*T2*C*T*F;
            end
            % visualize:
            if 0
                x=zeros(nemy,nemx); for i=1:max(nemx,nemy), x(1+round(rand*nemx*nemy-1))=1; end
                figure; 
                t=0*x; t(round(nemy/2),round(nemx/2))=1; t2=t; t2(:)=C'*T'*t(:);
                subplot(3,2,1); imagesc(t2); axis image; title('Tool');                
                subplot(3,2,2); imagesc(x); axis image; title('x');
                Tx=x; Tx(:)=T*x(:);
                subplot(3,2,3); imagesc(Tx); axis image; title('T*x');
                CTx=Tx; CTx(:)=C*Tx(:);
                subplot(3,2,4); imagesc(CTx); axis image; title('C*T*x');
                CTCTx=CTx; CTCTx(:)=C*T*CTx(:);
                subplot(3,2,5); imagesc(CTCTx); axis image; title('C*T*C*T*x');
                CTCTCTx=CTCTx; CTCTCTx(:)=C*T*CTCTx(:);
                subplot(3,2,6); imagesc(CTCTCTx); axis image; title('C*T*C*T*C*T*x');
                setappdata(gcf,'t',t);
                setappdata(gcf,'x',x);
                setappdata(gcf,'Tx',Tx);
                setappdata(gcf,'CTx',CTx);
                setappdata(gcf,'T',T);
                setappdata(gcf,'C',C);      
                % figure; A=C*T; figure; for i=1:3, subplot(1,3,i); spy(full(A)); A=A*A; end
            end
            
        case {'mill*3','mill*3P'} % further arguments: P, angles, Pmask, MaskOffsetsX, MaskOffsetsY
            % Pnorm-variant: P, angles, Pmask, MaskOffsetsX, MaskOffsetsY, epsilon, m
            % as *, but with tool shape (KSmask/Pmask)
            error(nargchk(8,10,nargin));
            if ((nargin==8) && (fd.type(end)~='3')) error('wrong nr inputs'); end
            if ((nargin==10) && (fd.type(end)~='P')) error('wrong nr inputs'); end
            fd.data.P = -varargin{4}; % note: taking the negative!
            fd.data.angles = -varargin{5} +90;  % adjusted to give logical angles: 0 = from right, 90 is from top
            Pmask=varargin{6};
            mox=varargin{7}(:); 
            moy=varargin{8}(:); 
            usePnorm=0;
            if nargin==10
                fd.data.eps=varargin{9};
                fd.data.m=varargin{10};
                usePnorm=1;
            end
            fd.data.Power=2;%3;
            %
            % determine ghost layer dimensions:
            % form milling matrices:
            fd.data.N=length(fd.data.angles); N=fd.data.N;
            fd.data.F=cell(N,1);
            fd.data.B=cell(N,1);
            fd.data.C=cell(N,1);
            fd.data.nemx=zeros(N,1);
            fd.data.nemy=zeros(N,1);            
            fd.data.KSmask=cell(N,1);
            %
            for i=1:N
                [F,B,nemy,nemx]=FieldRotationMappings(fd.nely,fd.nelx,fd.data.angles(i));
                fd.data.F{i}=F; fd.data.B{i}=B;
                fd.data.C{i} = CumSumMatrix(nemy,nemx);
                fd.data.nemy(i)=nemy;
                fd.data.nemx(i)=nemx;
                %fd.data.M{i}=B*C*F;
                if usePnorm
                    fd.data.KSmask{i}=anyfilter('Pmask',nemx,nemy,Pmask,fd.data.eps,fd.data.m,mox,moy); % store P-norm filter data
                else
                    fd.data.KSmask{i}=anyfilter('KSmask',nemx,nemy,Pmask,mox,moy); % store KS filter data
                end
            end
            
        case {'mill*4','mill*4P'} % further arguments: epsilon_angles, m_angles, angles, Pmask, MaskOffsetsX, MaskOffsetsY
            % Pnorm-variant:                           epsilon_angles, m_angles, angles, Pmask, MaskOffsetsX, MaskOffsetsY, epsilon, m
            % as *, but with tool shape (KSmask/Pmask)
            error(nargchk(9,11,nargin));
            if ((nargin==9) && (fd.type(end)~='4')) error('wrong nr inputs'); end
            if ((nargin==11) && (fd.type(end)~='P')) error('wrong nr inputs'); end
            fd.data.epsilon_angles = varargin{4}; 
            fd.data.m_angles = varargin{5}; 
            fd.data.angles = -varargin{6} +90;  % adjusted to give logical angles: 0 = from right, 90 is from top
            fd.data.N=length(fd.data.angles); N=fd.data.N;
            if isstr(fd.data.m_angles)
                if fd.data.m_angles=='N'
                    fd.data.m_angles=fd.data.N;
                else
                    error('wrong input: m_angles');
                end
            end
            Pmask=varargin{7};
            mox=varargin{8}(:); 
            moy=varargin{9}(:); 
            usePnorm=0;
            if nargin==11
                fd.data.eps=varargin{10};
                fd.data.m=varargin{11};
                if isstr(fd.data.m)
                    if fd.data.m=='N'
                        fd.data.m=fd.data.N;
                    else
                        error('wrong input: m');
                    end
                end                
                usePnorm=1;
            end
       fd.data.Power=3;%3;
            %
            % determine ghost layer dimensions:
            % form milling matrices:
            fd.data.N=length(fd.data.angles); N=fd.data.N;
            fd.data.F=cell(N,1);
            fd.data.B=cell(N,1);
            fd.data.C=cell(N,1);
            fd.data.nemx=zeros(N,1);
            fd.data.nemy=zeros(N,1);            
            fd.data.KSmask=cell(N,1);
            %
            for i=1:N
                [F,B,nemy,nemx]=FieldRotationMappings(fd.nely,fd.nelx,fd.data.angles(i));
                fd.data.F{i}=F; fd.data.B{i}=B;
                fd.data.C{i} = CumSumMatrix(nemy,nemx);
                fd.data.nemy(i)=nemy;
                fd.data.nemx(i)=nemx;
                %fd.data.M{i}=B*C*F;
                if usePnorm
                    fd.data.KSmask{i}=anyfilter('Pmask',nemx,nemy,Pmask,fd.data.eps,fd.data.m,mox,moy); % store P-norm filter data
                else
                    fd.data.KSmask{i}=anyfilter('KSmask',nemx,nemy,Pmask,mox,moy); % store KS filter data
                end
            end
            
        case 'millX' % combination of mill*'s TL (tool length) and mill*3
            % further arguments: P, angles, Pmask, 
            % ToolMaskOffsetsX, ToolMaskOffsetsY, 
            % ToolLength, PToolLength, 
            % Power
            error(nargchk(11,11,nargin));
            fd.data.P = -varargin{4}; % note: taking the negative!
            fd.data.angles = -varargin{5} +90;  % adjusted to give logical angles: 0 = from right, 90 is from top
            Pmask=varargin{6};
            mox=varargin{7}(:); 
            moy=varargin{8}(:); 
            fd.data.TL=varargin{9};
            fd.data.PTL=varargin{10};            
            fd.data.Power=varargin{11}; %3;
            %
            % determine ghost layer dimensions:
            % form milling matrices:
            fd.data.N=length(fd.data.angles); N=fd.data.N;
            fd.data.F=cell(N,1);
            fd.data.B=cell(N,1);
            fd.data.C=cell(N,1);
            fd.data.nemx=zeros(N,1);
            fd.data.nemy=zeros(N,1);            
            fd.data.KSmask=cell(N,1);
            fd.data.topi=cell(N,1);
            fd.data.topj=cell(N,1);
            %
            for i=1:N
                [F,B,nemy,nemx]=FieldRotationMappings(fd.nely,fd.nelx,fd.data.angles(i));
                fd.data.F{i}=F; fd.data.B{i}=B;
                fd.data.C{i} = CumSumMatrix(nemy,nemx);
                fd.data.nemy(i)=nemy;
                fd.data.nemx(i)=nemx;
                %fd.data.M{i}=B*C*F;
                if ~isempty(mox)
                    fd.data.KSmask{i}=anyfilter('KSmask',nemx,nemy,Pmask,mox,moy); % store KS filter data
                    % also determine indcies of the 'top elements' of the
                    % inner domain:
                    tmp=diff(reshape(sum(F,2),nemy,nemx));
                    [r,c]=find(tmp>0.9);
                    fd.data.topi{i}=r+1; fd.data.topj{i}=c; % i->y, j->x
                else
                    fd.data.KSmask{i}=anyfilter('none',nemx,nemy);
                end
            end
            
        case 'none'
            ;
            
        case 'oneminusx'
            ;
            
        case 'Pmask' % further arguments: P, epsilon, m, MaskOffsetsX, MaskOffsetsY
            error(nargchk(8,8,nargin));
            fd.data.P=varargin{4};
            fd.data.eps=varargin{5};
            fd.data.m=varargin{6};
            mox=varargin{7}(:); fd.data.mox=mox;
            moy=varargin{8}(:); fd.data.moy=moy;
            fd.data.N=length(mox);
            if isstr(fd.data.m)
                if fd.data.m=='N'
                    fd.data.m=fd.data.N;
                end
            end
            % determine ghost layer dimensions:
            if 0
                gx=[min(0,min(mox)), max(fd.nelx,fd.nelx+max(mox))]; Lx=diff(gx);
                gy=[min(0,min(moy)), max(fd.nely,fd.nely+max(moy))]; Ly=diff(gy);
            else
                % adjust for flipping: 
                gx=[min([0,min(mox),min(-mox)]), max([fd.nelx,fd.nelx+max(mox),fd.nelx+max(-mox)])]; Lx=diff(gx);
                gy=[min([0,min(moy),min(-moy)]), max([fd.nely,fd.nely+max(moy),fd.nely+max(-moy)])]; Ly=diff(gy);
            end
            if (fd.data.P<0) ghostval=1; else ghostval=0; end % adjust these bounds as required! depends on ranges of input and output data.
            fd.data.gx=gx; fd.data.Lx=Lx; fd.data.x0=max(0,-gx(1)); % 0-based, still '1' to be added for indexing
            fd.data.gy=gy; fd.data.Ly=Ly; fd.data.y0=max(0,-gy(1));
            fd.data.ghostval=ghostval;
            % prepare mask indexing on embedded domain:
            fd.data.mi=moy+mox*Ly;
            % for sensitivities:            
            % flip mask in [0,0]:
            %%fd.data.mi2=-moy-mox*Ly; %BS dit is gewoon -m
            
        case 'range' % further arguments: radius, beta
            error(nargchk(5,5,nargin));
            fd.data.rmin=varargin{4};
            fd.data.beta=varargin{5};
            fd.data.dilate = anyfilter('dilate',fd.nelx,fd.nely,fd.data.rmin,fd.data.beta);
            fd.data.erode  = anyfilter('erode',fd.nelx,fd.nely,fd.data.rmin,fd.data.beta);
            
        case 'stdev' % further arguments: radius
            error(nargchk(4,4,nargin));
            fd.data.rmin=varargin{4};
            rmin=fd.data.rmin;
            iH = ones(fd.nelx*fd.nely*(2*(ceil(rmin)-1)+1)^2,1);
            jH = ones(size(iH));
            sH = zeros(size(iH));
            k = 0;
            for i1 = 1:fd.nelx
                for j1 = 1:fd.nely
                    e1 = (i1-1)*fd.nely+j1;
                    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),fd.nelx)
                        for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),fd.nely)
                            e2 = (i2-1)*fd.nely+j2;
                            k = k+1;
                            iH(k) = e1;
                            jH(k) = e2;
                            %sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
                            sH(k) = 1;
                        end
                    end
                end
            end
            fd.data.H = sparse(iH,jH,sH);
            fd.data.Hs = sum(fd.data.H,2);           
            
        case 'symm' % further arguments: xy
            error(nargchk(4,4,nargin));
            fd.data.xy=varargin{4};
            
        otherwise
            error(['filter type ',fd.type,' not implemented']);
    end
    varargout{1}=fd;
else
    if isstr(varargin{2})
        % another filter definition
        % [filterdata] = anyfilter(filterdata, type, nelx, nely, parameters)
        fd=varargin{1};
        % check if domain size is equal - this is required
        if all(fd(1).nelx==varargin{3}, fd(1).nely==varargin{4})
            fd2=anyfilter(varargin{2:end});
            varargout={[fd,fd2]};
        else
            error('domain size mismatch');
        end
    else        
        % design or sensitivity filtering
        fd=varargin{1};
        x=varargin{2};        
        if nargout~=2
            error('incorrect number of output variables');
        end
        Nfd=length(fd);
            
        switch nargin
            case 2 % 2. FORWARD FILTERING
                % [filterdata, x_filtered] = anyfilter(filterdata, x)
                if Nfd>1
                    fd_in=fd;
                    fd=fd(1);
                end
                fd.x_in=x; % store for later use in sensitivity filtering
                
                switch fd.type
                    case 'HS'
                        beta=fd.data.beta; eta=fd.data.eta; num=fd.data.num;
                        xf =(tanh(beta*eta)+tanh(beta*(x - eta)))/num;

                    case 'HS01'
                        beta=fd.data.beta; eta=0.5;
                        xf =0.5 + 0.5*tanh(beta*(x - eta));

                    case 'HS01flipped'
                        beta=fd.data.beta; eta=0.5;
                        xf =0.5 - 0.5*tanh(beta*(x - eta));

                    case 'density'                        
                        xf=zeros(fd.nely,fd.nelx);
                        xf(:) = (fd.data.H*x(:))./fd.data.Hs;
                        
                    case 'dilate'
                        % pick up maximum using KS:
                        dd=fd.data.dd; beta=fd.data.beta; Ns=fd.data.Ns; ii=fd.data.ii;
                        % embed in zero padding:                        
                        x2=zeros(fd.nely+2*dd,fd.nelx+2*dd);
                        x2(dd+[1:fd.nely],dd+[1:fd.nelx])=x;
                        xf=zeros(size(x));
                        ny=fd.nely+2*dd;
                        for i=1:fd.nelx
                            for j=1:fd.nely
                                %KS version:
                                xf(j,i)= log( (1/Ns(j,i))*sum( exp(beta*x2(dd+j+ny*(dd+i-1)+ii))) )/beta;
                            end
                        end
                        
                    case {'dmean', 'dmeanline'}
                        beta=fd.data.beta; eta=fd.data.eta; a=fd.data.a; b=fd.data.b; 
                        num=tanh(beta*eta)+tanh(beta*(1-eta));
                        xf=zeros(fd.nely,fd.nelx); xmean=xf;
                        xmean(:) = (fd.data.H*x(:))./fd.data.Hs;
                        % apply knock-down function locally:
                        S=(tanh(beta*eta)+tanh(beta*(xmean-eta)))/num;
                        f=1-a*S-b*xmean.*S;
                        xf=x.*f;

                    case 'erode'
                        % pick up minimum using KS:
                        dd=fd.data.dd; beta=fd.data.beta; Ns=fd.data.Ns; ii=fd.data.ii;
                        % embed in unity padding:
                        x2=ones(fd.nely+2*dd,fd.nelx+2*dd); % ones: do not erode from outside of domain
                        x2(dd+[1:fd.nely],dd+[1:fd.nelx])=x;
                        xf=zeros(size(x));
                        ny=fd.nely+2*dd;
                        for i=1:fd.nelx
                            for j=1:fd.nely
                                %KS version:
                                xf(j,i)= 1 - log( (1/Ns(j,i))*sum( exp(beta*(1-x2(dd+j+ny*(dd+i-1)+ii)))) )/beta;
                            end
                        end
                        
                    case 'gmean'
                        P=fd.data.P; beta=fd.data.beta; eta=fd.data.eta; a=fd.data.a; b=fd.data.b; 
                        num=tanh(beta*eta)+tanh(beta*(1-eta));
                        xf=zeros(fd.nely,fd.nelx); xmean=xf;
                        xmean(:) = (fd.data.H*x(:))./fd.data.Hs;
                        % aggregate using P-norm:
                        maxmean=(mean(xmean(:).^P))^(1/P);
                        % apply knock-down function:
                        S=(tanh(beta*eta)+tanh(beta*(maxmean-eta)))/num;
                        f=1-a*S-b*maxmean*S;
                        xf=x*f;                        
                        % SIMPLER for testing:
                        if 0
                            xf=x*mean(xmean(:));
                        end         
                        
                    case 'grad'
                        xf=zeros(fd.nely,fd.nelx); 
                        gx=xf; gy=xf;
                        gx(:) = fd.data.Hx*x(:);
                        gy(:) = fd.data.Hy*x(:);
                        % fd.data.out=varargin{5}; % 'mag' (magnitude), 'theta' (angle), 'gx', 'gy'
                        switch fd.data.out
                            case 'mag'
                                xf=sqrt(gx.^2 + gy.^2);
                            case 'theta'
                                xf=atan2(gy,gx);
                            case 'gx'
                                xf=gx;
                            case 'gy'
                                xf=gy;
                            case 'dot'
                                xf=fd.data.vec(1)*gx + fd.data.vec(2)*gy;
                        end                        
                        
                    case 'KSmask'
                        %mox=fd.data.mox; moy=fd.data.moy; gx=fd.data.gx; gy=fd.data.gy; 
                        N=fd.data.N; P=fd.data.P;                        
                        Lx=fd.data.Lx; Ly=fd.data.Ly;
                        x0=fd.data.x0; y0=fd.data.y0; mi=fd.data.mi; % mask indices (single-vector)
                        nely=fd.nely; nelx=fd.nelx;                        
                        % embed x in padding:
                        x2=fd.data.ghostval*ones(Ly,Lx);
                        x2(y0+[1:nely],x0+[1:nelx])=x;
                        xf=zeros(size(x));
                        for i=1:nelx
                            for j=1:nely
                                %KS version:
                                xf(j,i)= log(realmin+ (1/N)*sum( exp(P*( x2( (x0+i-1)*Ly+(y0+j)+mi )) )) )/P;
                            end
                        end
                            
                    case 'mean'
                        % use weights to deal with ghost boundary
                        dd=fd.data.dd; Ns=fd.data.Ns; ii=fd.data.ii; w=fd.data.w;
                        % embed in padding:
                        x2=zeros(fd.nely+2*dd,fd.nelx+2*dd); 
                        x2(dd+[1:fd.nely],dd+[1:fd.nelx])=x;                                                
                        xf=zeros(size(x)); % this will hold the mean values
                        ny=fd.nely+2*dd;
                        for i=1:fd.nelx
                            for j=1:fd.nely                                
                                %%xf(j,i)= sum( w(dd+j+ny*(dd+i-1)+ii).*x2(dd+j+ny*(dd+i-1)+ii) )/Ns(j,i);
                                % in fact multiplication with w does
                                % nothing, since 'outside' is already zero
                                % in x2.
                                % this is equivalent:
                                xf(j,i)= sum( x2(dd+j+ny*(dd+i-1)+ii) )/Ns(j,i);
                            end
                        end         
                        
                    case 'mean2'
                        xf=zeros(fd.nely,fd.nelx);
                        xf(:) = (fd.data.H*x(:))./fd.data.Hs;
                        
                    case 'maxfs1'
                        % maximum feature size                        
                        % standard deviation within a filter radius.
                        xmean=zeros(size(x)); 
                        msd=zeros(size(x)); 
                        stdev=zeros(size(x)); 
                        xf=zeros(size(x)); 
                        xmean(:) = (fd.data.H*x(:))./fd.data.Hs;                        
                        dev=x-xmean;
                        dev2=dev.^2;                           
                        msd(:) = (fd.data.H*dev2(:))./fd.data.Hs; % mean of squared deviations
                        stdev=sqrt(msd); 
                        % apply knockdown function:
                        a=fd.data.a; b=fd.data.b;
                        xf=x.*(a+(1-a)*tanh(b*stdev));                        
                        
                    case 'maxfs2'
                        % maximum feature size                        
                        % msd within a filter radius.
                        xmean=zeros(size(x)); 
                        msd=zeros(size(x)); 
                        xf=zeros(size(x)); 
                        xmean(:) = (fd.data.H*x(:))./fd.data.Hs;                        
                        dev=x-xmean;
                        dev2=dev.^2;                           
                        msd(:) = (fd.data.H*dev2(:))./fd.data.Hs; % mean of squared deviations
                        % apply knockdown function:
                        a=fd.data.a; b=fd.data.b;
                        xf=x.*(a+(1-a)*tanh(b*msd));
                        
                    case 'maxfs2line'
                        % maximum feature size. Copy of maxfs2        
                        % msd within a line segment.
                        xmean=zeros(size(x)); 
                        msd=zeros(size(x)); 
                        xf=zeros(size(x)); 
                        xmean(:) = (fd.data.H*x(:))./fd.data.Hs;                        
                        dev=x-xmean;
                        dev2=dev.^2;                           
                        msd(:) = (fd.data.H*dev2(:))./fd.data.Hs; % mean of squared deviations
                        % apply knockdown function:
                        a=fd.data.a; b=fd.data.b;
                        xf=x.*(a+(1-a)*tanh(b*msd));
                        
                    case 'mill'
                        % milling from above:
                        xf=zeros(size(x));
                        xf(1,:)=x(1,:);
                        P=fd.data.P;
                        for i=2:fd.nely
                            xf(i,:)=((x(i,:).^P+xf(i-1,:).^P)/fd.data.N).^(1/P);
                        end
                        
                    case 'mill2'
                        % milling from one direction:
                        switch fd.data.dir
                            case 1 % from above
                                xf=cumsum(x,1);
                            case 2 % from left
                                xf=cumsum(x,2);
                            case 3 % from below
                                xf=cumsum(x,1,'reverse');
                            case 4 % from right
                                xf=cumsum(x,2,'reverse');
                            otherwise
                                error('wrong direction argument');
                        end
                        
                        
                    case 'mill3'
                        % milling from two directions: TEST
                        xf=cumsum(x);
                        xf2=cumsum(x,2);
                        % overflow protection:
                        xf(xf>1)=log(exp(1)*xf(xf>1));
                        xf2(xf2>1)=log(exp(1)*xf2(xf2>1));
                        P=fd.data.P;
                        N=fd.data.N;
                        xf=log( (exp(P*xf)+exp(P*xf2))/N )/P;
                        
                    case 'mill4'
                        % milling from various directions: 
                        P=fd.data.P;
                        N=fd.data.N;
                        KS=zeros(size(x));
                        for d=1:N
                            dr=fd.data.dir(d);
                            if dr<=4
                                xf=cumsum(x,fd.data.csi(dr),fd.data.css{fd.data.csd(dr)});
                            else
                                xf=cumsum45(x,dr-4,'forward');
                            end
%                             switch fd.data.dir(d)
%                                 case 1 % from above
%                                     xf=cumsum(x,1);
%                                 case 2 % from left
%                                     xf=cumsum(x,2);
%                                 case 3 % from below
%                                     xf=cumsum(x,1,'reverse');
%                                 case 4 % from right
%                                     xf=cumsum(x,2,'reverse');
%                                 otherwise
%                                     error('wrong direction argument');
%                             end
                            % overflow protection:
                            xf(xf>1)=log(exp(1)*xf(xf>1));
                            % accumulate
                            KS=KS+exp(P*xf);
                        end
                        xf=log(KS/N)/P;
                        
                    case 'mill5'
                        xf=zeros(fd.nely,fd.nelx);
                        xf(:)=fd.data.M*x(:);
                        
                    case 'mill5a'
                        xf=zeros(fd.nely,fd.nelx);
                        x=x(:);
                        P=fd.data.P;
                        N=fd.nely*fd.nelx;                        
                        XF=zeros(N, fd.data.ndir);
                        for i=1:fd.data.ndir
                            % apply separate transformations
                            XF(:,i)=fd.data.M{i}*x(N*(i-1) + [1:N]);
                        end
                        % apply P-norm on XF, combining rows
                        xf(:)=sum(XF.^P,2).^(1/P);
                        
                    case {'mill*', 'mill*2'}
%                         xf=zeros(fd.nely,fd.nelx);
%                         P=fd.data.P;
%                         N=fd.nely*fd.nelx;                        
%                         nd=fd.data.ndir;
%                         XF=zeros(N, nd);
%                         for i=1:nd
%                             % apply separate transformations
%                             XF(:,i)=fd.data.M{i}*x(:);
%                         end
%                         % apply P-norm on XF, combining rows
%                         xf(:)=sum((1/nd)*XF.^P,2).^(1/P);
%                         
                        % milling from various directions:                        
                        xf=zeros(fd.nely,fd.nelx);
                        P=fd.data.P;
                        N=fd.data.N;
                        KS=zeros(size(x));
              if isfield(fd.data,'Power')
              x=x.^fd.data.Power;
              end
                        for i=1:N
                            xf(:)=fd.data.M{i}*x(:);
                            % overflow protection:
                            xf(xf>1)=log(exp(1)*xf(xf>1));
                            % accumulate
                            KS=KS+exp(P*xf);
                        end
                        xf=log(KS/N)/P;
                        
                        
                    case {'mill*3','mill*3P'}
                        % milling from various directions, with tool shape:
                        xf=zeros(fd.nely,fd.nelx);
                        P=fd.data.P;
                        N=fd.data.N;
                        KS=zeros(size(x));
              if isfield(fd.data,'Power')
              x=x.^fd.data.Power;
              end
                        for i=1:N
                            %xf(:)=fd.data.M{i}*x(:);
                            xm=zeros(fd.data.nemy(i),fd.data.nemx(i));
                            
                            xm(:)=fd.data.F{i}*x(:); 
                            xm(:)=fd.data.C{i}*xm(:);                             
                            [ignore, tmp]=anyfilter(fd.data.KSmask{i},xm);
                            xf(:)=fd.data.B{i}*tmp(:);
                            
                            % overflow protection:
                            xf(xf>1)=log(exp(1)*xf(xf>1));
                            % accumulate
                            KS=KS+exp(P*xf);
                        end
                        xf=log(KS/N)/P;
                        
                    case {'mill*4','mill*4P'}
                        % milling from various directions, with tool shape:
                        xf=zeros(fd.nely,fd.nelx);
                        epa=fd.data.epsilon_angles;
                        ma=fd.data.m_angles;
                        Pa=-1;
                        N=fd.data.N;
                        PNormSum=zeros(size(x));
              if isfield(fd.data,'Power')
              x=x.^fd.data.Power;
              end
                        for i=1:N
                            %xf(:)=fd.data.M{i}*x(:);
                            xm=zeros(fd.data.nemy(i),fd.data.nemx(i));
                            
                            xm(:)=fd.data.F{i}*x(:); 
                            xm(:)=fd.data.C{i}*xm(:);                             
                            [ignore, tmp]=anyfilter(fd.data.KSmask{i},xm);
                            xf(:)=fd.data.B{i}*tmp(:);
                            
                            % overflow protection:
                            %xf(xf>1)=log(exp(1)*xf(xf>1));
                            % accumulate
                            %KS=KS+exp(P*xf);
                            PNormSum=PNormSum+(xf+epa).^Pa;
                        end
                        %xf=log(KS/N)/P;
                        xf= (PNormSum/ma).^(1/Pa) - epa;                        
                        
                    case 'millX'
                        % milling from various directions, with tool shape & length:
                        xf=zeros(fd.nely,fd.nelx);
                        P=fd.data.P;
                        N=fd.data.N;
                        PTL=fd.data.PTL; TL=fd.data.TL; 
                        KS=zeros(size(x));
                        x=x.^fd.data.Power;
                        for i=1:N
                            %xf(:)=fd.data.M{i}*x(:);
                            xm=zeros(fd.data.nemy(i),fd.data.nemx(i));
                            
                            xm(:)=fd.data.F{i}*x(:); 
                            if ~isempty(TL) % account for tool length
                                if (TL<fd.data.nemy(i))
                                    sumePx=sum(exp(PTL*xm),2);
                                    rowmax=(1/PTL)*log((1/fd.data.nemx(i))*sumePx);
                                    Crowmax=cumsum(rowmax);
                                    Crowmax((TL+1):end)=Crowmax(1:(end-TL)); Crowmax(1:TL)=0;
                                    fCrowmax=Crowmax.*exp(1-Crowmax);
                                    for j=1:fd.data.nemy(i)
                                        xm(j,:)=xm(j,:)+fCrowmax(j);
                                    end
                                end
                            end
                            xm(:)=fd.data.C{i}*xm(:);                             
                            % set outside region to 1, to prevent fake
                            % toolshape erosions
                            %x2=ones(size(x)); xm2=xm; xm2(:)=fd.data.F{i}*x2(:); xm=xm+10*(1-xm2);
                            for k=1:length(fd.data.topi{i})
                                bordervalue=xm(fd.data.topi{i}(k),fd.data.topj{i}(k));                                
                                xm(1:fd.data.topi{i}(k), fd.data.topj{i}(k) )=bordervalue;
                            end
                            [~, xm]=anyfilter(fd.data.KSmask{i},xm);
                            xf(:)=fd.data.B{i}*xm(:);
                            
                            % overflow protection:
                            xf(xf>1)=log(exp(1)*xf(xf>1));
                            % accumulate
                            KS=KS+exp(P*xf);
                        end
                        xf=log(KS/N)/P;
                    
                    case 'none'
                        xf=x;
                        
                    case 'oneminusx'                        
                        xf=1-x;
                        
                    case 'Pmask'
                        %mox=fd.data.mox; moy=fd.data.moy; gx=fd.data.gx; gy=fd.data.gy; 
                        N=fd.data.N; P=fd.data.P; ep=fd.data.eps; m=fd.data.m;
                        Lx=fd.data.Lx; Ly=fd.data.Ly;
                        x0=fd.data.x0; y0=fd.data.y0; mi=fd.data.mi; % mask indices (single-vector)
                        nely=fd.nely; nelx=fd.nelx;                        
                        % embed x in padding:
                        x2=fd.data.ghostval*ones(Ly,Lx);
                        x2(y0+[1:nely],x0+[1:nelx])=x;
                        xf=zeros(size(x));
                        for i=1:nelx
                            for j=1:nely
                                %P-norm version:
                                xf(j,i)= ( (1/m)*sum( (x2( (x0+i-1)*Ly+(y0+j)+mi) + ep).^P ) ).^(1.0/P) ;
                            end
                        end
                        xf=xf-ep; % backshift
                        
                    case 'range'
                        % maximum - minimum within filter radius.
                        [fd.data.dilate,xmax]=anyfilter(fd.data.dilate,x);
                        [fd.data.erode,xmin]=anyfilter(fd.data.erode,x);
                        xf=xmax-xmin;
                        
                    case 'stdev'
                        % standard deviation within a filter radius.
                        xmean=zeros(size(x)); 
                        xf=xmean;
                        xmean(:) = (fd.data.H*x(:))./fd.data.Hs;                        
                        dev=x-xmean;
                        dev2=dev.^2;                                  
                        xf(:) = (fd.data.H*dev2(:))./fd.data.Hs; % mean of squared deviations
                        xf=sqrt(xf);                                                
                        
                    case 'symm' % symmetry in x- or y-direction, halfway the domain
                        switch fd.data.xy
                            case 'x'
                                xf=0.5*(x+fliplr(x));
                            case 'y'
                                xf=0.5*(x+flipud(x));
                            otherwise
                                error('symm: wrong argument');
                        end
                end
                if Nfd>1
                    % filter with multiple steps: process next step(s)
                    [fd_part, xf] = anyfilter(fd_in(2:end), xf);
                    fd=[fd fd_part];
                end
                varargout={fd,xf};
                
            case 3 % 3. SENSITIVITY FILTERING
                % [filterdata, sens_filtered] = anyfilter(filterdata, x, sens)
                s= varargin{3};
                if Nfd>1
                    fd_in=fd;
                    fd=fd(end);                    
                    x=fd.x_in;
                end                

                switch fd.type
                    case 'HS'
                        dx = fd.data.beta*(1 - (tanh(fd.data.beta*(x-fd.data.eta))).^2)/fd.data.num;                        
                        sf=zeros(size(s));
                        sf(:)=s(:).*dx(:);
                    
                    case 'HS01'
                        beta=fd.data.beta; eta=0.5;
                        % xf =0.5 + 0.5*tanh(beta*(x - eta));
                        dx=0.5*beta*sech(beta*(x - eta)).^2;
                        sf=zeros(size(s));
                        sf(:)=s(:).*dx(:);                        
                        
                    case 'HS01flipped'
                        beta=fd.data.beta; eta=0.5;
                        % xf =0.5 - 0.5*tanh(beta*(x - eta));
                        dx=-0.5*beta*sech(beta*(x - eta)).^2;
                        sf=zeros(size(s));
                        sf(:)=s(:).*dx(:);                        
                        
                    case 'density'
                        sf=zeros(size(s));
                        sf(:) = fd.data.H*(s(:)./fd.data.Hs);
                        
                    case 'dilate'
                        % KS version:
                        dd=fd.data.dd; beta=fd.data.beta; Ns=fd.data.Ns; ii=fd.data.ii;
                        % embed x in zero padding:
                        x2=zeros(fd.nely+2*dd,fd.nelx+2*dd);
                        x2(dd+[1:fd.nely],dd+[1:fd.nelx])=x;
                        tmp=zeros(size(x));
                        ny=fd.nely+2*dd;
                        for i=1:fd.nelx
                            for j=1:fd.nely
                                tmp(j,i)= sum( exp(beta*x2(dd+j+ny*(dd+i-1)+ii)) );
                            end
                        end
                        tmp(:)=s(:)./tmp(:);
                        tmp2=zeros(fd.nely+2*dd,fd.nelx+2*dd);
                        tmp2(dd+[1:fd.nely],dd+[1:fd.nelx])=tmp;
                        sf=zeros(size(x));
                        for i=1:fd.nelx
                            for j=1:fd.nely
                                sf(j,i)= sum( tmp2(dd+j+ny*(dd+i-1)+ii) );
                            end
                        end
                        sf=sf.*exp(beta*x);
                        sf=sf(:);
                    
                    case {'dmean','dmeanline'}
                        % preparations: redo forward computations to have
                        % all quantities available
                        beta=fd.data.beta; eta=fd.data.eta; a=fd.data.a; b=fd.data.b; 
                        num=tanh(beta*eta)+tanh(beta*(1-eta));
                        xmean=zeros(fd.nely,fd.nelx); 
                        xmean(:) = (fd.data.H*x(:))./fd.data.Hs;
                        % apply knock-down function locally:
                        S=(tanh(beta*eta)+tanh(beta*(xmean-eta)))/num;
                        f=1-a*S-b*xmean.*S;
                        xf=x.*f;
                        % start of sensitivity part ---
                        sf=zeros(fd.nely,fd.nelx); sf1=sf; sf2=sf;
                        sf1(:)=s(:).*f(:);
                        dSdxmean=(beta*(sech(beta*(xmean-eta))).^2 )/num;
                        dydxmean=x.*( (-a-b*xmean).*dSdxmean -b*S );                        
                        tmp=s(:).*dydxmean(:);
                        sf2(:) = fd.data.H*(tmp(:)./fd.data.Hs);                        
                        sf=sf1+sf2;
                        sf=sf(:);
                        
                    case 'erode'
                        % KS version:
                        dd=fd.data.dd; beta=fd.data.beta; Ns=fd.data.Ns; ii=fd.data.ii;
                        % embed in unity padding:
                        x2=ones(fd.nely+2*dd,fd.nelx+2*dd); % ones: do not erode from outside of domain
                        x2(dd+[1:fd.nely],dd+[1:fd.nelx])=x;
                        tmp=zeros(size(x));
                        ny=fd.nely+2*dd;
                        for i=1:fd.nelx
                            for j=1:fd.nely
                                tmp(j,i)= sum( exp(beta*(1-x2(dd+j+ny*(dd+i-1)+ii))) );
                            end
                        end
                        tmp(:)=s(:)./tmp(:);
                        tmp2=zeros(fd.nely+2*dd,fd.nelx+2*dd);
                        tmp2(dd+[1:fd.nely],dd+[1:fd.nelx])=tmp;
                        sf=zeros(size(x));
                        for i=1:fd.nelx
                            for j=1:fd.nely
                                sf(j,i)= sum( tmp2(dd+j+ny*(dd+i-1)+ii) );
                            end
                        end
                        sf=sf.*exp(beta*(1-x));
                        sf=sf(:);
                        
                    case 'gmean'
                        % preparations: redo forward computations to have
                        % all quantities available
                        P=fd.data.P; beta=fd.data.beta; eta=fd.data.eta; a=fd.data.a; b=fd.data.b; 
                        num=tanh(beta*eta)+tanh(beta*(1-eta));
                        xmean=zeros(fd.nely,fd.nelx); 
                        xmean(:) = (fd.data.H*x(:))./fd.data.Hs;
                        % aggregate using P-norm:
                        maxmean=(mean(xmean(:).^P))^(1/P);
                        % apply knock-down function:
                        S=(tanh(beta*eta)+tanh(beta*(maxmean-eta)))/num;
                        f=1-a*S-b*maxmean*S;
                        % start of sensitivity part ---
                        N=length(s(:));
                        sf=zeros(fd.nely,fd.nelx); sf1=sf; sf2=sf;
                        sf1=s*f;
                        dMdmu=(1/N)*(xmean/maxmean).^(P-1);
                        dSdM=(beta*(sech(beta*(maxmean-eta)))^2 )/num;
                        dydM=x*( (-a-b*maxmean)*dSdM -b*S );                        
                        tmp=repmat(s(:)'*dydM(:),fd.nely,fd.nelx);
                        tmp=tmp.*dMdmu;
                        sf2(:) = fd.data.H*(tmp(:)./fd.data.Hs);                        
                        sf=sf1+sf2;
                        % SIMPLER for testing:
                        if 0
                            xf=x*mean(xmean(:));
                            sf1=s*mean(xmean(:));
                            tmp=repmat(s(:)'*x(:),fd.nely,fd.nelx)/N;
                            sf2(:) = fd.data.H*(tmp(:)./fd.data.Hs);
                            sf=sf1+sf2;
                        end
                        disp(num2str([maxmean, max(xmean(:)), f]));
                        
                    case 'grad'
                        sf=zeros(size(s));
                        sfx=sf; sfy=sf;
                        sfx(:) = fd.data.Hx'*s(:);
                        sfy(:) = fd.data.Hy'*s(:);
                        
                        xf=zeros(fd.nely,fd.nelx); 
                        gx=xf; gy=xf;
                        gx(:) = fd.data.Hx*x(:);
                        gy(:) = fd.data.Hy*x(:);                        
                        % fd.data.out=varargin{5}; % 'mag' (magnitude), 'theta' (angle), 'gx', 'gy', 'dot'
                        switch fd.data.out
                            case 'mag'
                                xf=sqrt(gx.^2 + gy.^2);
                                sf(:)=(gx(:).*sfx(:) + gy(:).*sfy(:))./(xf(:) + eps); % adding some eps to prevent Inf
                            case 'theta'
                                xf=atan2(gy,gx);
                                % https://en.wikipedia.org/wiki/Atan2
                                sf(:)=(-gy(:).*sfx(:)+gx(:).*sfy(:))./(gx(:).^2 + gy(:).^2 + eps); % adding some eps to prevent Inf
                            case 'gx'
                                sf=sfx;
                            case 'gy'
                                sf=sfy;
                            case 'dot'
                                sf=fd.data.vec(1)*sfx+fd.data.vec(2)*sfy;
                        end           
                        sf=sf(:);
                        
                    case 'KSmask'
                        N=fd.data.N; P=fd.data.P;                        
                        Lx=fd.data.Lx; Ly=fd.data.Ly;
                        x0=fd.data.x0; y0=fd.data.y0; mi=fd.data.mi; mi2=fd.data.mi2; % mask indices (single-vector)
                        nely=fd.nely; nelx=fd.nelx;                        
                        % embed x in padding:
                        x2=fd.data.ghostval*ones(Ly,Lx);
                        x2(y0+[1:nely],x0+[1:nelx])=x;
                        xf=zeros(size(x));
                        tmp=zeros(size(x));
                        for i=1:nelx
                            for j=1:nely
                                %KS version:
                                %xf(j,i)= log( (1/N)*sum( exp(P*( x2( (x0+i-1)*Ly+(y0+j)+mi )) )) )/P;
                                tmp(j,i)=  sum(realmin + exp(P*( x2( (x0+i-1)*Ly+(y0+j)+mi )) ));
                            end
                        end
                        %tmp(:)=s(:)./(realmin+tmp(:)); % DIV0 prevention - now fixed above in the sum term                  
                        tmp(:)=s(:)./tmp(:);
                        tmp2=zeros(Ly,Lx);
                        tmp2(y0+[1:nely],x0+[1:nelx])=tmp;
                        sf=zeros(size(x));
                        if 1
                            for i=1:nelx
                                for j=1:nely
                                    %xf(j,i)= log( (1/N)*sum( exp(P*( x2( (x0+i-1)*Ly+(y0+j)+mi )) )) )/P;
                                    sf(j,i)=sum( tmp2( (x0+i-1)*Ly+(y0+j)+mi2 ) ); % use flipped mask
                                end
                            end
                        else
                            % Alternative: flip the domain, not the mask!                            
                            % elegant but I must define a 2-way-flippable
                            % embedding domain. (this is needed anyway)
                            tmp2=tmp2(end:-1:1,:); tmp2=tmp2(:,end:-1:1);
                            for i=1:nelx
                                for j=1:nely
                                    %xf(j,i)= log( (1/N)*sum( exp(P*( x2( (x0+i-1)*Ly+(y0+j)+mi )) )) )/P;
                                    sf(j,i)=sum( tmp2( (x0+i-1)*Ly+(y0+j)+mi ) );
                                end
                            end
                            % and flip result:
                            sf=sf(end:-1:1,:); sf=sf(:,end:-1:1);
                        end
                        sf=sf.*exp(P*x);
                        sf=sf(:);
                        
                    case 'maxfs1'
                        % maximum feature size
                        % stdev:
                        xmean=zeros(size(x)); 
                        msd=zeros(size(x)); 
                        stdev=zeros(size(x)); 
                        xmean(:) = (fd.data.H*x(:))./fd.data.Hs;                        
                        dev=x-xmean;
                        dev2=dev.^2;
                        msd(:) = (fd.data.H*dev2(:))./fd.data.Hs; % mean of squared deviations = var
                        stdev=sqrt(msd+eps); % xf = stdev
                        a=fd.data.a; b=fd.data.b;
                        %xf=x.*(a+(1-a)*tanh(b*stdev));
                        %
                        s=reshape(s,size(x));
                        sf=zeros(size(s)); tmp=sf; tmp2=sf;
                        dyds=zeros(size(s)); dydx=zeros(size(s));
                        dyds=b*(1-a)*x./(cosh(b*stdev)).^2;
                        dydx=a+(1-a)*tanh(b*stdev);
                        tmp=.5*(dyds.*s)./(stdev+eps); % tmp = s*dy/d(stdev)*d(stdev)/d(var)
                        tmp2(:) = fd.data.H*(tmp(:)./fd.data.Hs);
                        tmp=2*dev.*tmp2;
                        sf(:) = fd.data.H*(tmp(:)./fd.data.Hs);  % mean part of dev                        
                        sf=tmp-sf;
                        % still missing: the df/dy*dy/dx part
                        sf=sf + s.*dydx;
                        sf=sf(:);

                    case 'maxfs2'
                        % maximum feature size
                        % msd:
                        xmean=zeros(size(x)); 
                        msd=zeros(size(x));                         
                        xmean(:) = (fd.data.H*x(:))./fd.data.Hs;                        
                        dev=x-xmean;
                        dev2=dev.^2;
                        msd(:) = (fd.data.H*dev2(:))./fd.data.Hs; % mean of squared deviations = var                        
                        a=fd.data.a; b=fd.data.b;
                        %xf=x.*(a+(1-a)*tanh(b*msd));
                        %
                        s=reshape(s,size(x));
                        sf=zeros(size(s)); tmp=sf; tmp2=sf;
                        dyds=zeros(size(s)); dydx=zeros(size(s));
                        dyds=b*(1-a)*x./(cosh(b*msd)).^2;
                        dydx=a+(1-a)*tanh(b*msd);
                        tmp=dyds.*s; % tmp = s*dy/d(msd)
                        tmp2(:) = fd.data.H*(tmp(:)./fd.data.Hs);
                        tmp=2*dev.*tmp2;
                        sf(:) = fd.data.H*(tmp(:)./fd.data.Hs);  % mean part of dev                        
                        sf=tmp-sf;
                        % still missing: the df/dy*dy/dx part
                        sf=sf + s.*dydx;
                        sf=sf(:);

                    case 'maxfs2line'
                        % maximum feature size. Copy of maxfs2
                        % msd:
                        xmean=zeros(size(x)); 
                        msd=zeros(size(x));                         
                        xmean(:) = (fd.data.H*x(:))./fd.data.Hs;                        
                        dev=x-xmean;
                        dev2=dev.^2;
                        msd(:) = (fd.data.H*dev2(:))./fd.data.Hs; % mean of squared deviations = var                        
                        a=fd.data.a; b=fd.data.b;
                        %xf=x.*(a+(1-a)*tanh(b*msd));
                        %
                        s=reshape(s,size(x));
                        sf=zeros(size(s)); tmp=sf; tmp2=sf;
                        dyds=zeros(size(s)); dydx=zeros(size(s));
                        dyds=b*(1-a)*x./(cosh(b*msd)).^2;
                        dydx=a+(1-a)*tanh(b*msd);
                        tmp=dyds.*s; % tmp = s*dy/d(msd)
                        tmp2(:) = fd.data.H*(tmp(:)./fd.data.Hs);
                        tmp=2*dev.*tmp2;
                        sf(:) = fd.data.H*(tmp(:)./fd.data.Hs);  % mean part of dev                        
                        sf=tmp-sf;
                        % still missing: the df/dy*dy/dx part
                        sf=sf + s.*dydx;
                        sf=sf(:);

                    case 'mean'
                        % use weights to deal with ghost boundary:
                        dd=fd.data.dd; Ns=fd.data.Ns; ii=fd.data.ii; w=fd.data.w;
                        % embed in padding:
                        tmp=zeros(fd.nely+2*dd,fd.nelx+2*dd); tmp(dd+[1:fd.nely],dd+[1:fd.nelx])=s./Ns;
                        sf=zeros(size(x)); 
                        ny=fd.nely+2*dd;
                        for i=1:fd.nelx
                            for j=1:fd.nely
                                sf(j,i)= sum( tmp(dd+j+ny*(dd+i-1)+ii) );
                            end
                        end
                        sf=sf(:);
                        
                    case 'mean2'
                        sf=zeros(size(s));
                        sf(:) = fd.data.H*(s(:)./fd.data.Hs);                        
                        
                    case 'mill'
                        % milling from above:
                        % first, compute filtered field:
                        xf=zeros(size(x));
                        xf(1,:)=x(1,:);
                        P=fd.data.P; nely=fd.nely;
                        for i=2:nely
                            xf(i,:)=((x(i,:).^P+xf(i-1,:).^P)/fd.data.N).^(1/P);
                        end
                        % next, compute multipliers:
                        lam=zeros(size(x));
                        rhorho=xf; rhorho(1:(nely-1),:)=rhorho(1:(nely-1),:)./rhorho(2:nely,:);
                        rhorho(find(xf(:)==0))=1;
                        brho=x./xf; brho(find(xf(:)==0))=1;
                        fac=1/fd.data.N;
                        lam(end,:)=-s(end,:);
                        for i=(nely-1):-1:1                            
                            lam(i,:)=fac*lam(i+1,:).*rhorho(i,:).^(P-1) -s(i,:);
                        end
                        % sensitivities:
                        sf=zeros(nely,fd.nelx);
                        s=reshape(s,nely,fd.nelx);
                        sf(1:(nely-1),:)= fac*(s(1:(nely-1),:) - fac*lam(2:nely,:).*rhorho(1:(nely-1),:).^(P-1) ).*brho((1:nely-1),:).^(P-1);
                        sf(nely,:)= fac*s(nely,:).*brho(nely,:).^(P-1);                        

                    case 'mill2'                        
                        % milling from one direction:
                        switch fd.data.dir
                            case 1 % from above
                                sf=cumsum(s,1,'reverse');
                            case 2 % from left
                                sf=cumsum(s,2,'reverse');
                            case 3 % from below
                                sf=cumsum(s,1);
                            case 4 % from right
                                sf=cumsum(s,2);
                            otherwise
                                error('wrong direction argument');
                        end
                        
                    case 'mill3'
                        % milling from 2 directions: TEST
                        xf=cumsum(x);
                        xf2=cumsum(x,2);
                        % overflow protection
                        xf0=xf; xf20=xf2;
                        xf(xf>1)=log(exp(1)*xf(xf>1));
                        xf2(xf2>1)=log(exp(1)*xf2(xf2>1));
                        tmp=1./xf0; tmp2=1./xf20;
                        dxf=ones(size(xf)); dxf(xf0>1)=tmp(xf0>1);
                        dxf2=ones(size(xf2)); dxf2(xf20>1)=tmp2(xf20>1);

                        P=fd.data.P;
                        N=fd.data.N;
                        xf3=log( (exp(P*xf)+exp(P*xf2))/N )/P;
                        s=reshape(s,fd.nely,fd.nelx);
                        
                        tmp1=cumsum(exp(P*(xf-xf3)).*s/N.*dxf,1,'reverse');
                        tmp2=cumsum(exp(P*(xf2-xf3)).*s/N.*dxf2,2,'reverse');
                        %tmp1=tmp1.*dxf;
                        %tmp2=tmp2.*dxf2;
                        sf=tmp1 + tmp2;
                        
                    case 'mill4'
                        % milling from various directions: 
                        P=fd.data.P;
                        N=fd.data.N;
                        KS=zeros(size(x));
                        xfc=cell(N,1); % xf, corrected for overflow
                        dxf=cell(N,1);
                        for d=1:N
                            dr=fd.data.dir(d);                            
                            if dr<=4
                                xf=cumsum(x,fd.data.csi(dr),fd.data.css{fd.data.csd(dr)});
                            else
                                xf=cumsum45(x,dr-4,'forward');
                            end
%                             switch fd.data.dir(d)
%                                 case 1 % from above
%                                     xf=cumsum(x,1);
%                                 case 2 % from left
%                                     xf=cumsum(x,2);
%                                 case 3 % from below
%                                     xf=cumsum(x,1,'reverse');
%                                 case 4 % from right
%                                     xf=cumsum(x,2,'reverse');
%                                 otherwise
%                                     error('wrong direction argument');
%                             end
                            % overflow protection:
                            xfc{d}=xf; xfc{d}(xf>1)=log(exp(1)*xf(xf>1));
                            dxf{d}=ones(size(xf)); tmp=1./xf; dxf{d}(xf>1)=tmp(xf>1);                            
                            % accumulate
                            KS=KS+exp(P*xfc{d});
                        end
                        xf=log(KS/N)/P;
                        s=reshape(s,fd.nely,fd.nelx);
                        sf=zeros(size(s));
                        for d=1:N
                            dr=fd.data.dir(d);
                            if dr<=4
                                sf=sf+cumsum(exp(P*(xfc{d}-xf)).*s/N.*dxf{d},fd.data.csi(dr),fd.data.css{3-fd.data.csd(dr)});
                            else
                                sf=sf+cumsum45(exp(P*(xfc{d}-xf)).*s/N.*dxf{d},dr-4,'reverse');
                            end
%                             switch fd.data.dir(d)
%                                 case 1 % from above
%                                     sf=sf+cumsum(exp(P*(xfc{d}-xf)).*s/N,1,'reverse').*dxf{d};
%                                 case 2 % from left
%                                     sf=sf+cumsum(exp(P*(xfc{d}-xf)).*s/N,2,'reverse').*dxf{d};
%                                 case 3 % from below
%                                     sf=sf+cumsum(exp(P*(xfc{d}-xf)).*s/N,1).*dxf{d};
%                                 case 4 % from right
%                                     sf=sf+cumsum(exp(P*(xfc{d}-xf)).*s/N,2).*dxf{d};
%                             end
                        end
                        
                    case 'mill5'
                        % sums alll orientation contributions, instead of
                        % taking the (P-norm) maximum
                        sf=fd.data.M'*s(:);
                        
                    case 'mill5a'
                        s=s(:);
                        xf=zeros(fd.nely,fd.nelx);
                        x=x(:);
                        P=fd.data.P;
                        N=fd.nely*fd.nelx;
                        SF=zeros(N, fd.data.ndir);
                        XF=zeros(N, fd.data.ndir);
                        for i=1:fd.data.ndir
                            % apply separate transformations
                            XF(:,i)=fd.data.M{i}*x(N*(i-1) + [1:N]);
                        end
                        % apply P-norm on XF, combining rows
                        %xf(:)=sum(XF.^P,2).^(1/P);
                        sp=sum(XF.^P,2).^(1/P -1);
                        for i=1:fd.data.ndir
                            % apply separate transformations
                            %sfi=fd.data.M{i}'*s(N*(i-1) + [1:N])';
                            %SF(:,i)=sfi.*sp.*XF(:,i).^(P-1);
                            SF(:,i)=fd.data.M{i}'*(s.*sp.*XF(:,i).^(P-1));
                        end
                        sf=SF(:);
                        
                    case {'mill*', 'mill*2'}
%                         s=s(:);
%                         xf=zeros(fd.nely,fd.nelx);
%                         P=fd.data.P;
%                         N=fd.nely*fd.nelx;
%                         nd=fd.data.ndir;
%                         SF=zeros(N, nd);
%                         XF=zeros(N, nd);
%                         for i=1:nd
%                             % apply separate transformations
%                             XF(:,i)=fd.data.M{i}*x(:);
%                         end
%                         % apply P-norm on XF, combining rows                        
%                         %xf(:)=sum(XF.^P,2).^(1/P);
%                         sp=sum((1/nd)*XF.^P,2).^(1/P -1);
%                         for i=1:nd
%                             % apply separate transformations
%                             %sfi=fd.data.M{i}'*s(N*(i-1) + [1:N])';
%                             %SF(:,i)=sfi.*sp.*XF(:,i).^(P-1);
%                             SF(:,i)=fd.data.M{i}'*(s.*sp.*XF(:,i).^(P-1))/nd;
%                         end
%                         sf=sum(SF,2);                        
                        % milling from various directions: 
                        s=s(:);
                        xf=zeros(fd.nely,fd.nelx);
                        P=fd.data.P;
                        N=fd.data.N;
                        KS=zeros(size(x(:)));
                        xfc=cell(N,1); % xf, corrected for overflow
                        dxf=cell(N,1);
              if isfield(fd.data,'Power')
                  x0=x;
                  x=x.^fd.data.Power;                        
              end
                        for i=1:N
                            xf=fd.data.M{i}*x(:);
                            % overflow protection:
                            xfc{i}=xf; xfc{i}(xf>1)=log(exp(1)*xf(xf>1));
                            dxf{i}=ones(size(xf)); tmp=1./xf; dxf{i}(xf>1)=tmp(xf>1);                            
                            % accumulate
                            KS=KS+exp(P*xfc{i}(:));
                        end
                        xf(:)=log(KS/N)/P;
                        s=reshape(s,fd.nely,fd.nelx);
                        sf=zeros(size(s));
                        for i=1:N
                            sf(:)=sf(:)+fd.data.M{i}'*(exp(P*(xfc{i}(:)-xf(:))).*s(:)/N.*dxf{i}(:));
                        end         
              if isfield(fd.data,'Power')
                   sf=fd.data.Power*sf.*x0.^(fd.data.Power-1);
              end
                        
              
                    case {'mill*3','mill*3P'}
                        % milling from various directions: 
                        s=s(:);
                        xf=zeros(fd.nely,fd.nelx);
                        P=fd.data.P;
                        N=fd.data.N;
                        KS=zeros(size(x(:)));
                        xfc=cell(N,1); % xf, corrected for overflow
                        dxf=cell(N,1);
                        KSxin=cell(N,1);
              if isfield(fd.data,'Power')
                  x0=x;
                  x=x.^fd.data.Power;                        
              end
                        for i=1:N
                            %xf=fd.data.M{i}*x(:);                            
                            xm=zeros(fd.data.nemy(i),fd.data.nemx(i));                            
                            xm(:)=fd.data.F{i}*x(:); 
                            xm(:)=fd.data.C{i}*xm(:);                             
                            KSxin{i}=xm;
                            [ignore, xm]=anyfilter(fd.data.KSmask{i},KSxin{i});
                            xf(:)=fd.data.B{i}*xm(:);

                            % overflow protection:
                            xfc{i}=xf; xfc{i}(xf>1)=log(exp(1)*xf(xf>1));
                            dxf{i}=ones(size(xf)); tmp=1./xf; dxf{i}(xf>1)=tmp(xf>1);                            
                            % accumulate
                            KS=KS+exp(P*xfc{i}(:));
                        end
                        xf(:)=log(KS/N)/P;
                        s=reshape(s,fd.nely,fd.nelx);
                        sf=zeros(size(s)); sfi=sf;
                        for i=1:N
                            %sf(:)=sf(:)+fd.data.M{i}'*(exp(P*(xfc{i}(:)-xf(:))).*s(:)/N.*dxf{i}(:));                            
                            sm=zeros(fd.data.nemy(i),fd.data.nemx(i));
                            sm(:) = fd.data.B{i}'*(exp(P*(xfc{i}(:)-xf(:))).*s(:)/N.*dxf{i}(:));
                            [ignore, sKS]=anyfilter(fd.data.KSmask{i},KSxin{i},sm);
                            sm(:) = fd.data.C{i}'*sKS(:);
                            sfi(:) = fd.data.F{i}'*sm(:);
                            
                            sf(:)=sf(:)+sfi(:);                            
                        end         
              if isfield(fd.data,'Power')
                   sf=fd.data.Power*sf.*x0.^(fd.data.Power-1);
              end
             
                    case {'mill*4','mill*4P'}
                        % milling from various directions: 
                        s=s(:);
                        xf=zeros(fd.nely,fd.nelx);
                        epa=fd.data.epsilon_angles;
                        ma=fd.data.m_angles;
                        Pa=-1;
                        N=fd.data.N;
                        PNormSum=zeros(size(x));
                        %xfc=cell(N,1); % xf, corrected for overflow
                        %dxf=cell(N,1);
                        KSxin=cell(N,1);
                        Pxin=cell(N,1);
              if isfield(fd.data,'Power')
                  x0=x;
                  x=x.^fd.data.Power;                        
              end
                        for i=1:N
                            %xf=fd.data.M{i}*x(:);                            
                            xm=zeros(fd.data.nemy(i),fd.data.nemx(i));                            
                            xm(:)=fd.data.F{i}*x(:); 
                            xm(:)=fd.data.C{i}*xm(:);                            
                            KSxin{i}=xm;
                            [ignore, xm]=anyfilter(fd.data.KSmask{i},KSxin{i});
                            xf(:)=fd.data.B{i}*xm(:);
                            Pxin{i}=xf;

                            % overflow protection:
                            %xfc{i}=xf; xfc{i}(xf>1)=log(exp(1)*xf(xf>1));
                            %dxf{i}=ones(size(xf)); tmp=1./xf; dxf{i}(xf>1)=tmp(xf>1);                            
                            % accumulate
                            %KS=KS+exp(P*xfc{i}(:));
                            PNormSum=PNormSum+(xf+epa).^Pa;
                        end
                        %xf(:)=log(KS/N)/P;
                        xf(:)= (PNormSum/ma).^(1/Pa) - epa; 
                        s=reshape(s,fd.nely,fd.nelx);
                        sf=zeros(size(s)); sfi=sf;
                        for i=1:N
                            %sf(:)=sf(:)+fd.data.M{i}'*(exp(P*(xfc{i}(:)-xf(:))).*s(:)/N.*dxf{i}(:));                            
                            sm=zeros(fd.data.nemy(i),fd.data.nemx(i));
                            %sm(:) = fd.data.B{i}'*(exp(P*(xfc{i}(:)-xf(:))).*s(:)/N.*dxf{i}(:));
                            sm(:) = fd.data.B{i}'*(  ( s(:).*(1/ma*PNormSum(:)).^(1/Pa -1) ).*( (Pxin{i}(:)+epa).^(Pa-1) )/ma );
                            [ignore, sKS]=anyfilter(fd.data.KSmask{i},KSxin{i},sm);
                            sm(:) = fd.data.C{i}'*sKS(:);
                            sfi(:) = fd.data.F{i}'*sm(:);  
                            
                            sf(:)=sf(:)+sfi(:);
                        end         
              if isfield(fd.data,'Power')
                   sf=fd.data.Power*sf.*x0.^(fd.data.Power-1);
              end             
              
                    case 'millX'
                        % milling from various directions, with toolshape and -length: 
                        s=s(:);
                        xf=zeros(fd.nely,fd.nelx);
                        P=fd.data.P;
                        N=fd.data.N;
                        TL=fd.data.TL; PTL=fd.data.PTL; 
                        KS=zeros(size(x(:)));
                        xfc=cell(N,1); % xf, corrected for overflow
                        dxf=cell(N,1);
                        KSxin=cell(N,1);              
                        x0=x;
                        x=x.^fd.data.Power;
                        %
                        for i=1:N
                            %xf=fd.data.M{i}*x(:);                            
                            xm=zeros(fd.data.nemy(i),fd.data.nemx(i));                            
                            xm(:)=fd.data.F{i}*x(:); 
                            if ~isempty(TL) % account for tool length
                                if (TL<fd.data.nemy(i))
                                    sumePx=sum(exp(PTL*xm),2);
                                    rowmax=(1/PTL)*log((1/fd.data.nemx(i))*sumePx);
                                    Crowmax=cumsum(rowmax);
                                    Crowmax((TL+1):end)=Crowmax(1:(end-TL)); Crowmax(1:TL)=0;
                                    fCrowmax=Crowmax.*exp(1-Crowmax);
                                    for j=1:fd.data.nemy(i)
                                        xm(j,:)=xm(j,:)+fCrowmax(j);
                                    end
                                end
                            end
                            xm(:)=fd.data.C{i}*xm(:);  
                            for k=1:length(fd.data.topi{i})
                                bordervalue=xm(fd.data.topi{i}(k),fd.data.topj{i}(k));                                
                                xm(1:fd.data.topi{i}(k), fd.data.topj{i}(k) )=bordervalue;
                            end
                            KSxin{i}=xm;
                            [~, xm]=anyfilter(fd.data.KSmask{i},KSxin{i});
                            xf(:)=fd.data.B{i}*xm(:);

                            % overflow protection:
                            xfc{i}=xf; xfc{i}(xf>1)=log(exp(1)*xf(xf>1));
                            dxf{i}=ones(size(xf)); tmp=1./xf; dxf{i}(xf>1)=tmp(xf>1);                            
                            % accumulate
                            KS=KS+exp(P*xfc{i}(:));
                        end
                        xf(:)=log(KS/N)/P;
                        s=reshape(s,fd.nely,fd.nelx);
                        sf=zeros(size(s)); sfi=sf;
                        for i=1:N
                            %sf(:)=sf(:)+fd.data.M{i}'*(exp(P*(xfc{i}(:)-xf(:))).*s(:)/N.*dxf{i}(:));                            
                            sm=zeros(fd.data.nemy(i),fd.data.nemx(i));
                            sm(:) = fd.data.B{i}'*(exp(P*(xfc{i}(:)-xf(:))).*s(:)/N.*dxf{i}(:));
                            [~, sKS]=anyfilter(fd.data.KSmask{i},KSxin{i},sm);
                            sKS=reshape(sKS,fd.data.nemy(i),fd.data.nemx(i));
                            for k=1:length(fd.data.topi{i})
                                sKS(fd.data.topi{i}(k),fd.data.topj{i}(k))=sum( sKS(1:fd.data.topi{i}(k),fd.data.topj{i}(k)) );
                            end
                            sm(:) = fd.data.C{i}'*sKS(:);
                            if ~isempty(TL) % account for tool length
                                if (TL<fd.data.nemy(i))
                                    % recreate input:
                                    xm=zeros(fd.data.nemy(i),fd.data.nemx(i));
                                    xm(:)=fd.data.F{i}*x(:);
                                    %
                                    sumePx=sum(exp(PTL*xm),2);
                                    rowmax=(1/PTL)*log((1/fd.data.nemx(i))*sumePx);
                                    Crowmax=cumsum(rowmax);
                                    Crowmax((TL+1):end)=Crowmax(1:(end-TL)); Crowmax(1:TL)=0;
                                    % start sensitivity part:
                                    dfdC=(1-Crowmax).*exp(1-Crowmax);
                                    sumrowsens=sum(sm,2);
                                    df_srs=dfdC.*sumrowsens;
                                    df_srs(1:(end-TL))=df_srs((TL+1):end); df_srs((end-TL+1):end)=0;
                                    Cdf_srs=cumsum(df_srs,'reverse');
                                    %for j=1:(fd.data.nemy(i)-TL)
                                    for j=1:fd.data.nemy(i)
                                        %sm(j,:)=sm(j,:)+sumrowsens(j+TL)*(exp(PTL*xm(j,:)))/sumePx(j);
                                        %sm(j,:)=sm(j,:)+Cdf_srs(j+TL)*(exp(PTL*xm(j,:)))/sumePx(j);
                                        sm(j,:)=sm(j,:)+Cdf_srs(j)*(exp(PTL*xm(j,:)))/sumePx(j);
                                    end
                                end
                            end                            
                            sfi(:) = fd.data.F{i}'*sm(:);
                            
                            sf(:)=sf(:)+sfi(:);                            
                        end              
                        sf=fd.data.Power*sf.*x0.^(fd.data.Power-1);
              
                    case 'none'
                        sf=s;
                        
                    case 'oneminusx'
                        sf=-s;
                        
                    case 'Pmask'
                        N=fd.data.N; P=fd.data.P; ep=fd.data.eps; m=fd.data.m;
                        Lx=fd.data.Lx; Ly=fd.data.Ly;
                        x0=fd.data.x0; y0=fd.data.y0; mi=fd.data.mi; % mask indices (single-vector)
                        nely=fd.nely; nelx=fd.nelx;                        
                        % embed x in padding:
                        x2=fd.data.ghostval*ones(Ly,Lx);
                        x2(y0+[1:nely],x0+[1:nelx])=x;
                        xf=zeros(size(x));
                        tmp=zeros(size(x));
                        for i=1:nelx
                            for j=1:nely
                                %P-norm version:
                                %xf(j,i)= ( (1/m)*sum( (x2( (x0+i-1)*Ly+(y0+j)+mi) + ep).^P ) ).^(1.0/P) ;
                                tmp(j,i)= ( (1/m)*sum( (x2( (x0+i-1)*Ly+(y0+j)+mi) + ep).^P ) ).^(1.0/P-1);                                
                            end
                        end
                        tmp(:)=s(:).*tmp(:);
                        tmp2=zeros(Ly,Lx);
                        tmp2(y0+[1:nely],x0+[1:nelx])=tmp;
                        sf=zeros(size(x));
                        for i=1:nelx
                            for j=1:nely                                
                                sf(j,i)=sum( tmp2( (x0+i-1)*Ly+(y0+j) - mi ) ); % use flipped mask
                            end
                        end
                        sf=(1/m)*sf.*(x + ep).^(P-1);
                        sf=sf(:);
                        
                    case 'range'
                        % maximum - minimum within filter radius.
                        sfmax=s; sfmin=s;
                        [fd.data.dilate,sfmax(:)]=anyfilter(fd.data.dilate,x,s);
                        [fd.data.erode,sfmin(:)]=anyfilter(fd.data.erode,x,s);
                        sf=sfmax-sfmin;                
                        
                    case 'stdev'
                        xmean=zeros(size(x)); 
                        xf=xmean;
                        xmean(:) = (fd.data.H*x(:))./fd.data.Hs;                        
                        dev=x-xmean;
                        dev2=dev.^2;
                        xf(:) = (fd.data.H*dev2(:))./fd.data.Hs; % mean of squared deviations = var
                        xf=sqrt(xf); % xf = stdev
                        
                        sf=zeros(size(s)); tmp=sf; tmp2=sf;
                        tmp=.5*s./xf; % tmp = s*d(stdev)/d(var)
                        tmp2(:) = fd.data.H*(tmp(:)./fd.data.Hs);
                        tmp=2*dev.*tmp2;
                        sf(:) = fd.data.H*(tmp(:)./fd.data.Hs);  % mean part of dev
                        sf=tmp-sf;
                        
                        
                    case 'symm' % symmetry in x- or y-direction, halfway the domain
                        switch fd.data.xy
                            case 'x'
                                sf=0.5*(s+fliplr(s));
                            case 'y'
                                sf=0.5*(s+flipud(s));
                            otherwise
                                error('symm: wrong argument');
                        end
                        
                end
                if Nfd>1
                    % filter with multiple steps: process previous step(s)
                    [fd_part, sf] = anyfilter(fd_in(1:(Nfd-1)), fd_in(Nfd-1).x_in, sf);
                    fd=[fd_part, fd_in(Nfd)];
                end
                varargout={fd,sf};
        end
    end
end

end


function fd=prep(fd_in)
% for fd_in having nelx, nely and data.rmin, this function generates:
% data.ii: index-stencil, relative indices of elements within the filter radius
% data.Ns: number of elements contained in the filter radius at each position
% data.dd: maximum padding
% usage: see implementation of dilate and erode filters.
%
% PREPARATIONS:
fd=fd_in;
% define index shifts in x and y direction, that cover the circle
% defined by rmin
fd.data.dd=ceil(fd.data.rmin);
dd=fd.data.dd; rmin=fd.data.rmin;
sx=zeros(dd*dd,1); sy=sx; k=1;
for i=-dd:dd
    for j=-dd:dd
        if sqrt(i^2+j^2)<=rmin
            sx(k)=i; sy(k)=j;
            k=k+1;
        end
    end
end
ns=k-1;
sx=sx(1:ns); sy=sy(1:ns);
ii=sx*(fd.nely+2*dd)+sy;
% sweep over 0-padded 1-field, to compute Ns:
Ns=zeros(fd.nely,fd.nelx);
xdummy=zeros(fd.nely+2*dd,fd.nelx+2*dd); xdummy(dd+[1:fd.nely],dd+[1:fd.nelx])=1;
for i=dd+[1:fd.nelx]
    for j=dd+[1:fd.nely]
        Ns(j-dd,i-dd)=sum(xdummy(j+(fd.nely+2*dd)*(i-1)+ii));
    end
end
fd.data.ii=ii;
fd.data.Ns=Ns;
end



function fd=prep2(fd_in)
% for fd_in having nelx, nely and data.rmin, this function generates:
% data.ii: index-stencil, relative indices of elements within the filter radius
% data.w:  weight field for entries (1) and ghost boundary (0)
% data.dd: maximum padding
% usage: see implementation of mean and stdev filters.
%
% PREPARATIONS:
fd=fd_in;
% define index shifts in x and y direction, that cover the circle
% defined by rmin
fd.data.dd=ceil(fd.data.rmin);
dd=fd.data.dd; rmin=fd.data.rmin;
sx=zeros(dd*dd,1); sy=sx; k=1;
for i=-dd:dd
    for j=-dd:dd
        if sqrt(i^2+j^2)<=rmin
            sx(k)=i; sy(k)=j;
            k=k+1;
        end
    end
end
ns=k-1;
sx=sx(1:ns); sy=sy(1:ns);
ii=sx*(fd.nely+2*dd)+sy;
fd.data.ns=ns;
% sweep over 0-padded 1-field, to compute Ns:
Ns=zeros(fd.nely,fd.nelx);
xdummy=zeros(fd.nely+2*dd,fd.nelx+2*dd); xdummy(dd+[1:fd.nely],dd+[1:fd.nelx])=1;
for i=dd+[1:fd.nelx]
    for j=dd+[1:fd.nely]
        Ns(j-dd,i-dd)=sum(xdummy(j+(fd.nely+2*dd)*(i-1)+ii));
    end
end
fd.data.ii=ii;
fd.data.Ns=Ns;
fd.data.w=xdummy;
end


function cs=cumsum45(x,dir,fw_rev)
% cumsum in 45-direction
% dir: 1 2 3 4
% fw_rev: 'forward' / 'reverse'
[nely,nelx]=size(x);
nelx2=nelx+nely-1;
x2=zeros(nely,nelx2);
cs=zeros(size(x));
fr={'forward','reverse'};
switch fw_rev
    case 'forward'
        fri=1;
    case 'reverse'
        fri=2;
    otherwise
        error('invalid input')
end
switch dir
    case 1
        % shift right, cumsum 1 forward
        frj=fri;
        for i=1:nely
            x2(i,[1:nelx]+i-1)=x(i,:);
        end
        x2=cumsum(x2,1,fr{frj});
        for i=1:nely
            cs(i,:)=x2(i,[1:nelx]+i-1);
        end
    case 2
        % shift left, cumsum 1 forward
        frj=fri;
        for i=1:nely
            x2(i,[1:nelx]+nely-i)=x(i,:);
        end
        x2=cumsum(x2,1,fr{frj});
        for i=1:nely
            cs(i,:)=x2(i,[1:nelx]+nely-i);
        end
    case 3 
        % shift right, cumsum 1 reverse
        frj=3-fri;
        for i=1:nely
            x2(i,[1:nelx]+i-1)=x(i,:);
        end
        x2=cumsum(x2,1,fr{frj});
        for i=1:nely
            cs(i,:)=x2(i,[1:nelx]+i-1);
        end
    case 4
        % shift left, cumsum 1 reverse
        frj=3-fri;
        for i=1:nely
            x2(i,[1:nelx]+nely-i)=x(i,:);
        end
        x2=cumsum(x2,1,fr{frj});
        for i=1:nely
            cs(i,:)=x2(i,[1:nelx]+nely-i);
        end
    otherwise
        error('invalid direction');
end
end



function T = ToolShapeMatrix(nely,nelx,Len,Diam)
    % form tool shape matrix
    N=nelx*nely;
    m=(nely-max(Len)+1)*(max(Diam)-1) + N; % number of sparse entries (upper bound);
    is=ones(m,1); js=is; vs=zeros(m,1);    
    k=1;
    for ii=1:length(Len)
        L=Len(ii); D=round(Diam(ii)/2); 
        if (ii==1), Din=0; else Din=round(Diam(ii-1)/2); end; % previous diameter
        for i=L:nely % this tool segment (y-index, in milling direction)
            for j=1:nelx % x-index
                for jj=-D:D
                    if all([abs(jj)>Din,j+jj>=1,j+jj<=nelx]) % skip jj=0 case, stay within domain
                        is(k)=nely*(j-1)+i;
                        js(k)=nely*(j-1+jj)+i-L+1;
                        vs(k)=1;
                        %figure(33); rho=zeros(nely,nelx); rho(i,j)=1; rho(is(k))=0.5; rho(js(k))=0.25; imagesc(rho); axis image; drawnow; pause(0.01);
                        k=k+1;
                    end
                end
            end
        end
    end
    T=sparse(is,js,vs,N,N);
    T=T+speye(N);    
    % testing code:
    if 1
        %T=T';
        rho=zeros(nely,nelx);
        rho(floor(nely/2),floor(nelx/2))=1;
        rho2=rho; rho2(:)=T*rho(:);
        figure; subplot(1,2,1); imagesc(rho); axis image; colorbar;
        subplot(1,2,2); imagesc(rho2); axis image; colorbar;
    end
end



function T2 = ToolShapeMatrix2(nely,nelx,Len,Diam)
    % form tool shape matrix
    N=nelx*nely;
    m=(nely-max(Len)+1)*(max(Diam)-1) + N; % number of sparse entries (upper bound);
    is=ones(m,1); js=is; vs=zeros(m,1);    
    k=1;
    for ii=1:length(Len)
        L=Len(ii); D=round(Diam(ii)/2); 
        if (ii==1), Din=0; else Din=round(Diam(ii-1)/2); end; % previous diameter
        for i=1:(nely+1-L) % this tool segment (y-index, in reversed milling direction)
            for j=1:nelx % x-index
                for jj=-D:D
                    if all([abs(jj)>Din,j+jj>=1,j+jj<=nelx]) % skip jj=0 case, stay within domain
                        is(k)=nely*(j-1)+i;
                        js(k)=nely*(j-1+jj)+i+L-1;
                        vs(k)=1;
                        %figure(33); rho=zeros(nely,nelx); rho(i,j)=1; rho(is(k))=0.5; rho(js(k))=0.25; imagesc(rho); axis image; drawnow; pause(0.01);
                        k=k+1;
                    end
                end
            end
        end
    end
    T2=sparse(is,js,vs,N,N);
    T2=T2+speye(N);    
    % testing code:
    if 1
        %T=T';
        rho=zeros(nely,nelx);
        rho(floor(nely/2),floor(nelx/2))=1;
        rho2=rho; rho2(:)=T2*rho(:);
        figure; subplot(1,2,1); imagesc(rho); axis image; colorbar;
        subplot(1,2,2); imagesc(rho2); axis image; colorbar;
    end
end


function T = ToolShapeMatrix3(nely,nelx,Len,Diam)
    % form tool shape matrix
    N=nelx*nely;
    m=(nely-max(Len)+1)*(max(Diam)-1) + N; % number of sparse entries (upper bound);
    is=ones(m,1); js=is; vs=zeros(m,1);    
    k=1;
    for ii=1:length(Len)
        L=Len(ii); D=round(Diam(ii)/2); 
        if (ii==1), Din=0; else Din=round(Diam(ii-1)/2); end; % previous diameter
        for i=1:(nely-L) % this tool segment (y-index, in milling direction)
            for j=1:nelx % x-index
                for jj=-D:D
                    if all([abs(jj)>Din,j+jj>=1,j+jj<=nelx]) % skip jj=0 case, stay within domain
                        % just sum stuff:
                        is(k)=i+(j-1)*nely;
                        js(k)=i+(j-1+jj)*nely;
                        % and shift output downward by L:
                        is(k)=is(k)+L;
                        
                        %is(k)=nely*(j-1)+i;
                        %js(k)=nely*(j-1+jj)+i-L+1;
                        vs(k)=1;
                        %figure(33); rho=zeros(nely,nelx); rho(i,j)=1; rho(is(k))=0.5; rho(js(k))=0.25; imagesc(rho); axis image; drawnow; pause(0.01);
                        k=k+1;
                    end
                end
            end
        end
    end
    T=sparse(is,js,vs,N,N);
    T=T+speye(N);    
    % testing code:
    if 1
        %T=T';
        rho=zeros(nely,nelx);
        rho(floor(nely/2),floor(nelx/2))=1;
        rho2=rho; rho2(:)=T*rho(:);
        figure; subplot(1,2,1); imagesc(rho); axis image; colorbar;
        subplot(1,2,2); imagesc(rho2); axis image; colorbar;        
    end
end


function TL = ToolLengthMatrix(nely,nelx,Len)
    nrow=nely;
    entries=nrow*nrow*(nely-Len);
    is=ones(entries,1); js=ones(entries,1); ss=zeros(entries,1);
    k=0;
    for i=(Len+1):nely
        for j=1:nelx
            for jj=1:nelx
                k=k+1;
                is(k)=i+(jj-1)*nely;
                js(k)=i-Len+(j-1)*nely;
                ss(k)=1;
            end
        end
    end
    N=nely*nelx;
    TL=sparse(is,js,ss,N,N);
end


function C = CumSumMatrix(nely,nelx)
    % form cumsum milling matrix (from above)
    n=nely*(nely+1)/2;
    i1=zeros(n,1); j1=i1;
    k=1;
    for i=1:nely
        i1(k:(k+i-1))=i;
        j1(k:(k+i-1))=1:i;
        k=k+i;
    end
    i2=zeros(nelx*n,1); j2=i2;
    for i=0:(nelx-1)
        i2(i*n+[1:n])=i1+i*nely;
        j2(i*n+[1:n])=j1+i*nely;
    end
    C=sparse(i2,j2,ones(size(i2)));
end


function [F,B,nemy,nemx]=FieldRotationMappings(nely,nelx,degrees)
    % constructs mapping matrices for rotating a field [nely,nelx] over 'degrees'.
    % F = forward mapping, B = return mapping. nemy, nemx are dimensions of the
    % rotated field.
    %
    % improved version, more memory-efficient & easier to turn into 3D
    % see newFieldRotationMappings1.m for comparison with original (below)
    if length(degrees)~=1, error('expecting single orientation'); end
    em_hx=1; em_hy=1; ha=1; hb=1; % no need to make this user-controlled here
    DEBUG_PLOTS=0; % careful, don't use too big meshes with this one
    % corner points of reference domain:
    cp=[0 0; nelx 0; nelx nely; 0 nely]; % corner vertices
    % rotated FE domain: Rcp
    [Rcx,Rcy]=rotfield(cp(:,1),cp(:,2),degrees,0,0); Rcp=[Rcx Rcy];
    x1=min(Rcx); x2=max(Rcx);
    y1=min(Rcy); y2=max(Rcy);
    % embedding domain: vbox
    vbox=[x1 y1; x2 y1; x2 y2; x1 y2];
    % embedding domain rotated back:
    [Rvbx,Rvby]=rotfield(vbox(:,1),vbox(:,2),-degrees,0,0); Rvbox=[Rvbx Rvby];
    if DEBUG_PLOTS
        figure; cla; subplot(2,2,1);
        xdesign = ones(nely,nelx);
        for i=1:nely, for j=1:nelx, xdesign(i,j)=.5+0.2*j/nelx+0.5*sin((i+j)/10); end; end
        f1p1=plotpatch(1-xdesign,em_hx,em_hy); set(f1p1,'edgecolor','k');
        % rotated FE domain: Rcp
        patch('faces',1:4,'vertices',Rcp,'facecolor','none','edgecolor','r');
        % embedding AM domain
        patch('faces',1:4,'vertices',vbox,'facecolor','none','edgecolor','b');
        % embedding domain rotated back:
        patch('faces',1:4,'vertices',Rvbox,'facecolor','none','edgecolor','m','linewidth',2);
        patch('faces',1:2,'vertices',Rvbox(1:2,:),'facecolor','none','edgecolor','g','linewidth',2);
        axis tight; axis equal;
    end
    % determine dimensions of embedding domain:
    RL=x2-x1; RH=y2-y1;
    em_nx=ceil(0.9999*RL/em_hx); em_ny=ceil(0.9999*RH/em_hy); em_xskip = (RL-em_nx*em_hx)/2; em_yskip=(RH-em_ny*em_hy)/2;
    [emx, emy]=meshgrid(0.5*em_hx+[em_xskip:em_hx:0.9999*RL],0.5*em_hy+[em_yskip:em_hy:0.9999*RH]); % element centers
    %patch('faces',[1:((em_nx+1)*(em_ny+1))]','vertices',[emx(:) emy(:)],'marker','.','edgecolor','k');
    [Remx, Remy]=rotfield(emx+x1,emy+y1,-degrees,0,0);
    if DEBUG_PLOTS
        % plot centroid nodes
        patch('faces',[1:em_nx*em_ny]','vertices',[Remx(:) Remy(:)],'marker','.','edgecolor','k');
    end
    if DEBUG_PLOTS
        % determine rotated cells:
        tmpx=repmat(Remx(:),1,4);ntmp=size(tmpx,1); tmpx=tmpx';
        tmpy=repmat(Remy(:),1,4);tmpy=tmpy';
        tmpv=[tmpx(:) tmpy(:)];
        vrts=[-1 -1; 1 -1; 1 1; -1 1]/2;
        vrts(:,1)=vrts(:,1)*em_hx;
        vrts(:,2)=vrts(:,2)*em_hy;
        [vrx,vry]=rotfield(vrts(:,1),vrts(:,2),-degrees,0,0);
        vrts=[vrx, vry];
        tmpf=4*(1:ntmp)'-3; tmpf=[tmpf tmpf+1 tmpf+2 tmpf+3];
        % flip y direction to match data ordering convention:
        k=reshape(1:size(tmpf,1),size(Remx,1),size(Remx,2));
        k=k(end:-1:1,:); tmpf=tmpf(k(:),:);
        %
        Remf=tmpf;
        Remv=tmpv+repmat(vrts,ntmp,1);
        % also plot 'cells':
        if 1
            patch('faces',Remf,'vertices',Remv,'facecolor','none','edgecolor','r','linestyle',':');
        end
    end
    % determine sampling points in [-.5 .5]^2 square:
    %https://pomax.github.io/bezierinfo/legendre-gauss.html
    GaussPts = 1/sqrt(3)*[1 1; -1 1; -1 -1; 1 -1]/2;
    %[xg,yg]=meshgrid(-0.495:0.01:0.495,-0.495:0.01:0.495); GaussPts0=[xg(:) yg(:)]; GaussPts=GaussPts0; % ULTRA-PRECISE!!
    % 10x10 scheme: very precise
    %[xg,yg]=meshgrid(-0.45:0.1:0.45,-0.45:0.1:0.45); GaussPts0=[xg(:) yg(:)]; GaussPts=GaussPts0;
    % 5x5 scheme: just fine
    [xg,yg]=meshgrid(-0.4:0.2:0.4,-0.4:0.2:0.4); GaussPts0=[xg(:) yg(:)]; GaussPts=GaussPts0;
    GaussPts(:,1)= GaussPts0(:,1)*em_hx; GaussPts(:,2)= GaussPts0(:,2)*em_hy;
    nGaussPts = length(GaussPts);
    [RGx, RGy]=rotfield(GaussPts(:,1),GaussPts(:,2),-degrees,0,0);
    if DEBUG_PLOTS
        % plot a few for illustration
        if 0
            for i=1:em_ny
                for j=1:em_nx
                    patch('faces',[1:nGaussPts],'vertices',[RGx RGy]+repmat([Remx(i,j) Remy(i,j)],nGaussPts,1),'edgecolor','r','marker','.','linestyle','none','facecolor','none');
                end
            end
            %             for k=-1:3
            %                 i=floor(em_nx/2)*em_ny-floor(em_ny/2)+k*em_ny;
            %                 for j=-2:2
            %                     patch('faces',[1 3; 2 4],'vertices',[RGx RGy]+repmat([Remx(i+j) Remy(i+j)],4,1),'edgecolor','r','marker','.');
            %                 end
            %             end
        end
    end
    % FORWARD MAPPING
    % x_extended(:) = F*x_FE(:)
    % determine mapping of points of the embedding mesh to FE elements (densities):
    [nemy,nemx]=size(Remx);
    elhx=1; elhy=1; % element dimensions
    % first, count BED elements that have an interaction with the x-domain:
    Remx=Remx(:); Remy=Remy(:);
    xin_sure = (ceil(Remx/elhx)>=0) & (ceil(Remx/elhx)<=(nelx+1));
    yin_sure = (ceil(Remy/elhy)>=0) & (ceil(Remy/elhy)<=(nely+1));
    in_sure = xin_sure & yin_sure;
    n_in = sum(in_sure); n_BED = length(in_sure);
    % alloc sufficient space for ijs data:
    patchsize=3*3;
    avghits=4; % expected nr of elements overlapped by rotated element
    sz=n_in*avghits;
    ivec=ones(sz,1); jvec=ones(sz,1); svec=zeros(sz,1);
    % loop over BED elements to add their contributions to ijs data:
    % init 3x3 patch offsets
    patch_offset_y=[-1 0 1 -1 0 1 -1 0 1]';
    patch_offset_x=[-1 -1 -1 0 0 0 1 1 1]';
    k = 1;
    for i=1:n_BED
        xi=ceil(Remx(i)); yi=ceil(Remy(i)); % index of BED element center in x-domain
        xp=patch_offset_x+xi; yp=patch_offset_y+yi;
        patch_i=yp+nely*(xp-1);
        patch_in_domain=all([xp>0, xp<=nelx, yp>0, yp<=nely],2);
        %if in_sure(i) % skip all 'outside' cases
        npin=sum(patch_in_domain);
        if npin>0 % skip all 'outside' cases
            Gx=RGx+Remx(i);
            Gy=RGy+Remy(i);
            xk=ceil(Gx/elhx); yk=ceil(Gy/elhy); % indices of Gausspoints in x-domain
            in_domain = all([xk>0, xk<=nelx, yk>0, yk<=nely],2);
            relxk=xk-xi+1+1; relyk=yk-yi+1+1; % +1 for 1-based indexing, +1 for offset from corner of 3x3 patch
            hitcount=accumarray(relyk+(relxk-1)*3, in_domain,[patchsize 1], @sum, 0); % Gauss point hits per element of 3x3 patch
            hitsum=sum(hitcount); % number of Gauss points in x-domain
            if hitsum
                % add entries to ijs vectors:
                hit_j=find(hitcount>0)';
                for j=hit_j
                    ivec(k)=i;
                    jvec(k)=patch_i(j);
                    svec(k)=hitcount(j)/hitsum;
                    k=k+1;
                end
            end
        end
    end
    disp(sprintf('%.1f%% of preallocated storage used\n',(k-1)/sz*100));
    if k>sz, warning('more hits than expected, adapt code'); end
    % flip y-axis, as is customary in top110 code:
    tmp=1:nely*nelx; tmp=reshape(tmp(:),nely,nelx); tmp=tmp(end:-1:1,:); tmp=tmp(:);
    jvec=tmp(jvec);
    tmp=1:nemy*nemx; tmp=reshape(tmp(:),nemy,nemx); tmp=tmp(end:-1:1,:); tmp=tmp(:);
    ivec=tmp(ivec);
    
    F=sparse(ivec,jvec,svec,nemx*nemy,nelx*nely);
    
    % BACKWARD MAPPING
    % x_FE(:) = B*x_extended(:)
    [fex,fey]=meshgrid(ha*([1:nelx]-.5),hb*([1:nely]-.5)); % FE mesh element centers
    [Rfex,Rfey]=rotfield(fex,fey,degrees,0,0);
    % align left and bottom boundary
    % emx/emy are the centroids of the embedding grid, before rotating.
    
    % shift to align with the embedding grid
    %Rfex=Rfex-x1-(emx(1)+0.5*em_hx);
    %Rfey=Rfey-y1-(emy(1)+0.5*em_hy);
    Rfex=Rfex-(emx(1)+x1)+em_hx/2;
    Rfey=Rfey-(emy(1)+y1)+em_hy/2;
    if DEBUG_PLOTS
        figure; set(plot(Rfex,Rfey,'x'),'linestyle','none'); axis equal
        hold on;
        % for plot, shift the em-points, as I have now positioned the Rfe*
        % points relative to the em_grid, so left points should be at
        % em_hx/2.
        set(plot(emx-emx(1)+em_hx/2,emy-emy(1)+em_hy/2,'.k'),'linestyle','none');
        [tmpx,tmpy]=rotfield(Remv(:,1),Remv(:,2),degrees,0,0);
        tmpx=tmpx-min(tmpx); tmpy=tmpy-min(tmpy);
        patch('faces',Remf,'vertices',[tmpx tmpy],'facecolor','none','edgecolor','k','linestyle',':');
    end
    if DEBUG_PLOTS
        % determine rotated cells:
        tmpx=repmat(Rfex(:),1,4);ntmp=size(tmpx,1); tmpx=tmpx';
        tmpy=repmat(Rfey(:),1,4);tmpy=tmpy';
        tmpv=[tmpx(:) tmpy(:)];
        vrts=[-1 -1; 1 -1; 1 1; -1 1]/2;
        vrts(:,1)=vrts(:,1)*ha;
        vrts(:,2)=vrts(:,2)*hb;
        [vrx,vry]=rotfield(vrts(:,1),vrts(:,2),degrees,0,0);
        vrts=[vrx, vry];
        tmpf=4*(1:ntmp)'-3; tmpf=[tmpf tmpf+1 tmpf+2 tmpf+3];
        % flip y direction to match data ordering convention:
        k=reshape(1:size(tmpf,1),size(Rfex,1),size(Rfex,2));
        k=k(end:-1:1,:); tmpf=tmpf(k(:),:);
        % also plot 'cells':
        if 1
            patch('faces',tmpf,'vertices',tmpv+repmat(vrts,ntmp,1),'facecolor','none','edgecolor','r','linestyle',':');
        end
    end
    % determine sampling points
    GaussPts=GaussPts0;
    GaussPts(:,1)= GaussPts0(:,1)*ha; GaussPts(:,2)= GaussPts0(:,2)*hb;
    [RGfex, RGfey]=rotfield(GaussPts(:,1),GaussPts(:,2),degrees,0,0);
    if DEBUG_PLOTS
        if 1
            for i=1:nely
                for j=1:nelx
                    patch('faces',[1:nGaussPts],'vertices',[RGfex RGfey]+repmat([Rfex(i,j) Rfey(i,j)],nGaussPts,1),'edgecolor','r','marker','.','linestyle','none','facecolor','none');
                end
            end
        end
    end
    % first, count x-elements that have an interaction with the BED-domain:
    Rfex=Rfex(:); Rfey=Rfey(:);
    xin_sure = (ceil(Rfex/em_hx)>=0) & (ceil(Rfex/em_hx)<=(nemx+1));
    yin_sure = (ceil(Rfey/em_hy)>=0) & (ceil(Rfey/em_hy)<=(nemy+1));
    in_sure = xin_sure & yin_sure;
    n_in = sum(in_sure); n_x = length(in_sure);
    % alloc sufficient space for ijs data:
    patchsize=3*3;
    avghits=4; % expected nr of elements overlapped by rotated element
    sz=n_in*avghits;
    ivec=ones(sz,1); jvec=ones(sz,1); svec=zeros(sz,1);
    % loop over x-elements to add their contributions to ijs data:
    % init 3x3 patch offsets
    patch_offset_y=[-1 0 1 -1 0 1 -1 0 1]';
    patch_offset_x=[-1 -1 -1 0 0 0 1 1 1]';
    k = 1;
    for i=1:n_x
        xi=ceil(Rfex(i)); yi=ceil(Rfey(i)); % index of x-element center in BED-domain
        xp=patch_offset_x+xi; yp=patch_offset_y+yi;
        patch_i=yp+nemy*(xp-1);
        patch_in_domain=all([xp>0, xp<=nemx, yp>0, yp<=nemy],2);
        %if in_sure(i) % skip all 'outside' cases
        npin=sum(patch_in_domain);
        if npin>0 % skip all 'outside' cases
            Gx=RGfex+Rfex(i);
            Gy=RGfey+Rfey(i);
            xk=ceil(Gx/em_hx); yk=ceil(Gy/em_hy); % indices of Gausspoints in BED-domain
            in_domain = all([xk>0, xk<=nemx, yk>0, yk<=nemy],2);
            relxk=xk-xi+1+1; relyk=yk-yi+1+1; % +1 for 1-based indexing, +1 for offset from corner of 3x3 patch
            hitcount=accumarray(relyk+(relxk-1)*3, in_domain,[patchsize 1], @sum, 0); % Gauss point hits per element of 3x3 patch
            hitsum=sum(hitcount); % number of Gauss points in x-domain
            if hitsum
                % add entries to ijs vectors:
                hit_j=find(hitcount>0)';
                for j=hit_j
                    ivec(k)=i;
                    jvec(k)=patch_i(j);
                    svec(k)=hitcount(j)/hitsum;
                    k=k+1;
                end
            end
        end
    end
    disp(sprintf('%.1f%% of preallocated storage used\n',(k-1)/sz*100));
    if k>sz, warning('more hits than expected, adapt code'); end
    % flip y-axis, as is customary in top110 code:
    tmp=1:nely*nelx; tmp=reshape(tmp(:),nely,nelx); tmp=tmp(end:-1:1,:); tmp=tmp(:);
    ivec=tmp(ivec);
    tmp=1:nemy*nemx; tmp=reshape(tmp(:),nemy,nemx); tmp=tmp(end:-1:1,:); tmp=tmp(:);
    jvec=tmp(jvec);
    
    B=sparse(ivec,jvec,svec,nelx*nely,nemx*nemy);
    %B=F'; NOT CORRECT
    if DEBUG_PLOTS
        % test mapping
        rho_test=zeros(nemy,nemx);
        rho_test(:)=F*(1-xdesign(:));
        backmapped=zeros(nely,nelx);
        backmapped(:)=B*rho_test(:);
        figure; subplot(2,2,1);
        if 1
            % using rotate:
            rt=plotpatch(backmapped,em_hx,em_hy); set(rt,'edgecolor','k');
            %rotate(rt,[0 0 1],-F(fi).degrees)
        else
            % using generated datastructures:
            %             patch('faces',F(fi).Remf,'vertices',F(fi).Remv,...
            %                 'facevertexcdata',rho_test(:),'edgecolor','k','facecolor','flat');
        end
        axis tight; axis equal;
        subplot(2,2,2); imagesc(1-xdesign); axis tight; axis equal;
        subplot(2,2,3); rt=plotpatch(rho_test,em_hx,em_hy); set(rt,'edgecolor','k'); axis tight; axis equal;
        subplot(2,2,4); imagesc(full(B*F)); axis tight; axis equal; colorbar;
        %disp('Press key to proceed ...'); pause
    end
end



function [F,B,nemy,nemx]=Original_FieldRotationMappings(nely,nelx,degrees)
% constructs mapping matrices for rotating a field [nely,nelx] over 'degrees'.
% F = forward mapping, B = return mapping. nemy, nemx are dimensions of the
% rotated field.
%
if length(degrees)~=1, error('expecting single orientation'); end
em_hx=1; em_hy=1; ha=1; hb=1; % no need to make this user-controlled here
DEBUG_PLOTS=0; % careful, don't use too big meshes with this one
% corner points of reference domain:
cp=[0 0; nelx 0; nelx nely; 0 nely]; % corner vertices
% rotated FE domain: Rcp
[Rcx,Rcy]=rotfield(cp(:,1),cp(:,2),degrees,0,0); Rcp=[Rcx Rcy];
x1=min(Rcx); x2=max(Rcx);
y1=min(Rcy); y2=max(Rcy);
% embedding domain: vbox
vbox=[x1 y1; x2 y1; x2 y2; x1 y2];
% embedding domain rotated back:
[Rvbx,Rvby]=rotfield(vbox(:,1),vbox(:,2),-degrees,0,0); Rvbox=[Rvbx Rvby];
if DEBUG_PLOTS
    figure; cla; subplot(2,2,1);
    xdesign = ones(nely,nelx);
    for i=1:nely, for j=1:nelx, xdesign(i,j)=.5+0.2*j/nelx+0.5*sin((i+j)/10); end; end
    f1p1=plotpatch(1-xdesign,em_hx,em_hy); set(f1p1,'edgecolor','k');
    % rotated FE domain: Rcp
    patch('faces',1:4,'vertices',Rcp,'facecolor','none','edgecolor','r');
    % embedding AM domain
    patch('faces',1:4,'vertices',vbox,'facecolor','none','edgecolor','b');
    % embedding domain rotated back:
    patch('faces',1:4,'vertices',Rvbox,'facecolor','none','edgecolor','m','linewidth',2);
    patch('faces',1:2,'vertices',Rvbox(1:2,:),'facecolor','none','edgecolor','g','linewidth',2);
    axis tight; axis equal;
end
% determine dimensions of embedding domain:
RL=x2-x1; RH=y2-y1;
em_nx=ceil(0.9999*RL/em_hx); em_ny=ceil(0.9999*RH/em_hy); em_xskip = (RL-em_nx*em_hx)/2; em_yskip=(RH-em_ny*em_hy)/2;
[emx, emy]=meshgrid(0.5*em_hx+[em_xskip:em_hx:0.9999*RL],0.5*em_hy+[em_yskip:em_hy:0.9999*RH]); % element centers
%patch('faces',[1:((em_nx+1)*(em_ny+1))]','vertices',[emx(:) emy(:)],'marker','.','edgecolor','k');
[Remx, Remy]=rotfield(emx+x1,emy+y1,-degrees,0,0);
if DEBUG_PLOTS
    % plot centroid nodes
    patch('faces',[1:em_nx*em_ny]','vertices',[Remx(:) Remy(:)],'marker','.','edgecolor','k');
end
if DEBUG_PLOTS
    % determine rotated cells:
    tmpx=repmat(Remx(:),1,4);ntmp=size(tmpx,1); tmpx=tmpx';
    tmpy=repmat(Remy(:),1,4);tmpy=tmpy';
    tmpv=[tmpx(:) tmpy(:)];
    vrts=[-1 -1; 1 -1; 1 1; -1 1]/2;
    vrts(:,1)=vrts(:,1)*em_hx;
    vrts(:,2)=vrts(:,2)*em_hy;
    [vrx,vry]=rotfield(vrts(:,1),vrts(:,2),-degrees,0,0);
    vrts=[vrx, vry];
    tmpf=4*(1:ntmp)'-3; tmpf=[tmpf tmpf+1 tmpf+2 tmpf+3];
    % flip y direction to match data ordering convention:
    k=reshape(1:size(tmpf,1),size(Remx,1),size(Remx,2));
    k=k(end:-1:1,:); tmpf=tmpf(k(:),:);
    %
    Remf=tmpf;
    Remv=tmpv+repmat(vrts,ntmp,1);
    % also plot 'cells':
    if 1
        patch('faces',Remf,'vertices',Remv,'facecolor','none','edgecolor','r','linestyle',':');
    end
end
% determine sampling points in [-.5 .5]^2 square:
%https://pomax.github.io/bezierinfo/legendre-gauss.html
GaussPts = 1/sqrt(3)*[1 1; -1 1; -1 -1; 1 -1]/2;
%[xg,yg]=meshgrid(-0.495:0.01:0.495,-0.495:0.01:0.495); GaussPts0=[xg(:) yg(:)]; GaussPts=GaussPts0; % ULTRA-PRECISE!!
% 10x10 scheme: very precise
%[xg,yg]=meshgrid(-0.45:0.1:0.45,-0.45:0.1:0.45); GaussPts0=[xg(:) yg(:)]; GaussPts=GaussPts0;
% 5x5 scheme: just fine
[xg,yg]=meshgrid(-0.4:0.2:0.4,-0.4:0.2:0.4); GaussPts0=[xg(:) yg(:)]; GaussPts=GaussPts0;
GaussPts(:,1)= GaussPts0(:,1)*em_hx; GaussPts(:,2)= GaussPts0(:,2)*em_hy;
nGaussPts = length(GaussPts);
[RGx, RGy]=rotfield(GaussPts(:,1),GaussPts(:,2),-degrees,0,0);
if DEBUG_PLOTS
    % plot a few for illustration
    if 0
        for i=1:em_ny
            for j=1:em_nx
                patch('faces',[1:nGaussPts],'vertices',[RGx RGy]+repmat([Remx(i,j) Remy(i,j)],nGaussPts,1),'edgecolor','r','marker','.','linestyle','none','facecolor','none');
            end
        end
        %             for k=-1:3
        %                 i=floor(em_nx/2)*em_ny-floor(em_ny/2)+k*em_ny;
        %                 for j=-2:2
        %                     patch('faces',[1 3; 2 4],'vertices',[RGx RGy]+repmat([Remx(i+j) Remy(i+j)],4,1),'edgecolor','r','marker','.');
        %                 end
        %             end
    end
end
% MAPPING
% determine mapping of points of the embedding mesh to FE elements (densities):
% nmap(i,j): number of mapped elements linked to a certain point
% imap(i,j,4): index 'i' into xdesign of elements linked to point (i,j)
% jmap(i,j,4): idem, index 'j'
% so if nmap(p,q)=1 --> link to xdesign(imap(p,q,1),jmap(p,q,1))
[nemy,nemx]=size(Remx);
nmap=zeros(nemy,nemx);
imap=zeros(nemy,nemx,nGaussPts)+nan; jmap=zeros(nemy,nemx,nGaussPts)+nan;
elhx=1; elhy=1; % element dimensions
for i=1:nemy
    for j=1:nemx
        % Gauss point locations around point (i,j) of embedding grid:
        Gx=RGx+Remx(i,j);
        Gy=RGy+Remy(i,j);
        m=0;
        for k=1:nGaussPts
            xk=ceil(Gx(k)/elhx); yk=ceil(Gy(k)/elhy);
            if all([xk>0, xk<=nelx, yk>0, yk<=nely])
                % in domain: add to list
                m=m+1;
                imap(i,j,m)=yk;
                jmap(i,j,m)=xk;
            end
        end
        nmap(i,j)=m;
    end
end
% flip y-axis, as is customary in top110 code:
nmap=nmap(end:-1:1,:);
imap=(nely+1)-imap(end:-1:1,:,:);
jmap=jmap(end:-1:1,:,:);
% next: put it in a matrix for easy mapping.
% suppose I have a field 'emb' on the embedding points, and a field 'x' on
% the FE mesh. I want to map values of 'x' on the 'emb' field.
% nicest is to do this in a vectorized way, i.e. forming the design fields
% into column vectors:
%   emb(ii)=M*x(:)
% where M is a (nemy*nemx)-by-(nely*nelx) sparse matrix.
sz=sum(nmap(:));
ivec=zeros(sz,1); jvec=zeros(sz,1); svec=zeros(sz,1);
k=1; % index in sz/ivec/jvec/svec
% (:) gives a column-wise numbering. This results in 'i + ny*(j-1)'
% statements below:
for i=1:nemy
    for j=1:nemx
        for ii=1:nmap(i,j)
            ivec(k)=i+nemy*(j-1);
            jvec(k)=imap(i,j,ii)+nely*(jmap(i,j,ii)-1);
            svec(k)=1/nmap(i,j);
            k=k+1;
        end
    end
end
F=sparse(ivec,jvec,svec,nemx*nemy,nelx*nely);
if DEBUG_PLOTS
    % test mapping
    rho_test=zeros(nemy,nemx);
    rho_test(:)=F*(1-xdesign(:));
    subplot(2,2,2);
    if 0
        % using rotate:
        rt=plotpatch(rho_test,em_hx,em_hy); set(rt,'edgecolor','k');
        rotate(rt,[0 0 1],-degrees)
    else
        % using generated datastructures:
        patch('faces',Remf,'vertices',Remv,...
            'facevertexcdata',rho_test(:),'edgecolor','k','facecolor','flat');
    end
    axis tight; axis equal;
    subplot(2,2,3); imagesc(1-xdesign); axis tight; axis equal;
    subplot(2,2,4); rt=plotpatch(rho_test,em_hx,em_hy); set(rt,'edgecolor','k'); axis tight; axis equal;
    %disp('Press key to proceed ...'); pause
end
% reduce memory footprint:
clear imap jmap ivec jvec svec 
% Reverse Mapping
% x_FE(:) = M2*x_extended(:)
[fex,fey]=meshgrid(ha*([1:nelx]-.5),hb*([1:nely]-.5)); % FE mesh element centers
[Rfex,Rfey]=rotfield(fex,fey,degrees,0,0);
% align left and bottom boundary
% emx/emy are the centroids of the embedding grid, before rotating.

% shift to align with the embedding grid
%Rfex=Rfex-x1-(emx(1)+0.5*em_hx);
%Rfey=Rfey-y1-(emy(1)+0.5*em_hy);
Rfex=Rfex-(emx(1)+x1)+em_hx/2;
Rfey=Rfey-(emy(1)+y1)+em_hy/2;
if DEBUG_PLOTS
    figure; set(plot(Rfex,Rfey,'x'),'linestyle','none'); axis equal
    hold on;
    % for plot, shift the em-points, as I have now positioned the Rfe*
    % points relative to the em_grid, so left points should be at
    % em_hx/2.
    set(plot(emx-emx(1)+em_hx/2,emy-emy(1)+em_hy/2,'.k'),'linestyle','none');
    [tmpx,tmpy]=rotfield(Remv(:,1),Remv(:,2),degrees,0,0);
    tmpx=tmpx-min(tmpx); tmpy=tmpy-min(tmpy);
    patch('faces',Remf,'vertices',[tmpx tmpy],'facecolor','none','edgecolor','k','linestyle',':');
end
if DEBUG_PLOTS
    % determine rotated cells:
    tmpx=repmat(Rfex(:),1,4);ntmp=size(tmpx,1); tmpx=tmpx';
    tmpy=repmat(Rfey(:),1,4);tmpy=tmpy';
    tmpv=[tmpx(:) tmpy(:)];
    vrts=[-1 -1; 1 -1; 1 1; -1 1]/2;
    vrts(:,1)=vrts(:,1)*ha;
    vrts(:,2)=vrts(:,2)*hb;
    [vrx,vry]=rotfield(vrts(:,1),vrts(:,2),degrees,0,0);
    vrts=[vrx, vry];
    tmpf=4*(1:ntmp)'-3; tmpf=[tmpf tmpf+1 tmpf+2 tmpf+3];
    % flip y direction to match data ordering convention:
    k=reshape(1:size(tmpf,1),size(Rfex,1),size(Rfex,2));
    k=k(end:-1:1,:); tmpf=tmpf(k(:),:);    
    % also plot 'cells':
    if 1
        patch('faces',tmpf,'vertices',tmpv+repmat(vrts,ntmp,1),'facecolor','none','edgecolor','r','linestyle',':');
    end
end
% determine sampling points
GaussPts=GaussPts0;
GaussPts(:,1)= GaussPts0(:,1)*ha; GaussPts(:,2)= GaussPts0(:,2)*hb;
[RGfex, RGfey]=rotfield(GaussPts(:,1),GaussPts(:,2),degrees,0,0);
if DEBUG_PLOTS
    if 1
        for i=1:nely
            for j=1:nelx
                patch('faces',[1:nGaussPts],'vertices',[RGfex RGfey]+repmat([Rfex(i,j) Rfey(i,j)],nGaussPts,1),'edgecolor','r','marker','.','linestyle','none','facecolor','none');
            end
        end
    end
end
% REVERSE MAPPING
% determine mapping of points of the embedding mesh to FE elements (densities):
nmap=zeros(nely,nelx);
imap=zeros(nely,nelx,nGaussPts)+nan; jmap=zeros(nely,nelx,nGaussPts)+nan;
%elhx=1; elhy=1; % element dimensions
for i=1:nely
    for j=1:nelx
        % Gauss point locations around point (i,j) of rotated FE grid:
        Gx=RGfex+Rfex(i,j);
        Gy=RGfey+Rfey(i,j);
        m=0;
        for k=1:nGaussPts
            xk=ceil(Gx(k)/em_hx); yk=ceil(Gy(k)/em_hy);
            if all([xk>0, xk<=nemx, yk>0, yk<=nemy])
                m=m+1;
                imap(i,j,m)=yk;
                jmap(i,j,m)=xk;
            end
        end
        nmap(i,j)=m;
    end
end
% flip y-axis, as is customary in top110 code:
nmap=nmap(end:-1:1,:);
imap=(nemy+1)-imap(end:-1:1,:,:);
jmap=jmap(end:-1:1,:,:);
% next: put it in a matrix for easy mapping.
% suppose I have a field 'emb' on the embedding points.
% I want to map values of 'emb' on the FE mesh.
% nicest is to do this in a vectorized way, i.e. forming the design fields
% into column vectors:
%   x(:)==RM*emb(:)
% where RM is a (nely*nelx)-by-(nemy*nemx) sparse matrix.
sz=sum(nmap(:));
ivec=zeros(sz,1); jvec=zeros(sz,1); svec=zeros(sz,1);
k=1; % index in sz/ivec/jvec/svec
% (:) gives a column-wise numbering. This results in 'i + ny*(j-1)'
% statements below:
for i=1:nely
    for j=1:nelx
        for ii=1:nmap(i,j)
            ivec(k)=i+nely*(j-1);
            jvec(k)=imap(i,j,ii)+nemy*(jmap(i,j,ii)-1);
            svec(k)=1/nmap(i,j);
            k=k+1;
        end
    end
end
B=sparse(ivec,jvec,svec,nelx*nely,nemx*nemy);
%B=F'; NOT CORRECT
if DEBUG_PLOTS
    % test mapping
    rho_test=zeros(nemy,nemx);
    rho_test(:)=F*(1-xdesign(:));
    backmapped=zeros(nely,nelx);
    backmapped(:)=B*rho_test(:);
    figure; subplot(2,2,1);
    if 1
        % using rotate:
        rt=plotpatch(backmapped,em_hx,em_hy); set(rt,'edgecolor','k');
        %rotate(rt,[0 0 1],-F(fi).degrees)
    else
        % using generated datastructures:
        %             patch('faces',F(fi).Remf,'vertices',F(fi).Remv,...
        %                 'facevertexcdata',rho_test(:),'edgecolor','k','facecolor','flat');
    end
    axis tight; axis equal;
    subplot(2,2,2); imagesc(1-xdesign); axis tight; axis equal;
    subplot(2,2,3); rt=plotpatch(rho_test,em_hx,em_hy); set(rt,'edgecolor','k'); axis tight; axis equal;
    subplot(2,2,4); imagesc(full(B*F)); axis tight; axis equal; colorbar;
    %disp('Press key to proceed ...'); pause
end
end


function h=plotpatch(field,ha,hb)
[nely,nelx]=size(field);
s=surface(zeros(nely+1,nelx+1)); pdat=surf2patch(s); delete(s);
pdat.facevertexcdata=field(:);
pdat.vertices(:,[1 2])=pdat.vertices(:,[1 2])-1;
pdat.vertices(:,1)=pdat.vertices(:,1)*ha;
pdat.vertices(:,2)=pdat.vertices(:,2)*hb;
pdat.vertices(:,2)=nely*hb-pdat.vertices(:,2);
h=patch(pdat,'edgecolor','none','facecolor','flat');
end

function [Rx,Ry]=rotfield(x,y,degrees,xC,yC)
% rotates fields x,y (arbitrary dimensions) by a certain angle (CCW) around
% a center point (xC, yC).
angle=degrees*pi/180; s=sin(angle); c=cos(angle);
%R=[c -s; s c];
dx=x-xC; dy=y-yC;
Rx=c*dx-s*dy + xC;
Ry=s*dx+c*dy + yC;
end    