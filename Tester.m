clc
clear
close all
% ------------ 0. PUT PolyTop.m (and PolyMesher.m) ON THE MATLAB PATH -----
% PolyMesher generates unstructured polygon meshes of rectangles and
% other 2-D regions.  Download it from the authors’ website or GitHub repo
% and keep both .m files in the same folder.

% ------------ 1. CREATE A POLYGONAL FEM MESH -----------------------------
Domain  = @(p)  dRectangle([0 120],[0  40],p);
[Node,Element] = PolyMesher(Domain,2400,0);   % 2 400 polygons ≈ 60×40 grid

% ------------ 2. BUILD THE FEM STRUCTURE ---------------------------------
fem.Node   = Node;                % NNode × 2
fem.Element= Element;             % {NElem × 1} cell
fem.NNode  = size(Node,1);
fem.NElem  = numel(Element);
fem.E0     = 1.0;                 % Young’s modulus (solid)
fem.Nu0    = 0.3;                 % Poisson’s ratio
fem.Reg    = 1;                  % all elements identical → reuse stiffness
% Loads: node#,  Fx,  Fy  (here: downward load at mid-span)
mid      = find(abs(Node(:,1)-60)<1e-6 & abs(Node(:,2)-40)<1e-6);
fem.Load = [mid  0  -1];
% Supports: node#,  fix-x?, fix-y? (roller, left edge; pin, right edge)
left     = find(abs(Node(:,1))<1e-6);
right    = find(abs(Node(:,1)-120)<1e-6);
fem.Supp = [left 1 1;   right 0 1];  % rollers on left, pin on right

% ------------ 3. SET OPTIMISATION PARAMETERS -----------------------------
opt.VolFrac = 0.4;
opt.Tol     = 1e-3;
opt.MaxIter = 50;
opt.zMin    = 0.0;
opt.zMax    = 1.0;
opt.zIni    = opt.VolFrac*ones(fem.NElem,1);
opt.OCMove  = 0.1;   % 10 % change allowed per iteration
opt.OCEta   = 0.5;   % damping exponent
% --- density filter matrix P  (identity → no filter) ---
opt.P       = speye(fem.NElem);
% --- SIMPLE material interpolation (SIMP with p = 3) ---
penalty     = 3;
opt.MatIntFnc = @(rho) deal( ...
    1e-3 + rho.^penalty.*(1-1e-3),              ... % E(rho)
    penalty*(1-1e-3)*rho.^(penalty-1),          ... % dE/drho
    rho,                                        ... % V(rho)
    ones(size(rho)));                               % dV/drho

% ------------ 4. RUN POLYTOP ---------------------------------------------
[z,V,fem] = PolyTop(fem,opt);
