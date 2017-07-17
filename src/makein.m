function [tmpfilename] = makein(inputfile)
% This function creates the input .mat file for the main function from the simpler
% text input file. The .mat file produced will contain all the needed variables to
% begin a new porous convection run. All units are SI, unless otherwise specified.
%
% Timothy Crone (tjcrone@gmail.com)

% read input file
readinput(inputfile);

% time stepping information
outputinterval = outputinterval*365*24*60*60; % convert output interval in years to seconds

% domain geometry
x = linspace(d/2,(nx-1)*d,nx);
z = linspace(d/2,(nz-1)*d,nz);
[~,Z] = meshgrid(x,z);

% some constants
rhom = 2900; % rock or grain density (basalt)
cm = 1000; % rock heat capacity (basalt)
lamdam = 2; % rock thermal conductivity (basalt)
alpham = 2e-5; % rock thermal expansion coefficient (basalt)
phi = ones(nz,nx)*phi; % porosity
g = 9.8; % gravitational constant

% cracked boolean
cracked = logical(Z*0+1);
cracked(Z>frontdepth) = 0;

% some poroelastic stuff (not used for now)
%Ks = 80e9; % rock bulk modulus (basalt/diabase)
%Vpu = ones(nz,nx)*7100; % undrained p-wave velocity
%nuu = ones(nz,nx)*0.2498511904; % undrained Poisson's ratio
%Tvp = 600; %temperature condition for p-wave value
%Pvp = 30e6; %pressure condition (Pa) for p-wave value
%[G,K,Ku,lamdal,nu,biot,gamma,Se] = elastconst(Vpu,nuu,phi,Ks,rhom,Tvp,Pvp);

% initial temperature conditions
rng('default');
%T = Z*0+Tcold+rand(nz,nx)*Thot/2;
T = repmat(-sin(linspace(0,1,nx)*2*pi*8.5)*(min(min(Z/(d*nz)*Thot))), nz,1) + ...
  Z/(d*nz)*Thot + rand(nz,nx);
T(~cracked) = Thot;

% restart temperature and uncracked fields if required
if restart==1
  R = load(restartfile, 'T2', 'cracked', 'P2');
  [m, n] = size(R.T2);
  if n~=nx
    error('Restarted domains must have the same number of columns as the original.');
  end
  cracked(1:m,:) = R.cracked;
  T(1:m,:) = R.T2;
  % map Tres onto current geometry if necessary
  %if sum(size(T)==size(Tres))~=2
  %  if steady==0
  %    error('Geometry resizing only allowed when restarting into another steady state run.');
  %  end

  %if sum(size(T)==size(Tres))~=2 % changing geometry
  %  [m, n] = size(Tres);
  %  if n~=nx
  %    error('Restarted temperature field must have the same number of columns.');
  %  end
  %  T(1:m,:) = Tres;
  %else % same size as restart
  %  T = Tres;
  %  cracked(1:m,:) = R.crackedout(:,:,end);
  %end 

  %[m, n] = size(Tres);
  %if n~=nx
  %  error('Restarted temperature field must have the same number of columns.');
  %end
end

%initial permeability
kx = ones(nz,nx)*kon;  % permeability in x-direction
kz = ones(nz,nx)*kon;  % permeability in z-direction
kx(~cracked) = koff;
kz(~cracked) = koff;

% temperature boundary conditions (0=Neumann 1=Dirichlet)
% first row/column is value, second is type
Tbt = [ones(1,nx)*Tcold; ones(1,nx)*1]; % Dirichlet cold
Tbb = [ones(1,nx)*Tbottomvalue; ones(1,nx)*Tbottomtype];
Tbr = [ones(nz,1)*0 ones(nz,1)*0]; % Neumann zero
Tbl = [ones(nz,1)*0 ones(nz,1)*0]; % Neumann zero
%Tbl = [ones(nz,1)*Thot ones(nz,1)*1]; % Dirichlet hot

% load or globalize thermodynamic tables
global TT PP RHO CP
if isempty(TT)
   load('../hydrotables/hydrotab8.mat');
end

% calculate starting pressure field uxing initp.m
[P,Pbound,dPdzbound,rhobound] = initp(nx,nz,T,Tbt,Tbb,Ptop,TT, ...
    PP,RHO,g,d);
if restart==1
  P(1:m,:) = R.P2;
  %if sum(size(P)==size(Pres))~=2
   % error('Geometry resizing not yet allowed');
  %else % same size as restart
  %  P = Pres;
  %end 
end

% pressure boundary conditions (0=Neumann 1=Dirichlet)
Pbt = [ones(1,nx).*Ptop;ones(1,nx)*1]; % open
Pbb = [ones(1,nx).*0;ones(1,nx)*0]; % closed
Pbr = [ones(nz,1).*0 ones(nz,1)*0]; % closed
Pbl = [ones(nz,1).*0 ones(nz,1)*0]; % closed

% save the output to a temporary file
inoutdir = '../in_out/';
[status, tmpfilename] = system(sprintf('mktemp %sinput.XXXXXX', inoutdir));

% save variables to an input .mat file
save(tmpfilename(1:end-1), 'nstep','outputinterval','nx','nz','d','cm', ...
  'lamdam','phi','rhom','kx','kz','g','T','P','Tbb','Tbl','Tbr','Tbt','Ptop','Pbt', ...
  'Pbb','Pbl','Pbr','alpham','rhobound','Pbound','Ttopconduction','cracked','Thot', ...
  'kon','koff','Z','steady', 'Apress','Atemp','KIc','Pw','Tve','rhoR', 'stopyear', ...
  'maxdT', '-v7.3');
