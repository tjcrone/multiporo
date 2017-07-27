function [input] = procinput(filename)
% This function processes FILENAME into a structure called INPUT that main.m
% uses for a model run.

% read input file into cells
fid = fopen(filename);
s = textscan(fid, '%s%s', 'CommentStyle', '%', 'Delimiter', {'=', ' ', '\t'}, ...
  'MultipleDelimsAsOne', 1);
fid = fopen(fid);

% assign cell parameters into an input structure
for i = 1:length(s{1})
  try
    eval(sprintf('input.%s = %s;', s{1}{i}, s{2}{i})); % read line as double
  catch
    eval(sprintf('input.%s = ''%s'';', s{1}{i}, s{2}{i})); % else read line as string
  end
  % assign into caller workspace (not used)
  %eval(sprintf('assignin(''caller'', ''%s'', %s)', s{1}{i}, s{1}{i})); 
end

% build out domain geometry
input.x = linspace(input.d/2, (input.nx-1)*input.d, input.nx);
input.z = linspace(input.d/2, (input.nz-1)*input.d, input.nz);
[input.X, input.Z] = meshgrid(input.x, input.z);

% build out porosity field
if strcmp(input.porosity_type, 'uniform')
  input.phi = ones(input.nz, input.nx)*input.phi;
end

% build out permeability field
if strcmp(input.permeability_type, 'uniform')
  input.kx = ones(input.nz, input.nx)*input.kx_1;
  input.kz = ones(input.nz, input.nx)*input.kz_1;
end

% build out temperature field
if strcmp(input.T_type, 'uniform')
  input.T = ones(input.nz, input.nx)*input.T;
end
%rng('default');
%T = Z*0+Tcold+rand(nz,nx)*Thot/2;
%T = repmat(-sin(linspace(0,1,nx)*2*pi*8.5)*(min(min(Z/(d*nz)*Thot))), nz,1) + ...
%  Z/(d*nz)*Thot + rand(nz,nx);
%T(~cracked) = Thot;

% temperature boundary conditions (0=Neumann 1=Dirichlet)
% first row/column is value, second is type
input.Tbt = [ones(1,input.nx)*input.Ttop_value; ones(1,input.nx)*input.Ttop_type];
input.Tbb = [ones(1,input.nx)*input.Tbottom_value; ones(1,input.nx)*input.Tbottom_type];
input.Tbr = [ones(input.nz,1)*0 ones(input.nz,1)*0]; % Neumann zero
input.Tbl = [ones(input.nz,1)*0 ones(input.nz,1)*0]; % Neumann zero

% pressure boundary conditions (0=Neumann 1=Dirichlet)
input.Pbt = [ones(1,input.nx).*input.Ptop;ones(1,input.nx)*1]; % open
input.Pbb = [ones(1,input.nx).*0;ones(1,input.nx)*0]; % closed
input.Pbr = [ones(input.nz,1).*0 ones(input.nz,1)*0]; % closed
input.Pbl = [ones(input.nz,1).*0 ones(input.nz,1)*0]; % closed

% starting pressure field
%[P,Pbound,dPdzbound,rhobound] = initp(nx,nz,T,Tbt,Tbb,Ptop,TT, PP,RHO,g,d);

% cracked stuff
%cracked = logical(Z*0+1);
%cracked(Z>frontdepth) = 0;
%input.kx(~cracked) = koff;
%input.kz(~cracked) = koff;

% poroelastic stuff
%Ks = 80e9; % rock bulk modulus (basalt/diabase)
%Vpu = ones(nz,nx)*7100; % undrained p-wave velocity
%nuu = ones(nz,nx)*0.2498511904; % undrained Poisson's ratio
%Tvp = 600; %temperature condition for p-wave value
%Pvp = 30e6; %pressure condition (Pa) for p-wave value
%[G,K,Ku,lamdal,nu,biot,gamma,Se] = elastconst(Vpu,nuu,phi,Ks,rhom,Tvp,Pvp);

% restart stuff
if isfield(input, 'restart_file')
  R = load(restartfile, 'T2', 'P2', 't');
  input.T = T2;
  input.P = P2;
  input.t = t;
end
%  [m, n] = size(R.T2);
%  if n~=nx
%    error('Restarted domains must have the same number of columns as the original.');
%  end
%  cracked(1:m,:) = R.cracked;
%  T(1:m,:) = R.T2;
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
%end

% load or globalize thermodynamic tables
%global TT PP RHO CP
%if isempty(TT)
%   load('../hydrotables/hydrotab8.mat');
%end

%if restart==1
%  P(1:m,:) = R.P2;
  %if sum(size(P)==size(Pres))~=2
   % error('Geometry resizing not yet allowed');
  %else % same size as restart
  %  P = Pres;
  %end 
%end
