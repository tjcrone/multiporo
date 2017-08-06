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
end

% fix output directory
if ~strcmp(input.output_dir(end),'/')
  input.output_dir = [input.output_dir,'/'];
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
  input.kx = ones(input.nz, input.nx)*input.kx;
  input.kz = ones(input.nz, input.nx)*input.kz;
elseif strcmp(input.permeability_type, 'fault')
  [input.kx, input.kz] = fault_permeability(input);
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

% timing/stepping information
if ~isfield(input, 'nstep')
  input.nstep = 0;
end

% starting pressure field
[input.P,input.Pbound,~,input.rhobound] = initp(input.nx,input.nz,input.T, ...
  input.Tbt,input.Tbb,input.Ptop,input.g,input.d,input.thermo_tables);

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
%[G,K,Ku,lamdal,nu,biot,gamma,Se] = elastconst(Vpu,nuu,phi,Ks,rhom,Tvp,Pvp, input.thermo_tables);

% restart options
if isfield(input, 'restart_file')
  R = load(input.restart_file, 'T2', 'P2', 't', 'd');
  [nz, nx] = size(R.T2);
  d = R.d;

  % restart time
  input.t = R.t;
  if isfield(input, 'stop_time')
    if input.stop_time <= input.t
      error('Input ''stop_time'' must be greater than the time associated with the restart file.');
    end
  end

  % resize if necessary
  if input.nx == nx && input.nz == nz && input.d == d % no resize
    input.T = R.T2;
    input.P = R.P2;
  else % resize
    if ~isfield(input, 'resize_type')
      error('When changing geometry the ''resize_type'' variable must be set.');
    end
    if strcmp(input.resize_type, 'resolve') % resolve-type resize
      if input.nx/input.nz ~= nx/nz || d/input.d ~= input.nx/nx
        error('New geometry incompatible with a resolve-type restart.');
      end
      input.T = R.T2(1+floor((0:input.nx-1)/(d/input.d)),1+floor((0:input.nx-1)/(d/input.d)));
      input.P = R.P2(1+floor((0:input.nx-1)/(d/input.d)),1+floor((0:input.nx-1)/(d/input.d)));
    elseif strcmp(input.resize_type, 'crack') % cracking-front-type vertical expansion
      if input.nx ~= nx
        error('Restarted temperature field must have the same number of columns.');
      end
      input.T(1:nz,:) = R.T2;
      input.P(1:nz,:) = R.P2;
      % input.cracked(1:m,:) = R.cracked;
    else
      error('Input ''resize_type'' not recognized');
    end
  end
end
