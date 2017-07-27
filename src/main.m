function [] = main(inputfile)
% This is the main multiporo Matlab function.

% process input file
input = procinput(inputfile);

% break variables out of input (this should be temporary)
nx = input.nx;
nz = input.nz;
d = input.d;
T1 = input.T;
Tbt = input.Tbt;
Tbb = input.Tbb;
Tbr = input.Tbr;
Tbl = input.Tbl;
kx = input.kx;
kz = input.kz;
Ptop = input.Ptop;
Pbt = input.Pbt;
Pbb = input.Pbb;
Pbr = input.Pbr;
Pbl = input.Pbl;
g = input.g;
nstep = input.nstep;
output_interval = input.output_interval;
rhom = input.rhom;
cm = input.cm;
phi = input.phi;
lambdam = input.lambdam;

% globalize thermodynamic tables
global TT PP RHO CP
if isempty(TT)
    load('../hydrotables/hydrotab8.mat');
end

% calculate starting pressure field if not defined
if ~exist('init.P')
  [P1,Pbound,dPdzbound,rhobound] = initp(nx,nz,T1,Tbt,Tbb,Ptop,g,d);
else
  P1 = input.P;
  Pbound = input.Pbound;
  dPdzbound = input.dPdzbound;
end

% compute T-P dependent fluid properties (t=1)
mu1 = dynvisc(T1); %fluid viscosity
rhof1 = interptim(PP,TT,RHO,P1./100000,T1); %fluid density
cf1 = interptim(PP,TT,CP,P1./100000,T1); %fluid heat capacity
Se1 = T1*0+1e-32;

% compute boundary fluid properties
rhobt = interptim(PP,TT,RHO,Tbt(1,:)*0+Ptop./100000,Tbt(1,:));
rhobb = rhobound;
rhobr = interptim(PP,TT,RHO,P1(:,end)./100000,Tbr(:,1));
rhobl = interptim(PP,TT,RHO,P1(:,1)  ./100000,Tbl(:,1));
cfbt = interptim(PP,TT,CP,Tbt(1,:)*0+Ptop./100000,Tbt(1,:));
cfbb = interptim(PP,TT,CP,Pbound./100000,Tbb(1,:));
cfbr = interptim(PP,TT,CP,P1(:,end)./100000,Tbr(:,1));
cfbl = interptim(PP,TT,CP,P1(:,1)  ./100000,Tbl(:,1));

% compute darcy velocities (t=1)
[qx1,qz1] = darcy(nx,nz,P1,rhof1,rhobb,kx,kz,mu1,g,d,Pbt,Pbb,Pbr,Pbl,T1);

% get output filename base
underloc = strfind(inputfile, '_');
outfilenamebase = inputfile(1:underloc(end)-1);

% create tentative values at t=2
rhof2 = rhof1;
cf2 = cf1;
qx2 = qx1;
qz2 = qz1;
Se2 = Se1;
mu2 = mu1;
P2 = P1;
T2 = T1;

% store t=1 output
t = 0;
outfilename = [outfilenamebase, sprintf('_out_%010.0f.mat', t)];

if strcmp(input.model_type,'cracking_convection')
  save(outfilename, '-v7.3', 'rhof2', 'cf2', 'T2', 'P2', 'qx2', 'qz2', 'cracked', 't');
elseif strcmp(input.model_type,'convection')
  save(outfilename, '-v7.3', 'rhof2', 'cf2', 'T2', 'P2', 'qx2', 'qz2', 't');
end
nout = 1;

% initialize Tmax and dTmax
Tmax = max(max(T1));
dTmax = 0;

% start timer
tic;
etime = toc;

% initialize other flags and counters
dtadjust = 0;
j = 0;
stepsdone = 0;

% time loop
for i = 1:nstep-1

  % crack domain if cracking front model
  if strcmp(input.model_type,'cracking_convection')
    KI=Atemp*(Tve-T1)-Apress*(Pw+rhoR*g*Z);
    cracked = ((KI>=KIc)+cracked)>0;
    kx(:,:) = koff;
    kz(:,:) = koff;
    kx(cracked) = kon;
    kz(cracked) = kon;
  end

  % dt based on stepper_type
  if strcmp(input.stepper_type,'dT')
    j = j + 1;
    if i == 1
      dt = input.initial_dt;
      fprintf('Starting dt: %0.6f h\n', dt/60/60);
    end

    if dtadjust == 1
      dt = dtlast;
      dtadjust = 0;
      fprintf('Adjusting dt to %0.4f hours at t = %0.2f years\n', dt/60/60, t/60/60/24/365);
    end

    if Tmax > input.Thot + 1 || dTmax > input.maxdT
      dt = dt*0.8;
      fprintf('Reducing dt to: %0.6f h\n', dt/60/60);
      fprintf('Tmax: %0.2f\n', Tmax);
    elseif j >= input.increase_interval
      dt = dt*1.1;
      if dt > output_interval
        dt = output_interval
      end
      j = 0;
      fprintf('Increasing dt to %0.4f hours at t = %0.2f years\n', dt/60/60, t/60/60/24/365);
    end
    t2 = t+dt;

    % if crossing over output_interval, adjust dt
    if t2 > output_interval*nout
      dtlast = dt;
      dtadjust = 1;
      t2 = output_interval*nout;
      dt = t2-t;
      fprintf('Adjusting dt to %0.4f hours at t = %0.2f years\n', dt/60/60, t/60/60/24/365);
    end
  end

  % compute beta for temperature equation
  beta2 = reshape(rhom.*cm.*(1-phi) + rhof2.*cf2.*phi,nx*nz,1);

  % create del-squared stiffness matrix for diffusion term in heat eq.
  [AimpT,CimpT] = tdiffcoeff(nx,nz,d,lambdam,Tbr,Tbl,Tbb,Tbt,input.Ttop_conduction);

  % compute T2 using implicit technique
  [BimpT,DimpT] = tadvectcoeff(nx,nz,d,qx2,qz2,rhof2,cf2,Tbb,Tbt,Tbr,Tbl, ...
    rhobt,rhobb,rhobr,rhobl,cfbt,cfbb,cfbr,cfbl);
        
  % single step implicit left hand side stiffness:
  Tstiff = spdiags(beta2,0,nx*nz,nx*nz) - dt*(AimpT + BimpT);
        
  % single step implicit right hand side
  RHS = reshape(T1,nx*nz,1).*beta2 + dt*(CimpT + DimpT); %see p. 464
        
  % and solve:
  T2 = Tstiff\RHS;
  T2 = reshape(T2,nz,nx);
  T2(T2<0) = 0; % kluge to prevent negative temperatures

  % if max of T2 is greater than Thot, start step over
  Tmax = max(max(T2));
  dTmax = max(max(T2-T1));
  if Tmax > input.Thot + 1 || dTmax > input.maxdT
    continue;
  end

  % update T-P dependent fluid properties and darcy velocities
  mu2 = dynvisc(T2);
  rhof2 = interptim(PP,TT,RHO,P2./100000,T2); %fluid density
  cf2 = interptim(PP,TT,CP,P2./100000,T2); %fluid heat capacity
  [qx2,qz2] = darcy(nx,nz,P2,rhof2,rhobb,kx,kz,mu2,g,d,Pbt,Pbb,Pbr,Pbl,T2);

  % compute P2 using implicit technique
  [AimpP,BimpP,CimpP] = pstiff(nx,nz,d,Se2,rhof2,rhobt,rhobb,rhobr,rhobl, ...
    qx2,qz2,kx,kz,mu2,g,T2,Pbt,Pbb,Pbr,Pbl);

  % single step implicit left hand side stiffness:
  Pstiffness = AimpP;

  % single step implicit right hand side
  PRHS = -BimpP-CimpP;

  % and solve:
  P2 = Pstiffness\PRHS;
  P2 = reshape(P2,nz,nx);

  % update T-P dependent fluid properties and darcy velocities
  mu2 = dynvisc(T2);
  rhof2 = interptim(PP,TT,RHO,P2./100000,T2); %fluid density
  cf2 = interptim(PP,TT,CP,P2./100000,T2); %fluid heat capacity
  [qx2,qz2] = darcy(nx,nz,P2,rhof2,rhobb,kx,kz,mu2,g,d,Pbt,Pbb,Pbr,Pbl,T2);

  % shift variables
  P1 = P2;
  T1 = T2;   
  t = t2;
  stepsdone = stepsdone + 1;

  % write outputs to file
  if t == output_interval*nout
    %t_years = output_interval*nout/60/60/24/365;
    outfilename = [outfilenamebase, sprintf('_out_%010.0f.mat', t)];

    if strcmp(input.model_type,'cracking_convection')
      save(outfilename, '-v7.3', 'rhof2', 'cf2', 'T2', 'P2', 'qx2', 'qz2', 'cracked', 't');
    elseif strcmp(input.model_type,'convection')
      save(outfilename, '-v7.3', 'rhof2', 'cf2', 'T2', 'P2', 'qx2', 'qz2', 't');
    end
    nout = nout + 1;

    % output information
    laptime = toc-etime;
    slashloc = strfind(outfilename, '/');
    fprintf('\nFile %s saved\n', outfilename(slashloc(end)+1:end));
    %fprintf('Year: %i\n', t_years);
    fprintf('Step: %i\n', stepsdone);
    %fprintf('Average steps/year: %.0f\n', stepsdone/t_years);
    %fprintf('Wall time per %i years: %0.f s\n\n', output_interval/60/60/24/365, laptime);
    etime = toc;
  end

  % stop at stop_time
  if t == input.stop_time
    break;
  end
end

% print timing info
etime = toc;
fprintf('\nSimulation complete\n');
fprintf('Total wall time\t\t\t%.1f s\n',etime);
fprintf('Total model time\t\t%.1f years\n', t_years);
fprintf('Number of model steps\t\t%i steps\n',stepsdone);
fprintf('Average wall time per step\t%.2f s\n',etime/(stepsdone));
fprintf('Average wall time per %i years\t%0.2f s\n', output_interval/60/60/24/365, ...
  etime/t*output_interval);
