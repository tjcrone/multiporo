function [] = main(inputfile)
% This is the main poroelastic convection model function. It creates and loads the
% input variables specified by INPUTFILE, then commences a convection run
% saving the outputs to a file with the same prefix as the input file.
%
% Timothy Crone (tjcrone@gmail.com)

% run makein
tmpfilename = makein(inputfile);

% load the result of makein and delete the temporary mat file
kx = 0;
load(tmpfilename(1:end-1), '-mat');
system(sprintf('rm %s', tmpfilename(1:end-1)));

% globalize thermodynamic tables
global TT PP RHO CP
if isempty(TT)
    load('../hydrotables/hydrotab8.mat');
end

% add step signifier to some variables
T1 = T;
P1 = P;

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
t_years = 0;
outfilename = [outfilenamebase, sprintf('_out_%07.0f.mat', t_years)];
tout = 0;
save(outfilename, '-v7.3', 'rhof2', 'cf2', 'T2', 'P2', 'qx2', 'qz2', 'cracked', 'tout');
nout = 1;

% initialize Tmax
Tmax = 0;

% start timer
tic;
etime = toc;

% time loop
for i = 1:nstep-1

  % if this is not a steady state run, crack
  if steady==0
    KI=Atemp*(Tve-T1)-Apress*(Pw+rhoR*g*Z);
    cracked = ((KI>=KIc)+cracked)>0;
    kx(:,:) = koff;
    kz(:,:) = koff;
    kx(cracked) = kon;
    kz(cracked) = kon;
  end

  % set dt using adaptive or predefined time stepping
  if adaptivetime==1
    % adaptive time stepping based on CFL condition
    maxV = max(max(max(qx2)),max(max(qz2)));
    if i==1
      dt = 24*3600;
    elseif i<10
      dt = min([0.001*d/maxV 0.001*d^2/1e-6]);
    else
      dt = min([stepfraction*d/maxV stepfraction*d^2/1e-6]);
    end
    t(i+1)=t(i)+dt;
  else
    % use t vector for dt
    dt = t(i+1)-t(i);
    error('don''t do this.');
  end

  % if crossing over outputinterval, adjust dt
  if t(i+1) > outputinterval*nout
    t(i+1) = outputinterval*nout;
    dt = t(i+1)-t(i);
  end

  % picard iterations
  for j = 1:maxpicard
    % compute beta for temperature equation
    beta2 = reshape(rhom.*cm.*(1-phi) + rhof2.*cf2.*phi,nx*nz,1);
        
    % create del-squared stiffness matrix for diffusion term in heat eq.
    [AimpT,CimpT] = tdiffcoeff(nx,nz,d,lamdam,Tbr,Tbl,Tbb,Tbt,Ttopconduction);

    % compute T2 using implicit technique
    [BimpT,DimpT] = tadvectcoeff(nx,nz,d,qx2,qz2,rhof2,cf2,Tbb,Tbt,Tbr,Tbl, ...
      rhobt,rhobb,rhobr,rhobl,cfbt,cfbb,cfbr,cfbl);
        
    % single step implicit left hand side stiffness:
    Tstiff = spdiags(beta2,0,nx*nz,nx*nz) - dt*(AimpT + BimpT);
        
    % single step implicit right hand side
    RHS = reshape(T1,nx*nz,1).*beta2 + dt*(CimpT + DimpT); %see p. 464
        
    % and solve:
    %T2int = Tstiff\RHS;
    T2 = Tstiff\RHS;
    %T2int = reshape(T2int,nz,nx);
    T2 = reshape(T2,nz,nx);
    %T2int(T2int<0) = 0; % kluge to prevent negative temperatures
    T2(T2<0) = 0; % kluge to prevent negative temperatures
    %T2(T2>Thot) = Thot; % kluge to prevent overshoots

    % damping
    %tdamp = min([j*0.04 0.4]);
    %tdamp = j*0.01;
    %tdamp = 0.3;
    %T2 = T2last+tdamp*(T2-T2last);
    %T2 = T1+tdamp*(T2-T1);

    % update T-P dependent fluid properties and darcy velocities
    mu2 = dynvisc(T2);
    rhof2 = interptim(PP,TT,RHO,P2./100000,T2); %fluid density
    cf2 = interptim(PP,TT,CP,P2./100000,T2); %fluid heat capacity
    [qx2,qz2] = darcy(nx,nz,P2,rhof2,rhobb,kx,kz,mu2,g,d,Pbt,Pbb,Pbr,Pbl,T2);

    % calculate Se2 here using Ku2
        
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

    % convergence check
    %picardthresh = 0.01;
    %maxdel = find(abs(T2-T1)==max(max(abs(T2-T1))),1,'first');
    %maxdel = find(abs(T2-T2last)==max(max(abs(T2-T2last))),1,'first');
    %if abs(T2(maxdel)-T2last(maxdel)) < 0.001 && j > 2 % picardthresh && j>10
    %if abs(T2(maxdel)-T2last(maxdel))/T2last(maxdel) < picardthresh && j>10
    %if max(max(((abs(T2-T2last))./(abs(T2-T1))))) < picardthresh
    %  disp(j)
    %  outfileobj.npicard(1,i) = j;
    %  break;
    %end

    % break after max picard iterations
    if j==maxpicard
      break;
    end
  end

  % shift variables
  P1 = P2;
  T1 = T2;   

  % track Tmax
  Tmax = max([Tmax max(max(T1))]);

  % write outputs to file
  if t(i+1) == outputinterval*nout
    t_years = outputinterval*nout/60/60/24/365;
    outfilename = [outfilenamebase, sprintf('_out_%07.0f.mat', t_years)];
    tout = t(i+1);
    save(outfilename, '-v7.3', 'rhof2', 'cf2', 'T2', 'P2', 'qx2', 'qz2', 'cracked', 'tout');
    nout = nout + 1;
    if Tmax > Thot + 0.01
      %error('Tmax is greater than Thot.');
      fprintf('\nAdjusting T. Tmax: %.2f\n', Tmax);
      T1(T1>Thot) = Thot;
      T2(T2>Thot) = Thot;
    end

    % output information
    laptime = toc-etime;
    slashloc = strfind(outfilename, '/');
    fprintf('\nFile %s saved\n', outfilename(slashloc(end)+1:end));
    fprintf('Year: %i\n', t_years);
    fprintf('Step: %i\n', i);
    fprintf('Average steps/year: %.0f\n', i/t_years);
    fprintf('Wall time per %i years: %0.f s\n', outputinterval/60/60/24/365, laptime);
    etime = toc;
  end

  % stop at stopyear
  if t_years == stopyear
    break;
  end
end

% print timing info
etime = toc;
fprintf('\nSimulation complete\n');
fprintf('Total wall time\t\t\t%.1f s\n',etime);
fprintf('Total model time\t\t%.1f years\n', tout/60/60/24/365);
fprintf('Number of model steps\t\t%i steps\n',i+1);
fprintf('Average wall time per step\t%.2f s\n',etime/(i+1));
fprintf('Average wall time per %i years\t%0.2f s\n', outputinterval/60/60/24/365, ...
  etime/tout*outputinterval);
