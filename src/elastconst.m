function [G,K,Ku,lamdal,nu,biot,gamma, Se] = elastconst(Vpu,nuu,phi,Ks,rhom,Tvp,Pvp)
%This function computes all of the needed drained moduli for the
%poroelastic convection code.  This function uses the undrained (in situ)
%seismic parameters Vpu and nuu, and assumes that the pore fluid to obtain
%these seismic values is cold and has cold thermodynamic and transport
%properties.  This function then uses the porosity, phi and the grain
%density rhom to compute the total saturated material density.  It uses
%a simple relationship to compute the undrained Vsu, and uses the velocity
%form of Gassmans's relations to determine the drained bulk modulus, K.  It
%then computes the shear modulus G, the biot-willis parameter biot, 
%and the lame parameter, lamdal.  The values output from this function 
%should not change during the course of a PEC2D run, and can be included
%in the input file.  The function uses Tvp and Pvp for the fluid conditions
%(temp and pressure) to compute the fluid density and bulk modulus, and 
%corresponding composite undrained properties.  A series of checks at the end
%of this function warns the user about possible problems with the parameter
%calculations.


%load or globalize thermodynamic tables
global TT PP RHO CP BETA
if isempty(TT)
  %load('../hydrotables/hydrotab8.mat');
  load('../../Hydro/hydrotab5.mat');
end

%determine the cold fluid properties
rhof = interptim(PP,TT,RHO,Pvp./100000,Tvp); %fluid density
betaf = interptim(PP,TT,BETA,Pvp./100000,Tvp); %fluid compressibility
Kf = 1./betaf; %fluid bulk modulus

%determine some composite, undrained properties
rhosat = (1-phi).*rhom + phi.*rhof; %saturated medium density
Vsu = -Vpu.*(2*(nuu-1).*(2*nuu-1)).^.5./(2*nuu-2);  %shear wave velocity (undrained)
Gu = rhosat.*Vsu.^2; %undrained shear modulus (same as drained)
Ku = rhosat.*(Vpu.^2 - (4/3).*Vsu.^2); %undrained bulk modulus
lamdalu = Ku - 2.*Gu./3; %undrained lame parameter lamda

%use gassmann's relation to determine drained bulk modulus
K = drainedbulk(Vpu,Vsu,rhosat,Kf,Ks,phi); %drained bulk modulus
%Kutest = pec_gassmann(K,Kf,Ks,phi);

%determine other drained moduli from K and Gu
G = Gu; %shear modulus should not change under drained conditions
lamdal = K - 2.*G./3; %drained lame parameter lamda
nu = lamdal ./ (3.*K - lamdal); %drained Poisson's raio
biot = 1 - K./Ks; %Biot-Willis parameter

%determine other important parameters
Se = biot.^2./(Ku - K); %compute Se using Wang 3.41
B = (1-K./Ku)./(1-K./Ks); %Skempton's coefficient (cold)
gamma = B.*(1+nuu)./3./(1-nuu); %Loading efficiency (cold)

%some wang & davis checks:
alpha_wd = biot;
beta_wd = alpha_wd./(alpha_wd + phi.*K.*(1./Kf - 1./Ks));
gamma_wd = beta_wd.*(1+nu)./(3.*(1-nu)-2.*alpha_wd.*beta_wd.*(1-2.*nu));
wronggamma = beta_wd.*(1+nuu)./(3.*(1-nuu)-2.*alpha_wd.*beta_wd.*(1-2.*nuu));

%parameter value checks
if min(min(K))<=0
    warning('Drained frame bulk modulus is negative.');
end
if max(max(K>Ks))>0
    warning('K is greater than Ks');
end
if min(min(lamdal)) <= 0
   warning('Drained Lame parameter is negative');
end
if min(min(nu))<0
    warning('Drained Poisson''s ratio is negative');
end
if max(max((nuu./nu)))>1.5
    warning('The ratio of undrained to drained Poisson''s ratio is over 1.5');
end
if max(max(B))>1
    warning('Skempton''s coefficient is over unity');
end
if min(min(B))<0
    warning('Skempton''s coefficient is negative');
end
