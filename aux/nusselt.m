function [Nu] = nusselt(nx,nz,d,T,P,lamdam,Tbt,Tbb,kx,kz,g, ...
   Pbt,Pbb,Pbr,Pbl)
% This function computes the nusselt number for horizontal
% layers in a 2-D porous convection computational domain that
% is heated from below.  The top and bottom boundaries should
% be maintained at a constant uniform temperature, and the box
% should be closed to fluid flow.  The sides should be insulated.
% Nusselt numbers are computed only for the boundaries that
% separate horizontal control volume bands in the model domain.

% load thermo tables
global TT PP RHO CP
if isempty(TT)
   load('../hydrotables/hydrotab8.mat');
end

% fluid properties
rhof = interptim(PP,TT,RHO,P./100000,T);
cf = interptim(PP,TT,CP,P./100000,T);
mu = dynvisc(T);
[qx,qz] = darcy(nx,nz,P,rhof,rhobb,kx,kz,mu,g,d,Pbt,Pbb,Pbr,Pbl,T);

% conductive transport
Qc = -lamdam.*(Tbb(1,1)-Tbt(1,1)).*((nx+1)/(nz+1));

% conductive heat transport across each boundary (with convection)
Qcc = -lamdam.*(T(2:end,:)-T(1:end-1,:));

% advective heat transport
Qan = T(2:end,:) .* rhof(2:end,:) .* cf(2:end,:) .* qz(2:end-1,:) .* ...
  double(qz(2:end-1,:)<0) * d;
Qas = T(1:end-1,:) .* rhof(1:end-1,:) .* cf(1:end-1,:) .* qz(2:end-1,:) .* ...
  double(qz(2:end-1,:)>=0) * d;
Qa = Qan + Qas;

% total Q
Q = sum(Qcc + Qa, 2);

% Nusselt number (for each band of horizontal boxes)
Nu = Q./Qc;
