function [Qout] = heatflux(d,T,rhof,cf,qz,cracked)

% conductive heat transport across each boundary (with convection)
Qcc = -lamdam.*(T(2:end,:)-T(1:end-1,:))./d;

% advective heat transport
Qan = T(2:end,:).*rhof(2:end,:).*cf(2:end,:).*qz(2:end-1,:).*double(qz(2:end-1,:)<0);
Qas = T(1:end-1,:).*rhof(1:end-1,:).*cf(1:end-1,:).*qz(2:end-1,:).*double(qz(2:end-1,:)>=0);
Qa = (Qan + Qas);

% total Q
Q = Qcc + Qa(cracked(2:end,:)==1);

%cracked_boundary = double(~any(cat(3, cracked(1:end-1,:), cracked(2:end,:)), 3));

%Qa = (Qan + Qas).*cracked_boundary;

Qout=mean(Q(cracked(2:end,:)==1));
