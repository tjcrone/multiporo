function K = drainedbulk(Vpu,Vsu,rhosat,Kf,Ks,phi)
%This function uses the velocity form of Gassmann's Relations
%to compute the drained bulk modulus K from the seismic 
%velocities Vpu and Vsu, the fluid saturated rock density rhosat,
%the bulk modulus of the rock grains Ks, and the bulk modulus
%of the fluid Kf.  See page 172 in "The Rock Physics Handbook"
%by Gary Mavko for more information on this form of the Gassmann
%equation.  The inputs may be scalars, vectors, or 2-dimensional
%or 3-dimensional matrices, but they all must be the same size.
%The units of the inputs should be SI, with the bulk moduli in Pa,
%the velocities in m/s, and the density in kg/m^3;


Gu = rhosat.*Vsu.^2; %Shear modulus from density and Vsu

K = Ks.*(3.*Vpu.^2.*Gu.*Ks.*phi + 3.*Vpu.^2.*Gu.*Kf - ...
   3.*Vpu.^2.*Gu.*Kf.*phi - 3.*Kf.*Vsu.^2.*Ks - ...
   4.*Vsu.^2.*Gu.*Ks.*phi - 4.*Vsu.^2.*Gu.*Kf + ...
   4.*Vsu.^2.*Gu.*Kf.*phi) ./ (3.*Vpu.^2.*Gu.*Kf - ...
   3.*Kf.*Vsu.^2.*Ks + 3.*Vsu.^2.*phi.*Ks.^2 - ...
   3.*Vsu.^2.*Kf.*Ks.*phi - 4.*Vsu.^2.*Gu.*Kf);
