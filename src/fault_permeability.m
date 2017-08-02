function [ kx, kz] = fault_permeability(input)
% [ kx, kz] = fault_permeability(input)
%   Calculates the permeability field including a fault
%   of permeability k_fault surrounded by background permeability k_bg
%   fault crosses top of model at fault_xs, has dip fault_dip (degrees)
%   and width fault_width, extends down to max depth of fault_depth

% initialize permeability array with background permeability
kx = ones(input.nz, input.nx)*input.kx;
kz = ones(input.nz, input.nx)*input.kz;

% construct an array of "distance to the fault line"
fault_distance = abs( tand(input.fault_dip).*input.X  + input.Z - ...
  tand(input.fault_dip)*input.fault_xs) ./ sqrt(1+tand(input.fault_dip).^2);

% set permeability wthiin fault
kx(fault_distance<=input.fault_width/2 & input.Z<=input.fault_depth) = input.kx_fault;
kz(fault_distance<=input.fault_width/2 & input.Z<=input.fault_depth) = input.kz_fault;
