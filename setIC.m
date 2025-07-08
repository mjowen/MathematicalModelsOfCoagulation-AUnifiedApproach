function [c0,odIdx] = setIC(IC)
% Generates a full initial condition vector from a set of donor specific
% initial conditions (IC)
% The initial conditions IC are inputted in the following order
% TF, II, V, VII, VIII, IX, X, XI, AT, TFPI, Fbg
c0 = zeros(79,1);
c0([1,8,10,2,12,16,6,14,23,20,70]) = IC;
% Set the other (non-donor specific) initial conditions
c0([54,58,60,62,64]) = [4e-5, 9.75e-7, 3e-6, 1.7e-6, 4.6e-10]; %a1AT, a2AP, a2M, C1-inh, PAI1
c0(3) = c0(2)/100; % VIIa
c0(51) = 5e-3; %Substrate - remove if measuring thrombin
odIdx = 52;
