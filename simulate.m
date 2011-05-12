function foobar()
	
	tend    = 30;
	
	omega_x = [1; 0; 0];
	omega_y = [0; 1; 0];
	omega_z = [0; 0; 1];
	
	[Tx, Ox, Txs, Oxs, Stx, Sx, Ssx] = calc(omega_x, tend);
	[Ty, Oy, Tys, Oys, Sty, Sy, Ssy] = calc(omega_y, tend);
	[Tz, Oz, Tzs, Ozs, Stz, Sz, Ssz] = calc(omega_z, tend);
	
	figure(1)
	display_speeds(Tx, Ox, Stx, Sx);
	
	figure(2)
	display_speeds(Txs, Oxs, Stx, Ssx);
	
	figure(3)
	display_speeds(Ty, Oy, Sty, Sy);
	
	figure(4)
	display_speeds(Tys, Oys, Sty, Ssy);
	
	figure(5)
	display_speeds(Tz, Oz, Stz, Sz);
	
	figure(6)
	display_speeds(Tzs, Ozs, Stz, Ssz);
	
	
function [] = display_speeds(t, omega, St, S)
	clf
	hold on
	plot(t, omega);
	plot(St, S, '*');
	hold off

function [T, omega, Ts, omegas, St, S, Ss] = calc(omega0, tend)
% CALC   Calculates the simulated rotation-speeds, also returns stability
%
%	SYNTAX
%		[T, omega, Ts, omegas, St, S, Ss] = calc(omegadot, tend)
%
%   ARGUMENTS
%		omega0 = A list of initial rotation speeds
%		tend   = The time to stop simulating
%		T      = Vector contianing time-points from ode-solver
%		omega  = Vector containing the speeds corresponding to T
%		Ts     = Vector containing time-points from the ode-solver for the perturbed system
%		omegas = Vector containing the speeds corresponding to Ts
%		St     = Vector containing the time points for the stability calculations
%		S      = Vector containing the stability values corresponding to St for the initial solution
%		Ss     = Vector containing the stability values corresponding to St for the perturbed solution
%
%	Created by Martin Wernståhl on 2011-05-12.
%	Copyright (C)  Martin Wernståhl. All rights reserved.
%
	SOL = ode45(@spin, [0, tend], omega0, odeset('RelTol', 1e-6));
	SOLs = ode45(@spin, [0, tend], omega0 + 0.1 * rand(3, 1), odeset('RelTol', 1e-6));
	
	% Solution without noise
	T      = SOL.x;
	omega  = SOL.y;
	% Solution with noise
	Ts     = SOLs.x;
	omegas = SOLs.y;
	% Time steps for calculating the stability
	St     = 0:1:tend;
	
	% Calculate the stability for the system with and without noise in St
	for i = 1:length(St)
		S(i)  = maxeig(deval(SOL, St(i)));
		Ss(i) = maxeig(deval(SOLs, St(i)));
	end

function [omegadot] = spin(t, omega)
	
	gamma    = [0.81778; -0.99681; 0.96923];
	
	omegadot = [omega(2) * omega(3) * gamma(1);
	            omega(1) * omega(3) * gamma(2);
	            omega(1) * omega(2) * gamma(3)];

function [maxeig] = maxeig(omega)
	
	gamma    = [0.81778; -0.99681; 0.96923];
	
	J = diag(gamma)*[0 omega(3) omega(2);
	                 omega(3) 0 omega(1);
	                 omega(2) omega(1) 0];
	
	maxeig = max(real(eig(J)));