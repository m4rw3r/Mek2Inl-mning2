function foobar()
	
	% Time and initial speeds
	tend   = 30;
	omega0 = [1 0 0;
	          0 1 0;
	          0 0 1];
	
	for i = 1:3
		[T, O, Ts, Os, St, S, Ss] = calc(omega0(:, i), tend);
		
		h = figure(i);
		display_speeds(T, O, St, S);
		saveas(h, sprintf('simulated/%d_std', i), 'epsc');
		
		h = figure(i+1);
		display_speeds(Ts, Os, St, Ss);
		saveas(h, sprintf('simulated/%d_per', i), 'epsc');
	end
	
function [] = display_speeds(t, omega, St, S)
	clf
	hold on
	plot(t, omega);
	plot(St, S, '*');
	hold off

function [T, omega, Ts, omegas, St, S, Ss] = calc(omega0, tend)
% CALC   Calculates the simulated rotation-speeds,
%        also returns stability
%
%	SYNTAX
%		[T, omega, Ts, omegas, St, S, Ss] = calc(omegadot, tend)
%
%   ARGUMENTS
%		omega0 = A list of initial rotation speeds
%		tend   = The time to stop simulating
%		T      = Vector contianing time-points from ode-solver
%		omega  = Vector containing the speeds corresponding to T
%		Ts     = Vector containing time-points from the ode-solver
%                for the perturbed system
%		omegas = Vector containing the speeds corresponding to Ts
%		St     = Vector containing the time points for the stability
%                calculations
%		S      = Vector containing the stability values corresponding
%                to St for the initial solution
%		Ss     = Vector containing the stability values corresponding
%                to St for the perturbed solution
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
	
	g        = [0.81778; -0.99681; 0.96923];
	
	omegadot = [omega(2) * omega(3) * g(1);
	            omega(1) * omega(3) * g(2);
	            omega(1) * omega(2) * g(3)];

function [maxeig] = maxeig(omega)
	
	g = [0.81778; -0.99681; 0.96923];
	
	J = diag(g)*[0 omega(3) omega(2);
	             omega(3) 0 omega(1);
	             omega(2) omega(1) 0];
	
	maxeig = max(real(eig(J)));