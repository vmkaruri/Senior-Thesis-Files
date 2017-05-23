% DIMELORENZ.M
%    Demonstration of the program  DIME for calculation
%    of the "corelation dimension" of a dataset.
%    Computes and plots the dimension of the classic
%    Lorenz attractor by two different ways

% Kirill Pankratov, Feb. 22, 1994

echo on
%   This example illustrate the calculation of the correlation dimension
% of chaotic timeseries, obtained from simulating the Lorenz attractor

% Enter parameters for Lorenz equation:
global SIGMA RHO BETA
SIGMA = 10.;
RHO = 28.;
BETA = 8./3.;

% Enter initial conditions and integration time:
y0 = [33.3 -12.5 -11.4]';
tfin = 15;

% Now compute the time series of the phase trajectory of the Lorenz equation

  [t,y] = ode23('lorenzeq',0,tfin,y0);

% Now calculate the correlation dimension of obtained timeseries.
% We shall use two different ways: 
%   1) Take 1 component of the timeseries  y and put it succesively into
%      different "embedding" dimensions.
%   1) Compute the dimension of the manifold  y in actual 3-dimensional
%      phase space.

%  (Please be patient, this will take a few minutes)

  dime(y(:,1),5);   % 1-d timeseries, "embedded" into 2,3,4,5-dimensional space
  drawnow

  dime(y);          % full 3-d series, actual phase space

%  Now look at two different plots:
%    What is the correlation dimension of the Lorenz attractor?
%    Are the results much different for two ways of calculations?

echo off
