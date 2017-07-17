% PlyNsFlr.m	: Polynomial_and_Noise_Floor_Ft : computes error in the predictions of a function defined by a polynomial at early times and a flat noise floor at later times. 

% Written by James Traer - jtraer@mit.edu - November 2014
% based on asym_error.m written Dec 6 2012 - Josh McDermott
% edited by James Traer Jan 2014

%% Inputs:
% xx	-- time-series data
% tt	-- time values
% prms 	-- an element vector of N-model parameters to fit an (N-1)th order polynomial
% 	   The 1st element of prms is the time of the inflection point and the rest 
%	   are polynomial coefficients as in POLY

function [total_error] = PlyNsFlrFt(xx,tt,prms)

%set(0,'DefaultFigureVisible','off');

% get index of inflection point
Tinf=prms(1); prms=prms(2:end);
[~,inf_ndx]=min(abs(tt-Tinf));
% compute the modeled data
xx_mdl=polyval(prms,tt(1:inf_ndx));
xx_mdl=[xx_mdl xx_mdl(end)*ones(1,length(tt)-length(xx_mdl))];
%--- this is unstable -----xx_mdl=10.^(xx_mdl/20);
% compute error
total_error = rms(xx-xx_mdl);
