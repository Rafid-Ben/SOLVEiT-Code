close all
% clear all
clc

% Alter filepath as needed to include full path to file
filepath = "EFIELD_DATABASE\1.9\";

state_filename = filepath + "state_323_1ps_2000ps_truncatedNew.bin";

events_filename = filepath + "events_323_1ps_2000ps.bin";

Sdata = load_particle_data(state_filename);

events_data = load_events_data(events_filename);

nSteps = length(Sdata);
[nBodies,~] = size(Sdata(nSteps).pos);

nShow = nBodies;
dim = 0;
p = 0;
view = 'x';
title_opts = ['nBodies = ', num2str(nBodies)];
dt = 1;

output = cell2mat(inputdlg('GIF (G) or video (V)?','Output format',1,{'V'})); 

nBody_animate_new(Sdata, events_data, nSteps, nShow, dim, p, view, title_opts, output, dt, nBodies, 1, "test")

