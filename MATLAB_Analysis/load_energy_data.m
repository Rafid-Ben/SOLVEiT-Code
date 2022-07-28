function Edata = load_energy_data(filename)
%LOAD_ENERGY_DATA Loads the energy data throughout the simulation
%   Description:
%       Given an input file, this function returns the energy conservation
%           data throughout the simultation. The output, Edata, contains an
%           array of structs, each containing the current kinetic,
%           potential, and field energy, as well as the kinetic, potential,
%           and field energy from injecting particles
%   Input:
%       filename - name of file containing energy data for simulation
%   Output:
%       Edata - array of structs containing energy data for simulation

% Open file
[fileID, errmsg] = fopen(filename);
if (fileID < 0)
    disp(errmsg);
    error("Invalid file name.");
end

Edata = struct('KE_inj',{},'PE_inj',{},'FE_inj',{},'KE',{},'PE',{},'FE',{});

iterNum = 1;

test_val = single(fread(fileID, 1, 'float'));
while ~feof(fileID)
    tempKE_inj = test_val;
    tempPE_inj = single(fread(fileID, 1, 'float'));
    tempFE_inj = single(fread(fileID, 1, 'float'));
    tempKE = single(fread(fileID, 1, 'float'));
    tempPE = single(fread(fileID, 1, 'float'));
    tempFE = single(fread(fileID, 1, 'float'));
    
    Edata(iterNum).KE_inj = tempKE_inj;
    Edata(iterNum).PE_inj = tempPE_inj;
    Edata(iterNum).FE_inj = tempFE_inj;
    Edata(iterNum).KE = tempKE;
    Edata(iterNum).PE = tempPE;
    Edata(iterNum).FE = tempFE;

    iterNum = iterNum + 1;
    
    test_val = single(fread(fileID, 1, 'float'));
end

fclose(fileID);
end

