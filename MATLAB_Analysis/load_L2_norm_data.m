function Ldata = load_L2_norm_data(filename)
%LOAD_PARTICLE_DATA Loads the L2-norm data for a given simulation
%   Description:
%       Given an input file, this function returns the L2-norm data of a
%           simulation across a sequence of timesteps. The output, Ldata,
%           contains an array of structs containing the number of
%           particles, the L2-norm of the system, the maximum L2-norm
%           value, and the maximum L2-norm particle index.
%   Input:
%       filename - name of file containing L2-norm data for simulation
%   Output:
%       Ldata - array of structs containing L2-norm data for simulation

% Open file
[fileID, errmsg] = fopen(filename);
if (fileID < 0)
    disp(errmsg);
    error("Invalid file name.");
end

Ldata = struct('numParticles',{},'L2norm',{},'maxL2norm',{},'maxL2normIndex',{});

iterNum = 1;

test_val = uint32(fread(fileID, 1, 'uint32'));
while ~feof(fileID)
    tempParticles = test_val;
    
    tempL2norm = single(fread(fileID, 1, 'float'));
    tempMaxL2norm = single(fread(fileID, 1, 'float'));

    tempMaxL2normIndex = uint32(fread(fileID, 1, 'uint32'));
    
    Ldata(iterNum).numParticles = tempParticles;
    Ldata(iterNum).L2norm = tempL2norm;
    Ldata(iterNum).maxL2norm = tempMaxL2norm;
    Ldata(iterNum).maxL2normIndex = tempMaxL2normIndex;
    
    iterNum = iterNum + 1;
    
    test_val = uint32(fread(fileID, 1, 'uint32'));
end

fclose(fileID);
end

