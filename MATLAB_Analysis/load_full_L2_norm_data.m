function Ldata = load_full_L2_norm_data(filename)
%LOAD_PARTICLE_DATA Loads the full L2-norm data for a given simulation
%   Description:
%       Given an input file, this function returns the full L2-norm data of
%           a simulation across a sequence of timesteps. The output, Ldata,
%           contains an array of structs containing the number of
%           particles and the L2-norm for each particle in every snapshot
%   Input:
%       filename - name of file containing full L2-norm data for simulation
%   Output:
%       Ldata - array of structs containing full L2-norm data for simulation

% Open file
[fileID, errmsg] = fopen(filename);
if (fileID < 0)
    disp(errmsg);
    error("Invalid file name.");
end

Ldata = struct('numParticles',{},'L2norm',{});

tempParticles = uint32(fread(fileID, 1, 'uint32'));
iterNum = 1;

while ~feof(fileID) && ~isempty(tempParticles)
    Ldata(iterNum).numParticles = tempParticles;
    Ldata(iterNum).L2norm = zeros(1, tempParticles);
    for i = 1:tempParticles
        tempL2norm = single(fread(fileID, 1, 'float'));
        Ldata(iterNum).L2norm(i) = tempL2norm;
    end
    
    tempParticles = uint32(fread(fileID, 1, 'uint32'));
    
    iterNum = iterNum + 1;
end

fclose(fileID);
end

