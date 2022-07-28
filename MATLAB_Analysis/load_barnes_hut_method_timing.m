function Tdata = load_barnes_hut_method_timing(filename)
%LOAD_BARNES_HUT_METHOD_TIMING Loads the timing data for a simulation using
%       the Barnes-Hut force method to compute accelerations
%   Description:
%       Given an input file, this function returns the timing data of a
%           simlation using the Barnes-Hut force method. The output, Tdata,
%           contains an array of structs containing the number of
%           particles, precalc time, acceleration time, and total time for
%           each snapshot
%   Input:
%       filename - name of file containing timing data for simulation
%   Output:
%       Tdata - array of structs containing timing data for simulation

% Open file
[fileID, errmsg] = fopen(filename);
if (fileID < 0)
    disp(errmsg);
    error("Invalid file name.");
end

Tdata = struct('numParticles',{},'accelerationTime',{},'precalcTime',{},'totalTime',{});

iterNum = 1;

test_val = uint32(fread(fileID, 1, 'uint32'));
while ~feof(fileID)
    tempParticles = test_val;
    tempPrecalc = single(fread(fileID, 1, 'float'));
    tempAcceleration = single(fread(fileID, 1, 'float'));
    
    Tdata(iterNum).numParticles = tempParticles;
    Tdata(iterNum).precalcTime = tempPrecalc;
    Tdata(iterNum).accelerationTime = tempAcceleration;
    Tdata(iterNum).totalTime = tempPrecalc + tempAcceleration;
    
    iterNum = iterNum + 1;
    
    test_val = uint32(fread(fileID, 1, 'uint32'));
end

fclose(fileID);
end

