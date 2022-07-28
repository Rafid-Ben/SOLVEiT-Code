function Tdata = load_direct_method_timing(filename)
%LOAD_DIRECT_METHOD_TIMING Loads the timing data for a simulation using the
%       direct force method to compute accelerations
%   Description:
%       Given an input file, this function returns the timing data of a
%           simlation using the direct force method. The output, Tdata,
%           contains an array of structs containing the number of
%           particles, acceleration time, and total time for each snapshot
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

Tdata = struct('numParticles',{},'accelerationTime',{},'totalTime',{});

iterNum = 1;

test_val = uint32(fread(fileID, 1, 'uint32'));
while ~feof(fileID)
    tempParticles = test_val;
    tempAcceleration = single(fread(fileID, 1, 'float'));
    
    Tdata(iterNum).numParticles = tempParticles;
    Tdata(iterNum).accelerationTime = tempAcceleration;
    Tdata(iterNum).totalTime = tempAcceleration;

    iterNum = iterNum + 1;
    
    test_val = uint32(fread(fileID, 1, 'uint32'));
end

fclose(fileID);
end
