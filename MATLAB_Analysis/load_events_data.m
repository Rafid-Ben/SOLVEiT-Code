function Edata = load_events_data(filename)
%LOAD_EVENTS_DATA Loads the relevant data for important events occuring
%       throughout the simulation
%   Description:
%       Given an input file, this function returns the data of important
%           events throughout the simultation. The output, Edata, contains
%           an array of structs, each containing the position, velocity,
%           species, index, event, and time of occurance.
%   Input:
%       filename - name of file containing event data for simulation
%   Output:
%       Edata - array of structs containing event data for simulation

% Open file
[fileID, errmsg] = fopen(filename);
if (fileID < 0)
    disp(errmsg);
    error("Invalid file name.");
end

Edata = struct('pos',{},'vel',{},'species',{},'ind',{},'event',{},'time',{});

i = 1;
test_val = single(fread(fileID, 1, 'float'));
temp_pos = zeros(1,3);
temp_vel = zeros(1,3);
while ~feof(fileID)
    Edata(i).time = test_val;
    
    for j = 1:3
        temp_pos(j) = single(fread(fileID, 1, 'float'));
    end
    Edata(i).pos = temp_pos;
    
    for j = 1:3
        temp_vel(j) = single(fread(fileID, 1, 'float'));
    end
    Edata(i).vel = temp_vel;
    
    Edata(i).ind = uint32(fread(fileID, 1, 'uint32'));
    Edata(i).species = uint32(fread(fileID, 1, 'uint32'));
    Edata(i).event = uint32(fread(fileID, 1, 'uint32'));
    i = i + 1;
    
    test_val = single(fread(fileID, 1, 'float'));
end

fclose(fileID);
end

