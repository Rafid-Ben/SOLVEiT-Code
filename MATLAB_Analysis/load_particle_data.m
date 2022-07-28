function Sdata = load_particle_data(filename)
%LOAD_PARTICLE_DATA Loads the state data of a system of particles
%   Description:
%       Given an input file, this function returns the state data of a
%           system of particles across a sequence of timesteps. The output,
%           Sdata, contains an array of structs containing the position,
%           velocity, Laplace electric field, Poisson electric field,
%           species, and index of each particle, as well as the time at
%           which each snapshot is captured.
%   Input:
%       filename - name of file containing state data of particles
%   Output:
%       Sdata - array of structs containing state data of particles

% Open file
[fileID, errmsg] = fopen(filename);
if (fileID < 0)
    disp(errmsg);
    error("Invalid file name.");
end

Sdata = struct('pos',{},'vel',{},'Efield1',{},'Efield2',{},'species',{},'ind',{},'time',{});

% Constant indicator used to separate data between snapshots
nextIterCheck = intmax('uint32');

snapshot_num = 0;

test_val = uint32(fread(fileID, 1, 'uint32'));
while ~feof(fileID)
    if (test_val == nextIterCheck)
        % Snapshot completed, move to new snapshot
        if (snapshot_num > 0)
            single_vec = single_vec(1:i,:);
            int_vec = int_vec(1:i,:);
            
            Sdata(snapshot_num).pos = single_vec(:, 1:3);
            Sdata(snapshot_num).vel = single_vec(:, 4:6);
            Sdata(snapshot_num).Efield1 = single_vec(:, 7:9);
            Sdata(snapshot_num).Efield2 = single_vec(:, 10:12);
            Sdata(snapshot_num).ind = int_vec(:, 1);
            Sdata(snapshot_num).species = int_vec(:, 2);
        end
        snapshot_num = snapshot_num + 1;
        i = 1;
        numParticlesEstimate = uint32(fread(fileID, 1, 'uint32'));
        single_vec = single(zeros(numParticlesEstimate, 12));
        int_vec = uint32(zeros(numParticlesEstimate, 2));
        Sdata(snapshot_num).time = single(fread(fileID, 1, 'float'));
        
        int_vec(i, 1) = uint32(fread(fileID, 1, 'uint32'));
        int_vec(i, 2) = uint32(fread(fileID, 1, 'uint32'));
    else
        % Snapshot not completed, continue reading
        i = i + 1;
        int_vec(i, 1) = test_val;
        int_vec(i, 2) = uint32(fread(fileID, 1, 'uint32'));
    end
    
    % Read in float data
    for j = 1:12
        single_vec(i, j) = fread(fileID, 1, 'float');
    end
    
    test_val = uint32(fread(fileID, 1, 'uint32'));
end

single_vec = single_vec(1:i,:);
int_vec = int_vec(1:i,:);

Sdata(snapshot_num).pos = single_vec(:, 1:3);
Sdata(snapshot_num).vel = single_vec(:, 4:6);
Sdata(snapshot_num).Efield1 = single_vec(:, 7:9);
Sdata(snapshot_num).Efield2 = single_vec(:, 10:12);
Sdata(snapshot_num).ind = int_vec(:, 1);
Sdata(snapshot_num).species = int_vec(:, 2);

fclose(fileID);
end
