function Tdata = load_multipole_method_timing(filename)
%LOAD_MULTIPOLE_METHOD_TIMING Loads the timing data for a simulation using
%       the direct force method to compute accelerations
%   Description:
%       Given an input file, this function returns the timing data of a
%           simlation using the multipole method. The output, Tdata,
%           contains an array of structs containing the number of
%           particles, precalc time, acceleration time, total time for each
%           snapshot, and a breakdown of timing for different steps in the
%           algorithm, as well as relevant parameters impacting the
%           performance
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

Tdata = struct('numParticles',{},'accelerationTime',{},'precalcTime',{},'totalTime',{},...
    'interactionListTime',{},'directTime',{},'p2p_time',{},'p2m_time',{},...
    'm2m_time',{},'m2l_time',{},'l2l_time',{},'l2p_time',{},'m2p_time',{},...
    'sm2m_time',{},'sm2l_time',{},'sl2l_time',{},'sm2p_time',{},'sp2p_time',{},...
    'sp2m_time',{},'sl2p_time',{},'outlier_p2p_time',{},'p2outlier_p_time',{},...
    'sp2outlier_p_time',{},'m2outlier_p_time',{},'sm2outlier_p_time',{},...
    'xdims',{},'ydims',{},'zdims',{},'maxLevel',{},'numOutliers',{});

iterNum = 1;

test_val = uint32(fread(fileID, 1, 'uint32'));
while ~feof(fileID)
    tempParticles = test_val;
    tempNumOutliers = uint32(fread(fileID, 1, 'uint32'));
    tempxdims = uint32(fread(fileID, 1, 'uint32'));
    tempydims = uint32(fread(fileID, 1, 'uint32'));
    tempzdims = uint32(fread(fileID, 1, 'uint32'));
    tempMaxLevel = uint32(fread(fileID, 1, 'uint32'));

    tempPrecalc = single(fread(fileID, 1, 'float'));
    tempInteractionListTime = single(fread(fileID, 1, 'float'));
    tempDirect = single(fread(fileID, 1, 'float'));
    temp_p2p = single(fread(fileID, 1, 'float'));
    temp_p2m = single(fread(fileID, 1, 'float'));
    temp_m2m = single(fread(fileID, 1, 'float'));
    temp_m2l = single(fread(fileID, 1, 'float'));
    temp_l2l = single(fread(fileID, 1, 'float'));
    temp_l2p = single(fread(fileID, 1, 'float'));
    temp_m2p = single(fread(fileID, 1, 'float'));
    temp_sm2m = single(fread(fileID, 1, 'float'));
    temp_sm2l = single(fread(fileID, 1, 'float'));
    temp_sl2l = single(fread(fileID, 1, 'float'));
    temp_sm2p = single(fread(fileID, 1, 'float'));
    temp_sp2p = single(fread(fileID, 1, 'float'));
    temp_sp2m = single(fread(fileID, 1, 'float'));
    temp_sl2p = single(fread(fileID, 1, 'float'));
    temp_outlier_p2p = single(fread(fileID, 1, 'float'));
    temp_p2outlier_p = single(fread(fileID, 1, 'float'));
    temp_sp2outlier_p = single(fread(fileID, 1, 'float'));
    temp_m2outlier_p = single(fread(fileID, 1, 'float'));
    temp_sm2outlier_p = single(fread(fileID, 1, 'float'));
    tempTotalTime = single(fread(fileID, 1, 'float'));
    
    Tdata(iterNum).numParticles = tempParticles;
    Tdata(iterNum).xdims = tempxdims;
    Tdata(iterNum).ydims = tempydims;
    Tdata(iterNum).zdims = tempzdims;
    Tdata(iterNum).maxLevel = tempMaxLevel;
    Tdata(iterNum).numOutliers = tempNumOutliers;
    
    Tdata(iterNum).precalcTime = tempPrecalc;
    Tdata(iterNum).directTime = tempDirect;
    Tdata(iterNum).interactionListTime = tempInteractionListTime;
    Tdata(iterNum).p2p_time = temp_p2p;
    Tdata(iterNum).p2m_time = temp_p2m;
    Tdata(iterNum).m2m_time = temp_m2m;
    Tdata(iterNum).m2l_time = temp_m2l;
    Tdata(iterNum).l2l_time = temp_l2l;
    Tdata(iterNum).l2p_time = temp_l2p;
    Tdata(iterNum).m2p_time = temp_m2p;
    Tdata(iterNum).sm2m_time = temp_sm2m;
    Tdata(iterNum).sm2l_time = temp_sm2l;
    Tdata(iterNum).sl2l_time = temp_sl2l;
    Tdata(iterNum).sm2p_time = temp_sm2p;
    Tdata(iterNum).sp2p_time = temp_sp2p;
    Tdata(iterNum).sp2m_time = temp_sp2m;
    Tdata(iterNum).sl2p_time = temp_sl2p;
    Tdata(iterNum).outlier_p2p_time = temp_outlier_p2p;
    Tdata(iterNum).p2outlier_p_time = temp_p2outlier_p;
    Tdata(iterNum).sp2outlier_p_time = temp_sp2outlier_p;
    Tdata(iterNum).m2outlier_p_time = temp_m2outlier_p;
    Tdata(iterNum).sm2outlier_p_time = temp_sm2outlier_p;
    
    Tdata(iterNum).accelerationTime = tempDirect + temp_p2p + temp_p2m + temp_m2m + ...
        temp_m2l + temp_l2l + temp_l2p + temp_m2p + temp_sm2m + temp_sm2l + temp_sl2l + ...
        temp_sm2p + temp_sp2p + temp_sp2m + temp_sl2p + temp_outlier_p2p + ...
        temp_p2outlier_p + temp_sp2outlier_p + temp_m2outlier_p + temp_sm2outlier_p;
    Tdata(iterNum).totalTime = tempTotalTime;
    
    iterNum = iterNum + 1;
    
    test_val = uint32(fread(fileID, 1, 'uint32'));
end

fclose(fileID);
end

