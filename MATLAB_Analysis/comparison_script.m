close all
% clear all
clc

filepath = "C:\Users\Marshall\Documents\N-Body_Simulations\" + ...
    "General_N-Body_Code\Updated_N-Body_Code_v2\Updated_N-Body_Code_v2\" + ...
    "EFIELD_DATABASE\1.9\";

timing_prefix = "Timing_Comparisons\timing_323_1ps_10000ps_";
% timing_filenames_short = ["direct","0theta_barnesHut",...
%     "100theta_barnesHut","500theta_barnesHut","1000theta_barnesHut"];
% timing_filenames = filepath + timing_prefix + timing_filenames_short + ".bin";

% timing_prefix = "Timing_Comparisons\timing_323_1ps_1000ps_";
timing_filenames_short = ["direct",...
    "3expansions_0m_fmm","3expansions_-1m_fmm",...
    "5expansions_0m_fmm","5expansions_-1m_fmm",...
    "7expansions_0m_fmm","7expansions_-1m_fmm"];
timing_filenames = filepath + timing_prefix + timing_filenames_short + ".bin";


L2norm_prefix = "L2-Norm_Comparisons\L2norm_323_1ps_10000ps_";
% L2norm_filenames_short = ["0theta_barnesHut","100theta_barnesHut",...
%     "500theta_barnesHut","1000theta_barnesHut"];
% L2norm_filenames = filepath + L2norm_prefix + L2norm_filenames_short + ".bin";

L2norm_filenames_short = ["3expansions_0m_fmm","3expansions_-1m_fmm",...
    "5expansions_0m_fmm","5expansions_-1m_fmm",...
    "7expansions_0m_fmm","7expansions_-1m_fmm"];
L2norm_filenames = filepath + L2norm_prefix + L2norm_filenames_short + ".bin";

timing_data = struct("Method",{},"Timing_Data",{});

for i = length(timing_filenames):-1:1
    if (contains(timing_filenames(i),'direct'))
        timing_data(i).Method = timing_filenames_short(i);
        timing_data(i).Timing_Data = load_direct_method_timing(timing_filenames(i));
    elseif (contains(timing_filenames(i),'barnesHut'))
        timing_data(i).Method = timing_filenames_short(i);
        timing_data(i).Timing_Data = load_barnes_hut_method_timing(timing_filenames(i));
    elseif (contains(timing_filenames(i),'fmm'))
        timing_data(i).Method = timing_filenames_short(i);
        timing_data(i).Timing_Data = load_multipole_method_timing(timing_filenames(i));
    end
end

figure
hold on
for i = 1:length(timing_data)
    plot([timing_data(i).Timing_Data.numParticles],[timing_data(i).Timing_Data.totalTime],'.');
end
legend(timing_filenames_short)
xlabel('Number of Particles')
ylabel('Time (ms)')
title('Timing Comparisons')
set(gca,'XScale','log','YScale','log')


L2norm_data = struct("Method",{},"L2Norm_Data",{});

for i = length(L2norm_filenames):-1:1
    L2norm_data(i).Method = L2norm_filenames_short(i);
    L2norm_data(i).L2Norm_Data = load_L2_norm_data(L2norm_filenames(i));
end

figure
hold on
for i = 1:length(L2norm_data)
    plot([L2norm_data(i).L2Norm_Data.numParticles],[L2norm_data(i).L2Norm_Data.L2norm],'.');
end
legend(L2norm_filenames_short)
xlabel('Number of Particles')
ylabel('L2-Norm')
title('L2-Norm Comparisons')
set(gca,'YScale','log')
