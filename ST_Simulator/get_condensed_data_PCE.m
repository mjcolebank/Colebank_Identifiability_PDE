%% Put together all the data files
clear; clc; close all;

n_ves = 3;
n_par = 11;%%8;
N_total = 100;
p_PCE = zeros(32,n_ves,N_total);
q_PCE = zeros(32,n_ves,N_total);
A_PCE = zeros(32,n_ves,N_total);

load new_test_data_A.mat
load new_test_data_p.mat
load new_test_data_q.mat

counter = 1;
leaveout = [];
for i=1:N_total
    if ~isempty(p_results{i})
        if length(p_results{i})==1024
        p_PCE(:,:,counter) = p_results{i}(1:32:end,4:6);
        q_PCE(:,:,counter) = q_results{i}(1:32:end,4:6);
        A_PCE(:,:,counter) = A_results{i}(1:32:end,4:6);
        
        else
            p_PCE(:,:,counter) = p_results{i}(1025:32:end,4:6);
        q_PCE(:,:,counter) = q_results{i}(1025:32:end,4:6);
        A_PCE(:,:,counter) = A_results{i}(1025:32:end,4:6);
        end
        counter = counter+1;
    else
        leaveout(end+1) = i;
    end

end


%% 
p_PCE = p_PCE(:,:,1:counter-1); 
q_PCE = q_PCE(:,:,1:counter-1); 
A_PCE = A_PCE(:,:,1:counter-1); 

load new_test_data_pars.mat LB param_sample UB
low = LB;
upp = UB;
param_sample(leaveout,:) = [];
param_sample = param_sample(:,1:7);

%% FIltering step
% max_p = max(p_PCE(:,1,:));
% min_p = min(p_PCE(:,1,:));
% pulse_P = max_p-min_p;
% id_pulse = find(pulse_P<6);
% id_max = find(max_p>100);
% id_remove = unique([id_pulse; id_max]);
% p_PCE(:,:,id_remove) = [];
% q_PCE(:,:,id_remove) = [];
% A_PCE(:,:,id_remove) = [];
% param_sample(id_remove,:) = [];
%%

% upp(1) = exp(upp(1));
% low(1) = exp(low(1));

save('pq_ST_NEW_testdata','p_PCE','q_PCE','A_PCE','param_sample','upp','low')
