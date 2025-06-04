%% Put together all the data files
clear; clc; close all;

n_ves = 3;
n_par = 7;%%8;
N_total = 100;
p_PCE = zeros(32,n_ves,N_total);
q_PCE = zeros(32,n_ves,N_total);
A_PCE = zeros(32,n_ves,N_total);

load testdata_WK_A.mat
load testdata_WK_p.mat
load testdata_WK_q.mat
counter = 1;
leaveout = [];
for i=1:N_total
    if ~isempty(p_results{i})
        p_PCE(:,:,counter) = p_results{i}(1:16:end,4:6);
        q_PCE(:,:,counter) = q_results{i}(1:16:end,4:6);
        A_PCE(:,:,counter) = A_results{i}(1:16:end,4:6);
        counter = counter+1;
    else
        leaveout(end+1) = i;
    end

end
%% 
p_PCE = p_PCE(:,:,1:counter-1); 
q_PCE = q_PCE(:,:,1:counter-1); 
A_PCE = A_PCE(:,:,1:counter-1); 

load testdata_WK_pars low param_sample upp
param_sample(leaveout,:) = [];
upp(1) = exp(upp(1));
low(1) = exp(low(1));

save('pq_test','p_PCE','q_PCE','A_PCE','param_sample','upp','low')
