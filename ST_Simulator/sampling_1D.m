
%% Define the parameters
% LB = [log(1e6) -30 log(3e4) 0.8 0.6   5 0.001 0.8 0.6   5 0.001];
% UB = [log(1e8) -10 log(5e5) 0.9 0.74 50 0.01  0.9 0.74 50 0.01];

% Updated 1/13/25
% LB = [log(1e6) -30 log(8e4) 0.84 0.64   5 0.001 0.84 0.64   5 0.001];
% UB = [log(1e8) -10 log(5e5) 0.9 0.74   50 0.01  0.9 0.74 50 0.01];


% Updated 1/22/25
LB = [log(1e6) -30 log(1e5) 0.75 0.60   15 0.001 0.75 0.60   15 0.001];
UB = [log(1e8) -10 log(8e5) 0.88 0.72   60 0.01  0.85 0.74   60 0.01];
% n_samp = 3000;
n_samp = 100;%
X = lhsdesign(n_samp,11);
param_sample = LB + X.*(UB-LB);
param_sample(:,1) = exp(param_sample(:,1));
param_sample(:,3) = exp(param_sample(:,3));
param_sample(:,8:11) = param_sample(:,4:7);


p_results = cell(n_samp,1);
q_results = cell(n_samp,1);
A_results = cell(n_samp,1);

% parpool(4);
parfor i=1:n_samp
    pars = [param_sample(i,:),i];
    pars_str = mat2str(pars);
    %% call the model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    out = unix(sprintf('sor06.exe  %s',pars_str(2:end-1)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 
    disp(i);
    if out~=1
        fname = strcat('output_',num2str(i),'.2d');
        data = load(fname);
        [time,x,p,q,A,C] = gnuplot(data);
        p_results{i} = p;
        q_results{i} = q;
        A_results{i} = A;
        
        delete(fname);

        % figure; plot(p);
    end
    
end
% save('simple_model_ST_1_25_p','p_results');
% save('simple_model_ST_1_25_q','q_results');
% save('simple_model_ST_1_25_A','A_results');
% save('simple_model_ST_1_25_pars','param_sample','LB','UB');

save('new_test_data_p','p_results');
save('new_test_data_q','q_results');
save('new_test_data_A','A_results');
save('new_test_data_pars','param_sample','LB','UB');
