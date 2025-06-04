% Model usines a simple linear elastic wall model and a three element Windkessl outflow
% boundary condition. The system is driven by imposing an inflow profile.

clc; clear all; close all;
format shortg;
%% Define the parameters
upp = [log(1e8) -10 log(5e5) 10 10 20 20 1.2 1.2];
low = [log(1e6) -30 log(3e4) 0.5 0.5 0.8 0.8 1e-5 1e-5];

% n_samp = 2500;
n_samp = 100;
X = lhsdesign(n_samp,9);
param_sample = low + X.*(upp-low);
param_sample(:,[1 3]) = exp(param_sample(:,[1 3]));

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
    end
    
end
save('testdata_WK_p','p_results');
save('testdata_WK_q','q_results');
save('testdata_WK_A','A_results');
save('testdata_WK_pars','param_sample','low','upp');


