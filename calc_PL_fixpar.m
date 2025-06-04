% New (simpler?) Profile Likelihood code

function [J_save,par_save,param_global_all] = calc_PL_fixpar(f,data,params,bounds,err_var,par_fix)
opt_options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt',...
    'Display','none','UseParallel',true,'MaxIterations',300,'FiniteDifferenceStepSize',1e-4);
% Fix one parameter, infer the others
num_par = length(params);
param0 = params;
N = 201;
Nhalf = (N-1)/2 + 1;
par_ids = 1:num_par;
par_ids(par_fix) = [];
% Upper and lower bounds for OPTIMIZATION
lower = bounds(:,1);
upper = bounds(:,2);


%% Run an optimization with all parameters to find the "global" minimum
% Run all parameters through optimization
% opt_options.UseParallel = false;
% tic
% [xall,Jall,res,~,~,~,~] = lsqnonlin(@(q) eval_PL(q, [],[],data,num_par,f,err_var,1:num_par,param0), param0, lower,upper, opt_options);
% toc
% opt_options.UseParallel = true;
% tic
[xall,Jall,res,~,~,~,~] = lsqnonlin(@(q) eval_PL(q, [],[],data,num_par,f,err_var,1:num_par,param0), param0, lower,upper, opt_options);
% toc
% Fix the fixed parameters at the true value
% [xall,Jall,res,~,~,~,~] = lsqnonlin(@(q) eval_PL(q, [],[],data,num_par,f,err_var,par_ids,param0), param0, lower,upper, opt_options);

param_global_all = xall;
params = xall(par_ids);

% Upper and lower bounds for OPTIMIZATION
lower = bounds(par_ids,1);
upper = bounds(par_ids,2);
%%
% %Range for PROFILING a single parameter
% temp_max = min(exp(1.2.*log(params)),upper');
% temp_min = max(exp(0.8.*log(params)),lower');
temp_max = min(1.5.*params,upper');
temp_min = max(0.5.*params,lower');
PL_max = max(temp_max,temp_min);
PL_min = min(temp_max,temp_min);
par_step = (PL_max-PL_min)./(N-1);

par_profile = 1:length(params);
num_profile = length(par_profile);
par_save = zeros(num_profile,num_profile,N);
J_save   = zeros(num_profile,N);
% Assume N points left and right of the values in params
for par_i=1:num_profile

    % Start at middle value, then move right
    par_space = PL_min(par_i):par_step(par_i):PL_max(par_i);
    qest      = params(par_profile~=par_i);
    lower_i   = lower(par_profile~=par_i);
    upper_i   = upper(par_profile~=par_i);
    qfix      = par_space(Nhalf);

    % Midpoint evaluation
    try
        [x,J,res,~,~,~,~] = lsqnonlin(@(q) eval_PL(q, qfix,par_i,data,num_par,f,err_var,par_ids,param_global_all), qest, lower_i,upper_i, opt_options);
        par_save(par_i,par_profile~=par_i,Nhalf) = x;
        par_save(par_i,par_i,Nhalf) = qfix;
        J_save(par_i,Nhalf)   = J;
        qest_half = x;
    catch
        par_save(par_i,par_profile~=par_i,Nhalf) = nan;
        par_save(par_i,par_i,Nhalf) = nan;
        J_save(par_i,Nhalf)   = nan;
    end
    % figure(100+par_i); clf; hold on;
    % plot(Nhalf,J_save(par_i,Nhalf),'*');
    % Values greater than half
    qest = qest_half;
    for N_i = Nhalf+1:N
        disp([par_i N_i])
        qfix = par_space(N_i);
        try
            [x,J,res] = lsqnonlin(@(q) eval_PL(q, qfix,par_i,data,num_par,f,err_var,par_ids,param_global_all), qest, lower_i,upper_i, opt_options);
            par_save(par_i,par_profile~=par_i,N_i) = x;
            par_save(par_i,par_i,N_i) = qfix;
            J_save(par_i,N_i)   = J;
            qest = x;

            % plot(N_i,J_save(par_i,N_i),'*');
        catch
            par_save(par_i,par_profile~=par_i,N_i) = nan;
            par_save(par_i,par_i,N_i) = nan;
            J_save(par_i,N_i)   = nan;
        end
    end

    % Values less than half
    qest = qest_half;
    for N_i = (Nhalf-1):-1:1
        disp([par_i N_i]);
        qfix = par_space(N_i);

        try
            [x,J,res] = lsqnonlin(@(q) eval_PL(q, qfix,par_i,data,num_par,f,err_var,par_ids,param_global_all), qest, lower_i,upper_i, opt_options);
            par_save(par_i,par_profile~=par_i,N_i) = x;
            par_save(par_i,par_i,N_i) = qfix;
            J_save(par_i,N_i)   = J;
            qest = x;


            % plot(N_i,J_save(par_i,N_i),'*');
        catch
            par_save(par_i,par_ids~=par_i,N_i) = nan;
            par_save(par_i,par_i,N_i) = nan;
            J_save(par_i,N_i)   = nan;
        end
    end

end

end

function out = eval_PL(qest,qfix,par_i,data,num_par,f,err_var,par_ids,param_global_all)
q_eval = param_global_all;%zeros(num_par,1);

if isempty(par_i)
    q_eval = qest;
else
    q_eval(par_ids(par_ids~=par_ids(par_i))) = qest;
    q_eval(par_ids(par_i)) = qfix;
end
model_out = feval(f,q_eval(:));

out = data(:) - model_out(:);
% out = sum((data(:) - model_out(:)).^2);
out = out./sqrt(err_var);
end

function out = run_model(q)
% PUT PCE HERE
x = linspace(0,1,20);
out = q(1) + q(2).*x + q(3).*x.^2;
end