% Run the Profile likelihood with a PCE surrogate
clear; clc; close all;
tic
%%
Names = {'$f1_{A}$','$f2_{A}$','$f3_{A}$',...
    '$\alpha$','$\beta$','$\ell_{rr}$','$r_{min}$'};
uqlab
% load('PCE_ST_PCA_5comp_NEW_deg51.mat','PCE_P','PCE_A','coeff_P','coeff_A','mu_P','mu_A');
load('PCE_ST_PCA_5comp_nofilter_51.mat','PCE_P','PCE_A','coeff_P','coeff_A','mu_P','mu_A');

pce_p1 = PCE_P;
pce_a1 = PCE_A;
pca_coeff_p1 = coeff_P;
pca_coeff_a1 = coeff_A;
pca_mu_p1    = mu_P;
pca_mu_a1    = mu_A;


load('PCE_ST_PCA_5comp_nofilter_52.mat','PCE_Q','PCE_A','coeff_Q','coeff_A','mu_Q','mu_A');
pce_q2 = PCE_Q;
pce_a2 = PCE_A;
pca_coeff_q2 = coeff_Q;
pca_coeff_a2 = coeff_A;
pca_mu_q2    = mu_Q;
pca_mu_a2    = mu_A;

load('PCE_ST_PCA_5comp_nofilter_53.mat','PCE_Q','PCE_A','coeff_Q','coeff_A','mu_Q','mu_A');
pce_q3 = PCE_Q;
pce_a3 = PCE_A;
pca_coeff_q3 = coeff_Q;
pca_coeff_a3 = coeff_A;
pca_mu_q3    = mu_Q;
pca_mu_a3    = mu_A;

% load simple_model_ST_1_25_pars.mat LB UB
load pq_ST_NEW_testdata.mat
p1 = squeeze(p_PCE(:,1,:));
a1 = squeeze(A_PCE(:,1,:));
q2 = squeeze(q_PCE(:,2,:));
a2 = squeeze(A_PCE(:,2,:));
q3 = squeeze(q_PCE(:,3,:));
a3 = squeeze(A_PCE(:,3,:));
param_test = param_sample;

upper = upp(1:7);%UB(1:7);
lower = low(1:7);%LB(1:7);
param_test(:,[1 3]) = log(param_test(:,[1 3]));

% param_sample(:,end) = param_sample(:,end)./1e-3;
% upper(end) = upper(end)./1e-3;
% lower(end) = lower(end)./1e-3;
param_fix_cell = {[],[1 2],[1 2 5], [1 2 7] , [1 2 5 7]};
fname_cell = {'','_fixf1f2','_fixf1f2beta','_fixf1f2rmin','_fixf1f2betarmin'};
for whichfix = 1:5
    for which_test_data = [47 33 3 43];%[31 37 100 15 71]
        %%
        for design = 3:-1:1
            %% for pressure only
            if design==1
                fname = strcat('NEW_V2_test_results/PL_PCA_data_',num2str(which_test_data),'_deg5_p1',fname_cell{whichfix},'_test');
                y_test =p1(:,which_test_data)';
                x_test =param_test(which_test_data,:);

                f = @(q) uq_eval_uq_metamodel(PCE_P,q(:)')*pca_coeff_p1(:,1:5)' + pca_mu_p1;

                y_pred = f(x_test);
                res0 = (y_test-y_pred).^2;
                sig_vec = ones(32*1,1).*mean(y_test);

            elseif design == 2
                %% MPA pressure, LPA/RPA areas
                fname = strcat('NEW_V2_test_results/PL_PCA_data_',num2str(which_test_data),'_deg5_p1a2a3',fname_cell{whichfix},'_test');

                y_test =[p1(:,which_test_data)'...
                    a2(:,which_test_data)'...
                    a3(:,which_test_data)'];
                x_test =param_test(which_test_data,:);
                % disp([lower; x_test; upper]')


                % Calculate error variance for multicomponent profile likelihood
                f = @(q) [uq_eval_uq_metamodel(pce_p1,q(:)')*pca_coeff_p1(:,1:5)' + pca_mu_p1 ...
                    uq_eval_uq_metamodel(pce_a2,q(:)')*pca_coeff_a2(:,1:5)' + pca_mu_a2 ...
                    uq_eval_uq_metamodel(pce_a3,q(:)')*pca_coeff_a3(:,1:5)' + pca_mu_a3];
                y_pred = f(x_test);
                res0 = (y_test-y_pred).^2;

                sig_vec = ones(32*3,1);
                for i=1:3
                    start_id = (i-1).*32+1;
                    end_id   = i.*32;
                    sig_vec(start_id:end_id)    = mean(y_test(start_id:end_id));
                end


            elseif design==3
                %% MPA pressure, LPA/RPA flows
                fname = strcat('NEW_V2_test_results/PL_PCA_data_',num2str(which_test_data),'_deg5_p1q2q3',fname_cell{whichfix},'_test');

                y_test =[p1(:,which_test_data)'...
                    q2(:,which_test_data)'...
                    q3(:,which_test_data)'];
                x_test =param_test(which_test_data,:);


                % Calculate error variance for multicomponent profile likelihood
                f = @(q) [uq_eval_uq_metamodel(pce_p1,q(:)')*pca_coeff_p1(:,1:5)' + pca_mu_p1 ...
                    uq_eval_uq_metamodel(pce_q2,q(:)')*pca_coeff_q2(:,1:5)' + pca_mu_q2...
                    uq_eval_uq_metamodel(pce_q3,q(:)')*pca_coeff_q3(:,1:5)' + pca_mu_q3];
                y_pred = f(x_test);
                res0 = (y_test-y_pred).^2;

                sig_vec = ones(32*3,1);

                sig_vec = ones(32*3,1);
                for i=1:3
                    start_id = (i-1).*32+1;
                    end_id   = i.*32;
                    sig_vec(start_id:end_id)    = mean(y_test(start_id:end_id));
                end

            else
                error('bad design')
            end

            %%
            true_par = x_test;
            data = y_test;


            params = true_par;

            bounds = [lower; upper]';
            PCA_flag = 1;
            % [J_save,par_save]  = calc_PL(f,data,params,bounds,options,sig_vec);
            par_fix = param_fix_cell{whichfix};
            [J_save,par_save,param_global_opt]  = calc_PL_fixpar(f,data,params,bounds,sig_vec,par_fix);

            save(fname,'x_test','y_test','J_save','par_save','data','sig_vec','f','param_global_opt');
        end
    end
end