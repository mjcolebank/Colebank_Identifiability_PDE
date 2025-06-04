Names = {'f1_{A}','f2_{A}','f3_{A}','Rp1','Rp2','Rd1','Rd2','CT1','CT2'};

uqlab
load('PCE_PCA_simple_deg5_1.mat','PCE_P','PCE_A','coeff_P','coeff_A','mu_P','mu_A');
pce_p1 = PCE_P;
pce_a1 = PCE_A;
pca_coeff_p1 = coeff_P;
pca_coeff_a1 = coeff_A;
pca_mu_p1    = mu_P;
pca_mu_a1    = mu_A;


load('PCE_PCA_simple_deg5_2.mat','PCE_Q','PCE_A','coeff_Q','coeff_A','mu_Q','mu_A');
pce_q2 = PCE_Q;
pce_a2 = PCE_A;
pca_coeff_q2 = coeff_Q;
pca_coeff_a2 = coeff_A;
pca_mu_q2    = mu_Q;
pca_mu_a2    = mu_A;

load('PCE_PCA_simple_deg5_3.mat','PCE_Q','PCE_A','coeff_Q','coeff_A','mu_Q','mu_A');
pce_q3 = PCE_Q;
pce_a3 = PCE_A;
pca_coeff_q3 = coeff_Q;
pca_coeff_a3 = coeff_A;
pca_mu_q3    = mu_Q;
pca_mu_a3    = mu_A;

% load simple_model_ST_1_25_pars.mat LB UB
load Cluster\pq_test.mat
p1 = squeeze(p_PCE(:,1,:));
a1 = squeeze(A_PCE(:,1,:));
q2 = squeeze(q_PCE(:,2,:));
a2 = squeeze(A_PCE(:,2,:));
q3 = squeeze(q_PCE(:,3,:));
a3 = squeeze(A_PCE(:,3,:));
param_test = param_sample;

upper = upp; upper(1) = log(upper(1));
lower = low; lower(1) = log(lower(1));
param_test(:,[1 3]) = log(param_test(:,[1 3]));

% param_sample(:,end) = param_sample(:,end)./1e-3;
% upper(end) = upper(end)./1e-3;
% lower(end) = lower(end)./1e-3;

for which_test_data = [74 22 14 72 34]
    %%
    for design = 3:-1:1
        %% for pressure only
        if design==1
            fname = strcat('test_results/PL_PCA_data_',num2str(which_test_data),'_deg5_p1_fixk1k2_CT_test');
            y_test =p1(:,which_test_data)';
            x_test =param_test(which_test_data,:);

            f = @(q) uq_eval_uq_metamodel(PCE_P,q(:)')*pca_coeff_p1(:,1:5)' + pca_mu_p1;

            y_pred = f(x_test);
            res0 = (y_test-y_pred).^2;
            sig_vec = ones(32*1,1).*mean(y_test);

        elseif design == 2
            %% MPA pressure, LPA/RPA areas
            fname = strcat('test_results/PL_PCA_data_',num2str(which_test_data),'_deg5_p1a2a3_fixk1k2_CT_test');
           
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
            fname = strcat('test_results/PL_PCA_data_',num2str(which_test_data),'_deg5_p1q2q3_fixk1k2_CT_test');

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
        par_fix = [1 2 8 9];
        [J_save,par_save,param_global_opt]  = calc_PL_fixpar(f,data,params,bounds,sig_vec,par_fix);

        save(fname,'x_test','y_test','J_save','par_save','data','sig_vec','f','param_global_opt');
    end
end