clear; clc; close all;
fname = 'PCE_PCA_simple_deg5_';
Names = {'$k_1$','$k_2$','$k_3$','$R_{p,1}$','$R_{p,2}$','$R_{d,1}$','$R_{d,2}$','$C_{T,1}$','$C_{T,2}$'};

% nsamp = 331;%2300;
nsamp = 2500;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

load Cluster\pq_simple.mat A_PCE p_PCE q_PCE

num_par = 9;
num_ves = 3;
n_pca = 5;
num_sob_terms = 5;
time_t = linspace(0,0.85,32);
plot_style = {'-.k','--b','.r',':c','m','-.g','c',':r','--m'};
color_type = turbo(num_par);

%%
for which_ves = 1:3

    fname_save = strcat(fname,num2str(which_ves));
    load(fname_save);

    % Try to calculate the true time dependent sobol indices following
    % https://doi.org/10.1016/j.ress.2019.106737
    % SiY_P_comp = zeros(num_ves,num_par,32); STY_P_comp = zeros(num_ves,num_par,32);
    % score(:,1:n_pca)*coeff_P(:,1:n_pca)'+mu_P;
    % Z is the score, phi is coeff
    p_time = squeeze(p_PCE(:,which_ves,:));
    q_time = squeeze(q_PCE(:,which_ves,:));
    A_time = squeeze(A_PCE(:,which_ves,:));

    var_p_t = var(p_time,[],2)';%*sum(p_explain(1:end)./100);
    var_q_t = var(q_time,[],2)';%*sum(q_explain(1:end)./100);
    var_A_t = var(A_time,[],2)';%*sum(A_explain(1:end)./100);

    % p_time = (PCE_P.ExpDesign.Y*coeff_P(:,1:n_pca)'+mu_P)';
    % q_time = (PCE_Q.ExpDesign.Y*coeff_Q(:,1:n_pca)'+mu_Q)';
    % A_time = (PCE_A.ExpDesign.Y*coeff_A(:,1:n_pca)'+mu_A)';


    PCE_ids = full(PCE_P.PCE(1).Basis.Indices);
    % term_si_p = zeros(num_par,32);
    % term_st_p = zeros(num_par,32);
    ST_p_time = zeros(num_par,32);
    Si_p_time = zeros(num_par,32);

    % term_si_q = zeros(num_par,32);
    % term_st_q = zeros(num_par,32);
    ST_q_time = zeros(num_par,32);
    Si_q_time = zeros(num_par,32);

    % term_si_A = zeros(num_par,32);
    % term_st_A = zeros(num_par,32);
    ST_A_time = zeros(num_par,32);
    Si_A_time = zeros(num_par,32);

    %% Calculate Sobol indices using PCE coefficients
    Zi_p = zeros(num_par,n_pca); Zi_q = zeros(num_par,n_pca); Zi_A = zeros(num_par,n_pca);
    ZTi_p = zeros(num_par,n_pca); ZTi_q = zeros(num_par,n_pca); ZTi_A = zeros(num_par,n_pca);
    weights = 1./(2.*sum(PCE_ids,2)+1);
    varZ_p = var(PCE_P.ExpDesign.Y);
    varZ_q = var(PCE_Q.ExpDesign.Y);
    varZ_A = var(PCE_A.ExpDesign.Y);
    S_i_index = zeros(num_par,num_sob_terms);
    ids_full = 1:length(PCE_P.PCE(1).Coefficients);
    for i=1:num_par
        for j=1:num_sob_terms
            ids_1 = ids_full(find(PCE_ids(1:end,i)==j));
            ids_1 = ids_1(find(sum(PCE_ids(ids_1,:),2)==j));
            S_i_index(i,j) = ids_1;
        end
    end
    for i=1:n_pca
        Zi_p(:,i) = sum((PCE_P.PCE(i).Coefficients(S_i_index).^2) ./ PCE_P.PCE(i).Moments.Var,2);
        Zi_q(:,i) = sum((PCE_Q.PCE(i).Coefficients(S_i_index).^2) ./ PCE_Q.PCE(i).Moments.Var,2);
        Zi_A(:,i) = sum((PCE_A.PCE(i).Coefficients(S_i_index).^2) ./ PCE_A.PCE(i).Moments.Var,2);

        % Zi_p(:,i) = sum((PCE_P.PCE(i).Coefficients(S_i_index).^2),2) ./ varZ_p(i);
        % Zi_q(:,i) = sum((PCE_Q.PCE(i).Coefficients(S_i_index).^2),2) ./ varZ_q(i);
        % Zi_A(:,i) = sum((PCE_A.PCE(i).Coefficients(S_i_index).^2),2) ./ varZ_A(i);
    end
    for i=1:num_par
        ids_i = find(PCE_ids(1:end,i)>0);
        % ids_i = find(PCE_ids(1:end,i)==0);
        % ids_i(1) = [];
        for j=1:n_pca
            ZTi_p(i,j) = sum(PCE_P.PCE(j).Coefficients(ids_i).^2) ./ PCE_P.PCE(j).Moments.Var;
            ZTi_q(i,j) = sum(PCE_Q.PCE(j).Coefficients(ids_i).^2) ./ PCE_Q.PCE(j).Moments.Var;
            ZTi_A(i,j) = sum(PCE_A.PCE(j).Coefficients(ids_i).^2) ./ PCE_A.PCE(j).Moments.Var;

            % ZTi_p(i,j) = sum(PCE_P.PCE(j).Coefficients(ids_i).^2) ./ varZ_p(j);
            % ZTi_q(i,j) = sum(PCE_Q.PCE(j).Coefficients(ids_i).^2) ./ varZ_q(j);
            % ZTi_A(i,j) = sum(PCE_A.PCE(j).Coefficients(ids_i).^2) ./ varZ_A(j);
        end
    end

    %%
    for i=1:num_par
        % ids_1 = find(sum(PCE_ids,2)==1);
        ids_1 = S_i_index(i,:);
        % ids_noti = find(PCE_ids(1:end,i)==0);
        % ids_noti(1) = [];
        ids_noti = find(PCE_ids(1:end,i)>0);
        ids_order = sum(PCE_ids(ids_noti,:),2);
        % ids_remove = ids_order>7;
        % ids_noti(ids_remove) = [];
        % ids_order = sum(PCE_ids(ids_noti,:),2);
        norm_fac_st = (1./(2.*ids_order + 1)).^(0);
        term_si_p = zeros(1,32);
        term_si_q = zeros(1,32);
        term_si_A = zeros(1,32);

        term_st_p = zeros(1,32);
        term_st_q = zeros(1,32);
        term_st_A = zeros(1,32);

        for j=1:n_pca

            for jj=j+1:n_pca
                for t=1:32
                    term_si_p(t) = term_si_p(t) + sum(PCE_P.PCE(j).Coefficients(ids_1).*PCE_P.PCE(jj).Coefficients(ids_1)) .*(coeff_P(t,j).*coeff_P(t,jj))';

                    term_si_q(t) = term_si_q(t) + sum(PCE_Q.PCE(j).Coefficients(ids_1).*PCE_Q.PCE(jj).Coefficients(ids_1)) .*(coeff_Q(t,j).*coeff_Q(t,jj))';

                    term_si_A(t) = term_si_A(t) + sum(PCE_A.PCE(j).Coefficients(ids_1).*PCE_A.PCE(jj).Coefficients(ids_1)) .*(coeff_A(t,j).*coeff_A(t,jj))';

                    term_st_p(t) = term_st_p(t) + (PCE_P.PCE(j).Coefficients(ids_noti)'*(PCE_P.PCE(jj).Coefficients(ids_noti).*norm_fac_st)).*(coeff_P(t,j).*coeff_P(t,jj))';
                    term_st_q(t) = term_st_q(t) + (PCE_Q.PCE(j).Coefficients(ids_noti)'*(PCE_Q.PCE(jj).Coefficients(ids_noti).*norm_fac_st)).*(coeff_Q(t,j).*coeff_Q(t,jj))';
                    term_st_A(t) = term_st_A(t) + (PCE_A.PCE(j).Coefficients(ids_noti)'*(PCE_A.PCE(jj).Coefficients(ids_noti).*norm_fac_st)).*(coeff_A(t,j).*coeff_A(t,jj))';
                end
                % counter = 1;
                % for jjj=ids_noti'
                %     term_st_p(i,:) = term_st_p(i,:) + ((PCE_P.PCE(j).Coefficients(jjj)'*(PCE_P.PCE(jj).Coefficients(jjj).*norm_fac_st(counter))).*coeff_P(:,j).*coeff_P(:,jj))';
                %     term_st_q(i,:) = term_st_q(i,:) + ((PCE_Q.PCE(j).Coefficients(jjj)'*(PCE_Q.PCE(jj).Coefficients(jjj).*norm_fac_st(counter))).*coeff_Q(:,j).*coeff_Q(:,jj))';
                %     term_st_A(i,:) = term_st_A(i,:) + ((PCE_A.PCE(j).Coefficients(jjj)'*(PCE_A.PCE(jj).Coefficients(jjj).*norm_fac_st(counter))).*coeff_A(:,j).*coeff_A(:,jj))';
                %     counter = counter+1;
                % end
            end
            % Si_p_time(i,:) = Si_p_time(i,:) + (Zi_p(i,j).*varZ_p(j))*((coeff_P(:,j)').^2);
            % Si_q_time(i,:) = Si_q_time(i,:) + (Zi_q(i,j).*varZ_q(j))*((coeff_Q(:,j)').^2);
            % Si_A_time(i,:) = Si_A_time(i,:) + (Zi_A(i,j).*varZ_A(j))*((coeff_A(:,j)').^2);
            %
            % ST_p_time(i,:) = ST_p_time(i,:) + ((1-ZTi_p(i,j)).*varZ_p(j))*((coeff_P(:,j)').^2);
            % ST_q_time(i,:) = ST_q_time(i,:) + ((1-ZTi_q(i,j)).*varZ_q(j))*((coeff_Q(:,j)').^2);
            % ST_A_time(i,:) = ST_A_time(i,:) + ((1-ZTi_A(i,j)).*varZ_A(j))*((coeff_A(:,j)').^2);


            % % If we assume the sums are nested
            % Si_p_time(i,:) = Si_p_time(i,:) + (Zi_p(i,j).*PCE_P.PCE(j).Moments.Var).*((coeff_P(:,j)').^2)./var(p_time,[],2)' ...
            %     +2.*term_si_p./var(p_time,[],2)';
            %
            % Si_q_time(i,:) = Si_q_time(i,:) + (Zi_q(i,j).*PCE_Q.PCE(j).Moments.Var).*((coeff_Q(:,j)').^2)./var(q_time,[],2)' ...
            %     +2.*term_si_q./var(q_time,[],2)';
            %
            % Si_A_time(i,:) = Si_A_time(i,:) + (Zi_A(i,j).*PCE_A.PCE(j).Moments.Var).*((coeff_A(:,j)').^2)./var(A_time,[],2)'...
            %     +2.*term_si_A./var(A_time,[],2)';
            %
            % ST_p_time(i,:) = ST_p_time(i,:) + ((1-ZTi_p(i,j)).*PCE_P.PCE(j).Moments.Var)*((coeff_P(:,j)').^2)./var(p_time,[],2)'...
            %     +2.*term_st_p./var(p_time,[],2)';
            %
            % ST_q_time(i,:) = ST_q_time(i,:) + ((1-ZTi_q(i,j)).*PCE_Q.PCE(j).Moments.Var)*((coeff_Q(:,j)').^2)./var(q_time,[],2)'...
            %     +2.*term_st_q./var(q_time,[],2)';
            %
            % ST_A_time(i,:) = ST_A_time(i,:) + ((1-ZTi_A(i,j)).*PCE_A.PCE(j).Moments.Var)*((coeff_A(:,j)').^2)./var(A_time,[],2)'...
            %     +2.*term_st_A./var(A_time,[],2)';

            % If we assume the sums are NOT nested
            %   Si_p_time(i,:) = Si_p_time(i,:) + (Zi_p(i,j).*PCE_P.PCE(j).Moments.Var)*((coeff_P(:,j)').^2)./var_p_t;
            %
            % Si_q_time(i,:) = Si_q_time(i,:) + (Zi_q(i,j).*PCE_Q.PCE(j).Moments.Var)*((coeff_Q(:,j)').^2)./var_q_t;
            %
            % Si_A_time(i,:) = Si_A_time(i,:) + (Zi_A(i,j).*PCE_A.PCE(j).Moments.Var)*((coeff_A(:,j)').^2)./var_A_t;
            %
            % ST_p_time(i,:) = ST_p_time(i,:) + ((1-ZTi_p(i,j)).*PCE_P.PCE(j).Moments.Var)*((coeff_P(:,j)').^2)./var_p_t;
            %
            % ST_q_time(i,:) = ST_q_time(i,:) + ((1-ZTi_q(i,j)).*PCE_Q.PCE(j).Moments.Var)*((coeff_Q(:,j)').^2)./var_q_t;
            %
            % ST_A_time(i,:) = ST_A_time(i,:) + ((1-ZTi_A(i,j)).*PCE_A.PCE(j).Moments.Var)*((coeff_A(:,j)').^2)./var_A_t;

            for t=1:32
                Si_p_time(i,t) = Si_p_time(i,t) + (Zi_p(i,j).*varZ_p(j))*((coeff_P(t,j)').^2)./var_p_t(t);

                Si_q_time(i,t) = Si_q_time(i,t) + (Zi_q(i,j).*varZ_q(j))*((coeff_Q(t,j)').^2)./var_q_t(t);

                Si_A_time(i,t) = Si_A_time(i,t) + (Zi_A(i,j).*varZ_A(j))*((coeff_A(t,j)').^2)./var_A_t(t);

                ST_p_time(i,t) = ST_p_time(i,t) + ((1-ZTi_p(i,j)).*varZ_p(j))*((coeff_P(t,j)').^2)./var_p_t(t);

                ST_q_time(i,t) = ST_q_time(i,t) + ((1-ZTi_q(i,j)).*varZ_q(j))*((coeff_Q(t,j)').^2)./var_q_t(t);

                ST_A_time(i,t) = ST_A_time(i,t) + ((1-ZTi_A(i,j)).*varZ_A(j))*((coeff_A(t,j)').^2)./var_A_t(t);
            end
        end
        % ST_p_time(i,:) = sum((Sobol_P.Results.Total(ids_Tnoti,:).*Sobol_P.Results.TotalVariance)*(coeff_P(:,1:n_pca)').^2,1);
        % Si_p_time(i,:) = (Sobol_P.Results.FirstOrder(i,:).*Sobol_P.Results.TotalVariance)*((coeff_P(:,1:n_pca)').^2);
        % Si_q_time(i,:) = (Sobol_Q.Results.FirstOrder(i,:).*Sobol_Q.Results.TotalVariance)*(coeff_Q(:,1:n_pca)').^2;
        % Si_A_time(i,:) = (Sobol_A.Results.FirstOrder(i,:).*Sobol_A.Results.TotalVariance)*(coeff_A(:,1:n_pca)').^2;
        %
        % ST_p_time(i,:) = ((1-Sobol_P.Results.Total(i,:)).*Sobol_P.Results.TotalVariance)*(coeff_P(:,1:n_pca)').^2;
        % ST_q_time(i,:) = ((1-Sobol_Q.Results.Total(i,:)).*Sobol_Q.Results.TotalVariance)*(coeff_Q(:,1:n_pca)').^2;
        % ST_A_time(i,:) = ((1-Sobol_A.Results.Total(i,:)).*Sobol_A.Results.TotalVariance)*(coeff_A(:,1:n_pca)').^2;

        %% For when the loops are NOT nested
        Si_p_time(i,:) = Si_p_time(i,:) + 2.*term_si_p./var_p_t;

        Si_q_time(i,:) = Si_q_time(i,:) + 2.*term_si_q./var_q_t;

        Si_A_time(i,:) = Si_A_time(i,:) + 2.*term_si_A./var_A_t;

        ST_p_time(i,:) = ST_p_time(i,:) + 2.*term_st_p./var_p_t;

        ST_q_time(i,:) = ST_q_time(i,:) + 2.*term_st_q./var_q_t;

        ST_A_time(i,:) = ST_A_time(i,:) + 2.*term_st_A./var_A_t;

    end


    Si_p_FINAL = Si_p_time;% + 2.*term_si_p./var(p_time,[],2)';
    Si_q_FINAL = Si_q_time;% + 2.*term_si_q./var(q_time,[],2)';
    Si_A_FINAL = Si_A_time;% + 2.*term_si_A./var(A_time,[],2)';

    ST_p_FINAL = 1 - (ST_p_time);% + 2.*term_st_p./var(p_time,[],2)');
    ST_q_FINAL = 1 - (ST_q_time);% + 2.*term_st_q./var(q_time,[],2)');
    ST_A_FINAL = 1 - (ST_A_time);% + 2.*term_st_A./var(A_time,[],2)');

    % figure(1);clf;
    % subplot(1,3,1); plot(ST_p_FINAL','LineWidth',3);
    % subplot(1,3,2); plot(Si_p_FINAL','LineWidth',3);
    % subplot(1,3,3); plot(ST_p_FINAL'-Si_p_FINAL','LineWidth',3);
    % legend(Names);
    % 
    % figure(2);clf;
    % subplot(1,3,1); plot(ST_q_FINAL','LineWidth',3);
    % subplot(1,3,2); plot(Si_q_FINAL','LineWidth',3);
    % subplot(1,3,3); plot(ST_q_FINAL'-Si_q_FINAL','LineWidth',3);
    % legend(Names);
    % 
    % figure(3);clf;
    % subplot(1,3,1); plot(ST_A_FINAL','LineWidth',3);
    % subplot(1,3,2); plot(Si_A_FINAL','LineWidth',3);
    % subplot(1,3,3); plot(ST_A_FINAL'-Si_A_FINAL','LineWidth',3);
    % 
    % legend(Names);

    figure(which_ves);

    subplot(2,3,1); hold on;
    for i=1:num_par
        plot(time_t,Si_p_FINAL(i,:)','LineWidth',3,'Color',color_type(i,:));
    end
    grid on; set(gca,'FontSize',20);
    ylim([0 1]); xlim([0 0.85]);
    ylabel('$S_i^{Y}(t)$');
    title('Pressure')

    subplot(2,3,2); hold on;
    for i=1:num_par
        plot(time_t,Si_q_FINAL(i,:)','LineWidth',3,'Color',color_type(i,:));
    end
    grid on; set(gca,'FontSize',20);
    ylim([0 1]); xlim([0 0.85]);
    title('Flow')

    subplot(2,3,3); hold on;
    for i=1:num_par
        plot(time_t,Si_A_FINAL(i,:)','LineWidth',3,'Color',color_type(i,:));
    end
    grid on; set(gca,'FontSize',20);
    ylim([0 1]); xlim([0 0.85]);
    title('Area')

    subplot(2,3,4); hold on;
    for i=1:num_par
        plot(time_t,ST_p_FINAL(i,:)','LineWidth',3,'Color',color_type(i,:));
    end
    grid on; set(gca,'FontSize',20);
    ylim([0 1]); xlim([0 0.85]);
    ylabel('$S_{T_i}^{Y}(t)$');
    xlabel('Time (s)')

    subplot(2,3,5); hold on;
    for i=1:num_par
        plot(time_t,ST_q_FINAL(i,:)','LineWidth',3,'Color',color_type(i,:));
    end
    grid on; set(gca,'FontSize',20);
    ylim([0 1]); xlim([0 0.85]);
    xlabel('Time (s)')

    subplot(2,3,6); hold on;
    for i=1:num_par
        plot(time_t,ST_A_FINAL(i,:)','LineWidth',3,'Color',color_type(i,:));
    end
    grid on; set(gca,'FontSize',20);
    ylim([0 1]); xlim([0 0.85]);
    xlabel('Time (s)')

    h = legend(Names);
    set(h,'Position',[0.93, 0.5, 0.05, 0.05])

   
    % pause

end

figure(1);
annotation('textbox', [0.05, 0.98, 0, 0], 'String', '(a)', ...
'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 26,'FontWeight','bold');

figure(2);
annotation('textbox', [0.05, 0.98, 0, 0], 'String', '(b)', ...
'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 26,'FontWeight','bold');

figure(3);
annotation('textbox', [0.05, 0.98, 0, 0], 'String', '(c)', ...
'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 26,'FontWeight','bold');


