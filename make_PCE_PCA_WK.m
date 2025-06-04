%% PCE METAMODELING: MULTIPLE OUTPUTS
%
% This example showcases an application of polynomial chaos expansion
% (PCE) to the metamodeling of a simply supported beam model
% with multiple outputs.
% The model computes the deflections at several points along the length
% of the beam subjected to a uniform random load.

%% 1 - INITIALIZE UQLAB
%
% Clear all variables from the workspace, set the random number generator
% for reproducible results, and initialize the UQLab framework:
% clearvars
rng(2200,'twister')
uqlab
pce_deg = 5;
n_pca = 5;
load Cluster\pq_simple.mat A_PCE p_PCE q_PCE param_sample low upp
fname = 'PCE_PCA_simple_deg5_'
vel_power = 9;
visc = 0.032; % 1 g/cm/s is 1 dyne s /cm^2
%% 2- Names
Names = {'f1_{A}','f2_{A}','f3_{A}','Rp1','Rp2','Rd1','Rd2','CT1','CT2'};

param_sample(:,[1 3]) = log(param_sample(:,[1 3]));
upp(1) = log(upp(1));
low(1) = log(low(1));
%% 3 - PROBABILISTIC INPUT MODEL
%
% The simply supported beam model has five inputs,
% modeled by independent lognormal random variables.
% The detailed model is given in the following table:
MetaOpts.Display = 'quiet';
[num_samp,num_par] = size(param_sample);

Input = [];
% Define an INPUT object with the following marginals:
for i=1:num_par
    Input.Marginals(i).Name = Names{i};  % beam width
    Input.Marginals(i).Type = 'Uniform';
    Input.Marginals(i).Parameters = [low(i) upp(i)];  % (m)
end
myInput = uq_createInput(Input);



% Initialize arrays
subsamp = 1;
num_pts_outs = 32;
QoI_P    = zeros(num_pts_outs,num_samp);
QoI_Q    = zeros(num_pts_outs,num_samp);
QoI_CS   = zeros(num_pts_outs,num_samp);
QoI_WSS  = zeros(num_pts_outs,num_samp);
x_interp = linspace(0,1,num_pts_outs);

%%

N_train = length(param_sample);%1900;%min(nchoosek(num_par+pce_deg,pce_deg).*2,length(par_sample));
samps = 1:N_train;
% Data
param_sample = param_sample(samps,:);
X_train   = param_sample;

MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'PCE';
MetaOpts.Degree = pce_deg;
MetaOpts.Method = 'OLS';
MetaOpts.ExpDesign.X = X_train;

PCESobol.Type = 'Sensitivity';
PCESobol.Method = 'Sobol';
PCESobol.Sobol.Order = 3;
t = linspace(0,0.85,num_pts_outs);

%%
for which_ves = 1:3
        %Look at the QoI
        
    QoI_P = squeeze(p_PCE(1:subsamp:end,which_ves,samps));
    QoI_Q = squeeze(q_PCE(1:subsamp:end,which_ves,samps));
    par_stiff = param_sample(:,1);

    if which_ves>15
        QoI_Q = -QoI_Q;
        par_stiff = param_sample(:,2);
    end
    QoI_A = squeeze(A_PCE(1:subsamp:end,which_ves,samps));
    

    %%
    % Create the PCE metamodels, but use PCA
    % score(10,1:4)*coeff(:,1:4)' + mu
    
    [coeff_P, score, ~,~,explained, mu_P] = pca(QoI_P');
    y_pca = score(:,1:n_pca);
    MetaOpts.ExpDesign.Y = y_pca;
    PCE_P = uq_createModel(MetaOpts);
    Sobol_P = uq_createAnalysis(PCESobol);
    sum(explained(1:n_pca))

    [coeff_Q, score, ~,~,explained, mu_Q] = pca(QoI_Q');
    y_pca = score(:,1:n_pca);
    MetaOpts.ExpDesign.Y = y_pca;
    PCE_Q = uq_createModel(MetaOpts);
    Sobol_Q = uq_createAnalysis(PCESobol);
    sum(explained(1:n_pca))

     [coeff_A, score, ~,~,explained, mu_A] = pca(QoI_A');
    y_pca = score(:,1:n_pca);
    MetaOpts.ExpDesign.Y = y_pca;
    PCE_A = uq_createModel(MetaOpts);
    Sobol_A = uq_createAnalysis(PCESobol);
    sum(explained(1:n_pca))
  

            fname_save = strcat(fname,num2str(which_ves));
            save(fname_save,'PCE_P','PCE_Q','PCE_A', ...
                  'Sobol_P','Sobol_Q','Sobol_A',...
                  'coeff_P','coeff_Q','coeff_A',...
                  'mu_P','mu_Q','mu_A');

end

