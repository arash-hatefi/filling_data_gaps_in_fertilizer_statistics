clear;
clc;
close all;
%% Step 1: Importing the Data

A_RCS_PATH = "./data/aRCS.csv";
M_RCS_PATH = "./data/mRCS.csv";
M_RC_PATH = "./data/mRC.csv";
M_CS_PATH = "./data/data/mCS.csv";
M_C_PATH = "./data/mC.csv";
M_R_PATH = "./data/mR.csv";

aRCS_table = readtable(A_RCS_PATH);
mRCS_table = readtable(M_RCS_PATH);
mRC_table = readtable(M_RC_PATH);
mCS_table = readtable(M_CS_PATH);
mC_table = readtable(M_C_PATH);
mR_table = readtable(M_R_PATH);
%% Step 2: Estimating the Value of m (Total Amount of Fertilizer)

mRCS_table_sum_m = sum(mRCS_table.m);
mRC_table_sum_m = sum(mRC_table.m);
mCS_table_sum_m = sum(mCS_table.m);
mC_table_sum_m = sum(mC_table.m);
mR_table_sum_m = sum(mR_table.m);

estimated_m = nanmean([mRCS_table_sum_m, mRC_table_sum_m, mCS_table_sum_m, mC_table_sum_m, mR_table_sum_m]);
%% Step 3: Setting Up the Constraint Equations

X = sym('x', [height(aRCS_table),1]);
%%
mRCS_table.p = mRCS_table.m/estimated_m;
mRCS_table.equation = X;
%%
mRC_table.p = mRC_table.m/estimated_m;
simbolic_values = [];

for row_idx=1:height(mRC_table)
    region = mRC_table.region(row_idx);
    crop = mRC_table.crop(row_idx);
    table = mRCS_table(strcmp(mRCS_table.region, region) & strcmp(mRCS_table.crop, crop), :);
    simbolic_values = [simbolic_values; sum(table.equation)];
end

mRC_table.equation = simbolic_values;
%%
mCS_table.p = mCS_table.m/estimated_m;
simbolic_values = [];

for row_idx=1:height(mCS_table)
    strategy = mCS_table.strategy(row_idx);
    crop = mCS_table.crop(row_idx);
    table = mRCS_table(strcmp(mRCS_table.strategy, strategy) & strcmp(mRCS_table.crop, crop), :);
    simbolic_values = [simbolic_values; sum(table.equation)];
end

mCS_table.equation = simbolic_values;
%%
mC_table.p = mC_table.m/estimated_m;
simbolic_values = [];

for row_idx=1:height(mC_table)
    crop = mC_table.crop(row_idx);
    table = mRCS_table(strcmp(mRCS_table.crop, crop), :);
    simbolic_values = [simbolic_values; sum(table.equation)];
end

mC_table.equation = simbolic_values;
%%
mR_table.p = mR_table.m/estimated_m;
simbolic_values = [];

for row_idx=1:height(mR_table)
    region = mR_table.region(row_idx);
    table = mRCS_table(strcmp(mRCS_table.region, region), :);
    simbolic_values = [simbolic_values; sum(table.equation)];
end

mR_table.equation = simbolic_values;
%% Step 4: Setting Up the Lagrange Formulation

% Parameters
ESTIMATION_METHOD = "MinCrossEntropy";  % Either "MinCrossEntropy" or "MaxEntropy"
CONSTRAINTS = "C2"; % Either "C1" or "C2" or "C3" or "C4" (For more info please see the report)
PRIOR_DIST = "Uniform";%"ProportionalToHarvest"  % Either "Uniform" or "ProportionalToHarvest"
%%
switch CONSTRAINTS
    case "C1"
        equations = [mC_table.equation; mR_table.equation];
        p_values = [mC_table.p; mR_table.p];  
    case "C2"
        equations = [mRC_table.equation; mCS_table.equation];
        p_values = [mRC_table.p; mCS_table.p];
    case "C3"
        equations = [mRC_table.equation; mCS_table.equation; mC_table.equation; mR_table.equation];
        p_values = [mRC_table.p; mCS_table.p; mC_table.p; mR_table.p];
    case "C4"
        equations = [mRC_table.equation; mCS_table.equation; mC_table.equation; mR_table.equation;mRCS_table.equation];
        p_values = [mRC_table.p; mCS_table.p; mC_table.p; mR_table.p;mRCS_table.p];
    otherwise
        assert(true, "Invalied CONSTRAINTS");
end

constraints = [];

for row_idx=1:length(p_values)
    equation = equations(row_idx);
    p = p_values(row_idx);
    if ~isnan(p)
        constraints = [constraints; (equation-p)];
    end
end
constraints = [constraints; sum(mRCS_table.equation)-1];

contraint_coeffs = sym('lambda', [length(constraints),1]);
contraints_term = 0;
for i=1:length(constraints)
     contraints_term = contraints_term + constraints(i)*contraint_coeffs(i);
end
%%
if (PRIOR_DIST=="Uniform")
    prior_dist = ones(height(mRCS_table),1)/height(mRCS_table);
elseif (PRIOR_DIST=="ProportionalToHarvest")
    prior_dist = aRCS_table.a;
    prior_dist = prior_dist / sum(prior_dist);   
else
    assert(true, "Invalied PRIOR_DIST");
end
%%
if (ESTIMATION_METHOD=="MaxEntropy")
    entropy_term = -sum(mRCS_table.equation.*log(mRCS_table.equation));
    main_term = entropy_term;
elseif (ESTIMATION_METHOD=="MinCrossEntropy")
    kl_divergence_term = (-sum(mRCS_table.equation.*log(prior_dist./mRCS_table.equation)));
    main_term = -kl_divergence_term;
else
    assert(true, "Invalied STIMATION_METHOD")
end
%%
lagrangean = main_term + contraints_term;
%%
all_variables = [X; contraint_coeffs];

lagrangean_diffs = [];
for i=1:length(all_variables)
    lagrangean_diffs = [lagrangean_diffs; diff(lagrangean, all_variables(i))];
end
%%
inputs = sym('input', [length(all_variables),1]);
inputs_array = [];
for i=1:length(all_variables)
    inputs_array = [inputs_array; inputs(i)];
end
%%
lagrangean_diffs_func = subs(lagrangean_diffs,all_variables,inputs_array);
lagrangean_diffs_func = matlabFunction(lagrangean_diffs_func,'vars',{inputs});
%% Step 5: Setting the Optimizer Parameters

lambda_min_value = -1000;
lambda_max_value = 1000;
lambda_init_value = 4;

lb = [zeros(length(prior_dist),1);ones(length(contraint_coeffs),1)*lambda_min_value];
ub = [ones(length(prior_dist),1);ones(length(contraint_coeffs),1)*lambda_max_value];
init_value = [prior_dist;ones(length(contraint_coeffs),1)*lambda_init_value];

options = optimoptions('lsqnonlin');
options.MaxFunctionEvaluations = 900000;
options.MaxIterations = 90000;
options.Display = 'iter';
%% Step 6: Solving the Optimization Problem

results = lsqnonlin(lagrangean_diffs_func,init_value,lb,ub);
%% Step 7: Post Process Results

predicted_probabilities = results(1:length(X));
predicted_m_values = predicted_probabilities .* estimated_m;
%%
rss = nanmean((predicted_probabilities-mRCS_table.p).^2);
tss = nanvar(mRCS_table.p);
r_squared = 1 - rss/tss;
%%
hold off;
scatter(mRCS_table.p, predicted_probabilities);
hold on;
limit = max([predicted_probabilities;mRCS_table.p]);
plot([0,limit*1.1],[0,limit*1.1]);
xlabel("Avalilabe Probabilities in Table",'FontSize',14);
ylabel("Predicted Probabilities",'FontSize',14);
text(limit*0.1,limit*0.95,"r^2="+string(r_squared),'FontSize',14);
grid on;
axis square;