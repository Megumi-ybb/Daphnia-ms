
#Set Daphnia-ms as the directory

#################################Target model########################################################
load("./Mixed-species/SIRJPF2/model/best_result.rda")
all_shared_ll = pf.loglik.of.mif.estimate
all_shared_ll_se = s.e.of.pf.loglik.of.mif.estimate
parameters = mif.estimate

# -------------------------------------Specific------------------------------------------------
load("./Mixed-species/SIRJPF2/specific_model/xi/specific_xi.RData")
xi_ll = pf.loglik.of.mif.estimate
xi_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Mixed-species/SIRJPF2/specific_model/theta_Si_theta_Sn/specific_theta_Si_theta_Sn.RData")
theta_Si_theta_Sn_ll = pf.loglik.of.mif.estimate
theta_Si_theta_Sn_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Mixed-species/SIRJPF2/specific_model/theta_P/specific_theta_P.RData")
theta_P_ll = pf.loglik.of.mif.estimate
theta_P_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Mixed-species/SIRJPF2/specific_model/theta_Ii_theta_In/specific_theta_Ii_theta_In.RData")
theta_Ii_theta_Ii_ll = pf.loglik.of.mif.estimate
theta_Ii_theta_Ii_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Mixed-species/SIRJPF2/specific_model/ri_rn/specific_ri_rn.RData")
ri_rn_ll = pf.loglik.of.mif.estimate
ri_rn_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Mixed-species/SIRJPF2/specific_model/probi_probn/specific_probi_probn.RData")
probi_probn_ll = pf.loglik.of.mif.estimate
probi_probn_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Mixed-species/SIRJPF2/specific_model/f_Si_f_Sn/specific_f_Si_f_Sn.RData")
f_Si_f_Sn_ll = pf.loglik.of.mif.estimate
f_Si_f_Sn_ll_se = s.e.of.pf.loglik.of.mif.estimate


# -------------------------------------Specific Block------------------------------------------------
load("./Mixed-species/SIRJPF2/specific_model_block/xi/specific_xi.RData")
xi_block_ll = pf.loglik.of.mif.estimate
xi_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Mixed-species/SIRJPF2/specific_model_block/theta_Si_theta_Sn/specific_theta_Si_theta_Sn.RData")
theta_Si_theta_Sn_block_ll = pf.loglik.of.mif.estimate
theta_Si_theta_Sn_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Mixed-species/SIRJPF2/specific_model_block/theta_P/specific_theta_P.RData")
theta_P_block_ll = pf.loglik.of.mif.estimate
theta_P_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Mixed-species/SIRJPF2/specific_model_block/theta_Ii_theta_In/specific_theta_Ii_theta_In.RData")
theta_Ii_theta_Ii_block_ll = pf.loglik.of.mif.estimate
theta_Ii_theta_Ii_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Mixed-species/SIRJPF2/specific_model_block/ri_rn/specific_ri_rn.RData")
ri_rn_block_ll = pf.loglik.of.mif.estimate
ri_rn_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Mixed-species/SIRJPF2/specific_model_block/probi_probn/specific_probi_probn.RData")
probi_probn_block_ll = pf.loglik.of.mif.estimate
probi_probn_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Mixed-species/SIRJPF2/specific_model_block/f_Si_f_Sn/specific_f_Si_f_Sn.RData")
f_Si_f_Sn_block_ll = pf.loglik.of.mif.estimate
f_Si_f_Sn_block_ll_se = s.e.of.pf.loglik.of.mif.estimate


K_all = length(parameters) - 2
K_two = length(parameters) - 4 + 2 * 8
K_one = length(parameters) - 3 + 1 * 8


target_para_ll_df =  list(
  all_shared = list(ll = all_shared_ll, ll_se = all_shared_ll_se, AIC = -2 * all_shared_ll + 2 * K_all),
  xi = list(ll = xi_ll, ll_se = xi_ll_se, AIC = -2 * xi_ll + 2 * K_one),
  theta_Si_theta_Sn = list(ll = theta_Si_theta_Sn_ll, ll_se = theta_Si_theta_Sn_ll_se, AIC = -2 * theta_Si_theta_Sn_ll + 2 * K_two),
  theta_P = list(ll = theta_P_ll, ll_se = theta_P_ll_se, AIC = -2 * theta_P_ll + 2 * K_one),
  theta_Ii_theta_Ii = list(ll = theta_Ii_theta_Ii_ll, ll_se = theta_Ii_theta_Ii_ll_se, AIC = -2 * theta_Ii_theta_Ii_ll + 2 * K_two),
  ri_rn = list(ll = ri_rn_ll, ll_se = ri_rn_ll_se, AIC = -2 * ri_rn_ll + 2 * K_two),
  probi_probn = list(ll = probi_probn_ll, ll_se = probi_probn_ll_se, AIC = -2 * probi_probn_ll + 2 * K_two),
  f_Si_f_Sn = list(ll = f_Si_f_Sn_ll, ll_se = f_Si_f_Sn_ll_se, AIC = -2 * f_Si_f_Sn_ll + 2 * K_two),
  xi_block = list(ll = xi_block_ll, ll_se = xi_block_ll_se, AIC = -2 * xi_block_ll + 2 * K_one),
  theta_Si_theta_Sn_block = list(ll = theta_Si_theta_Sn_block_ll, ll_se = theta_Si_theta_Sn_block_ll_se, AIC = -2 * theta_Si_theta_Sn_block_ll + 2 * K_two),
  theta_P_block = list(ll = theta_P_block_ll, ll_se = theta_P_block_ll_se, AIC = -2 * theta_P_block_ll + 2 * K_one),
  theta_Ii_theta_Ii_block = list(ll = theta_Ii_theta_Ii_block_ll, ll_se = theta_Ii_theta_Ii_block_ll_se, AIC = -2 * theta_Ii_theta_Ii_block_ll + 2 * K_two),
  ri_rn_block = list(ll = ri_rn_block_ll, ll_se = ri_rn_block_ll_se, AIC = -2 * ri_rn_block_ll + 2 * K_two),
  probi_probn_block = list(ll = probi_probn_block_ll, ll_se = probi_probn_block_ll_se, AIC = -2 * probi_probn_block_ll + 2 * K_two),
  f_Si_f_Sn_block = list(ll = f_Si_f_Sn_block_ll, ll_se = f_Si_f_Sn_block_ll_se, AIC = -2 * f_Si_f_Sn_block_ll + 2 * K_two)
)

parameter_table <- data.frame(
  ll = sapply(target_para_ll_df[!grepl("_block", names(target_para_ll_df))], function(x) x$ll),
  AIC = sapply(target_para_ll_df[!grepl("_block", names(target_para_ll_df))], function(x) x$AIC),
  block_ll = sapply(names(target_para_ll_df[!grepl("_block", names(target_para_ll_df))]), function(name) {
    block_name <- paste0(name, "_block")
    if (block_name %in% names(target_para_ll_df)) target_para_ll_df[[block_name]]$ll else NA
  }),
  block_AIC = sapply(names(target_para_ll_df[!grepl("_block", names(target_para_ll_df))]), function(name) {
    block_name <- paste0(name, "_block")
    if (block_name %in% names(target_para_ll_df)) target_para_ll_df[[block_name]]$AIC else NA
  }),
  row.names = names(target_para_ll_df[!grepl("_block", names(target_para_ll_df))])
)

parameter_table$max_AIC <- pmax(parameter_table$AIC, parameter_table$block_AIC, na.rm = TRUE)
parameter_table <- parameter_table[order(parameter_table$max_AIC), ]
parameter_table$max_AIC <- NULL
print(parameter_table)
target_para_parameter_table = parameter_table
save(target_para_ll_df = target_para_ll_df,target_para_parameter_table = target_para_parameter_table,file = './data/Target_dynamics/para/Target_para_loglik_df.rds')


#################################Target no para model########################################################
load("./Mixed-species/SRJF2/best_result.rda")
all_shared_ll = pf.loglik.of.mif.estimate
all_shared_ll_se = s.e.of.pf.loglik.of.mif.estimate
parameters = mif.estimate

# -------------------------------------Specific------------------------------------------------
load("./Mixed-species/SRJF2/specific_model/theta_Sn_theta_Si/specific_theta_Sn_theta_Si.RData")
theta_Si_theta_Sn_ll = pf.loglik.of.mif.estimate
theta_Si_theta_Sn_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Mixed-species/SRJF2/specific_model/rn_ri/specific_rn_ri.RData")
ri_rn_ll = pf.loglik.of.mif.estimate
ri_rn_ll_se = s.e.of.pf.loglik.of.mif.estimate


load("./Mixed-species/SRJF2/specific_model/f_Si_f_Sn/specific_f_Si_f_Sn.RData")
f_Si_f_Sn_ll = pf.loglik.of.mif.estimate
f_Si_f_Sn_ll_se = s.e.of.pf.loglik.of.mif.estimate


# -------------------------------------Specific Block------------------------------------------------
load("./Mixed-species/SRJF2/specific_model_block/theta_Sn_theta_Si/specific_theta_Sn_theta_Si.RData")
theta_Si_theta_Sn_block_ll = pf.loglik.of.mif.estimate
theta_Si_theta_Sn_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Mixed-species/SRJF2/specific_model_block/rn_ri/specific_rn_ri.RData")
ri_rn_block_ll = pf.loglik.of.mif.estimate
ri_rn_block_ll_se = s.e.of.pf.loglik.of.mif.estimate


load("./Mixed-species/SRJF2/specific_model_block/f_Si_f_Sn/specific_f_Si_f_Sn.RData")
f_Si_f_Sn_block_ll = pf.loglik.of.mif.estimate
f_Si_f_Sn_block_ll_se = s.e.of.pf.loglik.of.mif.estimate


K_all = length(parameters) - 2
K_two = length(parameters) - 4 + 2 * 9
K_one = length(parameters) - 3 + 1 * 9


target_no_para_ll_df =  list(
  all_shared = list(ll = all_shared_ll, ll_se = all_shared_ll_se, AIC = -2 * all_shared_ll + 2 * K_all),
  theta_Si_theta_Sn = list(ll = theta_Si_theta_Sn_ll, ll_se = theta_Si_theta_Sn_ll_se, AIC = -2 * theta_Si_theta_Sn_ll + 2 * K_two),
  ri_rn = list(ll = ri_rn_ll, ll_se = ri_rn_ll_se, AIC = -2 * ri_rn_ll + 2 * K_two),
  f_Si_f_Sn = list(ll = f_Si_f_Sn_ll, ll_se = f_Si_f_Sn_ll_se, AIC = -2 * f_Si_f_Sn_ll + 2 * K_two),
  theta_Si_theta_Sn_block = list(ll = theta_Si_theta_Sn_block_ll, ll_se = theta_Si_theta_Sn_block_ll_se, AIC = -2 * theta_Si_theta_Sn_block_ll + 2 * K_two),
  ri_rn_block = list(ll = ri_rn_block_ll, ll_se = ri_rn_block_ll_se, AIC = -2 * ri_rn_block_ll + 2 * K_two),
  f_Si_f_Sn_block = list(ll = f_Si_f_Sn_block_ll, ll_se = f_Si_f_Sn_block_ll_se, AIC = -2 * f_Si_f_Sn_block_ll + 2 * K_two)
)

parameter_table <- data.frame(
  ll = sapply(target_no_para_ll_df[!grepl("_block", names(target_no_para_ll_df))], function(x) x$ll),
  AIC = sapply(target_no_para_ll_df[!grepl("_block", names(target_no_para_ll_df))], function(x) x$AIC),
  block_ll = sapply(names(target_no_para_ll_df[!grepl("_block", names(target_no_para_ll_df))]), function(name) {
    block_name <- paste0(name, "_block")
    if (block_name %in% names(target_no_para_ll_df)) target_no_para_ll_df[[block_name]]$ll else NA
  }),
  block_AIC = sapply(names(target_no_para_ll_df[!grepl("_block", names(target_no_para_ll_df))]), function(name) {
    block_name <- paste0(name, "_block")
    if (block_name %in% names(target_no_para_ll_df)) target_no_para_ll_df[[block_name]]$AIC else NA
  }),
  row.names = names(target_no_para_ll_df[!grepl("_block", names(target_no_para_ll_df))])
)

parameter_table$max_AIC <- pmax(parameter_table$AIC, parameter_table$block_AIC, na.rm = TRUE)
parameter_table <- parameter_table[order(parameter_table$max_AIC), ]
parameter_table$max_AIC <- NULL
print(parameter_table)
target_no_para_parameter_table = parameter_table
save(target_no_para_ll_df = target_no_para_ll_df,target_no_para_parameter_table = target_no_para_parameter_table,file = './data/Target_dynamics/no_para/Target_no_para_loglik_df.rds')



#######################################################################Dent para model########################################################
load("./Single-species/Dent/SIRJPF/model/best_result.rda")
all_shared_ll = pf.loglik.of.mif.estimate
all_shared_ll_se = s.e.of.pf.loglik.of.mif.estimate
parameters = mif.estimate

# -------------------------------------Specific------------------------------------------------
load("./Single-species/Dent/SIRJPF/specific_model/xi/specific_xi.RData")
xi_ll = pf.loglik.of.mif.estimate
xi_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Dent/SIRJPF/specific_model/theta_Sn/specific_theta_Sn.RData")
theta_Sn_ll = pf.loglik.of.mif.estimate
theta_Sn_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Dent/SIRJPF/specific_model/theta_P/specific_theta_P.RData")
theta_P_ll = pf.loglik.of.mif.estimate
theta_P_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Dent/SIRJPF/specific_model/theta_In/specific_theta_In.RData")
theta_In_ll = pf.loglik.of.mif.estimate
theta_In_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Dent/SIRJPF/specific_model/rn/specific_rn.RData")
ri_rn_ll = pf.loglik.of.mif.estimate
ri_rn_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Dent/SIRJPF/specific_model/probn/specific_probn.RData")
probn_ll = pf.loglik.of.mif.estimate
probn_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Dent/SIRJPF/specific_model/f_Sn/specific_f_Sn.RData")
f_Sn_ll = pf.loglik.of.mif.estimate
f_Sn_ll_se = s.e.of.pf.loglik.of.mif.estimate


# -------------------------------------Specific Block------------------------------------------------
load("./Single-species/Dent/SIRJPF/specific_model_block/xi/specific_xi.RData")
xi_block_ll = pf.loglik.of.mif.estimate
xi_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Dent/SIRJPF/specific_model_block/theta_Sn/specific_theta_Sn.RData")
theta_Sn_block_ll = pf.loglik.of.mif.estimate
theta_Sn_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Dent/SIRJPF/specific_model_block/theta_P/specific_theta_P.RData")
theta_P_block_ll = pf.loglik.of.mif.estimate
theta_P_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Dent/SIRJPF/specific_model_block/theta_In/specific_theta_In.RData")
theta_In_block_ll = pf.loglik.of.mif.estimate
theta_In_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Dent/SIRJPF/specific_model_block/rn/specific_rn.RData")
rn_block_ll = pf.loglik.of.mif.estimate
rn_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Dent/SIRJPF/specific_model_block/probn/specific_probn.RData")
probn_block_ll = pf.loglik.of.mif.estimate
probn_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Dent/SIRJPF/specific_model_block/f_Sn/specific_f_Sn.RData")
f_Sn_block_ll = pf.loglik.of.mif.estimate
f_Sn_block_ll_se = s.e.of.pf.loglik.of.mif.estimate


K_all = length(parameters) - 1
K_two = length(parameters) - 3 + 2 * 8
K_one = length(parameters) - 2 + 1 * 8


dent_para_ll_df <- list(
  all_shared = list(ll = all_shared_ll, ll_se = all_shared_ll_se, AIC = -2 * all_shared_ll + 2 * K_all),
  xi = list(ll = xi_ll, ll_se = xi_ll_se, AIC = -2 * xi_ll + 2 * K_one),
  theta_Sn = list(ll = theta_Sn_ll, ll_se = theta_Sn_ll_se, AIC = -2 * theta_Sn_ll + 2 * K_one),
  theta_P = list(ll = theta_P_ll, ll_se = theta_P_ll_se, AIC = -2 * theta_P_ll + 2 * K_one),
  theta_In = list(ll = theta_In_ll, ll_se = theta_In_ll_se, AIC = -2 * theta_In_ll + 2 * K_one),
  rn = list(ll = ri_rn_ll, ll_se = ri_rn_ll_se, AIC = -2 * ri_rn_ll + 2 * K_one),
  probn = list(ll = probn_ll, ll_se = probn_ll_se, AIC = -2 * probn_ll + 2 * K_one),
  f_Sn = list(ll = f_Sn_ll, ll_se = f_Sn_ll_se, AIC = -2 * f_Sn_ll + 2 * K_one),
  xi_block = list(ll = xi_block_ll, ll_se = xi_block_ll_se, AIC = -2 * xi_block_ll + 2 * K_one),
  theta_Sn_block = list(ll = theta_Sn_block_ll, ll_se = theta_Sn_block_ll_se, AIC = -2 * theta_Sn_block_ll + 2 * K_one),
  theta_P_block = list(ll = theta_P_block_ll, ll_se = theta_P_block_ll_se, AIC = -2 * theta_P_block_ll + 2 * K_one),
  theta_In_block = list(ll = theta_In_block_ll, ll_se = theta_In_block_ll_se, AIC = -2 * theta_In_block_ll + 2 * K_one),
  rn_block = list(ll = rn_block_ll, ll_se = rn_block_ll_se, AIC = -2 * rn_block_ll + 2 * K_one),
  probn_block = list(ll = probn_block_ll, ll_se = probn_block_ll_se, AIC = -2 * probn_block_ll + 2 * K_one),
  f_Sn_block = list(ll = f_Sn_block_ll, ll_se = f_Sn_block_ll_se, AIC = -2 * f_Sn_block_ll + 2 * K_one)
)

dent_para_parameter_table <- data.frame(
  ll = sapply(dent_para_ll_df[!grepl("_block", names(dent_para_ll_df))], function(x) x$ll),
  AIC = sapply(dent_para_ll_df[!grepl("_block", names(dent_para_ll_df))], function(x) x$AIC),
  block_ll = sapply(names(dent_para_ll_df[!grepl("_block", names(dent_para_ll_df))]), function(name) {
    block_name <- paste0(name, "_block")
    if (block_name %in% names(dent_para_ll_df)) dent_para_ll_df[[block_name]]$ll else NA
  }),
  block_AIC = sapply(names(dent_para_ll_df[!grepl("_block", names(dent_para_ll_df))]), function(name) {
    block_name <- paste0(name, "_block")
    if (block_name %in% names(dent_para_ll_df)) dent_para_ll_df[[block_name]]$AIC else NA
  }),
  row.names = names(dent_para_ll_df[!grepl("_block", names(dent_para_ll_df))])
)

dent_para_parameter_table$max_AIC <- pmax(dent_para_parameter_table$AIC, dent_para_parameter_table$block_AIC, na.rm = TRUE)
dent_para_parameter_table <- dent_para_parameter_table[order(dent_para_parameter_table$max_AIC), ]
dent_para_parameter_table$max_AIC <- NULL
print(dent_para_parameter_table)

save(dent_para_ll_df = dent_para_ll_df,dent_para_parameter_table = dent_para_parameter_table,file = './data/Single-species/Dent/SRJF/para/Dent_para_loglik_df.rds')






#######################################################################Dent no para model########################################################
load("./Single-species/Dent/SRJF/model/best_result.rda")
all_shared_ll = pf.loglik.of.mif.estimate
all_shared_ll_se = s.e.of.pf.loglik.of.mif.estimate
parameters = mif.estimate
# -------------------------------------Specific------------------------------------------------
load("./Single-species/Dent/SRJF/specific_model/theta_Sn/specific_theta_Sn.RData")
theta_Sn_ll = pf.loglik.of.mif.estimate
theta_Sn_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Dent/SRJF/specific_model/rn/specific_rn.RData")
ri_rn_ll = pf.loglik.of.mif.estimate
ri_rn_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Dent/SRJF/specific_model/f_Sn/specific_f_Sn.RData")
f_Sn_ll = pf.loglik.of.mif.estimate
f_Sn_ll_se = s.e.of.pf.loglik.of.mif.estimate


# -------------------------------------Specific Block------------------------------------------------
load("./Single-species/Dent/SRJF/specific_model_block/theta_Sn/specific_theta_Sn.RData")
theta_Sn_block_ll = pf.loglik.of.mif.estimate
theta_Sn_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Dent/SRJF/specific_model_block/rn/specific_rn.RData")
rn_block_ll = pf.loglik.of.mif.estimate
rn_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Dent/SRJF/specific_model_block/f_Sn/specific_f_Sn.RData")
f_Sn_block_ll = pf.loglik.of.mif.estimate
f_Sn_block_ll_se = s.e.of.pf.loglik.of.mif.estimate


K_all = length(parameters) - 1
K_two = length(parameters) - 3 + 2 * 10
K_one = length(parameters) - 2 + 1 * 10


dent_no_para_ll_df <- list(
  all_shared = list(ll = all_shared_ll, ll_se = all_shared_ll_se, AIC = -2 * all_shared_ll + 2 * K_all),
  theta_Sn = list(ll = theta_Sn_ll, ll_se = theta_Sn_ll_se, AIC = -2 * theta_Sn_ll + 2 * K_one),
  rn = list(ll = ri_rn_ll, ll_se = ri_rn_ll_se, AIC = -2 * ri_rn_ll + 2 * K_one),
  f_Sn = list(ll = f_Sn_ll, ll_se = f_Sn_ll_se, AIC = -2 * f_Sn_ll + 2 * K_one),
  theta_Sn_block = list(ll = theta_Sn_block_ll, ll_se = theta_Sn_block_ll_se, AIC = -2 * theta_Sn_block_ll + 2 * K_one),
  rn_block = list(ll = rn_block_ll, ll_se = rn_block_ll_se, AIC = -2 * rn_block_ll + 2 * K_one),
  f_Sn_block = list(ll = f_Sn_block_ll, ll_se = f_Sn_block_ll_se, AIC = -2 * f_Sn_block_ll + 2 * K_one)
)

dent_no_para_parameter_table <- data.frame(
  ll = sapply(dent_no_para_ll_df[!grepl("_block", names(dent_no_para_ll_df))], function(x) x$ll),
  AIC = sapply(dent_no_para_ll_df[!grepl("_block", names(dent_no_para_ll_df))], function(x) x$AIC),
  block_ll = sapply(names(dent_no_para_ll_df[!grepl("_block", names(dent_no_para_ll_df))]), function(name) {
    block_name <- paste0(name, "_block")
    if (block_name %in% names(dent_no_para_ll_df)) dent_no_para_ll_df[[block_name]]$ll else NA
  }),
  block_AIC = sapply(names(dent_no_para_ll_df[!grepl("_block", names(dent_no_para_ll_df))]), function(name) {
    block_name <- paste0(name, "_block")
    if (block_name %in% names(dent_no_para_ll_df)) dent_no_para_ll_df[[block_name]]$AIC else NA
  }),
  row.names = names(dent_no_para_ll_df[!grepl("_block", names(dent_no_para_ll_df))])
)

dent_no_para_parameter_table$max_AIC <- pmax(dent_no_para_parameter_table$AIC, dent_no_para_parameter_table$block_AIC, na.rm = TRUE)
dent_no_para_parameter_table <- dent_no_para_parameter_table[order(dent_no_para_parameter_table$max_AIC), ]
dent_no_para_parameter_table$max_AIC <- NULL
print(dent_no_para_parameter_table)

save(dent_no_para_ll_df = dent_no_para_ll_df,dent_no_para_parameter_table = dent_no_para_parameter_table,file = './data/Simple_dynamics/Dent/no_para/Dent_no_para_loglik_df.rds')





#######################################################################Lum para model########################################################
load("./Single-species/Lum/SIRJPF/model/best_result.rda")
all_shared_ll = pf.loglik.of.mif.estimate
all_shared_ll_se = s.e.of.pf.loglik.of.mif.estimate
parameters = mif.estimate

# -------------------------------------Specific------------------------------------------------
load("./Single-species/Lum/SIRJPF/specific_model/xi/specific_xi.RData")
xi_ll = pf.loglik.of.mif.estimate
xi_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Lum/SIRJPF/specific_model/theta_Si/specific_theta_Si.RData")
theta_Si_ll = pf.loglik.of.mif.estimate
theta_Si_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Lum/SIRJPF/specific_model/theta_P/specific_theta_P.RData")
theta_P_ll = pf.loglik.of.mif.estimate
theta_P_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Lum/SIRJPF/specific_model/theta_Ii/specific_theta_Ii.RData")
theta_Ii_ll = pf.loglik.of.mif.estimate
theta_Ii_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Lum/SIRJPF/specific_model/ri/specific_ri.RData")
ri_rn_ll = pf.loglik.of.mif.estimate
ri_rn_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Lum/SIRJPF/specific_model/probi/specific_probi.RData")
probi_ll = pf.loglik.of.mif.estimate
probi_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Lum/SIRJPF/specific_model/f_Si/specific_f_Si.RData")
f_Si_ll = pf.loglik.of.mif.estimate
f_Si_ll_se = s.e.of.pf.loglik.of.mif.estimate



# -------------------------------------Specific Block------------------------------------------------
load("./Single-species/Lum/SIRJPF/specific_model_block/xi/specific_xi.RData")
xi_block_ll = pf.loglik.of.mif.estimate
xi_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Lum/SIRJPF/specific_model_block/theta_Si/specific_theta_Si.RData")
theta_Si_block_ll = pf.loglik.of.mif.estimate
theta_Si_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Lum/SIRJPF/specific_model_block/theta_P/specific_theta_P.RData")
theta_P_block_ll = pf.loglik.of.mif.estimate
theta_P_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Lum/SIRJPF/specific_model_block/theta_Ii/specific_theta_Ii.RData")
theta_Ii_block_ll = pf.loglik.of.mif.estimate
theta_Ii_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Lum/SIRJPF/specific_model_block/ri/specific_ri.RData")
ri_block_ll = pf.loglik.of.mif.estimate
ri_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Lum/SIRJPF/specific_model_block/probi/specific_probi.RData")
probi_block_ll = pf.loglik.of.mif.estimate
probi_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Lum/SIRJPF/specific_model_block/f_Si/specific_f_Si.RData")
f_Si_block_ll = pf.loglik.of.mif.estimate
f_Si_block_ll_se = s.e.of.pf.loglik.of.mif.estimate


K_all = length(parameters) - 1
K_two = length(parameters) - 3 + 2 * 9
K_one = length(parameters) - 2 + 1 * 9


lum_para_ll_df <- list(
  all_shared = list(ll = all_shared_ll, ll_se = all_shared_ll_se, AIC = -2 * all_shared_ll + 2 * K_all),
  xi = list(ll = xi_ll, ll_se = xi_ll_se, AIC = -2 * xi_ll + 2 * K_one),
  theta_Si = list(ll = theta_Si_ll, ll_se = theta_Si_ll_se, AIC = -2 * theta_Si_ll + 2 * K_one),
  theta_P = list(ll = theta_P_ll, ll_se = theta_P_ll_se, AIC = -2 * theta_P_ll + 2 * K_one),
  theta_Ii = list(ll = theta_Ii_ll, ll_se = theta_Ii_ll_se, AIC = -2 * theta_Ii_ll + 2 * K_one),
  ri = list(ll = ri_rn_ll, ll_se = ri_rn_ll_se, AIC = -2 * ri_rn_ll + 2 * K_one),
  probi = list(ll = probi_ll, ll_se = probi_ll_se, AIC = -2 * probi_ll + 2 * K_one),
  f_Si = list(ll = f_Si_ll, ll_se = f_Si_ll_se, AIC = -2 * f_Si_ll + 2 * K_one),
  xi_block = list(ll = xi_block_ll, ll_se = xi_block_ll_se, AIC = -2 * xi_block_ll + 2 * K_one),
  theta_Si_block = list(ll = theta_Si_block_ll, ll_se = theta_Si_block_ll_se, AIC = -2 * theta_Si_block_ll + 2 * K_one),
  theta_P_block = list(ll = theta_P_block_ll, ll_se = theta_P_block_ll_se, AIC = -2 * theta_P_block_ll + 2 * K_one),
  theta_Ii_block = list(ll = theta_Ii_block_ll, ll_se = theta_Ii_block_ll_se, AIC = -2 * theta_Ii_block_ll + 2 * K_one),
  ri_block = list(ll = ri_block_ll, ll_se = ri_block_ll_se, AIC = -2 * ri_block_ll + 2 * K_one),
  probi_block = list(ll = probi_block_ll, ll_se = probi_block_ll_se, AIC = -2 * probi_block_ll + 2 * K_one),
  f_Si_block = list(ll = f_Si_block_ll, ll_se = f_Si_block_ll_se, AIC = -2 * f_Si_block_ll + 2 * K_one)
)

lum_para_parameter_table <- data.frame(
  ll = sapply(lum_para_ll_df[!grepl("_block", names(lum_para_ll_df))], function(x) x$ll),
  AIC = sapply(lum_para_ll_df[!grepl("_block", names(lum_para_ll_df))], function(x) x$AIC),
  block_ll = sapply(names(lum_para_ll_df[!grepl("_block", names(lum_para_ll_df))]), function(name) {
    block_name <- paste0(name, "_block")
    if (block_name %in% names(lum_para_ll_df)) lum_para_ll_df[[block_name]]$ll else NA
  }),
  block_AIC = sapply(names(lum_para_ll_df[!grepl("_block", names(lum_para_ll_df))]), function(name) {
    block_name <- paste0(name, "_block")
    if (block_name %in% names(lum_para_ll_df)) lum_para_ll_df[[block_name]]$AIC else NA
  }),
  row.names = names(lum_para_ll_df[!grepl("_block", names(lum_para_ll_df))])
)

lum_para_parameter_table$max_AIC <- pmax(lum_para_parameter_table$AIC, lum_para_parameter_table$block_AIC, na.rm = TRUE)
lum_para_parameter_table <- lum_para_parameter_table[order(lum_para_parameter_table$max_AIC), ]
lum_para_parameter_table$max_AIC <- NULL
print(lum_para_parameter_table)

save(lum_para_ll_df = lum_para_ll_df,lum_para_parameter_table = lum_para_parameter_table,file = './data/Simple_dynamics/Lum/para/Lum_para_loglik_df.rds')






#######################################################################Lum no para model########################################################
load("./Single-species/Lum/SRJF/model/best_result.rda")
all_shared_ll = pf.loglik.of.mif.estimate
all_shared_ll_se = s.e.of.pf.loglik.of.mif.estimate
parameters = mif.estimate

# -------------------------------------Specific------------------------------------------------
load("./Single-species/Lum/SRJF/specific_model/theta_Si/specific_theta_Si.RData")
theta_Si_ll = pf.loglik.of.mif.estimate
theta_Si_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Lum/SRJF/specific_model/ri/specific_ri.RData")
ri_rn_ll = pf.loglik.of.mif.estimate
ri_rn_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Lum/SRJF/specific_model/f_Si/specific_f_Si.RData")
f_Si_ll = pf.loglik.of.mif.estimate
f_Si_ll_se = s.e.of.pf.loglik.of.mif.estimate



# -------------------------------------Specific Block------------------------------------------------
load("./Single-species/Lum/SRJF/specific_model_block/theta_Si/specific_theta_Si.RData")
theta_Si_block_ll = pf.loglik.of.mif.estimate
theta_Si_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Lum/SRJF/specific_model_block/ri/specific_ri.RData")
ri_block_ll = pf.loglik.of.mif.estimate
ri_block_ll_se = s.e.of.pf.loglik.of.mif.estimate

load("./Single-species/Lum/SRJF/specific_model_block/f_Si/specific_f_Si.RData")
f_Si_block_ll = pf.loglik.of.mif.estimate
f_Si_block_ll_se = s.e.of.pf.loglik.of.mif.estimate


K_all = length(parameters) - 1
K_two = length(parameters) - 3 + 2 * 10
K_one = length(parameters) - 2 + 1 * 10


lum_no_para_ll_df <- list(
  all_shared = list(ll = all_shared_ll, ll_se = all_shared_ll_se, AIC = -2 * all_shared_ll + 2 * K_all),
  theta_Si = list(ll = theta_Si_ll, ll_se = theta_Si_ll_se, AIC = -2 * theta_Si_ll + 2 * K_one),
  ri = list(ll = ri_rn_ll, ll_se = ri_rn_ll_se, AIC = -2 * ri_rn_ll + 2 * K_one),
  f_Si = list(ll = f_Si_ll, ll_se = f_Si_ll_se, AIC = -2 * f_Si_ll + 2 * K_one),
  theta_Si_block = list(ll = theta_Si_block_ll, ll_se = theta_Si_block_ll_se, AIC = -2 * theta_Si_block_ll + 2 * K_one),
  ri_block = list(ll = ri_block_ll, ll_se = ri_block_ll_se, AIC = -2 * ri_block_ll + 2 * K_one),
  f_Si_block = list(ll = f_Si_block_ll, ll_se = f_Si_block_ll_se, AIC = -2 * f_Si_block_ll + 2 * K_one)
)

lum_no_para_parameter_table <- data.frame(
  ll = sapply(lum_no_para_ll_df[!grepl("_block", names(lum_no_para_ll_df))], function(x) x$ll),
  AIC = sapply(lum_no_para_ll_df[!grepl("_block", names(lum_no_para_ll_df))], function(x) x$AIC),
  block_ll = sapply(names(lum_no_para_ll_df[!grepl("_block", names(lum_no_para_ll_df))]), function(name) {
    block_name <- paste0(name, "_block")
    if (block_name %in% names(lum_no_para_ll_df)) lum_no_para_ll_df[[block_name]]$ll else NA
  }),
  block_AIC = sapply(names(lum_no_para_ll_df[!grepl("_block", names(lum_no_para_ll_df))]), function(name) {
    block_name <- paste0(name, "_block")
    if (block_name %in% names(lum_no_para_ll_df)) lum_no_para_ll_df[[block_name]]$AIC else NA
  }),
  row.names = names(lum_no_para_ll_df[!grepl("_block", names(lum_no_para_ll_df))])
)

lum_no_para_parameter_table$max_AIC <- pmax(lum_no_para_parameter_table$AIC, lum_no_para_parameter_table$block_AIC, na.rm = TRUE)
lum_no_para_parameter_table <- lum_no_para_parameter_table[order(lum_no_para_parameter_table$max_AIC), ]
lum_no_para_parameter_table$max_AIC <- NULL
print(lum_no_para_parameter_table)

save(lum_no_para_ll_df = lum_no_para_ll_df,lum_no_para_parameter_table = lum_no_para_parameter_table,file = './data/Simple_dynamics/Lum/no_para/Lum_no_para_loglik_df.rds')

