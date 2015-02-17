source('Scripts/aDDM_Simulation/aDDM_optimization_fake_drift_1000.R')
source('Scripts/aDDM_Simulation/aDDM_optimization_fake_drift_2000.R')
source('Scripts/aDDM_Simulation/aDDM_optimization_fake_drift_3000.R')
source('Scripts/aDDM_Simulation/aDDM_optimization_fake_var_1000.R')
source('Scripts/aDDM_Simulation/aDDM_optimization_fake_var_2000.R')
source('Scripts/aDDM_Simulation/aDDM_optimization_fake_var_3000.R')
source('Scripts/aDDM_Simulation/aDDM_optimization_fake_theta_1000.R')
source('Scripts/aDDM_Simulation/aDDM_optimization_fake_theta_2000.R')
source('Scripts/aDDM_Simulation/aDDM_optimization_fake_theta_3000.R')
source('Scripts/aDDM_Simulation/aDDM_optimization_fake_drift_var.R')
source('Scripts/aDDM_Simulation/aDDM_optimization_fake_drift_theta.R')
source('Scripts/aDDM_Simulation/aDDM_optimization_fake_var_theta.R')

core.aDDM.fake.opti.drift.1000 = aDDM.optimization.fake.drift.1000(c(4,6,8),4)
core.aDDM.fake.opti.drift.2000 = aDDM.optimization.fake.drift.2000(c(4,6,8),4)
core.aDDM.fake.opti.drift.3000 = aDDM.optimization.fake.drift.3000(c(4,6,8),4)
core.aDDM.fake.opti.var.1000 = aDDM.optimization.fake.var.1000(c(4,6,8),4)
core.aDDM.fake.opti.var.2000 = aDDM.optimization.fake.var.2000(c(4,6,8),4)
core.aDDM.fake.opti.var.3000 = aDDM.optimization.fake.var.3000(c(4,6,8),4)
core.aDDM.fake.opti.theta.1000 = aDDM.optimization.fake.theta.1000(c(4,6,8),4)
core.aDDM.fake.opti.theta.2000 = aDDM.optimization.fake.theta.2000(c(4,6,8),4)
core.aDDM.fake.opti.theta.3000 = aDDM.optimization.fake.theta.3000(c(4,6,8),4)
core.aDDM.fake.opti.drift.var.1000 = aDDM.optimization.fake.drift.var(c(4,6,8),4)
core.aDDM.fake.opti.drift.theta.1000 = aDDM.optimization.fake.drift.theta(c(4,6,8),4)
core.aDDM.fake.opti.var.theta.1000 = aDDM.optimization.fake.var.theta(c(4,6,8),4)

# Full optimization

# aDDM fake
source('Scripts/aDDM_Simulation/aDDM_fake_optimization.R')
core.aDDM.fake.opti = aDDM.optimization.fake(c(4,6,8),4)

# aDDM
source('Scripts/aDDM_Simulation/aDDM_optimization.R')
core.aDDM.opti = aDDM.optimization(c(4,6,8),4)

# DDM
source('Scripts/DDM_Simulation/DDM_optimization.R')
core.DDM.opti = DDM.optimization(c(4,6,8),4)




