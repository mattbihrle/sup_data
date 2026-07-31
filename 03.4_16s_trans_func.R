# 03.4_16s_trans func


# Load data
load("output/data/mt_16s.RData")
load("output/data/full_16s.RData")
load("output/data/contam_16s.RData")
# Load Libraries 
packages <- c("tidyverse", "microeco")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
# Remove the list of packages and installed packages from the environment
rm(installed_packages, packages)
# Start Analysis ------------------------------------------------

# Create trans func obj

f_mt <- clone(mt_16s)

f_mt$tax_table <- f_mt$tax_table
f1 <- trans_func$new(f_mt)
f1$for_what <- "prok"
f1$cal_func(prok_database = "FAPROTAX")

save(f1, file = "output/data/16s_f1.RData")

for(i in 1:length(colnames(f1$res_func))){
  print(colnames(f1$res_func)[i])
  print(table(f1$res_func[,i]))
}

f1$cal_func_FR(abundance_weighted = TRUE)
f1$cal_func_FR_comm()
f1$plot_func_FR() +
  facet_wrap("deployment")

f1$show_prok_func(use_func = "intracellular_parasites")
