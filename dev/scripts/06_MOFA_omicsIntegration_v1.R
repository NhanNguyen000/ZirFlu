BiocManager::install("MOFA2")
library(MOFA2)

data <- make_example_data(
  n_views = 2, 
  n_samples = 200, 
  n_features = 1000, 
  n_factors = 10
)[[1]]

lapply(data,dim)

MOFAobject <- create_mofa(data)

N = ncol(data[[1]])
groups = c(rep("A",N/2), rep("B",N/2))

MOFAobject <- create_mofa(data, groups=groups)
print(MOFAobject)
plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)
head(data_opts)

model_opts <- get_default_model_options(MOFAobject)
head(model_opts)

train_opts <- get_default_training_options(MOFAobject)
head(train_opts)

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

outfile = "temp/model.hdf5"

MOFAobject.trained <- run_mofa(MOFAobject, outfile)



BiocManager::install("basilisk")

# NEXT TRY -----------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MOFA2")
