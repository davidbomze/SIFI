This is the code to calculate the SIFI using the 'myeloid' and 'retinopathy' sample datasets from the 'survival' package

```R
library(survival)
source("../SIFI.R")

# Example 1: 'myeloid' dataset
sifi(myeloid[, c("futime","death","trt")], plot_iteration = T, file_iteration = "myeloid_sifi_default.pdf")    # Default SIFI strategy
sifi_all(myeloid[ , c("futime","death","trt")], plot_iteration = T, file_prefix = "myeloid_sifi")      # All SIFI strategies

# Example 2: 'retinopathy' dataset
sifi(retinopathy[ , c("futime","status","laser")], treatment_arm = "argon", plot_iteration = T, file_iteration = "retinopathy_sifi_default_argon.pdf")    # Default SIFI strategy, define experimental arm
sifi(retinopathy[ , c("futime","status","laser")], treatment_arm = "xenon", plot_iteration = T, file_iteration = "retinopathy_sifi_default_xenon.pdf")    # Default SIFI strategy, define experimental arm
sifi(retinopathy[ , c("futime","status","laser")], plot_iteration = T, file_iteration = "retinopathy_sifi_default_agnostic.pdf")    # Default SIFI strategy, without defining experimental arm ('agnostically' choose 'xenon' as experimental in this case)
sifi_all(retinopathy[ , c("futime","status","laser")], plot_iteration = T, file_prefix = "retinopathy_sifi")    # All SIFI strategies
```

The figures generated can be found in this folder
