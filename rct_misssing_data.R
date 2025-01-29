# Multiple imputation in R using the mice library

# Install packages:----
# Only need to do this once on your machine
# install.packages("tidyverse")
# install.packages("mice")

# Load packages:----
# load package (needs to be done for every new R session)
library(tidyverse)
# library(purrr)
# library(mice)

# Read and view the data:----
data_ <- read_csv("~/R/rct_data.csv")
View(data_) # open up data viewer to check dataset looks reasonable

# Split the data based on treatment group and PROM:----
## Treatment-split data:----
data_t1 <- data_ %>%
  dplyr::as_tibble() %>%
  dplyr::filter(treatment == 1)
data_t2 <- data_ %>%
  dplyr::as_tibble() %>%
  dplyr::filter(treatment == 2)
## EQ5D data:----
EQ5D_df <- data_ %>%
  dplyr::as_tibble() %>%
  dplyr::select(-dplyr::contains("SF6D"))
EQ5D_df_t1 <- data_t1 %>%
  dplyr::as_tibble() %>%
  dplyr::select(-dplyr::contains("SF6D"))
EQ5D_df_t2 <- data_t2 %>%
  dplyr::as_tibble() %>%
  dplyr::select(-dplyr::contains("SF6D"))
## SF6D data:----
SF6D_df <- data_ %>%
  dplyr::as_tibble() %>%
  dplyr::select(-dplyr::contains("EQ5D"))
SF6D_df_t1 <- data_t1 %>%
  dplyr::as_tibble() %>%
  dplyr::select(-dplyr::contains("EQ5D"))
SF6D_df_t2 <- data_t2 %>%
  dplyr::as_tibble() %>%
  dplyr::select(-dplyr::contains("EQ5D"))

# Investigate missingness pattern:----
# it is great to add these plots to the assignment report:----
## Treatment 1 plot:----
mice::md.pattern(
  x = data_t1 %>%
    dplyr::select(-studyid, -treatment),
  rotate.names = TRUE
)
## Treatment 2 plot:----
mice::md.pattern(
  x = data_t2 %>%
    dplyr::select(-studyid, -treatment),
  rotate.names = TRUE
)

# MICE:----
## Set the number of imputations:----
imp_number_ <- 2
## Set the seed number:----
seed_no_ <- 1234
## Impute missing data using multiple imputation with chain equation:----
### Treatment 1 dataset:----
#### Run MICE:----
imp_data_ls_t1 <- mice::mice(
  data = data_t1 %>%
    dplyr::select(-treatment, -ethnict), # these variables are constants (do not change)
  m = imp_number_,
  seed = seed_no_,
  printFlag = FALSE
)
plot(imp_data_ls_t1)
#### Complete the missing values and arrange the resulting data sets in wide format:----
imp_data_t1 <- mice::complete(
  data = imp_data_ls_t1,
  action = "broad", # mice will put the imputed data sets next to each other
  include = FALSE
)
#### Extract columns with missing values
imp_cols_data_t1 <- data_t1 %>% # look back into the first subset
  dplyr::select(
    tidyselect::vars_select_helpers$where(
      fn = function(x) any(is.na(x))
    )
  ) %>% # select every variable with NA
  colnames() %>% # get the column (variable) names
  purrr::map_dfc( # loop through the column names
    .x = .,
    .f = function(.col_name_) { # for each column
      imp_data_t1 %>% # look into the imputed data set
        dplyr::select(
          dplyr::contains(.col_name_), # keep columns with that name
          dplyr::contains("studyid")
        ) %>% # keep patient id to avoid matching issues
        dplyr::rowwise() %>% # let dplyr know we want the next operations per row
        dplyr::mutate(
          studyid = mean( # all studyid columns have the same values, mean will
            dplyr::c_across( # only keep one column will the exact study ids.
              dplyr::contains("studyid")
            )
          ),
          {{ .col_name_ }} := mean( # {{.col_name :=}} allows R to create column names
            dplyr::c_across( # c.across() collects values from the named columns
              dplyr::contains(.col_name_)
            )
          )
        ) %>%
        {
          if (.col_name_ == data_t1 %>% # look back into the first subset
            dplyr::select(
              tidyselect::vars_select_helpers$where(
                fn = function(.x_) any(is.na(.x_))
              )
            ) %>% # select every variable with NA
            colnames() %>%
            .[1]) { # If .col_name_ was the first column with NA in the data set
            dplyr::select( # Keep both studyid and that column
              .data = ., # In if statements preceeded with %>%, we must pass "." to function
              studyid, {{ .col_name_ }}
            )
          } else { # otherwise
            dplyr::select( # Only keep that column as studyid was selected first column
              .data = .,
              {{ .col_name_ }}
            )
          }
        } %>%
        dplyr::ungroup() # disable row-wise
    }
  )
#### Finally, put imputed data together with complete data:----
completed_data_t1 <- data_t1 %>%
  dplyr::select( # select columns with complete data
    tidyselect::vars_select_helpers$where(fn = function(.x_) !any(is.na(.x_)))
  ) %>%
  dplyr::left_join( # join completed columns with the ones with complete data
    x = .,
    y = imp_cols_data_t1,
    by = "studyid"
  ) %>%
  dplyr::select(
    colnames(data_t1)
  ) # get the same order of column as the original data set
View(completed_data_t1) # open up data viewer to check dataset looks reasonable
#### Split completed data set by utility:----
completed_SF6D_df_t1 <- completed_data_t1 %>%
  dplyr::as_tibble() %>%
  dplyr::select(-dplyr::contains("EQ5D"))
completed_EQ5D_df_t1 <- completed_data_t1 %>%
  dplyr::as_tibble() %>%
  dplyr::select(-dplyr::contains("SF6D"))

### Treatment 2 dataset:----
#### Run MICE:----
imp_data_ls_t2 <- mice::mice(
  data = data_t2 %>%
    dplyr::select(-treatment, -ethnict), # these variables are constants (do not change)
  m = imp_number_,
  seed = seed_no_,
  printFlag = FALSE
)
plot(imp_data_ls_t2)
#### Complete the missing values and arrange the resulting data sets in wide format:----
imp_data_t2 <- mice::complete(
  data = imp_data_ls_t2,
  action = "broad", # mice will put the imputed data sets next to each other
  include = FALSE
)
#### Extract columns with missing values
imp_cols_data_t2 <- data_t2 %>% # look back into the first subset
  dplyr::select(
    tidyselect::vars_select_helpers$where(
      fn = function(x) any(is.na(x))
    )
  ) %>% # select every variable with NA
  colnames() %>% # get the column (variable) names
  purrr::map_dfc( # loop through the column names
    .x = .,
    .f = function(.col_name_) { # for each column
      imp_data_t2 %>% # look into the imputed data set
        dplyr::select(
          dplyr::contains(.col_name_), # keep columns with that name
          dplyr::contains("studyid")
        ) %>% # keep patient id to avoid matching issues
        dplyr::rowwise() %>% # let dplyr know we want the next operations per row
        dplyr::mutate(
          studyid = mean( # all studyid columns have the same values, mean will
            dplyr::c_across( # only keep one column will the exact study ids.
              dplyr::contains("studyid")
            )
          ),
          {{ .col_name_ }} := mean( # {{.col_name :=}} allows R to create column names
            dplyr::c_across( # c.across() collects values from the named columns
              dplyr::contains(.col_name_)
            )
          )
        ) %>%
        {
          if (.col_name_ == data_t2 %>% # look back into the first subset
            dplyr::select(
              tidyselect::vars_select_helpers$where(
                fn = function(.x_) any(is.na(.x_))
              )
            ) %>% # select every variable with NA
            colnames() %>%
            .[1]) { # If .col_name_ was the first column with NA in the data set
            dplyr::select( # Keep both studyid and that column
              .data = ., # In if statements preceeded with %>%, we must pass "." to function
              studyid, {{ .col_name_ }}
            )
          } else { # otherwise
            dplyr::select( # Only keep that column as studyid was selected first column
              .data = .,
              {{ .col_name_ }}
            )
          }
        } %>%
        dplyr::ungroup() # disable row-wise
    }
  )
#### Finally, put imputed data together with complete data:----
completed_data_t2 <- data_t2 %>%
  dplyr::select( # select columns with complete data
    tidyselect::vars_select_helpers$where(fn = function(.x_) !any(is.na(.x_)))
  ) %>%
  dplyr::left_join( # join completed columns with the ones with complete data
    x = .,
    y = imp_cols_data_t2,
    by = "studyid"
  ) %>%
  dplyr::select(
    colnames(data_t2)
  ) # get the same order of column as the original data set
View(completed_data_t2) # open up data viewer to check dataset looks reasonable
#### Split completed data set by utility:----
completed_SF6D_df_t2 <- completed_data_t2 %>%
  dplyr::as_tibble() %>%
  dplyr::select(-dplyr::contains("EQ5D"))
completed_EQ5D_df_t2 <- completed_data_t2 %>%
  dplyr::as_tibble() %>%
  dplyr::select(-dplyr::contains("SF6D"))

# #Average imputed data set for bootstrapping
# impute_mean <- rep(0, length(complete(impute.extra, 1)))
# for(i in 1:80){impute_mean <- impute_mean + complete(impute.extra, i)}
# impute_mean <- impute_mean/80
#
# #Bootstrap CI
# library(boot) # load the boot library
#
# # define function to be used for bootstrapping (taking the mean of TxSurv and NonTxSurv)
# boot_function1 <- function(dat, d) {
#   E <- dat[d,] # allows boot to select sample
#   return(c(mean(E$QALY)))
# }
#
# Bootstrap_impute <- boot(impute_mean, boot_function1, R=5000)
# Bootstrap_impute
#
#
# #Bootstrap CI
# boot.ci(Bootstrap_impute, type=c("norm", "perc", "bca"), index=1) #
