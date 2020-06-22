getOddsRatio = function(M, levels=c())
{

  # M = data_table

  if (is_empty(levels))
  {
    # print("empty levels")
    var_levels = colnames(M)
  } else
  {
    var_levels = levels
  }

  # var_levels=c("Black", "Other", "White")
  # print(var_levels)

  results = data.frame(
    comparison = character(),
    OR = numeric(),
    SE = numeric(),
    CI_low = numeric(),
    CI_high = numeric(),
    cvlog = numeric(),
    p_value = numeric()
  )

  reference = var_levels[1]
  comparisons = var_levels %>% setdiff(reference)

  for (i in 2:length(var_levels))
  {
    # print(i)
    var_to_compare = var_levels[i]
    # print(var_to_compare)
    M2x2 = M[,var_levels[c(1, i)]]
    # print(M2x2)

    odds_ratio = (M2x2[1,1]*M2x2[2,2])/(M2x2[1,2]*M2x2[2,1])

    se_log_or = sqrt(sum(1/M2x2))
    SE = exp(se_log_or)
    CI_high = odds_ratio*exp(1.96*se_log_or)
    CI_low = odds_ratio*exp(-1.96*se_log_or)
    cvlog = abs(se_log_or/log(odds_ratio))

    comparison = paste0(var_to_compare, '_vs_', reference)

    results %<>% add_row(
      comparison = comparison,
      OR = odds_ratio,
      SE = SE,
      CI_high = CI_high,
      CI_low = CI_low,
      cvlog = cvlog,
      p_value = chisq.test(M2x2, simulate.p.value = T)$p.value
    )
  }

  return(results)
}



doUniverariateTests = function(
  DV,
  IVs,
  data,
  pvalue_cutoff
)
{
  # DV = "Vital_Status"
  # IVs = colnames(data) %>% setdiff(DV) %>% unname()
  # print(IVs)

  # DV=DV
  # IVs=IVs
  # data=myopia

  univariate_stats = data.frame(
    IV = character(),
    Comparison = character(),
    Test = character(),
    Pvalue = numeric(),
    OR = numeric(),
    Abs_Diff = numeric(),
    SE = numeric(),
    CV = numeric(),
    has_zero_cell = logical(),
    has_exp_lt_5 = logical()
  )


  for (IV in IVs)
  {

    # print(IV)

    iv_data = data %>% pull(IV)

    is_factor = is.factor(iv_data)

    if(is_factor)
    {
      # print(sprintf("%s -- Factor", IV))
      ### Chi-square, Fishers, check for 0 cells, etc.

      data_table = data %>% select(DV, IV) %>% table()
      # print(data_table)
      cstest = chisq.test(data_table, simulate.p.value = T)
      or = getOddsRatio(data_table)

      ### Checks
      has_zero_cell = any(data_table == 0)
      has_exp_lt_5 = any(cstest$expected < 5)

      comparison = paste0(levels(DV)[2], '_vs_', levels(DV)[1])

      univariate_stats %<>%
        add_row(
          IV = IV,
          Comparison = or$comparison,
          Test = "Chi-Square-Wald",
          Pvalue = or$p_value,
          OR = or$OR,
          has_exp_lt_5 = has_exp_lt_5,
          has_zero_cell = has_zero_cell,
          SE = or$SE,
          # CI_low = or$CI_low,
          # CI_high = or$CI_high,
          CV = or$cvlog
        )

    } else
    {
      # print(sprintf("%s -- Continuous", IV))
      ### t=test vs DV

      formla = as.formula(paste0(IV, " ~ ", DV))
      t_test = t.test(formla, data)
      abs_diff = abs(t_test$estimate[[2]] - t_test$estimate[[1]])
      conf_int = t_test$conf.int %>% as.numeric()
      conf_width = abs(diff(conf_int))
      t_95 = qt(0.975, t_test$parameter)
      SE = conf_width/2/t_95

      univariate_stats %<>%
        add_row(
          IV = IV,
          Test = "T-Test",
          Pvalue = t_test$p.value,
          Abs_Diff = abs_diff,
          SE = t_test$stderr,
          # CI_low = t_test$conf.int[[1]],
          # CI_high = t_test$conf.int[[2]],
          CV = SE/abs_diff
        )
    }
  }

  print(univariate_stats)

  univariate_eliminated =
    univariate_stats %>%
    filter(Pvalue >= pvalue_cutoff) %>%
    mutate(as.character(IV)) %>%
    pull(IV) %>%
    droplevels() %>%
    unname() %>%
    as.character() %>%
    unique()

  # print(univariate_eliminated)

  univariate_kept = setdiff(IVs, univariate_eliminated) %>% unique()

  # print(univariate_kept)

  print(sprintf("All Variables: %s", IVs))
  print(sprintf("Eliminated Variables: %s", univariate_eliminated))
  print(sprintf("Kept Variables: %s", univariate_kept))

  return(list(
    All=IVs,
    Eliminated=univariate_eliminated,
    Kept=univariate_kept,
    UnivariateStats=univariate_stats
  ))
}


testPctDeltaBeta = function(
  full_model,
  reduced_model
)
{

  full_model_vars = full_model$terms %>% attr("term.labels")
  print(full_model_vars)

  full_model_results =
    tidy(full_model) %>%
    data.frame() %>%
    filter(term != "(Intercept)") %>%
    mutate(IV = full_model_vars) %>%
    rename(beta = estimate) %>%
    select(IV, everything()) %>%
    arrange(-p.value)


  reduced_model_vars = reduced_model$terms %>% attr("term.labels")
  print(reduced_model_vars)

  reduced_model_results =
    tidy(reduced_model) %>%
    data.frame() %>%
    filter(term != "(Intercept)") %>%
    mutate(IV = reduced_model_vars) %>%
    rename(beta = estimate) %>%
    select(IV, everything()) %>%
    arrange(-p.value)

  test_var = setdiff(full_model_vars, reduced_model_vars)
  print(test_var)

  print("*****    Calculating delta beta pct    *****")
  pct_delta_beta =
    inner_join(
      reduced_model_results,
      full_model_results,
      by="IV",
      suffix=c(".reduced", ".full")
    ) %>%
    mutate(pct_delta_beta = (abs((beta.full - beta.reduced)/beta.full))) %>%
    select(IV, beta.full, beta.reduced, pct_delta_beta)

  print(pct_delta_beta)

  any_pct_diff_beta_gt_cutoff = any(pct_delta_beta$pct_delta_beta > 0.2)
  if (any_pct_diff_beta_gt_cutoff)
  {
    print(sprintf("*****    Removing the variable \"%s\" changes the beta for at least one variable by more than 20%%", test_var))
    pct_delta_betas_gt_20 = F
  } else
  {
    print(sprintf("*****    Removing the variable \"%s\" does not change  the beta for at least one variable by more than 20%%", test_var))
    pct_delta_betas_gt_20 = T
  }

  return(any_pct_diff_beta_gt_cutoff)
}


doLRT = function(
  full_model,
  reduced_model,
  lrt_pval_cutoff = 0.1,
  data
)
{

  full_model_vars = full_model$terms %>% attr("term.labels")
  reduced_model_vars = reduced_model$terms %>% attr("term.labels")

  test_var = setdiff(full_model_vars, reduced_model_vars)

  print("*****    Performing LRT on full and reduced models    *****")
  print("*****    Testing H_0 that the two models are not significantly different    *****")

  lrt = anova(full_model, reduced_model, test="LRT")
  lrt_pvalue = lrt[2, "Pr(>Chi)"]

  print(sprintf("*****    LRT pvalue=%0.2f    *****", lrt_pvalue))

  if (lrt_pvalue > lrt_pval_cutoff)
  {
    print("*****    LRT is not significant    *****")
    print("Accept H_0")
    reject_H_0 = F
  } else
  {
    print("*****    LRT is significant    *****")
    print("Reject H_0")
    reject_H_0 = T
  }

  return(reject_H_0)
}



reduceModel = function(
  DV,
  IVs,
  data,
  model_pvalue_cutoff=0.05,
  model_comparison_pval_cutoff=0.1
)
{

  # DV = "Myopia_in_5_Years"
  # IVs = step_1_kept
  # data = myopia
  # model_pvalue_cutoff=0.05
  # model_comparison_pval_cutoff=0.1
  # model_pvalue_cutoff = 0.05
  # lrt_pval_cutoff = 0.1



  print("*****    Starting model reduction    *****")
  print(sprintf("starting variables: %s", paste(IVs, collapse=", ")))


  ### Build Full Model
  full_model_IVs = IVs

  print("*****    Building Full Model    *****")
  print(sprintf("Current full model variables: %s", paste(full_model_IVs, collapse=", ")))

  full_model_formula_string =
    paste(IVs, collapse=" + ") %>%
    paste0(DV, " ~ ", .)

  full_model_formula = as.formula(full_model_formula_string)
  print("*****    Full Model Formula    *****")
  print(full_model_formula)

  full_model =
    glm(formula=full_model_formula, data, family=binomial)

  print("*****    Full model results    *****")
  full_model_results =
    tidy(full_model) %>%
    data.frame() %>%
    filter(term != "(Intercept)") %>%
    mutate(IV = full_model_IVs) %>%
    rename(beta = estimate) %>%
    select(IV, everything()) %>%
    arrange(-p.value)

  print("*****    Full model results    *****")
  print(full_model_results)

  ### Determine which variables are not significant
  vars_to_test =
    full_model_results %>%
    filter(p.value > model_pvalue_cutoff) %>%
    arrange(-p.value) %>%
    pull(IV)

  print("*****    The following variables will be tested:    *****")
  print(vars_to_test)

  ### Repository for variables that will not be kept
  vars_to_remove = c()

  ### Test the variables 1 by 1
  for (var_to_test in vars_to_test)
  {
    print(sprintf("*****    Testing: %s", var_to_test))
    print("*****    Building Reduced Model    *****")

    ### Determine reduced model IVs
    reduced_model_IVs = setdiff(full_model_IVs, var_to_test)
    print(sprintf("Reduced model variables: %s", paste(reduced_model_IVs, collapse=", ")))

    ### Build reduced model
    reduced_model_formula_string =
      paste(reduced_model_IVs, collapse=" + ") %>%
      paste0(DV, " ~ ", .)

    reduced_model_formula = as.formula(reduced_model_formula_string)
    print("*****    Reduced Model Formula    *****")
    print(reduced_model_formula)

    reduced_model =
      glm(formula=reduced_model_formula, data, family=binomial)

    print("*****    Reduced model results    *****")
    reduced_model_results =
      tidy(reduced_model) %>%
      data.frame() %>%
      filter(term != "(Intercept)") %>%
      mutate(IV = reduced_model_IVs) %>%
      rename(beta = estimate) %>%
      select(IV, everything()) %>%
      arrange(-p.value)

    print(reduced_model_results)

    ### Do LRT. If TRUE, then reject H0 and keep the variable because removing
    ### it changes the model
    keep_var = doLRT(
      full_model,
      reduced_model,
      lrt_pval_cutoff = model_comparison_pval_cutoff,
      data=data
    )

    ### If accept H0, then still may want to keep if changes betas a lot
    if (!keep_var)
    {
      keep_var = testPctDeltaBeta(full_model, reduced_model)
    }

    ### If the variable is still going to be removed, then add it to the
    ### list of variables to remove
    if (!keep_var)
    {
      print(sprintf("*****    REMOVING VAR: %s    *****", var_to_test))
      vars_to_remove = c(vars_to_remove, var_to_test)
    } else
    {
      print(sprintf("*****    KEEPING VAR: %s    *****", var_to_test))
    }
  }


  vars_to_keep = setdiff(IVs, vars_to_remove)

  reduced_model_formula_string =
    paste(vars_to_keep, collapse=" + ") %>%
    paste0(DV, " ~ ", .)

  reduced_model_formula = as.formula(reduced_model_formula_string)

  reduced_model =
    glm(formula=reduced_model_formula, data, family=binomial)

  return(list(
    all_vars = IVs,
    vars_to_keep = vars_to_keep,
    vars_to_remove = vars_to_remove,
    reduced_model = reduced_model
  ))

}



testExcludedVariables = function(
  DV,
  IVs,
  test_IVs,
  data,
  pvalue_cutoff
)
{

  # DV = "Myopia_in_5_Years"
  # IVs = reduced_results$vars_to_keep
  # test_IVs = step_1_hold
  # data = myopia

  base_formula_string =
    paste0(IVs, collapse=" + ") %>%
    paste(DV, "~", .)
  print(base_formula_string)

  individual_IV_test_stats = data.frame(
    test_IV = character(),
    Comparison = character(),
    Pvalue = numeric(),
    OR = numeric(),
    SE = numeric(),
    has_zero_cell = logical(),
    has_exp_lt_5 = logical()
  )

  for (test_IV in test_IVs)
  {

    print(test_IV)

    test_IV_data = data %>% pull(test_IV)

    is_factor = is.factor(test_IV_data)

    formla = as.formula(paste(base_formula_string, "+", test_IV))
    test_model = glm(formula=formla, data=data, family=binomial)

    test_model_results = tidy(test_model)


    if(is_factor)
    {
      print(sprintf("%s -- Factor", test_IV))
      ### Chi-square, Fishers, check for 0 cells, etc.

      data_table = data %>% select(DV, test_IV) %>% table()
      # print(data_table)

      ### Checks
      cstest = chisq.test(data_table, simulate.p.value = T)
      has_zero_cell = any(data_table == 0)
      has_exp_lt_5 = any(cstest$expected < 5)

      all_levels = levels(test_IV_data)
      num_levels = length(all_levels)
      ref_level = all_levels[1]

      for (level in all_levels[2:num_levels])
      {
        print(sprintf("level %s of %s", level, toString(all_levels)))

        comparison = paste0(level, '_vs_', ref_level)
        term = paste0(test_IV, level)

        regression_data = test_model_results %>% filter(term==!!term)

        individual_IV_test_stats %<>%
          add_row(
            test_IV = test_IV,
            Comparison = comparison,
            Pvalue = regression_data$p.value,
            OR = exp(regression_data$estimate),
            has_exp_lt_5 = has_exp_lt_5,
            has_zero_cell = has_zero_cell,
            SE = exp(regression_data$std.error)
            # CI_low = or$CI_low,
            # CI_high = or$CI_high,
          )

      }
    } else
    {
      # print(sprintf("%s -- Continuous", IV))
      ### t=test vs DV

      term = test_IV
      regression_data = test_model_results %>% filter(term==!!term)

      individual_IV_test_stats %<>%
        add_row(
          test_IV = test_IV,
          Comparison = term,
          Pvalue = regression_data$p.value,
          OR = exp(regression_data$estimate),
          SE = exp(regression_data$std.error)
          # CI_low = or$CI_low,
          # CI_high = or$CI_high,
        )


    }
  }

  print(individual_IV_test_stats)

  vars_to_add =
    individual_IV_test_stats %>%
    filter(Pvalue < pvalue_cutoff) %>%
    pull(test_IV) %>%
    as.character()

  return(list(
    vars_to_add = vars_to_add,
    vars_kept_out = setdiff(test_IVs, vars_to_add),
    all_vars = union(IVs, vars_to_add)
  ))
}


