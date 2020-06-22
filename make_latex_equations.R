require(tidyverse)
require(magrittr)

makeLatexRegressionEquation = function(
  model,
  max_var_width = 4,
  true_coefficients = F,
  sig_digits = 3
)
{
  # model = testfit
  # true_coefficients = T

  IVs = model$terms %>% attr("term.labels")
  DV = model$terms %>% attr("variables") %>% as.character() %>% .[[2]]

  num_terms = length(IVs) + 1
  num_equal_length_terms = floor(num_terms/max_var_width)
  num_extra_terms = num_terms%%max_var_width

  # term_template = "\\beta_{%s} \\cdot {\\small{%s}}"

  term_template = "\\small{%s}"

  if (true_coefficients)
  {
    betas =
      model$coefficients[2:num_terms] %>%
      as.numeric() %>%
      signif(3) %>%
      sapply(function(x) ifelse(x>0, paste0('\\; + \\; ', x), x), USE.NAMES = F) %>%
      paste0(" \\cdot ") %>%
      c(signif(model$coefficients[[1]], 3), .)

  } else
  {
    betas =
      sprintf(' + \\; \\beta_{%s} \\cdot ', 1:(num_terms-1)) %>%
      c('\\beta_0', .)
  }

  eq_terms = sprintf(
    term_template,
    c(IVs)
  ) %>%
    gsub(':', ' \\\\times ', .) %>%
    c("", .) %>%
    paste0(betas, .) %>%
    gsub('-', ' \\\\; - \\\\; ', .)


  LHS = sprintf('\\log \\Big( \\frac{\\pi_{\\scriptsize{%s}}}{1 - \\pi_{\\scriptsize{%s}}} \\Big) & = \\; ', DV, DV)

  RHS = ""

  if (num_equal_length_terms > 0)
  {
    for (i in 1:num_equal_length_terms)
    {
      # i = 1
      low = max_var_width*(i-1) + 1
      high = max_var_width*i

      # eq_line = paste0(eq_terms[low:high], collapse=" + ")
      eq_line = paste0(eq_terms[low:high], collapse="")

      if(i > 1)
      {
        eq_line = paste0(" & ", eq_line)
      }

      if (num_equal_length_terms > 1 | num_extra_terms > 0)
      {
        eq_line =
          eq_line %>%
          paste0(., " \\\\ \n  ")
      }

      RHS = paste0(RHS, eq_line)
    }
  }

  if (num_extra_terms > 0)
  {
    low = num_terms - num_extra_terms + 1
    high = num_terms

    eq_line =
      paste0(eq_terms[low:high], collapse=" + ") %>%
      paste0(' & ', .)
    #   paste0(' & \\; ', .)

    RHS = paste0(RHS, eq_line)
  }


  eq_string = paste0(
    "\\begin{aligned}\n", LHS, RHS, "\n\\end{aligned}")

  return(eq_string)
}
