# 03-selecting_actual_hyperparameters_not_auto.R


if (!interactive()){
  stop("Stop!!! Please run with:
       \"R\" then \"source('03-selecting_actual_hyperparameters_not_auto.R', echo = TRUE)\",
       not: 
       \"Rscript 03-...\". We need interactive decisions.")
}


myreadline <- function(prompt){
  readline(prompt=paste0("---------------------------------------\n",
                         prompt))
}

parameters_as_string <- function(x, digits = 4){
  full_string <- ""
  for (l_idx in names(x)){
    if (inherits(x[[l_idx]], "numeric")){
    inner_string <- paste0(l_idx, ": ", signif(x[[l_idx]], digits = digits), "\n")
  } else {
    inner_string <- paste0(l_idx, ": ", x[[l_idx]], "\n")
  }
    full_string <- paste(full_string, inner_string)
  }
  return(full_string)
}

myreadline(paste0("Make sure you did \nsource('03-selecting_actual_hyperparameters_not_auto.R',",
              "echo = TRUE). \n\nOtherwise you won't be able to use the ",
              "interactive plotly code. \n(Enter for Next) "))


# initalization ------------------

suppressMessages(library(tidyverse))
suppressMessages(library(gridExtra))

input_args <- eval(parse(text = myreadline("Input Parameters - e.g.: c(1000, \"default\"): ")))
                         #as.numeric(commandArgs(trailingOnly=TRUE))
num_sims <- as.numeric(input_args[1]) # 1000



input_sigma_info <- input_args[2]
if (input_sigma_info == "default"){
  range_sigma <- paste0(c(1:5, seq(10, 85, by = 5)), "%")
  input_sigma_info_str <- input_sigma_info
} else {
  sigma_range_info <- input_sigma_info %>% stringr::str_split(":{1,2}", simplify = T) %>% 
    as.numeric()
  range_sigma <- paste0(seq(sigma_range_info[1], 
                            sigma_range_info[2], 
                            by = sigma_range_info[3]),
                        "%")
  input_sigma_info_str <- paste0(sigma_range_info, collapse = "-")
}

# select sigma (interactive) ----------------


load(paste0("data/ggvis_",num_sims,"_", 
            input_sigma_info_str, "_3.Rdata"))

## examine visuals:
myreadline(prompt=paste("Please examine the following 3 visuals to get a good\n",
                      "sense of the data and what you'll be using to select your parameters.\n",
                      "\nReady? "))

system2("open", args = paste0("data/vis_all_", num_sims,"_", 
            input_sigma_info_str, "_3.pdf"))
myreadline(prompt="Next? ")
system2("open", args = paste0("data/vis_all_paths_", num_sims,"_", 
                              input_sigma_info_str, "_3.pdf"))
myreadline(prompt="Next? ")
system2("open", args = paste0("data/vis_combo_", num_sims,"_", 
                              input_sigma_info_str, "_3.pdf"))

myreadline(prompt=paste0("You will now use a plotly interactive visual",
                       "\nto select optimal parameters?\nReady? "))


suppressWarnings(plotly::ggplotly((data_combo %>%
                    ggplot() +
                    geom_tile(aes(x = -log10(eps), y = -log10(diff_eps),
                                  fill = factor(correct_num_across),
                                  text = paste0("(log10(eps): ",log10(eps), 
                                               ", log10(diff_eps): ",
                                               log10(diff_eps), ")"))) +
                    labs(x = "-log10(eps)", y = "-log10(diff_eps)",
                         title = "Parameters with correct amount of modes",
                         subtitle = "(sum across X values, NA = False)",
                         fill = "# of correct modes")  +
                    facet_wrap(~.sigma_string) +
                    scale_fill_manual(values = c("grey70",'#fee8c8',
                                                 '#fdbb84','#e34a33'))),
                 tooltip = "text"))


# recommend
.sigma_string <- myreadline("Select sigma (without quotes, e.g. 25%): ")
.my_sigma_string <- .sigma_string
eps <- eval(parse(text = myreadline("Select eps (in form 10^-8): ")))
diff_eps <- eval(parse(text = myreadline("Select diff_eps (in form 10^-8): ")))
maxT <- 200 # should converge earlier

data_combo %>%
  filter(.sigma_string == .my_sigma_string) %>%
  ggplot() +
  geom_tile(aes(x = -log10(eps), y = -log10(diff_eps),
                fill = factor(correct_num_across))) +
  labs(x = "-log10(eps)", y = "-log10(diff_eps)",
       title = "Parameters with correct amount of modes",
       subtitle = paste("(sum across X values, NA = False) Sigma =", 
                        .my_sigma_string),
       fill = "# of correct modes")  +
  scale_fill_manual(values = c("0"= "grey70",
                               "1" = '#fee8c8',
                               "2" = '#fdbb84',
                               "3" = '#e34a33')) +
  geom_hline(data = data.frame(eps = eps,
                               diff_eps = diff_eps),
             aes(yintercept = -log10(diff_eps)),
             color = 'black') +
  geom_vline(data = data.frame(eps = eps,
                               diff_eps = diff_eps),
             aes(xintercept = -log10(eps)),
             color = 'black')




selected_parameters <- list(.sigma_string = .sigma_string,
                            .eps = eps, 
                            .diff_eps = diff_eps,
                            .maxT = maxT)

myreadline(prompt = paste0("Please look over the plot and examine if you selected",
                          " the correct parameters.\n",
                          "**The parameters you selected were:**\n",
                          parameters_as_string(selected_parameters), "Complete (Enter)? "))


save(selected_parameters, file = paste0("data/selected_parameters_", 
                                        num_sims, "_", input_sigma_info_str,
     "_3.Rdata"))
