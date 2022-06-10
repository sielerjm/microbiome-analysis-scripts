# Functions

###########################################################################
# Considerations for the future:
#   - Split these functions into separate scripts based on category:
#     - alpha-diversity functions, beta-diversity, export tables/figure, etc.
###########################################################################

# CLR Matrix -------------------------------------------------------
#   Description: 
#   Input: 
#   Output: CLR matrix

gen.clr.matrix <- function(
  asv.mat,
  min_reads,
  min_prop = 0.001,
  min_occur = 0,
  smpls_by_row = TRUE,
  method = "CZM",
  lab = 0
) {
  require(CoDaSeq)
  require(zCompositions)
  asv.mat.f <- codaSeq.filter(
    asv.mat,
    min.reads = min_reads,
    min.prop = min_prop,
    min.occurrence = min_occur,
    samples.by.row = smpls_by_row
  )
  # replace 0 values with an estimate
  asv.mat.fn0 <- cmultRepl(t(asv.mat.f), method = method, label = lab)
  return(codaSeq.clr(asv.mat.fn0))
}


# Sub-char in variable ----------------------------------------------------


sub_char_var <- function(var, old.char, new.char){
  
  tmp.var <- var
  if(old.char == "."){
    
    print(paste0("Replacing period with ", new.char))  # TEST
    tmp.var.new <- gsub("\\.", new.char, # replaces periods 9.) with space
                        tmp.var
    )
  } else{
    
    print(paste0("Replacing ", old.var," with ", new.char))  # TEST
    tmp.var.new <- gsub(old.char, new.char, 
                        tmp.var
    )
  }
  
  return(tmp.var.new) 
}


# P-value Format ----------------------------------------------------------
#   Description: Formats P-values for summary statistic tables
#     * P-values below a certain threshold will appear as "<0.001"
#   Input: 
#   Output: 

p_val_format <- function(x){
  z <- scales::pvalue_format()(x)
  z[!is.finite(x)] <- ""
  z
}


# CountTable_Cat ----------------------------------------------------------
#   Description: Takes in categorical variables and then produces count tables
#   Input: dataframe, list of variables, variable to arrange/sort
#   Output: table of counts

CountTable_Cat <- function(data, vars, anchor_var, subset_var = NA){
  # data = dataframe, var = list of variables, anchor_var = arrange by X
  # [] How to manage more than two variables dynamically?
  
  var <- get_Variables(vars, "categorical")
  
  if(is.na(subset_var)){
      var <- var
  } else {
      var <- var[!var %in% subset_var]
  }
  
  print(var)
  
  data %>%
    dplyr::count(.dots = var,  sort = T) %>%  # .dots allows dynamic variable use
    dplyr::arrange(anchor_var)
}


# HeatMap -----------------------------------------------------------------
#   Description: Creates a heatmap from tabular data
#   Input: table of counts
#   Output: heatmap
#   Todo: Needs to be updated to be dynamic

HeatMap_cat <- function(table){
  
  heat.map <- table %>%
    ggplot(aes(table[[1]], table[[2]])) +
    geom_tile(aes(fill = n), colour = "white", show.legend = F) + 
    geom_label(aes(label = n), size = 6, alpha = .75) +
    facet_wrap(table[[3]] ~ .) +  # Needs to be updated to be dynamic
    labs(
      # title = paste0("Number of fish by ", 
      #                 colnames(table[1]),
      #                 ", ", 
      #                 colnames(table[2]), " and ", 
      #                 colnames(table[3])), 
         
         x = colnames(table[1]), 
         y = colnames(table[2]), 
         caption = "") +
    theme(axis.text = element_text(size = 10)
    )
  return(heat.map)             
  
}



# Get variable names ------------------------------------------------------
# Description: Extracts a list of variables you want from a dataframe
# Input: dataframe 
# Output: list of variables from column 1 if subset.type condition met

# Extracts variables with more broadness

get_Variables <- function(data.vars, subset.type = NULL, col.1 = 1, col.2 = 2){
  # Expects a two column df/dt, but can handle any size if col's specified
  
  # returns a list of variables
  list.vars <- data.vars %>%
    
    # If subset.type is not specified, returns all variables in col.1
    filter(if (is.null(subset.type)) TRUE else across(col.2) == subset.type) %>% 
    pull(col.1)
  return(list.vars)
}



# get_Variables <- function(data.vars, subset.type = NULL){
#   list.vars <- data.vars %>%
#     # If subset.type is not specified, returns all variables in column 1
#     filter(if (is.null(subset.type)) TRUE else across(2) == subset.type) %>% 
#     pull(1)
#   return(list.vars)
# }


# # Retired functions (replaced by above)
# get_Variables_OLD <- function(vars, data.type = c("numerical", "categorical", "interaction")){
#   var <- vars %>% 
#     filter(type %in% data.type) %>% 
#     pull(variable)
#   return(var)
# }
# 
# # Same in principle as get_Variable(), but with microbiome diversity metrics
# get_Metrics_OLD <- function(vars, metric.type = c("alpha", "beta", "interaction")){
#   var <- vars %>% 
#     filter(type %in% metric.type) %>% 
#     pull(metric)
#   return(var)
# }


# Calculate Alpha Scores --------------------------------------------------
#   Description: Generates alpha diversity scores from a list of alpha methods
#   Input: phyloseq object, list of alpha div. methods, 
#   Output: dataframe of alpha-diversity scores

alpha_base <- function(phyloseq.obj, alpha.methods, smpl.col.name, phylo.div = F){  
  tax.dt <- estimate_richness(
    physeq = phyloseq.obj,
    measures = alpha.methods[1:3]
  ) %>% as.data.table(keep.rownames = smpl.col.name) %>% setkeyv(smpl.col.name)
  tax.dt[, se.chao1 := NULL]
  
  if(isTRUE(phylo.div)){
      ## Need phylogeny data for this to work
      phy.dt <- pd(samp = otu.matrix(phyloseq.obj), tree = phy_tree(phyloseq.obj)) %>%
          as.data.table(keep.rownames = smpl.col.name) %>%
          setkeyv(smpl.col.name)
      
      names(phy.dt)[2:3] <- alpha.methods[4:5]
      tax.dt[phy.dt]
      
      return (tax.dt[phy.dt])
  }

  return (tax.dt)
}



# Normalize Alpha Scores --------------------------------------------------
#   Description: normalizes alpha diversity scores based on their distributions
#   Input: dataframe of alpha diversity scores, metadata table
#   Output: normalized datatable of alpha scores (0 to 1)

norm_alpha_score <- function(alpha.base, sample.df, methods){
  model.data.base <- copy(alpha.base[alpha.base$Sample %in% 
                                       row.names(sample.df)])
  
  for (alpha in methods) {
    
    progress_Bar_Start(which(methods == alpha), length(methods))  # Progress bar
    
    if (ad.test(model.data.base[[alpha]])$p.value <= 0.05) {
      trans <- transformTukey(model.data.base[[alpha]], plotit = F, quiet = F, statistic = 2)
      trans <- (trans-min(trans))/(max(trans)-min(trans))   # Fixes normalization 0 to 1
      
      if (ad.test(trans)$p.value > 0.05) {
        model.data.base[[alpha]] <- trans
        # print(paste0("Transforming: ", alpha))  # TEST
        print(paste0("Finished: ", alpha))
      } else {
        model.data.base[[alpha]] <- (model.data.base[[alpha]] - min(model.data.base[[alpha]] ))/(max(model.data.base[[alpha]] )-min(model.data.base[[alpha]] ))  # Fixes normalization 0 to 1
        # print(paste0("dividing by max: ", alpha ))  # TEST
        print(paste0("Finished: ", alpha))
      } 
    } else {
      model.data.base[[alpha]] <- (model.data.base[[alpha]] - min(model.data.base[[alpha]] ))/(max(model.data.base[[alpha]] )-min(model.data.base[[alpha]] ))  # Fixes normalization 0 to 1
      print(paste0("Finished: ", alpha))
    }
  }
  
  return(model.data.base)
}


# Alpha GLM ---------------------------------------------------------------
#   Description: Runs a generalized linear model on alpha scores and variables
#   Input: sample datatable, normalized alpha scores, variables, terms, glm family
#   Output: glm model



alpha_glm <- function(data, alpha.base.norm, alpha.methods, vars, terms, glm.fam, rand.eff = NA) {
  
  cl <- makeCluster(num.cores, type = "FORK", outfile = "")
  registerDoParallel(cl, num.cores)
  
  var <- get_Variables(vars)  # Extract a list of variables
  
  model.data <- data[alpha.base.norm, on = "Sample"] %>% setkeyv("Sample")
  
  # print(paste0(" ~ (", paste0(var, collapse = "+"), ")^", terms, " + ( 1 |", rand.eff, ")" ) )
    
  # return()
  # Anova
  alpha.glm <- foreach(
    alpha = alpha.methods,
    .verbose = TRUE
  ) %dopar% {
    
    if(is.na(rand.eff)){
        
        # form <- if(length(vars) == 1){
        #   paste0("dist.mat ~ ", vars)
        # } else{paste0("dist.mat ~ (", paste0(vars, collapse = "+") ,")^", terms)}
        full.frm <- paste0(alpha, " ~ (", paste0(var, collapse = "+"), ")^", terms)
        full.mod <- glm(as.formula(full.frm), 
                                   data = model.data,
                                   family = glm.fam)
        
        return(Anova(full.mod, type = 2) %>%
                   tidy() %>%
                   mutate(sig = ifelse(p.value <= 0.05, "*", "")) %>%
                   mutate(metric = alpha, .before = 1) %>%
                   filter((p.value < 0.15 | is.na(p.value)) & df > 0) %>%
                   arrange(desc(statistic))  # highest to lowest effect size
        )
            
    } else {
        
        full.frm <- paste0(alpha, " ~ (", paste0(var, collapse = "+"), ")^", terms, " + ( 1 |", rand.eff, ")" )
        full.mod <- lme4::glmer(as.formula(full.frm),
                        data = model.data,
                        family = glm.fam)
        return(Anova(full.mod, type = 2) %>%
                   tidy() %>%
                   mutate(sig = ifelse(p.value <= 0.05, "*", "")) %>%
                   mutate(metric = alpha, .before = 1) %>%
                   # filter((p.value < 0.15 | is.na(p.value)) & df > 0) %>%
                   arrange(desc(statistic))  # highest to lowest effect size
        )
    }
    
      
    # return(Anova(full.mod, type = 2) %>%
    #   tidy() %>%
    #   mutate(sig = ifelse(p.value <= 0.05, "*", "")) %>%
    #   mutate(metric = alpha, .before = 1) %>%
    #   filter((p.value < 0.15 | is.na(p.value)) & df > 0) %>%
    #   arrange(desc(statistic))  # highest to lowest effect size
    # )
  } %>%
    bind_rows()
  
  stopCluster(cl)
  return(alpha.glm)
}

###

alpha_glm_2 <- function(data, alpha.base.norm, alpha.methods, vars, terms, glm.fam, rand.eff = NA) {
  
  cl <- makeCluster(num.cores, type = "FORK", outfile = "")
  registerDoParallel(cl, num.cores)
  
  # var <- get_Variables(vars)  # Extract a list of variables
  
  model.data <- data[alpha.base.norm, on = "Sample"] %>% setkeyv("Sample")
  
  # print(paste0(" ~ (", paste0(var, collapse = "+"), ")^", terms, " + ( 1 |", rand.eff, ")" ) )
  
  # return()
  # Anova
  alpha.glm <- foreach(
    alpha = alpha.methods,
    .verbose = TRUE
  ) %dopar% {
    
    if(is.na(rand.eff)){
      
      form <- if(length(vars) == 1){
        paste0(alpha, " ~ ", vars)
      } else{paste0(alpha, " ~ (", paste0(vars, collapse = "+") ,")^", terms)}
      # full.frm <- paste0(alpha, " ~ (", paste0(var, collapse = "+"), ")^", terms)
      full.mod <- glm(as.formula(form), 
                      data = model.data,
                      family = glm.fam)
    } else {
      
      form <- if(length(vars) == 1){
        paste0("alpha ~ ", vars)
      } else{paste0(paste0("alpha ~ (", paste0(vars, collapse = "+") ,")^", terms), " + ( 1 |", rand.eff, ")" )}
      # full.frm <- paste0(alpha, " ~ (", paste0(var, collapse = "+"), ")^", terms, " + ( 1 |", rand.eff, ")" )
      full.mod <- lme4::glmer(as.formula(form),
                              data = model.data,
                              family = glm.fam)

    }
    
    return(Anova(full.mod, type = 2) %>%
             tidy() %>%
             mutate(sig = ifelse(p.value <= 0.05, "*", "")) %>%
             mutate(metric = alpha, .before = 1) %>%
             # filter((p.value < 0.15 | is.na(p.value)) & df > 0) %>%
             arrange(desc(statistic))  # highest to lowest effect size
    )

  } %>%
    bind_rows()
  
  stopCluster(cl)
  return(alpha.glm)
}

###


# general GLM (updated) ---------------------------------------------------
#   Description: Runs a generalized linear model on user inputted x and y variables
#   Input: datatable, variables, terms, glm family
#   Output: glm model


gen_glm <- function(data, x.var = NULL, y.var = NULL, loop = NA, loop.var = NULL, 
                    extra.data = NA, samp.name = "Sample", terms = 1, glm.fam = "gaussian", 
                    rand.eff = NA, fact.y = F, run.anv = T) {
  
  cl <- makeCluster(num.cores, type = "FORK", outfile = "")
  registerDoParallel(cl, num.cores)
  
  # if(!isTRUE(run.anv)){
  #   glm.mod <- vector("list", length(loop.var))
  # }
  
  # Factorizes y variable if fact.y == T
  if(isTRUE(fact.y)){
    y.var <- factor(y.var)
  }
  
  # Checks to see if more than one datatable exists, then merges them if T
  if(!is.na(extra.data)){
    model.data <- data[extra.data, on = samp.name] %>% setkeyv(samp.name)
  } else {
    model.data <- data
  }
  print(paste0(y.var, " ~ ", paste0( c(loop.var, x.var), collapse = " + ") ))  #TEST
  # IF LOOPING ON AN X VAR
  if(!is.na(loop)){
    if(loop == "x" | loop == "X"){
      
      print(paste0("here: loop = ", loop, " .loop.var = ", loop.var)) # TEST
      
      glm.mod <- foreach(
        loop.v = loop.var,
        .verbose = TRUE
      ) %dopar% {
        # Build glm formula
        form <- glm_form_x(x.var, y.var, loop.v, terms, glm.fam)
        
        # Build glm() model
        tmp.mod <- glm(as.formula(form), data = model.data, family = glm.fam)
        
        # Run Anova on glm model
        if(isTRUE(run.anv)){
          # Returns ANOVA results
          return(gen_glm_anova(tmp.mod, y.var))
        } else {
          # Returns the glm model so you can run anova with more detail separately
          return(tmp.mod)
        }
        
                
      } %>%
        bind_rows()  # Combines multiple rows if looping
      
    # IF LOOPING ON Y VAR
    } else if(loop == "y" | loop == "Y"){
      
      print(paste0("here: loop = ", loop, " .loop.var = ", loop.var)) # TEST
      
      glm.mod <- foreach(
        loop.v = loop.var,
        .verbose = TRUE
      ) %dopar% {
        # Build glm formula
        form <- glm_form_y(x.var, loop.v, terms, glm.fam)
        
        # Build glm() model
        tmp.mod <- glm(as.formula(form), data = model.data, family = glm.fam)
        
        # Run Anova on glm model
        if(isTRUE(run.anv)){
          # Returns ANOVA results
          return(gen_glm_anova(tmp.mod, loop.v))
        } else {
          # Returns the glm model so you can run anova with more detail separately
          return(tmp.mod)
        }
        
      } %>% 
        bind_rows()
        # {if(isTRUE(run.anv)) bind_rows()}  # Combines multiple rows if looping

    } else if(is.null(loop.var)){
      warning("Error in loop: Cannot leave loop.var empty")
    }
  } else {

    form <- glm_form(x.var, y.var, terms, glm.fam)
    
    # Build glm() model
    tmp.mod <- glm(as.formula(form), data = model.data, family = glm.fam)
    
    # Run Anova on glm model
    if(isTRUE(run.anv)){
      # Returns ANOVA results
      return(gen_glm_anova(tmp.mod, y.var))
    } else {
      # Returns the glm model so you can run anova with more detail separately
      return(tmp.mod)
    }
    
  }
  
  print(glm.mod)
  
  
  stopCluster(cl)
  return(glm.mod)
}

# Builds formulas when looped variable is an x var
glm_form <- function(x.var, y.var, terms, glm.fam){
  if(terms == 1){
    tmp.form <- paste0(y.var, " ~ ", paste0( c(x.var), collapse = " + ") )
  } else {
    tmp.form <- paste0(y.var, " ~ (", paste0( c(x.var), collapse = " + ") ,")^", terms)
  }
  
  print(paste0("FORM: ", tmp.form))  # TEST
  return(tmp.form)
}


# Builds formulas when looped variable is an x var
glm_form_x <- function(x.var, y.var, loop.v, terms, glm.fam){
  # if(is.na(x.var)){
  #   if(terms == 1){
  #     tmp.form <- paste0(y.var, " ~ ", loop.v)
  #   } else { 
  #     tmp.form <- paste0(y.var, " ~ (", paste0( c(loop.v, x.var), collapse = " + ") ,")^", terms)
  #   }
  # } else {
  #   if(terms == 1){
  #     tmp.form <- paste0(y.var, " ~ ", loop.v)
  #   } else {
  #     tmp.form <- paste0(y.var, " ~ (", paste0( c(loop.v, x.var), collapse = " + ") ,")^", terms)
  #   }
  # }
  
  if(terms == 1){
    tmp.form <- paste0(y.var, " ~ ", paste0( c(loop.v, x.var), collapse = " + ") )
  } else {
    tmp.form <- paste0(y.var, " ~ (", paste0( c(loop.v, x.var), collapse = " + ") ,")^", terms)
  }
  
  print(paste0("FORM: ", tmp.form))  # TEST
  return(tmp.form)
}

# Builds formulas when looped variable is an y var
glm_form_y <- function(x.var, loop.v, terms, glm.fam){
  if(terms == 1){
    tmp.form <- paste0(loop.v, " ~ ", paste0( c(x.var), collapse = " + ") )
  } else {
    tmp.form <- paste0(loop.v, " ~ (", paste0(x.var, collapse = " + ") ,")^", terms)
  }

  print(paste0("FORM: ", tmp.form))  # TEST
  return(tmp.form)
}

gen_glm_anova <- function(tmp.mod, tmp.metric, filt.pval = 1){
  return(Anova(tmp.mod, type = 2) %>%
           tidy() %>%
           mutate(sig = ifelse(p.value <= 0.05, "*", "")) %>%
           mutate(metric = tmp.metric, .before = 1) %>%
           filter((p.value < filt.pval | is.na(p.value)) & df > 0) %>%
           arrange(desc(statistic))  # highest to lowest effect size
  )
}

gen_glm_anova_2 <- function(tmp.mod, tmp.metric, term = NULL){
  return(Anova(tmp.mod, type = 2) %>%
           tidy() %>%
           mutate(sig = ifelse(p.value <= 0.05, "*", "")) %>%
           mutate(term = replace(term, term != "Residuals", tmp.nm)) %>%
           mutate(metric = tmp.metric, .before = 1) %>%
           arrange(desc(statistic))  # highest to lowest effect size
  )
}

# Check Response Vars -----------------------------------------------------------
#   Description: 
#   Input: list of dataframes/datatables 
#   Output: list of of dataframes

check_Exp_Vars <- function(lists, excl){
  if(is.na(excl)){
    return(lists)
  } else{
    return(lists[excl])
  }
}


# Excluded Vars -----------------------------------------------------------
#   Description: 
#   Input: list of dataframes/datatables 
#   Output: list of of dataframes

exclude_vars <- function(data){
  
  tmp.list<- if(deparse(substitute(data)) == "metrics"){
                unlist(lapply(get_Variables(data), function(x){
                  lapply(get_Variables(data), function(y){
                    paste0(c(x,y), collapse = ":")
                  })
                }))
              } else if(deparse(substitute(data)) == "variables"){
                unlist(lapply(get_Variables(data), function(x){
                  lapply(get_Variables(data), function(y){
                    paste0(c(x,y), collapse = ":")
                  })
                }))
              } else {
                "Unknown data type"
              }
  print(paste0("inside of exclude_var() > tmp.list: ", tmp.list))  # TEST 
  return(tmp.list)
}



# Alt Alpha Analysis ---------------------------------------------------------------
#   Description: runs statistical analysis on physiology and alpha metrics
#   Input: 
#   Output: glm model



alt_alpha_analysis <- function(data, methods, vars, r.var){
    
    tmp.obj <- list()
    
    tmp.obj[["GLM"]] <- lapply(methods, function(alpha){
        glm_func(data = data, 
                 resp.vars = r.var,  # response variable
                 exp.vars = list(variables = vars,  # list of dataframes of response variables
                                 metrics = metrics[which(metrics$variable == alpha),]),  # list(<set_list_name> = <list>)
                 glm.fam = "gaussian")
    })
    
    print("Completed GLM")
    
    tmp.obj[["stats.table"]] <- lapply(methods, function(alpha){
        stats_table(tmp.obj[["GLM"]][[alpha]],  # import glm data
                    vars = vars, # imports variables of interest, numerical
                    terms = 2, # import factors
                    methods = methods,
                    formula = alpha
        ) 
    })
    
    print("Completed Stats Table")
    
    tmp.obj[["sig.vars"]] <- lapply(methods, function(alpha){
        print(sig_terms(tmp.obj[["GLM"]][[alpha]], 
                  merge(vars, metrics, all = T)))
    })
    
    print("Completed Sig Vars")
    
    return(tmp.obj)
    
}



# Generic GLM ---------------------------------------------------------------
#   Description: Runs a generalized linear model on explanatory and response variables 
#   Input: 
#   Output: glm model

# GLM Families
#   binomial(link = "logit") #  for continuous decimal data with normal distribution, like weight, length, et al.
#   gaussian(link = "identity")
#   Gamma(link = "inverse")
#   inverse.gaussian(link = "1/mu^2")
#   poisson(link = "log")
#   quasi(link = "identity", variance = "constant")
#   quasibinomial(link = "logit")
#   quasipoisson(link = "log")


# glm_func <- function(data, alpha.base.norm, alpha.methods, vars, terms, glm.fam) {
glm_func <- function(data, resp.vars, exp.vars, terms = 2, glm.fam = "quasibinomial") {  
  
  # cl <- makeCluster(num.cores, type = "FORK", outfile = "")
  # registerDoParallel(cl, num.cores)
  
  tmp.test <- lapply(exp.vars, get_Variables)
  # print(tmp.test)  # TEST
  exp.var <- unlist(lapply(exp.vars, get_Variables), use.names = F)
  # print(exp.var)  # TEST
  
  glm <- {

    lapply(resp.vars, function(resp) {

      progress_Bar_Start(which(resp.vars == resp), length(resp.vars))  # Progress bar

      full.frm <- paste0(resp, " ~ (", paste0(exp.var[exp.var != resp], collapse = "+"), ")^", terms)
      print(full.frm)  # TEST

      full.mod <- glm(as.formula(full.frm),
                      data = data,
                      family = glm.fam)
      return(
        Anova(full.mod, type = 2) %>%
          tidy() %>%
          mutate(sig = ifelse(p.value <= 0.05, "*", "")) %>%
          mutate(metric = resp, .before = 1) %>%
          filter(p.value < 0.15 | is.na(p.value)) %>%
          arrange(desc(statistic))  # highest to lowest effect size
      )
    }) %>%
      bind_rows()

  }
  
  # stopCluster(cl)
  return(glm)
  
}

# cl <- makeCluster(num.cores, type = "FORK", outfile = "")
# registerDoParallel(cl, num.cores)
# 
# # Anova
# res <- foreach(
#   beta = methods,
#   .final = names,
#   .verbose = TRUE
# ) %dopar% {
#   mod <- beta.model[[beta]]
#   anova(mod, by = "term") %>% 
#     tidy() %>% 
#     mutate(sig = ifelse(p.value <= 0.05, "*", "")) %>%
#     mutate(metric = beta, .before = 1) %>%
#     arrange(desc(statistic))  # highest to lowest effect size 
# } %>%
#   bind_rows()
# 
# stopCluster(cl)
# return(res)

# Variable Name -----------------------------------------------------------
#   Description: converts a variable's name to a string
#   Input: variable
#   Output: string of variable's name

var_name <- function(var.name){
  return (deparse(substitute(var.name)) )
}

## A better version of the above
var_name_sub <- function(var.name){
    
    new.var.name <- gsub("\\.", "_", # replaces periods (.) with dashes (-)
                         deparse(substitute(var.name)))
    # print(new.var.name)  # TEST
    return (new.var.name)
} 


# Export Table ------------------------------------------------------------
#   Description: exports an image of table
#   Input: table (from flextable)
#   Output: exports a ".png" image of table

export_table_img <- function(dataframe, var.name, output, type = "Tables", ID = Sys.Date()) {
  # export_table(<name_of_table>, var_name(<name_of_table>))
  #   Change "<name_of_table>" to the variable storing your table
  
  # If obj is actually a list of dataframes
  if(class(dataframe) == "list"){ 
    
    # Will write a separate csv for each dataframe
    for (x in seq_along(dataframe)) {
      print(dataframe[[x]])  # TEST
      df_name <- paste0(output,  # Output path
                        "/", type, "/",  # sub-directory
                        gsub("\\.", "_", # replaces periods 9.) with dashes (-)
                             var.name), "-", # Variable name
                        names(dataframe[x]),  # name of list
                        "_", ID,  # Analysis ID
                        ".png")  # image output type
      # Save table as image
      save_as_image(dataframe[[x]], df_name)
    }
    
  } else {
    df_name <- paste0(output,  # Output path
                      "/", type, "/",  # sub-directory
                      gsub("\\.", "_", # replaces periods 9.) with dashes (-)
                           var.name),  # Variable name
                      "_", ID,  # Analysis ID
                      ".png")  # image output type
    
    # Save table as image
    save_as_image(dataframe, df_name)
  }
  
}



# Export CSV --------------------------------------------------------------
#   Description: exports a ".CSV" 
#   Input: tabular data
#   Output: ".CSV" file

export_table_csv <- function(dataframe, var.name, output, ID = Sys.Date()){

  # If obj is actually a list of dataframes
  if(class(dataframe) == "list"){ 
    
    # Will write a separate csv for each dataframe
    for (x in seq_along(dataframe)) {
      # Write dataframe to CSV
      write.csv(dataframe[[x]],
                file = paste0(output,  # Output path
                              "/Tables/",  # Tables sub-directory
                              gsub("\\.", "_", # replaces periods 9.) with dashes (-)
                             var.name), "-", # Variable name
                              names(dataframe[x]),  # name of list
                             "_", ID,  # Analysis ID,  # Analysis ID
                             ".csv"),  # image output type,
                row.names=TRUE)
    }

  } else {
    write.csv((dataframe), 
              file = paste0(output,  # Output path
                            "/Tables/",  # Tables sub-directory
                            gsub("\\.", "_", # replaces periods 9.) with dashes (-)
                                 var.name),  # Variable name
                            "_", ID,  # Analysis ID,  # Analysis ID
                            ".csv"),  # image output type,
              row.names=TRUE) 
  }
  
}



# Melt Data Table ---------------------------------------------------------
#   Description: Melts sample data table for easy plotting
#   Input: normalized alpha scores, alpha methods, variable names
#   Output: a melted datatable

melt_to_datatable <- function(data.table, methods, var.name){
  model.data <- sample.dt[data.table, on = ("Sample")] %>% setkeyv("Sample")
  return( melt(
            model.data,
            measure.vars = methods, 
            variable.name = var.name, 
            value.name = "Score"
    )
  )
}

# More flexible version of above

melt_to_datatable_2 <- function(datatable1, datatable2, vars, var.name, samp.name = "Sample", val.name = "Score"){
  
  comb.data <- datatable1[datatable2, on = (samp.name)] %>% setkeyv(samp.name)
  
  return( melt(
    comb.data,  # combined datatables from above
    measure.vars = vars, # Measure variables for melting. Can be missing, vector, list, or pattern-based.
    variable.name = var.name, # Name for the measured variable names column.
    value.name = val.name  # Name for the molten data values column(s).
  )
  )
}


# Significant Terms -------------------------------------------------------
#   Description: Returns a list of (distinct) significant terms
#   Input: GLM model
#   Output: list of significant terms

sig_terms <- function(model, vars){
  # Extract significant variables from GLM model
  terms <- model %>%  # dataframe
    filter(p.value <0.05) %>%  # term must have p.val below 0.05
      distinct(term) %>%  # unique, aka no duplicates added
      pull() # pulls value and adds to terms

  # Dataframe containing significant terms/variables with empty datatype column
  tmp.df <- data.frame(variable = terms,
                        type = NA)
  
  # Assigns datatype to "type" column in tmp.df
  sig.vars <- tmp.df %>% 
    mutate(type = ifelse(grepl(":", variable),  # checks if var/term is in variables dataframe
                         "interaction",  # assigns var type from variables dataframe
                         vars$type[match(tmp.df$variable, vars$variable)]))
  
  # Returns a list of signficant variables/terms and their datatype
  return(sig.vars)
}



# sig_terms_multi <- function(model, data){
#     
#     # test <- list()
#     test <- lapply(model, function(x){
#                 # print(x)
#                 terms <- x %>%  # dataframe
#                     filter(p.value <0.05) %>%  # term must have p.val below 0.05
#                     distinct(term) %>%  # unique, aka no duplicates added
#                     pull() # pulls value and adds to terms
# 
#                 # Dataframe containing significant terms/variables with empty datatype column
#                 tmp.df <- data.frame(variable = terms,
#                                      type = NA)
#                 
#                 for (y in 1:length(data)){
#                     print(data[[3]][[1]]$type[match(tmp.df[[3]][[1]]$variable, data[[3]][[1]]$variable)])
#                     
#                     # Assigns datatype to "type" column in tmp.df
#                     # tmp.sig.vars <- tmp.df %>% 
#                     #     mutate(type = ifelse(grepl(":", variable),  # checks if var/term is in variables dataframe
#                     #                          "interaction",  # assigns var type from variables dataframe
#                     #                          data[[3]][[1]]$type[match(tmp.df[[3]][[1]]$variable, data[[3]][[1]]$variable)]))
#                 }
#                 
#         return(sig.vars)
#     })
#     
#     
#     # print(paste0("test: ", test))
#     
#     # print(paste0("tmp.df: ", tmp.df))
#     
#     # sig.vars <- lapply(test, function(x){
#     #     print(tmp.vars$type[match(x$variable, tmp.vars$variable)])
#     #     sig.vars <- x %>%
#     #         print(ifelse(grepl(":", variable),  # checks if var/term is in variables dataframe
#     #                      "interaction",  # assigns var type from variables dataframe
#     #                      tmp.vars$type[match(x$variable, tmp.vars$variable)]))
#     #     mutate(type = ifelse(grepl(":", variable),  # checks if var/term is in variables dataframe
#     #                          "interaction",  # assigns var type from variables dataframe
#     #                          tmp.vars$type[match(x$variable, tmp.vars$variable)]))
#     #     
#     # })
#     # # Assigns datatype to "type" column in tmp.df
#     # sig.vars <- tmp.df %>%
#     #     mutate(type = ifelse(grepl(":", variable),  # checks if var/term is in variables dataframe
#     #                          "interaction",  # assigns var type from variables dataframe
#     #                          vars$type[match(tmp.df$variable, vars$variable)]))
#         
#     
#     
#     # terms <- model %>%  # dataframe
#     #     filter(p.value <0.05) %>%  # term must have p.val below 0.05
#     #     distinct(term) %>%  # unique, aka no duplicates added
#     #     pull() # pulls value and adds to terms
#     # 
#     # # Dataframe containing significant terms/variables with empty datatype column
#     # tmp.df <- data.frame(variable = terms,
#     #                      type = NA)
#     # 
#     # # Assigns datatype to "type" column in tmp.df
#     # sig.vars <- tmp.df %>% 
#     #     mutate(type = ifelse(grepl(":", variable),  # checks if var/term is in variables dataframe
#     #                          "interaction",  # assigns var type from variables dataframe
#     #                          vars$type[match(tmp.df$variable, vars$variable)]))
#     
#     # Returns a list of signficant variables/terms and their datatype
#     return(test)
# }



# Box Plot ----------------------------------------------------------------
#   Description: Creates a box plot
#   Input: data, x variable, y variable, facet wrap conditions
#   Output: a box plot
#   To Do: 
#     [] Import dynamic color settings from RcolorBrewer

box_plot <- function(data, x.var, y.var, facet.var = NULL) {
  #print(data$facet.var)
  p <- ggplot(data,
              aes_string(x = x.var, y= y.var)) +
    geom_boxplot(aes(fill = eval(parse(text = x.var)) ) # converts string to var name
                 , alpha = 0.50  # changes opacity
    ) + 
    geom_violin(aes(fill = eval(parse(text = x.var)) ) # converts string to var name
                , alpha = 0.50  # changes opacity
    ) +
    #facet_grid(as.formula(paste(data$facet.var, " ~ .")), scales = "free_y")
    facet_grid(paste(facet.var, " ~ ."), scales = "free_y")
  
  return(p)
}


# Scatter Plot ------------------------------------------------------------
#   Description: Creates a scarter plot
#   Input: data, x variable, y variable, facet wrap conditions
#   Output: a scatter plot
#   To Do: 
#     [] Import dynamic color settings from RcolorBrewer

scatter_plot <- function(data, x.var, y.var, facet.var = NULL) {

  p <- ggplot(data,  
              aes_string(x = x.var, y= y.var)) +
    geom_point(aes(color = x.var)) + 
    geom_smooth(aes(fill = x.var), 
                #color = x.var, 
                method="glm", 
                size = 1, 
                alpha = .25) + 
    facet_grid(paste(facet.var, " ~ ."), scales = "free_y")
  
  return(p)
  # return(ggMarginal(p, fill = "grey", type = "density"))
}



# Generate Multiple Box Plots (Alpha) -----------------------------------------
#   Description: Generates alpha diversity plots of significant variables 
#   Input: data, significant variables, alpha score, facet conditions
#   Output: multiple plots displaying alpha diversity

gen_multi_alpha_box_plots <- function(data, vars, y.var, facet.var = NULL){
  
  print(paste0(vars))
  
  var <- get_Variables(vars, "categorical")
  
  print(var)
  
  tmp.list <- vector("list", length(var))
  for (x in 1:length(tmp.list)) {
    
    
    
    tmp.list[[x]] <- box_plot(data,  # data
                              var[x],  # x variable, significant variables 
                              y.var,  # y variable
                              facet.var  # facet grid variable
    )
  }
  
  return(tmp.list)
}


# Generate Multiple Scatter Plots (Alpha) -----------------------------------------
#   Description: Generates alpha diversity plots of significant variables 
#   Input: data, significant variables, alpha score, facet conditions
#   Output: multiple plots displaying alpha diversity

gen_multi_alpha_scatter_plots <- function(data, vars, y.var, facet.var = NULL){
  
  var <- get_Variables(vars, "numerical")
  
  tmp.list <- vector("list", length(var))
  for (x in 1:length(tmp.list)) {
    tmp.list[[x]] <- scatter_plot(data,  # data
                              var[x],  # x variable, significant variables 
                              y.var,  # y variable
                              facet.var  # facet grid variable
    )
  }
  
  return(tmp.list)
}

# Generate Multiple Plots -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

gen_multi_plots <- function(data, sig.vars, resp.var, all.vars, an.type = "", axis = NA, stats = NA, subset.arg = NA){
  
  if(!is.na(subset.arg)){
    data <- subset(data, eval(parse(text = subset.arg)))  # subset data based on this argument
  } else {
    data <- data  # Data isn't subsetted
  }
  
  print(data)
  
  # Creats a list of the significant variables
  vars <- get_Variables(sig.vars)
  
  # Creates an empty list to store plots 
  tmp.list <- vector("list", length(get_Variables(sig.vars)))
  # tmp.list <- vector("list")
  
  # Creates an additional empty list for storing plots with int terms
  tmp.list.int <- vector("list", length(get_Variables(sig.vars, subset.type = "interaction")))
  
  # Temporary work around in the case there are no significant interaction terms 
  if(length(get_Variables(sig.vars, subset.type = "interaction")) == 0){
      tmp.list.int <- NULL
  }
  
  for (x in 1:length(tmp.list)) {
    
    tmp.exp.var <- vars[x]
      
    print(paste0("tmp.exp.var: ", tmp.exp.var))  # TEST

    # progress_Bar_Start(iteration = x, n_iter = length(tmp.list)+length(tmp.list.int)) # Progress Bar
    

# alpha -------------------------------------------------------------------


    
    
    if(an.type == "alpha"){
      # print(paste0("gen_multi_plots > if: "))  # TEST
      
      if(check_Var_Type(sig.vars[x,]) == "numerical"){
        
        # Send to gen_scatter_plot(data, x.var, y.var, facet.var = NULL)
        tmp.list[[paste0(tmp.exp.var, "_", x, "_")]] <- gen_scatter_plot(data, tmp.exp.var, resp.var, facet.var.y = "alpha.methods")
        
      } else if(check_Var_Type(sig.vars[x,]) == "categorical") {
        
        # Send to gen_box_plot(data, x.var, y.var, facet.var = NULL)
        tmp.list[[paste0(tmp.exp.var, "_", x, "_")]] <- gen_box_plot(data, tmp.exp.var, resp.var, facet.var.y = "alpha.methods")
        
      } else if(check_Var_Type(sig.vars[x,]) == "interaction"){
        
        # Extract variables from interaction terms
        tmp.ex.int <- extract_Int_Terms(tmp.exp.var, all.vars)
        
        # Check if there is a numerical variable in int terms
        # if(tmp.ex.int$type %in% "categorical"){
        if("categorical" %in% tmp.ex.int$type ){
          
          # Check if both are categorical, either var can be facet
          if(all(tmp.ex.int$type %in% "categorical")){
            
            # Create a list of variables
            tmp.ex.int.vars <- get_Variables(tmp.ex.int)
            
            # Make a plot
            tmp.list[[paste0(tmp.ex.int.vars[2], ":", tmp.ex.int.vars[1], "_", x, "_")]] <- gen_box_plot(data, 
                                          tmp.ex.int.vars[2], 
                                          resp.var, 
                                          facet.var.x = tmp.ex.int.vars[1],
                                          facet.var.y = "alpha.methods")
            
            # Make additional plot for other way of visualization
            tmp.list.int[[paste0(tmp.ex.int.vars[1], ":", tmp.ex.int.vars[2], "_", x, "_")]] <- gen_box_plot(data, 
                                                                                            tmp.ex.int.vars[1], 
                                                                                            resp.var, 
                                                                                            facet.var.x = tmp.ex.int.vars[2],
                                                                                            facet.var.y = "alpha.methods")
            
            # If one var is numerical, then cat var is facet var
          } else {
            # Assign vars
            tmp.ex.int.var.num <- get_Variables(tmp.ex.int, subset.type = "numerical")
            tmp.ex.int.var.cat <- get_Variables(tmp.ex.int, subset.type = "categorical")
            
            # Make a plot
            tmp.list[[paste0(tmp.ex.int.var.num, ":", tmp.ex.int.var.cat, "_", x, "_")]] <- gen_scatter_plot(data, 
                                              tmp.ex.int.var.num, 
                                              resp.var, 
                                              facet.var.x = tmp.ex.int.var.cat,
                                              facet.var.y = "alpha.methods")
          }
          
          # Otherwise both int terms are numerical variables
        } else { 
          print("I dunno how to plot this. Two numerical interaction terms")  # FIX
            
            ### Two Numerical Interaction Variables ###
            
            # Create a list of variables
            tmp.ex.int.vars <- get_Variables(tmp.ex.int)
            
            # Make a plot
            tmp.list[[paste0(tmp.ex.int.vars[2], ":", tmp.ex.int.vars[1], "_", x, "_")]] <- gen_scatter_plot(data, 
                                                                                                         tmp.ex.int.vars[2], 
                                                                                                         resp.var, 
                                                                                                         facet.var.x = tmp.ex.int.vars[1],
                                                                                                         facet.var.y = "alpha.methods")
            
            # Make additional plot for other way of visualization
            tmp.list.int[[paste0(tmp.ex.int.vars[1], ":", tmp.ex.int.vars[2], "_", x, "_")]] <- gen_scatter_plot(data, 
                                                                                                             tmp.ex.int.vars[1], 
                                                                                                             resp.var, 
                                                                                                             facet.var.x = tmp.ex.int.vars[2],
                                                                                                             facet.var.y = "alpha.methods")
        }
        
      } else {
        print(print(paste0("Error, unknown variable type (", tmp.ex.int, ") in alpha")))
      }
      


# beta --------------------------------------------------------------------


        
        
    } else if(an.type == "beta"){
      
      # # Under construction: organize each section into sep functs
      # tmp.list[[x]] <- gen_multi_beta_plots(data = data, 
      #                      sig.vars = sig.vars, 
      #                      resp.var = resp.var, 
      #                      all.vars = all.vars, 
      #                      an.type = an.type, 
      #                      axis = axis, 
      #                      stats = stats)
      
      print("Inside of gen_multi_plots > beta")  # TEST

      if(check_Var_Type(sig.vars[x,]) == "numerical"){

        tmp.list[[paste0(tmp.exp.var, "_", x, "_")]] <- gen_PCoA_Num_Plot(data,
                                       tmp.exp.var,
                                       resp.var,
                                       axis,
                                       stats = stats)

      } else if(check_Var_Type(sig.vars[x,]) == "categorical") {

        tmp.list[[paste0(tmp.exp.var, "_", x, "_")]] <- gen_PCoA_Cat_Plot(data,
                                       tmp.exp.var,
                                       resp.var,
                                       axis,
                                       stats = stats)

      } else if(check_Var_Type(sig.vars[x,]) == "interaction"){

        # print("Inside of gen_multi_plots > beta > interaction")  # TEST

          # Extract variables from interaction terms
          tmp.ex.int <- extract_Int_Terms(tmp.exp.var, all.vars)
          # print(paste0("tmp.ex.int: ", tmp.ex.int))  #TEST

          # Create a list of variables
          tmp.ex.int.vars <- get_Variables(tmp.ex.int)

          # print(paste0("tmp.ex.int.vars: ", tmp.ex.int.vars))  #TEST

          # Check if there is a numerical variable in int terms
          # if(tmp.ex.int$type %in% "categorical"){
          if("categorical" %in% tmp.ex.int$type){

            # print("Inside of gen_multi_plots > beta > interaction > one cat")  # TEST

            # Check if both are categorical, either var can be facet
            if(all(tmp.ex.int$type %in% "categorical")){

              # print("Inside of gen_multi_plots > beta > interaction > two cat")  # TEST

              # # Create a list of variables
              # tmp.ex.int.vars <- get_Variables(tmp.ex.int)

              # Make a plot
              tmp.list[[paste0(tmp.ex.int.vars[1], ":", tmp.ex.int.vars[2], "_", x, "_")]] <- gen_PCoA_Int_Cat_Plot(data,
                                   tmp.ex.int.vars[1],
                                   x.var.2 = tmp.ex.int.vars[2],
                                   resp.var,
                                   axis,
                                   stats = stats)

              # Make an additional plot to store the reversed facets
              #   - Often one cat variable makes more sense to be the center of focus
              tmp.list.int[[paste0(tmp.ex.int.vars[2], ":", tmp.ex.int.vars[1], "_", x, "_")]] <- gen_PCoA_Int_Cat_Plot(data,
                                       tmp.ex.int.vars[2],
                                       x.var.2 = tmp.ex.int.vars[1],
                                       resp.var,
                                       axis,
                                       stats = stats)

              # If one var is numerical, then cat var is facet var
            } else {
              # print("Inside of gen_multi_plots > beta > interaction > one num")  # TEST

              # Assign cat var
              tmp.ex.int.var.cat <- get_Variables(tmp.ex.int, subset.type = "categorical")

              # Assign num var:

              # Check if second int term is "numerical" type
              if("numerical" %in% tmp.ex.int$type){
                tmp.ex.int.var.num <- get_Variables(tmp.ex.int, subset.type = "numerical")

                # Check if second int term is "alpha" type
              } else if("alpha" %in% tmp.ex.int$type){
                tmp.ex.int.var.num <- get_Variables(tmp.ex.int, subset.type = "alpha")

                # Check if second int term is "beta" type
              } else if("beta" %in% tmp.ex.int$type){
                tmp.ex.int.var.num <- get_Variables(tmp.ex.int, subset.type = "beta")

                # Error
              } else { print(paste0("Error, unknown variable type (", tmp.ex.int, ") in beta"))}


              # Make a plot
              tmp.list[[paste0(tmp.ex.int.var.cat, ":", tmp.ex.int.var.num, "_", x, "_")]] <- gen_PCoA_Int_Num_Plot(data,
                                             tmp.ex.int.var.cat,
                                             resp.var,
                                             axis,
                                             stats = stats,
                                             x.var.2 = tmp.ex.int.var.num)
            }

            # Otherwise both int terms are numerical variables
          } else {
            print("I dunno how to plot this. Two numerical interaction terms")  # FIX
              
              ### Two Interaction Variables ###
              
              
        }

      } else {
        print("Error, unknown data type")
      }

        

# other -------------------------------------------------------------------



    } else {
      # print(paste0("gen_multi_plots > else: "))  # TEST
      
      if(check_Var_Type(sig.vars[x,]) == "numerical"){
        
        # Send to gen_scatter_plot(data, x.var, y.var, facet.var = NULL)
        tmp.list[[x]] <- gen_scatter_plot(data, tmp.exp.var, resp.var)
        
      } else if(check_Var_Type(sig.vars[x,]) == "categorical") {
        
        # Send to gen_box_plot(data, x.var, y.var, facet.var = NULL)
        tmp.list[[x]] <- gen_box_plot(data, tmp.exp.var, resp.var)
        
      } else if(check_Var_Type(sig.vars[x,]) == "interaction"){
        
        # Extract variables from interaction terms
        tmp.ex.int <- extract_Int_Terms(tmp.exp.var, all.vars)
        # print(paste0("tmp.ex.int: ", tmp.ex.int))  #TEST
        
        # Create a list of variables
        tmp.ex.int.vars <- get_Variables(tmp.ex.int)
        
        # print(paste0("tmp.ex.int.vars: ", tmp.ex.int.vars))  #TEST
        
        # Check if there is a numerical variable in int terms
        # if(tmp.ex.int$type %in% "categorical"){
        if("categorical" %in% tmp.ex.int$type){
          
          # Check if both are categorical, either var can be facet
          if(all(tmp.ex.int$type %in% "categorical")){
              
              # print(paste0("Inside of Other > two cat vars "))  # TEST
            
            # Make a plot
            tmp.list[[x]] <- gen_box_plot(data, 
                                          tmp.ex.int.vars[1], 
                                          resp.var, 
                                          facet.var.x = tmp.ex.int.vars[2],)
            
            
            # Make an additional plot to store the reversed facets
            #   - Often one cat variable makes more sense to be the center of focus
            tmp.list.int[[x]] <- gen_box_plot(data, 
                                          tmp.ex.int.vars[2], 
                                          resp.var, 
                                          facet.var.x = tmp.ex.int.vars[1])
          
            # If one var is numerical, then cat var is facet var
          } else {
            
            # Assign cat var
            tmp.ex.int.var.cat <- get_Variables(tmp.ex.int, subset.type = "categorical")
            
            # Assign num var:
            
            # Check if second int term is "numerical" type
            if("numerical" %in% tmp.ex.int$type){
              tmp.ex.int.var.num <- get_Variables(tmp.ex.int, subset.type = "numerical")
            
            # Check if second int term is "alpha" type  
            } else if("alpha" %in% tmp.ex.int$type){
              tmp.ex.int.var.num <- get_Variables(tmp.ex.int, subset.type = "alpha")
              
            # Check if second int term is "beta" type  
            } else if("beta" %in% tmp.ex.int$type){
            tmp.ex.int.var.num <- get_Variables(tmp.ex.int, subset.type = "beta")
            
            # Error
            } else { print(paste0("Error, unknown variable type (", tmp.ex.int, ") in other")) }
            
            
            # Make a plot
            # print(paste0("Inside of Other > one num var + one cat var "))  # TEST
            tmp.list[[x]] <- gen_scatter_plot(data, 
                                              tmp.ex.int.var.num, 
                                              resp.var, 
                                              facet.var.x = tmp.ex.int.var.cat)
          }

        # Otherwise both int terms are numerical variables
        } else { 
          # print("I dunno how to plot this. Two numerical interaction terms")  # FIX
            
            ### Two Numerical Interaction Variables ###
            
            # Create a list of variables
            tmp.ex.int.vars <- get_Variables(tmp.ex.int)
            
            # print(paste0("int vars: ", tmp.ex.int.vars[1], " and ", tmp.ex.int.vars[2]))  # TEST
            # print(paste0("Inside of Other > two num vars"))  # TEST
            
            # Make a plot
            tmp.list[[paste0(tmp.ex.int.vars[2], ":", tmp.ex.int.vars[1], "_", x, "_")]] <- gen_scatter_plot(data, 
                                                                                                             tmp.ex.int.vars[2], 
                                                                                                             resp.var, 
                                                                                                             # facet.var.x = NULL,
                                                                                                             # facet.var.y = NULL,
                                                                                                             size = tmp.ex.int.vars[1])
            
            # Make additional plot for other way of visualization
            tmp.list.int[[paste0(tmp.ex.int.vars[1], ":", tmp.ex.int.vars[2], "_", x, "_")]] <- gen_scatter_plot(data, 
                                                                                                                 tmp.ex.int.vars[1], 
                                                                                                                 resp.var, 
                                                                                                                 # facet.var.x = NULL,
                                                                                                                 # facet.var.y = NULL,
                                                                                                                 size = tmp.ex.int.vars[2])
            
            
        }
        
      } else {
        print(paste0("Error, unknown variable type (", tmp.ex.int, ") in ?"))
      }
    }
  }
  # return(tmp.list)
  tmp.list[sapply(tmp.list, is.null)] <- NULL # Remove empty lists, temp workaround
  print(paste0("tmp.list.int: ", tmp.list.int))
  tmp.list.int[sapply(tmp.list.int, is.null)] <- NULL # Remove empty lists, temp workaround
  return(append(tmp.list, tmp.list.int))
}


# Generate alpha plots -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

gen_multi_alpha_plots <- function(data, an.type, sig.vars){
  
  if(an.type == "numerical"){
    
  } 
  else if(an.type == "categorical"){
    
  } 
  else if(an.type == "interaction"){
    
  }
  else{
    
  }
  
}




# Generate beta plots -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

gen_multi_beta_plots <- function(data, sig.vars, resp.var, all.vars, an.type = "", axis = NA, stats = NA){
  
  print("Inside of gen_multi_plots > beta > gen_multi_beta_plots")  # TEST
  
    tmp.beta.plot <- { if(check_Var_Type(sig.vars[x,]) == "numerical"){
      
      tmp.list[[tmp.exp.var]] <- gen_PCoA_Num_Plot(data, 
                                         tmp.exp.var, 
                                         resp.var, 
                                         axis,
                                         stats = stats)
      
    } else if(check_Var_Type(sig.vars[x,]) == "categorical") {
      
      tmp.list[[tmp.exp.var]] <- gen_PCoA_Cat_Plot(data, 
                                         tmp.exp.var, 
                                         resp.var, 
                                         axis,
                                         stats = stats)
      
    } else if(check_Var_Type(sig.vars[x,]) == "interaction"){
      
      # print("Inside of gen_multi_plots > beta > interaction")  # TEST
      
      # Extract variables from interaction terms
      tmp.ex.int <- extract_Int_Terms(tmp.exp.var, all.vars)
      # print(paste0("tmp.ex.int: ", tmp.ex.int))  #TEST
      
      # Create a list of variables
      tmp.ex.int.vars <- get_Variables(tmp.ex.int)
      
      # print(paste0("tmp.ex.int.vars: ", tmp.ex.int.vars))  #TEST
      
      # Check if there is a numerical variable in int terms
      # if(tmp.ex.int$type %in% "categorical"){
      if("categorical" %in% tmp.ex.int$type){
        
        # print("Inside of gen_multi_plots > beta > interaction > one cat")  # TEST
        
        # Check if both are categorical, either var can be facet
        if(all(tmp.ex.int$type %in% "categorical")){
          
          # print("Inside of gen_multi_plots > beta > interaction > two cat")  # TEST
          
          # # Create a list of variables
          # tmp.ex.int.vars <- get_Variables(tmp.ex.int)
          
          # Make a plot
          tmp.list[[x]] <- gen_PCoA_Int_Cat_Plot(data, 
                                                 tmp.ex.int.vars[1],
                                                 x.var.2 = tmp.ex.int.vars[2], 
                                                 resp.var, 
                                                 axis,
                                                 stats = stats)
          
          # Make an additional plot to store the reversed facets
          #   - Often one cat variable makes more sense to be the center of focus
          tmp.list.int[[x]] <- gen_PCoA_Int_Cat_Plot(data, 
                                                     tmp.ex.int.vars[2],
                                                     x.var.2 = tmp.ex.int.vars[1], 
                                                     resp.var, 
                                                     axis,
                                                     stats = stats)
          
          # If one var is numerical, then cat var is facet var
        } else {
          # print("Inside of gen_multi_plots > beta > interaction > one num")  # TEST
          
          # Assign cat var
          tmp.ex.int.var.cat <- get_Variables(tmp.ex.int, subset.type = "categorical")
          
          # Assign num var:
          
          # Check if second int term is "numerical" type
          if("numerical" %in% tmp.ex.int$type){
            tmp.ex.int.var.num <- get_Variables(tmp.ex.int, subset.type = "numerical")
            
            # Check if second int term is "alpha" type  
          } else if("alpha" %in% tmp.ex.int$type){
            tmp.ex.int.var.num <- get_Variables(tmp.ex.int, subset.type = "alpha")
            
            # Check if second int term is "beta" type  
          } else if("beta" %in% tmp.ex.int$type){
            tmp.ex.int.var.num <- get_Variables(tmp.ex.int, subset.type = "beta")
            
            # Error
          } else { print("Error, unknown variable type") }
          
          
          # Make a plot
          tmp.list[[x]] <- gen_PCoA_Int_Num_Plot(data, 
                                                 tmp.ex.int.var.cat, 
                                                 resp.var, 
                                                 axis,
                                                 stats = stats,
                                                 x.var.2 = tmp.ex.int.var.num)
        }
        
        # Otherwise both int terms are numerical variables
      } else { 
        print("I dunno how to plot this. Two numerical interaction terms")  # FIX
      }
      
    } else {
      print("Error, unknown data type")
    }
  }
  
  return(unlist(tmp.beta.plot, recursive = F))
}


# Generic scatter plot -----------------------------------------------------------
#   Description: Generates scatter plot
#   Input: data, significant explanatory variables, resp variable, facet conditions
#   Output: scatter plot


gen_scatter_plot <- function(data, x.var, y.var, facet.var.x = ".", facet.var.y = ".", size = NULL){
    
    p <- if(is.null(size)){
        ggplot(data, aes_string(x = x.var, y = y.var)) +
            geom_point(aes(fill = x.var, color = x.var)) + 
            geom_smooth(aes(fill = x.var), 
                        # color = x.var,
                        method="glm", 
                        size = 1, 
                        alpha = .20) + 
            scale_size_area() +
            facet_grid(paste0(facet.var.y, " ~ ", facet.var.x), scales = "free_y") + 
            labs(title = paste0(y.var, " ~ ", x.var)) + 
            scale_fill_brewer(palette = "Dark2") +
            scale_color_brewer(palette = "Dark2") +
            # scale_color_distiller(palette = "Spectral") +
            theme(legend.position = "bottom") + 
            guides(size = guide_legend(title=size), alpha = F, fill=guide_legend(title=NULL), color=guide_legend(title=NULL))
    } else {
        # If 
        ggplot(data, aes_string(x = x.var, y = y.var, size = size
                                , color = size
                                )) +
            geom_point(aes(fill = x.var), alpha = 0.8) + 
            geom_smooth(aes(fill = x.var), 
                        #color = x.var, 
                        method="glm", 
                        size = 1, 
                        alpha = .25) + 
            scale_size_area() +
            facet_grid(paste0(facet.var.y, " ~ ", facet.var.x), scales = "free_y") + 
            labs(title = paste0(y.var, " ~ ", x.var)) + 
            scale_fill_brewer(palette = "Dark2") +
            # scale_color_brewer(palette = "Dark2") +
            # scale_color_distiller(palette = "Spectral") +
            scale_colour_viridis_c() +
            theme(legend.position = "bottom" ) + 
            guides(size = guide_legend(title=size), 
                   alpha = F, 
                   fill=guide_legend(title=NULL), 
                   color=guide_legend(title=size),
                   )
        
    }
 
  return(p)
}

# An improved version of the above for single plot making

gen_scatter_plot_2 <- function(data, x.var, y.var, facet.var.x = ".", facet.var.y = ".", 
                               y.lab = y.var, colors = NA, title = "", caption = ""){
  
  p <- ggplot(data,
              aes_string(x = x.var, y = y.var)) +
    geom_point(aes(fill = x.var ) # converts string to var name
    ) +
    geom_smooth(color = "blue",
                method="glm", 
                size = 1, 
                alpha = .20) + 
    facet_grid(paste0(facet.var.y, " ~ ", facet.var.x), scales = "free_y") + 
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    guides(fill=guide_legend(title=x.var)) +
    # theme(legend.position = "non") +
    labs(
      title = ifelse(title != "", sub_char_var(title, ".", " "), ""),
      caption = ifelse(caption != "", caption, ""),
      y = sub_char_var(y.lab, ".", " "),
      x = sub_char_var(x.var, ".", " ")
    )
  
  return(p)
}



# Generic Box Plots -----------------------------------------
#   Description: Generates box plot 
#   Input: data, significant explanatory variables, resp variable, facet conditions
#   Output: box plot

gen_box_plot <- function(data, x.var, y.var, facet.var.x = ".", facet.var.y = "."){
    print(paste0(x.var, y.var, facet.var.x, facet.var.y))
  p <- ggplot(data,
              aes_string(x = x.var, y = y.var)) +
    geom_violin(aes(fill = eval(parse(text = x.var)) ) # converts string to var name
                , alpha = 0.70  # changes opacity
    ) +
    geom_boxplot(color = "black", fill = NA, width = 0.25) +
    # geom_boxplot(aes(fill = eval(parse(text = x.var)) ) # converts string to var name
    #              , alpha = 0.50  # changes opacity
    # ) +  
    #facet_grid(as.formula(paste(data$facet.var, " ~ .")), scales = "free_y")
    facet_grid(paste0(facet.var.y, " ~ ", facet.var.x), scales = "free_y") + 
    labs(title = paste0(y.var, " ~ ", x.var)) + 
    scale_color_brewer(palette = "Dark2") + 
    scale_fill_brewer(palette = "Dark2") + 
    theme(legend.position = "bottom") + 
      guides(size = F, alpha = F, fill=guide_legend(title=NULL, override.aes = list(alpha = 1)), color=guide_legend(title=NULL, override.aes = list(alpha = 1))) 
  
  return(p)
}

# An improved version of the above for single boxplot making

gen_box_plot_2 <- function(data, x.var, y.var, facet.var.x = ".", facet.var.y = ".", 
                           y.lab = y.var, colors = pal.dark2, title = "", caption = "",
                           comparisons = ""){
  
  p <- ggplot(data,
              aes_string(x = x.var, y = y.var)) +
    geom_boxplot(aes(fill = eval(parse(text = x.var)) ) # converts string to var name
    ) +
    # geom_jitter(color="black", size=0.4, alpha=0.9) +
    geom_quasirandom() + 
    facet_grid(as.formula(paste0(facet.var.y, " ~ ", facet.var.x)), scales = "free_y") + 
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    guides(fill=guide_legend(title=x.var)) +
    # theme(legend.position = "bottom") + 
    labs(
      title = ifelse(title != "", sub_char_var(title, ".", " "), ""),
      caption = ifelse(caption != "", caption, ""),
      y = sub_char_var(y.lab, ".", " "),
      x = sub_char_var(x.var, ".", " ")
    ) #+
    # ggpubr::stat_compare_means(comparisons = comparisons, 
    #                            label = "p.signif" ) + 
    # scale_y_continuous(breaks = scales::breaks_pretty(n = 8), limits = c(0, 1.1)) 

  return(p)
}


# Bar Plot ----------------------------------------------------------------



gen_bar_plot <- function(data, x.var, y.var = "Count", facet.var.x = ".", facet.var.y = ".",
                         y.lab = y.var, colors = pal.dark2, title = "", caption = ""){
  
  p <- ggplot(data,
              aes_string(x = x.var)) +
    geom_bar(aes(fill = eval(parse(text = x.var)) ), # converts string to var name
              position = "dodge") +
    facet_grid(paste0(facet.var.y, " ~ ", facet.var.x), scales = "fixed") + 
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    guides(fill=guide_legend(title=x.var)) +
    # theme(legend.position = "bottom") + 
    labs(
      title = ifelse(title != "", sub_char_var(title, ".", " "), ""),
      caption = ifelse(caption != "", caption, ""),
      y = sub_char_var(y.lab, ".", " "),
      x = sub_char_var(x.var, ".", " ")
    ) 
  
  return(p)
}

# Generic PCoA Plots -----------------------------------------
#   Description: Generates box plot 
#   Input: data, significant explanatory variables, resp variable, facet conditions
#   Output: box plot

gen_PCoA_Cat_Plot <- function(data, x.var, y.var, axis, stats){
  
  print("Inside: gen_PCoA_Cat_Plot()")  # TEST
  
  centroid <- calc_Centroid(data, x.var)
  
  p <- {lapply(names(y.var), function(beta){
    
      # Temp list for storing ggplots as they are created
      tmp.p <- list()
          
      # Extract statistic
      statistic = round(stats[which(stats[2] == x.var &
                                      stats[1] == beta),]$statistic, 3)
      # Extract p.value
      p.val = stats[which(stats[2] == x.var &
                            stats[1] == beta),]$p.value
      
      if(length(p.val != 0)){
        p.val <- p.val
        statistic <- statistic
      } else {
        p.val <- ">> 0.05"
        statistic <- "..."
      }

      tmp.p[[beta]] <- ggplot(data[Dist == beta], aes(x = data[Dist == beta][[2]], y = data[Dist == beta][[3]])) +
        geom_point(aes(color = eval(parse(text = x.var)),
                       fill = eval(parse(text = x.var)),
                       size = .5), 
                        alpha = 0.5) +
        scale_size_area() +
        stat_ellipse(aes(color = eval(parse(text = x.var)) )) +
        geom_label_repel(data = centroid[Dist == beta],
                         aes(label = eval(parse(text = x.var)) ), 
                         alpha = 0.80) +
        labs(title = paste0(beta, ": Beta diversity ~ ", x.var),
          caption = paste0("Statistic = ", statistic ,", P-value = ", p.val),
          x = axis[Dist == beta]$X.lab,
          y = axis[Dist == beta]$Y.lab) + 
        theme(legend.position="bottom") + guides(size = F, alpha = F, fill=guide_legend(title=x.var), color=guide_legend(title=x.var)) + 
        scale_color_brewer(palette = "Dark2") + 
        scale_fill_brewer(palette = "Dark2")
      
      
      ggMarginal(tmp.p[[beta]], size = 7, colour = NA, groupFill = T, groupColour = T, type = "density")
    })
  }
    
  return(p)
}

gen_PCoA_Cat_Plot_2 <- function(
  data, 
  x.var, 
  y.var, 
  axis, 
  stats, 
  x.var.2 = "", 
  elps.var = x.var, 
  title = "", 
  caption = "", 
  disp.stats = F,
  colors = pal.dark2
  ){
  
  print("Inside: gen_PCoA_Cat_Plot_2()")  # TEST
  
  # centroid <- calc_Centroid(data, x.var)
  
  p <- {lapply(names(y.var), function(beta){
    
    # If disp.stats == T:
    if(isTRUE(disp.stats)){
      # Extract statistic
      statistic = round(stats[which(stats[2] == x.var &
                                      stats[1] == beta),]$statistic, 3)
      # Extract p.value
      p.val = stats[which(stats[2] == x.var &
                            stats[1] == beta),]$p.value
      
      if(length(p.val != 0)){
        p.val <- p.val
        statistic <- statistic
      } else {
        p.val <- ">> 0.05"
        statistic <- NULL
      }
    }
    
    # Temp list for storing ggplots as they are created
    tmp.p <- list()
    
    tmp.p[[beta]] <- ggplot(data[Dist == beta], aes(x = CAP1, y = data[Dist == beta][[3]])) +
      # geom_point(aes(#color = eval(parse(text = x.var)),
      #   fill = eval(parse(text = x.var)),
      #   shape = eval(parse(text = x.var.2))#,
      #   #size = 1, 
      #   #alpha = 1
      #   ), color = "black", alpha = .75) +
      geom_point(aes(color = eval(parse(text = x.var)),
                     fill = eval(parse(text = x.var)),
                     shape = eval(parse(text = x.var.2)),
                     size = 1, 
                     alpha = .9)
                 ) +
      # geom_text(aes(label = eval(parse(text = label.var)) )) +
      stat_ellipse(aes(color = eval(parse(text = x.var)) )) +
      # stat_ellipse(aes(linetype = eval(parse(text = elps.var)) )) +
      # geom_label_repel(data = centroid[Dist == beta],
      #                  aes(label = eval(parse(text = x.var)),
      #                      size = 2), alpha = 0.80) +
      labs(
        # title = ifelse(title == "", paste0(beta, " ~ ", paste0(c(x.var, x.var.2), collapse = " + ")), title),
        # title = ifelse(title == "", paste0("Beta Score ~ ", paste0(c(sub_char_var(x.var, ".", " "), sub_char_var(x.var.2, ".", " ")), collapse = " + ")), sub_char_var(title, ".", " ")),
        title = ifelse(title == "", paste0("Beta Score ~ ", ifelse(x.var.2 == "", sub_char_var(x.var, ".", " "), paste0(c(sub_char_var(x.var, ".", " "), sub_char_var(x.var.2, ".", " ")), collapse = " + ")) ), sub_char_var(title, ".", " ")),
        caption = ifelse(isTRUE(disp.stats), paste0("Statistic = ", statistic ,", P-value = ", p.val, ", capscale(), ", beta), caption),
        x = axis[Dist == beta]$X.lab,
        y = axis[Dist == beta]$Y.lab) + 
      guides(fill=guide_legend(title=x.var), color = guide_legend(title=x.var, override.aes = list(size=3)), 
             # shape = guide_legend(title=x.var.2,
             #                      override.aes = list(color = c(col.exposure[2], col.exposure[1]),
             #                                          size = 3)),
             shape = guide_legend(title=x.var.2, 
                                  override.aes = list(alpha = 1, size = 3)),
            linetype=guide_legend(title="Ellipses (95% CI)"),
             size = F,
             alpha = F) +
      theme(legend.position = "bottom") + 
      # scale_fill_brewer(palette = colors) +
      # scale_color_brewer(palette = colors)
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) 
      # scale_shape_manual(values = c(16, 15), labels = c("exposed", "unexposed"))
    
    
    
    # ggMarginal(tmp.p[[beta]], size = 7, colour = NA, groupFill = T, groupColour = T, type = "density")
  })
  }
  
  return(p)
}

## Adds vectors

# gen_PCoA_Cat_Plot_vector <- function(
#   data,
#   mod,
#   x.var, 
#   y.var, 
#   axis, 
#   stats, 
#   x.var.2 = "", 
#   elps.var = x.var, 
#   title = "", 
#   caption = "", 
#   disp.stats = F,
#   colors = pal.dark2
# ){
#   
#   print("Inside: gen_PCoA_Cat_Plot_vector()")  # TEST
#   
#   vars <- c(x.var, x.var.2)
#   
#   # centroid <- calc_Centroid(data, x.var)
#   
#   p <- {lapply(names(y.var), function(beta){
#     
#     # If disp.stats == T:
#     if(isTRUE(disp.stats)){
#       # Extract statistic
#       statistic = round(stats[which(stats[2] == x.var &
#                                       stats[1] == beta),]$statistic, 3)
#       # Extract p.value
#       p.val = stats[which(stats[2] == x.var &
#                             stats[1] == beta),]$p.value
#       
#       if(length(p.val != 0)){
#         p.val <- p.val
#         statistic <- statistic
#       } else {
#         p.val <- ">> 0.05"
#         statistic <- NULL
#       }
#     }
#     
#     # VECTORS
#     smpl.ord <- data.table()
#     arws.ord <- data.table()
#     cntr.ord <- data.table()
#     
#     permanova <- subset(stats, metric == beta)
#     
#     # Coords, Vectors, Centriods
#     dbrda.data <- phyloseqCompanion::get.biplot.data(
#       ps = ps.obj,
#       ord = mod[[beta]]
#     )
#     
#     dbrda.data$sample.coords[, Dist := beta]
#     dbrda.data$vector.coords[, Dist := beta]
#     dbrda.data$centroid.coords[, Dist := beta]
#     
#     smpl.ord <- rbind(smpl.ord, dbrda.data$sample.coords)[Dist == beta]
#     arws.ord <- rbind(arws.ord, dbrda.data$vector.coords)[Dist == beta]
#     cntr.ord <- rbind(cntr.ord, dbrda.data$centroid.coords)[Dist == beta]
# 
#     
#     # Extract Significant Variables
#     perm.dt <- as.data.table(permanova)
#     sig.intrxn <- perm.dt[str_detect(term, ":") & p.value <= 0.05]$term
#     sig.intrxn.patterns <- sapply(sig.intrxn, function(term) { {str_split(term, ":")[[1]]} %>% paste(collapse = ".+:") })
#     sig.main <- perm.dt[!str_detect(term, ":") & p.value <= 0.05]$term
#     sig.vector.coords <- if(length(sig.intrxn.patterns) == 0){
#       arws.ord[!str_detect(Variable, ":") & str_detect(Variable, paste(sig.main, collapse = "|"))]
#     } else {
#       rbind(
#         arws.ord[!str_detect(Variable, ":") & str_detect(Variable, paste(sig.main, collapse = "|"))],
#         arws.ord[str_detect(Variable, paste(sig.intrxn.patterns, collapse = "|"))]
#       )
#     }
#     
#     ## Remove variable prefixes
#     lapply(vars, function(tmp.var){
#       sig.vector.coords[, Variable := sub(tmp.var, "", Variable)]
#       # arws.ord[, Variable := sub(tmp.var, "", Variable)]
#     })
#     
#     ## Relevel factors
#     levels(cntr.ord$Variable) <- factor(levels(cntr.ord$Variable), levels = cntr.ord$Variable)
#     
#     
#     # Axis labels
#     # var.exp <- {
#     #   pc.exp <- mod[[beta]]$CCA$eig / sum(mod[[beta]]$CCA$eig)
#     #   unname(paste0(round(pc.exp[1:2] * 100, 1), "%")) 
#     # }
#     # 
#     # xlab <- paste0("CAP1 (", var.exp[1], ")")
#     # ylab <- paste0("CAP2 (", var.exp[2], ")")
#     
#     # Temp list for storing ggplots as they are created
#     tmp.p <- list()
#     
#     tmp.p[[beta]] <- ggplot(smpl.ord[Dist == beta], aes(x = CAP1, y = smpl.ord[Dist == beta][[3]])) +
#       geom_point(aes(color = eval(parse(text = x.var)),
#                      fill = eval(parse(text = x.var)),
#                      shape = eval(parse(text = x.var.2)) 
#                      # alpha = .975,
#                      # size = 1
#                      )
#       ) +
#       # geom_point(data = cntr.ord, aes(color = Variable), shape = 10, size = 5) +  # Doesn't work well with ggplot
#       stat_ellipse(aes(color = eval(parse(text = x.var)) )) +
#       geom_segment(
#         data = sig.vector.coords,
#         x = 0, y = 0,
#         aes(xend = CAP1, yend = CAP2),
#         arrow = arrow(length = unit(0.05, "npc")),
#       ) +
#       geom_label_repel(
#         data = sig.vector.coords,
#         aes(x = CAP1, y = CAP2, label = Variable),
#         alpha = .75,
#         size = 3
#       ) +
#       labs(
#         title = ifelse(title == "", paste0("Beta Score ~ ", ifelse(x.var.2 == "", sub_char_var(x.var, ".", " "), paste0(c(sub_char_var(x.var, ".", " "), sub_char_var(x.var.2, ".", " ")), collapse = " + ")) ), sub_char_var(title, ".", " ")),
#         caption = ifelse(isTRUE(disp.stats), paste0("Statistic = ", statistic ,", P-value = ", p.val, ", capscale(), ", beta), caption),
#         x = axis[Dist == beta]$X.lab,
#         y = axis[Dist == beta]$Y.lab) +
#       guides(fill=guide_legend(title=x.var), color = guide_legend(title=x.var, override.aes = list(size=2)),
#              shape = guide_legend(title=x.var.2,
#                                   override.aes = list(alpha = 1, size = 3)),
#              linetype=guide_legend(title="Ellipses (95% CI)"),
#              size = F,
#              alpha = F) +
#       theme(legend.position = "bottom") +
#       scale_fill_manual(values = colors) +
#       scale_color_manual(values = colors)
#   })
#   }
#   
#   return(p)
# }


# Vectors 2

gen_PCoA_Cat_Plot_vector <- function(
  data,
  mod,
  x.var, 
  y.var, 
  axis, 
  stats, 
  x.var.2 = "", 
  shape.var = "",
  size.var = 3.5,
  elps.var = x.var, 
  title = "", 
  caption = "", 
  disp.stats = F,
  colors = pal.dark2
){
  
  print("Inside: gen_PCoA_Cat_Plot_vector()")  # TEST
  
  vars <- c(x.var, x.var.2)
  
  # centroid <- calc_Centroid(data, x.var)
  
  p <- {lapply(names(y.var), function(beta){
    
    # If disp.stats == T:
    if(isTRUE(disp.stats)){
      # Extract statistic
      statistic = round(stats[which(stats[2] == x.var &
                                      stats[1] == beta),]$statistic, 3)
      # Extract p.value
      p.val = stats[which(stats[2] == x.var &
                            stats[1] == beta),]$p.value
      
      if(length(p.val != 0)){
        p.val <- p.val
        statistic <- statistic
      } else {
        p.val <- ">> 0.05"
        statistic <- NULL
      }
    }
    
    # VECTORS
    smpl.ord <- data.table()
    arws.ord <- data.table()
    cntr.ord <- data.table()
    
    permanova <- subset(stats, metric == beta)
    
    # Coords, Vectors, Centriods
    dbrda.data <- phyloseqCompanion::get.biplot.data(
      ps = ps.obj,
      ord = mod[[beta]]
    )
    
    dbrda.data$sample.coords[, Dist := beta]
    dbrda.data$vector.coords[, Dist := beta]
    dbrda.data$centroid.coords[, Dist := beta]
    
    smpl.ord <- rbind(smpl.ord, dbrda.data$sample.coords)[Dist == beta]
    arws.ord <- rbind(arws.ord, dbrda.data$vector.coords)[Dist == beta]
    cntr.ord <- rbind(cntr.ord, dbrda.data$centroid.coords)[Dist == beta]
    
    
    # Extract Significant Variables
    perm.dt <- as.data.table(permanova)
    sig.intrxn <- perm.dt[str_detect(term, ":") & p.value <= 0.05]$term
    sig.intrxn.patterns <- sapply(sig.intrxn, function(term) { {str_split(term, ":")[[1]]} %>% paste(collapse = ".+:") })
    sig.main <- perm.dt[!str_detect(term, ":") & p.value <= 0.05]$term
    sig.vector.coords <- if(length(sig.intrxn.patterns) == 0){
      arws.ord[!str_detect(Variable, ":") & str_detect(Variable, paste(sig.main, collapse = "|"))]
    } else {
      rbind(
        arws.ord[!str_detect(Variable, ":") & str_detect(Variable, paste(sig.main, collapse = "|"))],
        arws.ord[str_detect(Variable, paste(sig.intrxn.patterns, collapse = "|"))]
      )
    }
    
    ## Remove variable prefixes
    lapply(vars, function(tmp.var){
      sig.vector.coords[, Variable := sub(tmp.var, "", Variable)]
      cntr.ord[, Variable := sub(tmp.var, "", Variable)]
    })
    
    ## Relevel factors
    levels(cntr.ord$Variable) <- factor(levels(cntr.ord$Variable), levels = cntr.ord$Variable)
    
    
    # Temp list for storing ggplots as they are created
    tmp.p <- list()
    
    tmp.p[[beta]] <- ggplot(smpl.ord[Dist == beta], aes(x = CAP1, y = CAP2)) + 
      geom_point(aes(color = eval(parse(text = x.var)), 
                     shape = eval(parse(text = ifelse(shape.var == "", x.var.2, shape.var))),
                     size = eval(parse(text = size.var))
                     ),
                 alpha = 0.75) +
      geom_point(data = cntr.ord, aes(color = Variable), stroke = 1.5, shape = 13, size = 4, alpha = .95, show.legend = F) +
      geom_segment(
        data = sig.vector.coords, 
        x = 0, y = 0,
        aes(xend = CAP1, yend = CAP2),
        arrow = arrow(length = unit(0.05, "npc")),
      ) +
      stat_ellipse(aes(color = eval(parse(text = x.var)) )) +
      geom_label_repel(
        data = sig.vector.coords, 
        aes(x = CAP1, y = CAP2, label = Variable),
        alpha = .75,
        size = 3
      ) +
      scale_color_manual(name = "Variable", values = colors) +
      theme(
        legend.position = "bottom", 
        legend.box = "vertical",
        legend.box.just = "center"
      ) +
      guides(
        color = guide_legend(order = 1),
        shape = guide_legend(title=x.var.2),
        size = F
        ) + 
      labs(
          title = ifelse(title == "", paste0("Beta Score ~ ", ifelse(x.var.2 == "", sub_char_var(x.var, ".", " "), paste0(c(sub_char_var(x.var, ".", " "), sub_char_var(x.var.2, ".", " ")), collapse = " + ")) ), sub_char_var(title, ".", " ")),
          caption = ifelse(isTRUE(disp.stats), paste0("Statistic = ", statistic ,", P-value = ", p.val, ", capscale(), ", beta), caption),
          x = axis[Dist == beta]$X.lab,
          y = axis[Dist == beta]$Y.lab
          )
    
  })
  }
  
  return(p)
}

# PCOA Num Plot -----------------------------------------------------------



gen_PCoA_Num_Plot <- function(data, x.var, y.var, axis, stats){
  
  print("Inside: gen_PCoA_Num_Plot()")  # TEST
  print(paste0("gen_PCoA_Num_Plot: x.var: ", x.var, " y.var: ", y.var))  # TEST
  
  centroid <- calc_Centroid(data, x.var)
  
  p <- {lapply(names(y.var), function(beta){
    
    # Temp list for storing ggplots as they are created
    tmp.p <- list()
    
    # Extract statistic
    statistic = round(stats[which(stats[2] == x.var &
                                    stats[1] == beta),]$statistic, 3)
    # Extract p.value
    p.val = stats[which(stats[2] == x.var &
                          stats[1] == beta),]$p.value
    
    # exp.var[exp.var != resp]  # Needs work. Select the highest significant variable that's not the current x.var
    high.stat = stats[which(stats$statistic == max(stats$statistic, na.rm = T)),]$term
    
    if(length(p.val != 0)){
      p.val <- p.val
      statistic <- statistic
    } else {
      p.val <- ">> 0.05"
      statistic <- "..."
    }
    
    tmp.p[[beta]] <- ggplot(data[Dist == beta], aes(x = data[Dist == beta][[2]], y = data[Dist == beta][[3]])) +
                geom_point(aes(color = eval(parse(text = x.var)),
                               fill = eval(parse(text = x.var)),
                               size = eval(parse(text = x.var))), 
                               alpha = 0.5) +
                scale_size_area() +
                stat_ellipse(aes(linetype = eval(parse(text = high.stat)) )) +
                geom_label_repel(data = calc_Centroid(data, high.stat)[Dist == beta],
                                aes(label = eval(parse(text = high.stat)) ), 
                                alpha = 0.80) +
                labs(title = paste0(beta, ": Beta diversity ~ ", x.var),
                     caption = paste0("Statistic = ", statistic ,", P-value = ", p.val),
                     x = axis[Dist == beta]$X.lab,
                     y = axis[Dist == beta]$Y.lab) + 
      theme(legend.position="bottom") + 
      guides(size = F, 
             alpha = F, 
             fill=guide_legend(title=x.var), 
             color=guide_legend(title=x.var), 
             linetype=guide_legend(title=high.stat)) + 
      scale_color_distiller(palette = "Spectral") + 
      scale_fill_distiller(palette = "Spectral")
    
    ggMarginal(tmp.p[[beta]], fill = "grey", type = "density")
    
  })
  }
  
  return(p)
}


gen_PCoA_Num_Plot_2 <- function(data, x.var, y.var, axis, stats, x.var.2 = "", 
                                elps.var = x.var, title = "", disp.stats = F, 
                                colors = "RdYlBu"){
  
  print("Inside: gen_PCoA_Num_Plot_2()")  # TEST
  
  # centroid <- calc_Centroid(data, x.var)
  
  
  
  p <- {lapply(names(y.var), function(beta){
    
    # If disp.stats == T:
    if(isTRUE(disp.stats)){
      # Extract statistic
      statistic = round(stats[which(stats[2] == x.var &
                                      stats[1] == beta),]$statistic, 3)
      # Extract p.value
      p.val = stats[which(stats[2] == x.var &
                            stats[1] == beta),]$p.value
      
      if(length(p.val != 0)){
        p.val <- p.val
        statistic <- statistic
      } else {
        p.val <- ">> 0.05"
        statistic <- NULL
      }
    }
    
    # Temp list for storing ggplots as they are created
    tmp.p <- list()
    
    tmp.p[[beta]] <- ggplot(data[Dist == beta], aes(x = CAP1, y = data[Dist == beta][[3]])) +
      # geom_point(aes(#color = eval(parse(text = x.var)),
      #                # fill = eval(parse(text = x.var)),
      #                shape = eval(parse(text = x.var.2)),
      #                size = (eval(parse(text = x.var)) +.25)), color = "black",
      #            alpha = .9) +
      geom_point(aes(color = eval(parse(text = x.var)),
                     # fill = eval(parse(text = x.var)),
                     shape = eval(parse(text = x.var.2)),
                     size = eval(parse(text = x.var)) ),
                 alpha = .9) +
      scale_size_area() +
      stat_ellipse(aes(linetype = eval(parse(text = elps.var)) )) +
      # geom_label_repel(data = centroid[Dist == beta],
      #                  aes(label = eval(parse(text = x.var)),
      #                      size = 2), alpha = 0.80) +
      labs(
        # title = ifelse(title == "", paste0(beta, " ~ ", paste0(c(x.var, x.var.2), collapse = " + ")), title),
        # title = ifelse(title == "", paste0("Beta Score ~ ", paste0(c(sub_char_var(x.var, ".", " "), sub_char_var(x.var.2, ".", " ")), collapse = " + ")), sub_char_var(title, ".", " ")),
        title = ifelse(title == "", paste0("Beta Score ~ ", ifelse(x.var.2 == "", sub_char_var(x.var, ".", " "), paste0(c(sub_char_var(x.var, ".", " "), sub_char_var(x.var.2, ".", " ")), collapse = " + ")) ), sub_char_var(title, ".", " ")),
        # caption = ifelse(isTRUE(disp.stats), paste0("Statistic = ", statistic ,", P-value = ", p.val),""),
        caption = ifelse(isTRUE(disp.stats), paste0("Statistic = ", statistic ,", P-value = ", p.val, ", capscale(), ", beta), caption),
        x = axis[Dist == beta]$X.lab,
        y = axis[Dist == beta]$Y.lab) +  
      guides(size = F, 
             alpha = F, 
             # fill=guide_legend(title=x.var), 
             color=guide_legend(title=x.var, 
                                override.aes= list(alpha = 1, size = 3)
                                ),
             shape=guide_legend(title=x.var.2, 
                                override.aes= list(alpha = 1, size = 3)
                                ), 
             linetype=guide_legend(title="Ellipses (95% CI)")
             ) +
      theme(legend.position = "bottom") + 
      # scale_linetype_manual(name = "Ellipses", values = c([1:n])) +
      # scale_linetype_manual(name = "Ellipses", values = c([1:n])) +
      # scale_fill_brewer(palette = "Dark2") +
      # scale_color_brewer(palette = "Dark2") 
      # scale_color_distiller(palette = colors ) #scales::rescale(c(0,.5,1)))  
      scale_colour_gradient(
        low = pal.RdPu[3],
        high = pal.paired[10],
        space = "Lab",
        na.value = "grey50",
        guide = "colourbar",
        aesthetics = "colour"
      )
    
      # scale_fill_distiller(palette = "PuOr")
    
    
    # ggMarginal(tmp.p[[beta]], size = 7, colour = NA, groupFill = T, groupColour = T, type = "density")
  })
  }
  
  return(p)
}

# adds vectors


gen_PCoA_Num_Plot_vector <- function(data, 
                                     mod,
                                     x.var, 
                                     y.var, 
                                     axis, 
                                     stats, 
                                     x.var.2 = "", 
                                     elps.var = x.var, 
                                     title = "", 
                                     disp.stats = F, 
                                     colors = "RdYlBu"){
  
  print("Inside: gen_PCoA_Num_Plot_vector()")  # TEST
  
  # centroid <- calc_Centroid(data, x.var)
  
  vars <- c(x.var.2)
  
  p <- {lapply(names(y.var), function(beta){
    
    # If disp.stats == T:
    if(isTRUE(disp.stats)){
      # Extract statistic
      statistic = round(stats[which(stats[2] == x.var &
                                      stats[1] == beta),]$statistic, 3)
      # Extract p.value
      p.val = stats[which(stats[2] == x.var &
                            stats[1] == beta),]$p.value
      
      if(length(p.val != 0)){
        p.val <- p.val
        statistic <- statistic
      } else {
        p.val <- ">> 0.05"
        statistic <- NULL
      }
    }
    
    # VECTORS
    smpl.ord <- data.table()
    arws.ord <- data.table()
    cntr.ord <- data.table()
    
    permanova <- subset(stats, metric == beta)
    
    # Coords, Vectors, Centriods
    dbrda.data <- phyloseqCompanion::get.biplot.data(
      ps = ps.obj,
      ord = mod[[beta]]
    )
    
    dbrda.data$sample.coords[, Dist := beta]
    dbrda.data$vector.coords[, Dist := beta]
    dbrda.data$centroid.coords[, Dist := beta]
    
    smpl.ord <- rbind(smpl.ord, dbrda.data$sample.coords)[Dist == beta]
    arws.ord <- rbind(arws.ord, dbrda.data$vector.coords)[Dist == beta]
    cntr.ord <- rbind(cntr.ord, dbrda.data$centroid.coords)[Dist == beta]
    
    
    # Extract Significant Variables
    perm.dt <- as.data.table(permanova)
    sig.intrxn <- perm.dt[str_detect(term, ":") & p.value <= 0.05]$term
    sig.intrxn.patterns <- sapply(sig.intrxn, function(term) { {str_split(term, ":")[[1]]} %>% paste(collapse = ".+:") })
    sig.main <- perm.dt[!str_detect(term, ":") & p.value <= 0.05]$term
    sig.vector.coords <- if(length(sig.intrxn.patterns) == 0){
      arws.ord[!str_detect(Variable, ":") & str_detect(Variable, paste(sig.main, collapse = "|"))]
    } else {
      rbind(
        arws.ord[!str_detect(Variable, ":") & str_detect(Variable, paste(sig.main, collapse = "|"))],
        arws.ord[str_detect(Variable, paste(sig.intrxn.patterns, collapse = "|"))]
      )
    }
    
    ## Remove variable prefixes
    lapply(vars, function(tmp.var){
      sig.vector.coords[, Variable := sub(tmp.var, "", Variable)]
      cntr.ord[, Variable := sub(tmp.var, "", Variable)]
    })
    
    ## Relevel factors
    levels(cntr.ord$Variable) <- factor(levels(cntr.ord$Variable), levels = cntr.ord$Variable)
    
    
    # Temp list for storing ggplots as they are created
    tmp.p <- list()
      
    tmp.p[[beta]] <- ggplot(smpl.ord[Dist == beta], aes(x = CAP1, y = CAP2))  +
      geom_point(aes(color = eval(parse(text = x.var)), 
                     shape = eval(parse(text = x.var.2)),
                     size = eval(parse(text = x.var)) 
                     ),
                 alpha = 0.75) +
      geom_point(data = cntr.ord, color = "black", stroke = 1.5, shape = 13, size = 4, alpha = .95, show.legend = F) +
      stat_ellipse(aes(linetype = eval(parse(text = elps.var)) )) +
      geom_segment(
        data = sig.vector.coords, 
        x = 0, y = 0,
        aes(xend = CAP1, yend = CAP2),
        arrow = arrow(length = unit(0.05, "npc")),
      ) +
      geom_label_repel(
        data = sig.vector.coords, 
        aes(x = CAP1, y = CAP2, label = Variable),
        alpha = .75,
        size = 3
      ) +
      # scale_color_manual(name = "Variable", values = colors) +
      theme(
        legend.position = "bottom", 
        legend.box = "vertical",
        legend.box.just = "center"
      ) +
      guides(size = F, 
             alpha = F, 
             color=guide_legend(title=x.var, 
                                override.aes= list(alpha = 1, size = 3)
             ),
             shape=guide_legend(title=x.var.2, 
                                override.aes= list(alpha = 1, size = 3)
             ), 
             linetype=guide_legend(title="Ellipses (95% CI)")
      ) +
      labs(
        title = ifelse(title == "", paste0("Beta Score ~ ", ifelse(x.var.2 == "", sub_char_var(x.var, ".", " "), paste0(c(sub_char_var(x.var, ".", " "), sub_char_var(x.var.2, ".", " ")), collapse = " + ")) ), sub_char_var(title, ".", " ")),
        caption = ifelse(isTRUE(disp.stats), paste0("Statistic = ", statistic ,", P-value = ", p.val, ", capscale(), ", beta), caption),
        x = axis[Dist == beta]$X.lab,
        y = axis[Dist == beta]$Y.lab
      ) +
      scale_colour_gradient(
        low = pal.RdPu[3],
        high = pal.paired[10],
        space = "Lab",
        na.value = "grey50",
        guide = "colourbar",
        aesthetics = "colour"
      )
    
  })
  }
  
  return(p)
}


#

gen_PCoA_Num_Plot_3 <- function(data, x.var, y.var, axis, stats, x.var.2 = "", 
                                elps.var = x.var, title = "", disp.stats = F, 
                                colors = "RdYlBu"){
  
  print("Inside: gen_PCoA_Num_Plot_2()")  # TEST
  
  # centroid <- calc_Centroid(data, x.var)
  
  
  
  p <- {lapply(names(y.var), function(beta){
    
    # If disp.stats == T:
    if(isTRUE(disp.stats)){
      # Extract statistic
      statistic = round(stats[which(stats[2] == x.var &
                                      stats[1] == beta),]$statistic, 3)
      # Extract p.value
      p.val = stats[which(stats[2] == x.var &
                            stats[1] == beta),]$p.value
      
      if(length(p.val != 0)){
        p.val <- p.val
        statistic <- statistic
      } else {
        p.val <- ">> 0.05"
        statistic <- NULL
      }
    }
    
    # Temp list for storing ggplots as they are created
    tmp.p <- list()
    
    tmp.p[[beta]] <- ggplot(data[Dist == beta], aes(x = CAP1, y = data[Dist == beta][[3]])) +
      geom_point(aes(#color = eval(parse(text = x.var)),
        # fill = eval(parse(text = x.var)),
        # shape = eval(parse(text = x.var.2)),
        size = (eval(parse(text = x.var.2)) +.25)), color = "black",
        alpha = .9) +
      geom_point(aes(color = eval(parse(text = x.var)),
                     # fill = eval(parse(text = x.var)),
                     # shape = eval(parse(text = x.var.2)),
                     size = eval(parse(text = x.var.2)) ),
                 alpha = 1) +
      scale_size_area() +
      stat_ellipse(aes(linetype = eval(parse(text = elps.var)) )) +
      # geom_label_repel(data = centroid[Dist == beta],
      #                  aes(label = eval(parse(text = x.var)),
      #                      size = 2), alpha = 0.80) +
      labs(
        # title = ifelse(title == "", paste0(beta, " ~ ", paste0(c(x.var, x.var.2), collapse = " + ")), title),
        # title = ifelse(title == "", paste0("Beta Score ~ ", paste0(c(sub_char_var(x.var, ".", " "), sub_char_var(x.var.2, ".", " ")), collapse = " + ")), sub_char_var(title, ".", " ")),
        title = ifelse(title == "", paste0("Beta Score ~ ", ifelse(x.var.2 == "", sub_char_var(x.var, ".", " "), paste0(c(sub_char_var(x.var, ".", " "), sub_char_var(x.var.2, ".", " ")), collapse = " + ")) ), sub_char_var(title, ".", " ")),
        # caption = ifelse(isTRUE(disp.stats), paste0("Statistic = ", statistic ,", P-value = ", p.val),""),
        caption = ifelse(isTRUE(disp.stats), paste0("Statistic = ", statistic ,", P-value = ", p.val, ", capscale(), ", beta), caption),
        x = axis[Dist == beta]$X.lab,
        y = axis[Dist == beta]$Y.lab) +  
      guides(size = guide_legend(title=x.var.2), #size = F, 
             alpha = F, 
             # fill=guide_legend(title=x.var), 
             color=guide_legend(title=x.var, 
                                override.aes= list(alpha = 1, size = 3)
             ),
             shape=guide_legend(title=x.var.2, 
                                override.aes= list(alpha = 1, size = 3)
             ), 
             linetype=guide_legend(title="Ellipses (95% CI)")
      ) +
      theme(legend.position = "bottom") + 
      # scale_linetype_manual(name = "Ellipses", values = c([1:n])) +
      # scale_linetype_manual(name = "Ellipses", values = c([1:n])) +
      # scale_fill_brewer(palette = "Dark2") +
      # scale_color_brewer(palette = "Dark2") 
      # scale_color_distiller(palette = colors ) #scales::rescale(c(0,.5,1)))  
      scale_colour_gradient(
        low = pal.RdPu[3],
        high = pal.paired[10],
        space = "Lab",
        na.value = "grey50",
        guide = "colourbar",
        aesthetics = "colour"
      )
    
    # scale_fill_distiller(palette = "PuOr")
    
    
    # ggMarginal(tmp.p[[beta]], size = 7, colour = NA, groupFill = T, groupColour = T, type = "density")
  })
  }
  
  return(p)
}


# PCOA Int Cat Plot -------------------------------------------------------



gen_PCoA_Int_Cat_Plot <- function(data, x.var, x.var.2, y.var, axis, stats){
  
  print("Inside: gen_PCoA_Int_Cat_Plot()")  # TEST
  
  centroid <- calc_Centroid(data, x.var, x.var.2 = x.var.2)
  
  p <- {lapply(names(y.var), function(beta){
    
    # Temp list for storing ggplots as they are created
    tmp.p <- list()
    
    # Extract statistic
    statistic = round(stats[which(stats[2] == paste0(x.var, ":", x.var.2) &
                                    stats[1] == beta),]$statistic, 3)
    # Extract p.value
    p.val = stats[which(stats[2] == paste0(x.var, ":", x.var.2) &
                          stats[1] == beta),]$p.value
    
    high.stat = stats[which(stats$statistic == max(stats$statistic, na.rm = T)),]$term
    
    if(length(p.val != 0)){
      p.val <- p.val
      statistic <- statistic
    } else {
      p.val <- ">> 0.05"
      statistic <- "..."
    }
    
    tmp.p[[beta]] <- ggplot(data[Dist == beta], aes(x = data[Dist == beta][[2]], y = data[Dist == beta][[3]])) +
      geom_point(aes(
                    color = eval(parse(text = x.var)),
                     fill = eval(parse(text = x.var)),
                     shape = eval(parse(text = x.var.2)),
                     size = .5), 
                 alpha = 0.5) +
      scale_size_area() +
      stat_ellipse(aes(color = eval(parse(text = x.var.2)) )) +
      geom_label_repel(data = centroid[Dist == beta],
                          aes(label = paste0(eval(parse(text = x.var)), ":",
                                            eval(parse(text = x.var.2)) )), 
                       alpha = 0.80) +
      labs(title = paste0(beta, ": Beta diversity ~ ", x.var, ":", x.var.2),
          caption = paste0("Statistic = ", statistic ,", P-value = ", p.val),
           x = axis[Dist == beta]$X.lab,
           y = axis[Dist == beta]$Y.lab) + 
      theme(legend.position="bottom") + 
      guides(size = F, alpha = F, 
             fill=guide_legend(title=NULL), 
             color=guide_legend(title=NULL), 
             linetype=guide_legend(title=x.var.2),
             shape=guide_legend(title=x.var.2),) + 
      scale_color_brewer(palette = "Dark2") + 
      scale_fill_brewer(palette = "Dark2")
    
    ggMarginal(tmp.p[[beta]], size = 7, colour = NA, groupFill = T, groupColour = T, type = "density")
  })
  }
  
  return(p)
}


gen_PCoA_Int_Num_Plot <- function(data, x.var, x.var.2, y.var, axis, stats){
  
  print("Inside: gen_PCoA_Int_Num_Plot()")  # TEST
  
  centroid <- calc_Centroid(data, x.var)
  # print(paste0(x.var, ":", x.var.2))  # TEST
  
  p <- {lapply(names(y.var), function(beta){
    
    # Temp list for storing ggplots as they are created
    tmp.p <- list()
    
    # Extract statistic
    statistic = round(stats[which(stats[2] == paste0(x.var, ":", x.var.2) &
                                    stats[1] == beta),]$statistic, 3)
    # Extract p.value
    p.val = stats[which(stats[2] == paste0(x.var, ":", x.var.2) &
                          stats[1] == beta),]$p.value
    
    high.stat = stats[which(stats$statistic == max(stats$statistic, na.rm = T)),]$term
    
    if(length(p.val != 0)){
      p.val <- p.val
      statistic <- statistic
    } else {
      p.val <- ">> 0.05"
      statistic <- "..."
    }
    
    tmp.p[[beta]] <- ggplot(data[Dist == beta], aes(x = data[Dist == beta][[2]], y = data[Dist == beta][[3]])) +
      geom_point(aes(color = eval(parse(text = x.var)),
                     fill = eval(parse(text = x.var)),
                     size = .5), 
                 alpha = 0.5) +
      scale_size_area() +
      stat_ellipse(aes(linetype = eval(parse(text = high.stat)) )) +
      geom_label_repel(data = calc_Centroid(data, high.stat)[Dist == beta],
                       aes(label = eval(parse(text = high.stat)) ), 
                       alpha = 0.80) +
      labs(title = paste0(beta, ": Beta diversity ~ ", x.var, ":", x.var.2),
           caption = paste0("Statistic = ", statistic ,", P-value = ", p.val),
           x = axis[Dist == beta]$X.lab,
           y = axis[Dist == beta]$Y.lab) + 
      theme(legend.position="bottom") + guides(size = F,
                                               alpha = F, fill=guide_legend(title=NULL), color=guide_legend(title=NULL), linetype=guide_legend(title=high.stat)) + 
      scale_color_brewer(palette = "Dark2") + 
      scale_fill_brewer(palette = "Dark2")
    
    ggMarginal(tmp.p[[beta]], size = 7, colour = NA, groupFill = T, groupColour = T, type = "density")
  })
  }
  
  return(p)
}

 # Calculate Centroid -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

calc_Centroid <- function(data, x.var, y.var = "Dist", x.var.2 = NA){
  #x = data[Dist == beta][[2]], y = data[Dist == beta][[3]]
  if(!is.na(x.var.2)){
    tmp.centroid <- data[,.(AXIS1 = mean(data[Dist == beta][[2]]),
                            AXIS2 = mean(data[Dist == beta][[3]])),
                         by = c(x.var, x.var.2, y.var)]
  } else {
    tmp.centroid <- data[,.(AXIS1 = mean(data[Dist == beta][[2]]),
                            AXIS2 = mean(data[Dist == beta][[3]])),
                         by = c(x.var, y.var)]
  }
  
  return(tmp.centroid)
}


gen_PCoA_Combo_Plot_1 <- function(data, p.method = "PCoA", physeq, dist = c("bray", "canberra", "sor"), 
                                 facet.var.x = ".", facet.var.y = ".", 
                                 x.var = NULL, x.var.2 = "", title = "", colors = pal.dark2){
  
  tmp.ellipse <- ifelse(x.var.2 == "", x.var, x.var.2)
  
  # Set Seed
  set.seed(42)
  
  # If no x.var.2, temp workaround for an endless loop
  if(is.null(x.var.2)){
    x.var.2 <- x.var
  }
  
  p <- {lapply(dist, function(beta){
    
    # Prints which beta metric is being ran
    print(beta)
    
    tmp.p <- list() # temp list to store plots
    
    # Temp list to store ordination data
    tmp.ord <- ordinate(
      physeq = physeq, 
      method = p.method,  # c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA")
      distance = beta  # c("bray", "canberra", "sor")
    )
    
    # Temp list for storing ggplots as they are created
    tmp.p <- list()
    
    # Plot function
    tmp.p[[beta]] <- plot_ordination(
      physeq = physeq,
      ordination = tmp.ord,
    ) +
      geom_point(aes(fill = eval(parse(text = x.var)),
                     color = eval(parse(text = x.var)),
                     shape = eval(parse(text = x.var.2))),
                 alpha = 0.7, size = 4) +
      stat_ellipse(aes(linetype = eval(parse(text = tmp.ellipse)) )) +
      # geom_label_repel(data = calc_Centroid_2(tmp.ord, x.var = x.var),
      #                  aes(label = eval(parse(text = x.var)) ),
      #                  alpha = 0.80) +
      facet_grid(paste0(facet.var.y, " ~ ", facet.var.x), scales = "free_y") +
      # scale_color_brewer(palette = "Dark2") +
      # scale_fill_brewer(palette = "Dark2") +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      labs(
        # title = ifelse(title != "", title, paste0(c(p.method, beta), collapse = ", "))
        title = ifelse(title == "", paste0("Beta.Score ~ ", ifelse(x.var.2 == "", x.var, paste0(c(x.var, x.var.2), collapse = " + ")) ), title),
        caption = paste0(c(p.method, beta), collapse = ", ")
      ) +
      guides(fill = guide_legend(title=x.var),
             color = guide_legend(title=x.var),#, override.aes = list(size=3)),
             shape = guide_legend(title=x.var.2),
             #override.aes = list(alpha = 1, size = 3) ),
             linetype=guide_legend(title="Ellipses"),
             size = F,
             alpha = F)
    # 
  })}
  
  # return(p)
}



# Check Variable Type -----------------------------------------------------------
#   Description: Checks what type of data a variable is
#   Input: variable
#   Output: data type (categorical, numerical, interaction, etc.)


check_Var_Type <- function(vars){
  
  tmp.var <- get_Variables(vars)  # returns variable name
  tmp.type <- get_Variables(vars, col.1 = 2, col.2 = 1)  # returns category name
  # print(paste0("tmp.var: ", tmp.var, ". tmp.type: ", tmp.type))  # TEST
  
  # Testing
  # if(tmp.type == "numerical"){
  #   # print(paste0("Variable (", tmp.var,") is numerical"))  # TEST
  #   
  # } else if(tmp.type == "categorical") {
  #   # print(paste0("Variable (", tmp.var,") is categorical"))  # TEST
  #   
  # } else if(tmp.type == "interaction"){
  #   # print(paste0("Variable (", tmp.var,") is an interaction"))  # TEST
  #   
  # } else{
  #   # print(paste0("Variable (", var,") is of unknown type"))  # TEST
  #   }
  return(tmp.type)  # returns category type
}

# Extract Interaction Terms -----------------------------------------------------------
#   Description: extracts the interaction terms and their types
#   Input: an interaction term (<var1:var2>)
#   Output: a dataframe containing variable and type columns, with each variable and type of the interaction terms


extract_Int_Terms <- function(int.var, all.vars){
  
  tmp.var <- unlist(strsplit(int.var, split = ":"))  # splits int terms
  var.1 <- tmp.var[1]  # term on the left of the ":"
  var.2 <- tmp.var[2]  # term on the right
    
  int.var.df <- all.vars %>%
    filter(variable %in% c(var.1, var.2)) # grabs the rows that contain these vars
  
  return(int.var.df)  # returns a dataframe 
}




# Print Plots -----------------------------------------------------------
#   Description: 
#   Input: list of plots
#   Output: printed plots

print_Plots <- function(plots, num.cores, print.blank.page = FALSE){
  cl <- makeCluster(num.cores, type = "FORK", outfile = "")
  registerDoParallel(cl, num.cores)
  
  tmp.plot <- foreach(
                plot = plots,
                .verbose = TRUE
              ) %dopar% {
                print(plot, newpage = print.blank.page)
              }
  
  stopCluster(cl)
  
  return (tmp.plot)
}

# Function Name -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 


# Save Plots --------------------------------------------------------------
#   Description: Saves a list of plots as individual
#   Input: list of plots
#   Output: figures as Rdata object

save_plots <- function(list.plots, path, var.name, ID = Sys.Date()){
  save(list.plots, 
       file = ifelse(grepl(".",var.name), paste0(path,"/", gsub("\\.", "_", # replaces periods 9.) with dashes (-)
                                                                var.name), "_", ID, ".RData"),
              paste0(path,"/", var.name, "_", ID, ".RData")) 
       )
       # file = paste0(path,"/", var.name, ".RData"))
}

save_Rdata <- function(list.obj, path, var.name, ID = Sys.Date()){
  save(list.obj, 
       file = ifelse(grepl(".",var.name), paste0(path,"/", gsub("\\.", "_", # replaces periods 9.) with dashes (-)
                                                                var.name), "_", ID, ".RData"),
                     paste0(path,"/", var.name, "_", ID, ".RData")) 
  )
  # file = paste0(path,"/", var.name, ".RData"))
}



# Beta Full dbRDA  --------------------------------------------------------------
#   Description: Distance based redundancy analysis beta diversity metrics
#   Input: variables, distance list, beta methods, method names, phyloseq obj
#   Output: statistical model

full_dbrda_model <- function(vars, distance, methods, names, physeq, data, num.cores = num.cores){

  # num.cores = 4  # Need to add this as a function arg
  cl <- makeCluster(num.cores, type = "FORK", outfile = "")
  registerDoParallel(cl, num.cores)
  
  var <- get_Variables(vars)  # returns a list of variables as strings
  # data <- na_to_zero(data, vars)  # temp swaps NAs with 0s
  
  beta.model <- foreach(
    beta = methods,
    .final = names,
    .verbose = TRUE
  ) %dopar% {
    # progress_Bar_Start(which(methods == beta), length(methods))  # Progress bar
    dist.mat <- distance[[beta]]
    formula <- paste0("dist.mat ~(", paste0(var, collapse = "+") ,")^2")
    print(formula)
    capscale(as.formula(formula),
             data = data,
             comm = otu.matrix(physeq),
             na.action = na.omit # na.omit only non-missing site scores are shown
    )
  }
  
  # beta.model <- {for (beta in methods) {
  #     dist.mat <- distance[[beta]]
  #     # formula <- paste0("dist.mat ~(", paste0(var, collapse = "+") ,")^2")
  #     formula <- "dist.mat ~(Strain+Experiment+Pathology.Results)^2"
  #     print(formula)
  #     capscale(as.formula(formula),
  #              data = data,
  #              comm = otu.matrix(physeq)  
  #     )}}
  
  stopCluster(cl)
  return(beta.model)

}


# full_dbrda_model_2 <- function(vars, distance, methods, names, physeq, data, num.cores = num.cores){
full_dbrda_model_2  <- function(vars, distance, methods, names, physeq, data, terms = 1, num.cores = num.cores){  
  
  # Starts parallel computing
  cl <- makeCluster(num.cores, type = "FORK", outfile = "")
  registerDoParallel(cl, num.cores)
  
  # Build full model
  beta.model.full <- foreach(
    beta = methods,
    .final = names,
    .verbose = TRUE
  ) %dopar% {
    # progress_Bar_Start(which(methods == beta), length(methods))  # Progress bar
    dist.mat <- distance[[beta]]
    form <- if(length(vars) == 1){
      paste0("dist.mat ~ ", vars)
    } else{paste0("dist.mat ~ (", paste0(vars, collapse = "+") ,")^", terms)}
    print(form)  #  TEST
    capscale(as.formula(form),
             data = data,
             na.action = na.omit # na.omit only non-missing site scores are shown
             #comm = otu.matrix(physeq)  # Error might originate here, I think
    )
  }
  
  stopCluster(cl)
  return(beta.model.full)
  
}

# Beta Select dbRDA  --------------------------------------------------------------
#   Description: Selects best dbRDA model from full model
#   Input: variables, distance list, beta methods, method names, phyloseq obj
#   Output: statistical model
#   Note: Error in { : task 1 failed - "otu_table slot is empty."

select_dbrda_model  <- function(vars, distance, methods, names, physeq, data, num.cores = num.cores){
  
  # num.cores = 4  # Need to add this as a function arg 
  var <- get_Variables(vars)  # returns a list of variables as strings
  # data <- na_to_zero(data, vars)  # temp swaps NAs with 0s
  
  # Starts parallel computing
  cl <- makeCluster(num.cores, type = "FORK", outfile = "")
  registerDoParallel(cl, num.cores)

  # Build full model
  beta.model.full <- foreach(
    beta = methods,
    .final = names,
    .verbose = TRUE
  ) %dopar% {
    # progress_Bar_Start(which(methods == beta), length(methods))  # Progress bar
    dist.mat <- distance[[beta]]
    formula <- paste0("dist.mat ~(", paste0(var, collapse = "+") ,")^2")
    print(formula)
    capscale(as.formula(formula),
             data = data,
             na.action = na.omit # na.omit only non-missing site scores are shown
             #comm = otu.matrix(physeq)  # Error might originate here, I think
             )
  }

  # Build select model
  beta.model.select <- foreach(
    beta = methods,
    .final = names,
    .verbose = TRUE
  ) %dopar% {
    progress_Bar_Start(which(methods == beta), length(methods))  # Progress bar
    model <- beta.model.full[[beta]]
    select.mod <- ordistep(model, direction = "both") # Error triggered here, I think
    return(select.mod)
  }

  # Stops parallel computing
  stopCluster(cl)
  return(beta.model.select)
  
  
}



# Select dbRDA 2 ----------------------------------------------------------

select_dbrda_model_2  <- function(vars, distance, methods, names, physeq, data, terms = 1, num.cores = num.cores){
  
  
  
  # Starts parallel computing
  cl <- makeCluster(num.cores, type = "FORK", outfile = "")
  registerDoParallel(cl, num.cores)
  
  # Build full model
  beta.model.full <- foreach(
    beta = methods,
    .final = names,
    .verbose = TRUE
  ) %dopar% {
    # progress_Bar_Start(which(methods == beta), length(methods))  # Progress bar
    dist.mat <- distance[[beta]]
    form <- if(length(vars) == 1){
      paste0("dist.mat ~ ", vars)
    } else{paste0("dist.mat ~ (", paste0(vars, collapse = "+") ,")^", terms)}
    print(form)  #  TEST
    capscale(as.formula(form),
             data = data,
             na.action = na.omit # na.omit only non-missing site scores are shown
             #comm = otu.matrix(physeq)  # Error might originate here, I think
    )
  }
  
  # Build select model
  beta.model.select <- foreach(
    beta = methods,
    .final = names,
    .verbose = TRUE
  ) %dopar% {
    progress_Bar_Start(which(methods == beta), length(methods))  # Progress bar
    model <- beta.model.full[[beta]]
    select.mod <- ordistep(model, direction = "both") # Error triggered here, I think
    return(select.mod)
  }
  
  # Stops parallel computing
  stopCluster(cl)
  return(beta.model.select)
  
  
}



# Beta Anova --------------------------------------------------------------

beta_anova <- function(beta.model, methods, names, num.cores = num.cores){
  
  # Assigns number of cores/threads for parallel computing (aka %dopar%) 
  # num.cores = 4
  cl <- makeCluster(num.cores, type = "FORK", outfile = "")
  registerDoParallel(cl, num.cores)
  
  # Anova
  res <- foreach(
    beta = methods,
    .final = names,
    .verbose = TRUE
  ) %dopar% {
    mod <- beta.model[[beta]]
    anova(mod, by = "term") %>% 
      tidy() %>% 
      mutate(sig = ifelse(p.value <= 0.05, "*", "")) %>%
      mutate(metric = beta, .before = 1) %>%
      arrange(desc(statistic))  # highest to lowest effect size 
  } %>%
    bind_rows()
  
  stopCluster(cl)
  return(res)
}

beta_anova_list <- function(beta.model, methods, names, num.cores = num.cores){
  
  # Assigns number of cores/threads for parallel computing (aka %dopar%) 
  # num.cores = 4
  cl <- makeCluster(num.cores, type = "FORK", outfile = "")
  registerDoParallel(cl, num.cores)
  
  # Anova
  res <- foreach(
    beta = methods,
    .final = names,
    .verbose = TRUE
  ) %dopar% {
    mod <- beta.model[[beta]]
    anova(mod, by = "term") %>% 
      tidy() %>% 
      mutate(sig = ifelse(p.value <= 0.05, "*", "")) %>%
      mutate(metric = beta, .before = 1) %>%
      arrange(desc(statistic))  # highest to lowest effect size 
  }
  
  stopCluster(cl)
  return(res)
}

# anova  -------------------------------------------------------
#   Description: produces an anova stats dataframe
#   Input: linear model
#   Output: anova dataframe

anova_df <- function(model){
  return(Anova(model, type = 2) %>%
    tidy() %>%
    mutate(sig = ifelse(p.value <= 0.05, "*", "")) %>%
    # filter(p.value < 0.15 | is.na(p.value)) %>%
    arrange(desc(statistic))  # highest to lowest effect size
  )
}


# Stats Table -------------------------------------------------------
#   Description: produces a statistical table from anova stats
#   Input: anova results, variables, terms
#   Output: statistical table of anova results 

stats_table <- function(dataframe, vars, terms, methods, formula = NA){
  
  var <- get_Variables(vars)
  
  test <- dataframe %>%
      dplyr::count(.dots = c("metric"),  sort = T)
  
  test <- test[match(methods, test$metric),]

  print(paste0("test: ", test))    # TEST
  
  test.v <- unlist(test$n)
  test.tmp.v <- 0
  hline.num <- c()

  # test.tmp.list <- vector()
  for (i in test.v[1:(length(test.v)-1)]) {
      print(paste0("test.tmp.v: ", test.tmp.v))    # TEST
      test.tmp.v <- test.tmp.v + i
      hline.num <- append(hline.num, test.tmp.v)
  }
  
  if(is.na(hline.num)){
    hline.num <- nrow(dataframe)
  }
  
  print(paste0("hline.num: ", hline.num))  # TEST
  
  # caption <- paste0("beta score ~ (", paste0(var, collapse = "+"), ")^", terms)
  return (dataframe %>%
            flextable() %>%
            set_caption(caption = ifelse(is.na(formula), "", formula)) %>%
            #set_caption(caption = table.caption(caption)) %>%  # Uses autoNumCaptions
            align(j = c(1, 6), align = "left") %>%
            align(j = 2:5, align = "right") %>%
            colformat_double(j = 3, digits = 2) %>%
            colformat_double(j = 5, digits = 3) %>%
            merge_v(j = 1) %>%
            hline(i = hline.num, j = NULL, border = NULL, part = "body") %>%
            set_formatter(values = list("p.value" = p_val_format) )%>%
            autofit()
  )
}

stats_table_2 <- function(dataframe, vars, terms, methods, formula = NA){
  
  var <- get_Variables(vars)
  
  test <- dataframe %>%
    dplyr::count(.dots = c("metric"),  sort = T)
  
  test <- test[match(methods, test$metric),]
  
  print(paste0("test: ", test))    # TEST
  
  test.v <- unlist(test$n)
  test.tmp.v <- 0
  hline.num <- c()
  
  # test.tmp.list <- vector()
  for (i in test.v[1:(length(test.v)-1)]) {
    print(paste0("test.tmp.v: ", test.tmp.v))    # TEST
    test.tmp.v <- test.tmp.v + i
    hline.num <- append(hline.num, test.tmp.v)
  }
  
  if(is.na(hline.num)){
    hline.num <- nrow(dataframe)
  }
  
  print(paste0("hline.num: ", hline.num))  # TEST
  
  # caption <- paste0("beta score ~ (", paste0(var, collapse = "+"), ")^", terms)
  return (dataframe %>%
            flextable() %>%
            set_caption(caption = ifelse(is.na(formula), "", formula)) %>%
            #set_caption(caption = table.caption(caption)) %>%  # Uses autoNumCaptions
            align(j = c(1, 6), align = "left") %>%
            align(j = 2:5, align = "right") %>%
            colformat_double(j = 3, digits = 2) %>%
            colformat_double(j = 5, digits = 3) %>%
            merge_v(j = 1) %>%
            hline(i = hline.num, j = NULL, border = NULL, part = "body") %>%
            set_formatter(values = list("p.value" = p_val_format) ) %>%
            autofit()
  )
}

stats_table_3 <- function(dataframe, terms = NA, hline.num = NA, formula = NA, stat.desc = F){
  
  # Arrange dataframe by most to least significant
  if(!"sig" %in% colnames(dataframe)){
    dataframe<- dataframe %>%
      tidy() %>%
      mutate(sig = ifelse(p.value <= 0.05, "*", "")) %>%
      # arrange(desc(statistic))  # highest to lowest effect size
      arrange(if(!isTRUE(stat.desc)) TRUE else desc(statistic)) 
  }
  
  if(is.na(hline.num)){
    hline.num = seq(1, nrow(dataframe) -1)
  } 

  
  # caption <- paste0("beta score ~ (", paste0(var, collapse = "+"), ")^", terms)
  return (dataframe %>%
            flextable() %>%
            set_caption(caption = ifelse(is.na(formula), "", formula)) %>%
            #set_caption(caption = table.caption(caption)) %>%  # Uses autoNumCaptions
            # align(j = c(1, ncol(dataframe)), align = "left") %>%
            align(j = 2:5, align = "right") %>%
            colformat_double(j = 3, digits = 2) %>%
            colformat_double(j = 5, digits = 3) %>%
            merge_v(j = 1) %>%
            hline(i = hline.num, j = NULL, border = NULL, part = "body") %>%
            set_formatter(values = list("p.value" = p_val_format) ) %>%
            autofit()
  )
}

# Homogeneity of Dispersion  -------------------------------------------------
#   Description: 
#   Input: 
#   Output: Statistical results

beta_hom_disper <- function(data, physeq, vars, methods = "bray", plot = F, factor = F){
  
  tmp.list <- list()
  
  # Temp workaround in the case that data needs to be "factorized"
  if(isTRUE(factor)){
    data$vars <- factor(data[[vars]])
  }
  
  for (m in methods) {
    # print(m)
    tmp.dist <- phyloseq::distance(physeq, method = m)
    # print(tmp.dist)
    tmp.mod <- betadisper(tmp.dist, data[[vars]])
    tmp.anova <- anova(tmp.mod)
    tmp.permTest <- permutest(tmp.mod, pairwise = TRUE, permutations = 999) 
    tmp.mod.HSD <- TukeyHSD(tmp.mod)
    
    cat("\n")
    print(paste0("Method: ", m))
    # print(tmp.anova)
    print(tmp.permTest) 
    # tmp.permTest %>% pander()
    # tmp.list[m] <- c(tmp.permTest)
    
    if(isTRUE(plot)){
      # plot(tmp.mod)
      boxplot(tmp.mod)
      # plot(tmp.mod.HSD)
    }
    
  }
  # tmp.dist <- phyloseq::distance(physeq, method = methods)
  # tmp.mod <- betadisper(tmp.dist, data[[vars]])
  # tmp.permTest <- permutest(tmp.mod, pairwise = TRUE, permutations = 99)
  
  # return(tmp.mod)
}

# Unconstrained Ordination Plot -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

betaUncon_plot <- function(data, p.method = "PCoA", physeq, dist = c("bray", "canberra", "sor"), 
                           facet.var.x = ".", facet.var.y = ".", tmp.ellipse = NULL,
                           x.var = NULL, x.var.2 = "", title = "", colors = pal.dark2){
  
  # tmp.ellipse <- ifelse(x.var.2 == "", x.var, x.var.2)
  # if(!is.null(tmp.ellipse)){
  #   tmp.ellipse <- tmp.ellipse
  # } else{
  #   tmp.ellipse <- ifelse(x.var.2 == "", x.var, x.var.2)
  # }
  tmp.ellipse <- x.var
  
  # Set Seed
  set.seed(42)
  
  # If no x.var.2, temp workaround for an endless loop
  if(is.null(x.var.2)){
    x.var.2 <- x.var
  }
  
  p <- {lapply(dist, function(beta){
    
    # Prints which beta metric is being ran
    print(beta)

    tmp.p <- list() # temp list to store plots
    
    # Temp list to store ordination data
    tmp.ord <- ordinate(
      physeq = physeq, 
      method = p.method,  # c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA")
      distance = beta  # c("bray", "canberra", "sor")
    )
    
    # Temp list for storing ggplots as they are created
    tmp.p <- list()
    
    # Plot function
    tmp.p[[beta]] <- plot_ordination(
      physeq = physeq,
      ordination = tmp.ord,
      # color = x.var,12333
      justDF = T
      
    ) #+
    
    ggplot(data = tmp.p[[beta]], aes(x=Axis.1, y=Axis.2)) +
      geom_point(aes(fill = eval(parse(text = x.var)),
                     color = eval(parse(text = x.var)),
                     shape = eval(parse(text = x.var.2))),
                 alpha = 0.9, size = 2) +
      stat_ellipse(aes(color = eval(parse(text = tmp.ellipse)) )) +
      # geom_label_repel(data = calc_Centroid_2(tmp.ord, x.var = x.var),
      #                  aes(label = eval(parse(text = x.var)) ),
      #                  alpha = 0.80) +
      facet_grid(paste0(facet.var.y, " ~ ", facet.var.x), scales = "free_y") +
      # scale_color_brewer(palette = "Dark2") +
      # scale_fill_brewer(palette = "Dark2") +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      labs(
        # title = ifelse(title != "", title, paste0(c(p.method, beta), collapse = ", "))
        title = ifelse(title == "", paste0("Beta Score ~ ", ifelse(x.var.2 == "", sub_char_var(x.var, ".", " "), paste0(c(sub_char_var(x.var, ".", " "), sub_char_var(x.var.2, ".", " ")), collapse = " + ")) ), sub_char_var(title, ".", " ")),
        caption = paste0(c(p.method, beta), collapse = ", "),
        x = "Axis 1", # x = axis[Dist == beta]$X.lab,
        y = "Axis 2" # y = axis[Dist == beta]$Y.lab
      ) +
      theme(legend.position = "bottom") +
      guides(fill = guide_legend(title=x.var),
             color = guide_legend(title=x.var),#, override.aes = list(size=3)),
             shape = guide_legend(title=x.var.2),
                                  #override.aes = list(alpha = 1, size = 3) ),
             # linetype=guide_legend(title="Ellipses"),
             size = F,
             alpha = F)
    # 
  })}
  
  # return(p)
}


# sub_char_var(title, ".", " ")
# sub_char_var(x.var, ".", " ")
# sub_char_var(y.var, ".", " ")

betaUncon_plot_num <- function(data, p.method = "PCoA", physeq, dist = c("bray", "canberra", "sor"), 
                           facet.var.x = ".", facet.var.y = ".", tmp.ellipse = NULL,
                           x.var = NULL, x.var.2 = "", title = "", colors = "RdYlBu"){
  
  # tmp.ellipse <- ifelse(x.var.2 == "", x.var, x.var.2)
  # if(!is.null(tmp.ellipse)){
  #   tmp.ellipse <- tmp.ellipse
  # } else{
  #   tmp.ellipse <- ifelse(x.var.2 == "", x.var, x.var.2)
  # }
  tmp.ellipse <- x.var
  
  # Set Seed
  set.seed(42)
  
  # If no x.var.2, temp workaround for an endless loop
  if(is.null(x.var.2)){
    x.var.2 <- x.var
  }
  
  p <- {lapply(dist, function(beta){
    
    # Prints which beta metric is being ran
    print(beta)
    
    tmp.p <- list() # temp list to store plots
    
    # Temp list to store ordination data
    tmp.ord <- ordinate(
      physeq = physeq, 
      method = p.method,  # c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA")
      distance = beta  # c("bray", "canberra", "sor")
    )
    
    # Temp list for storing ggplots as they are created
    tmp.p <- list()
    
    # Plot function
    tmp.p[[beta]] <- plot_ordination(
      physeq = physeq,
      ordination = tmp.ord,
      # color = x.var,
      justDF = T
      
    ) #+
    
    ggplot(data = tmp.p[[beta]], aes(x=Axis.1, y=Axis.2)) +
      # geom_point(aes(fill = eval(parse(text = x.var)),
      #                color = eval(parse(text = x.var)),
      #                shape = eval(parse(text = x.var.2))),
      #            alpha = 0.7, size = 4) +
      geom_point(aes(color = eval(parse(text = x.var.2)),
                     # fill = eval(parse(text = x.var)),
                     shape = eval(parse(text = x.var)),
                     size = eval(parse(text = x.var.2)) ),
                 alpha = .9, size = 2) +
      scale_size_area() +
      stat_ellipse(aes(linetype = eval(parse(text = tmp.ellipse)) )) +
      # geom_label_repel(data = calc_Centroid_2(tmp.ord, x.var = x.var),
      #                  aes(label = eval(parse(text = x.var)) ),
      #                  alpha = 0.80) +
      facet_grid(paste0(facet.var.y, " ~ ", facet.var.x), scales = "free_y") +
      # scale_color_brewer(palette = "Dark2") +
      # scale_fill_brewer(palette = "Dark2") +
      scale_colour_gradient(
        low = pal.RdPu[3],
        high = pal.paired[10],
        space = "Lab",
        na.value = "grey50",
        guide = "colourbar",
        aesthetics = "colour"
      ) +
      labs(
        # title = ifelse(title != "", title, paste0(c(p.method, beta), collapse = ", "))
        title = ifelse(title == "", paste0("Beta Score ~ ", ifelse(x.var.2 == "", sub_char_var(x.var, ".", " "), paste0(c(sub_char_var(x.var, ".", " "), sub_char_var(x.var.2, ".", " ")), collapse = " + ")) ), sub_char_var(title, ".", " ")),
        caption = paste0(c(p.method, beta), collapse = ", ")
      ) + 
      theme(legend.position = "bottom") +
      guides(fill = guide_legend(title=x.var.2),
             color = guide_legend(title=x.var.2, override.aes = list(size=3)),
             shape = guide_legend(title=x.var, override.aes = list(size=3)),
             #override.aes = list(alpha = 1, size = 3) ),
             linetype=guide_legend(title="Ellipses"),
             size = F,
             alpha = F)
    # 
  })}
  
  # return(p)
}



calc_Centroid_2 <- function(data, x.var, x.var.2 = NA){
  #x = data[Dist == beta][[2]], y = data[Dist == beta][[3]]
  if(!is.na(x.var.2)){
    tmp.centroid <- data[,.(AXIS1 = mean(data[[2]]),
                            AXIS2 = mean(data[[3]])),
                         by = c(x.var, x.var.2)]
  } else {
    tmp.centroid <- data[,.(AXIS1 = mean(data[[2]]),
                            AXIS2 = mean(data[[3]])),
                         by = c(x.var)]
  }
  
  return(tmp.centroid)
}

calc_Centroid_CAP <- function(data, x.var, x.var.2 = NA){
  #x = data[Dist == beta][[2]], y = data[Dist == beta][[3]]
  if(!is.na(x.var.2)){
    tmp.centroid <- data[,.(CAP1 = mean(data[[2]]),
                            CAP2 = mean(data[[3]])),
                         by = c(x.var, x.var.2)]
  } else {
    tmp.centroid <- data[,.(CAP1 = mean(data[[2]]),
                            CAP2 = mean(data[[3]])),
                         by = c(x.var)]
  }
  
  return(tmp.centroid)
}


# Check Project Structure -------------------------------------------------
#   Description: Checks to make sure the important folders are present
#   Input: Project (working) directory
#   Output: True/False if project structure is good

check_proj_struct <- function(dir.path, sub.dir.path){
  # Checks if directory exists, creates folder if not
  ifelse(!dir.exists(file.path(dir.path, sub.dir.path)), 
         dir.create(file.path(dir.path, sub.dir.path)), FALSE)
  return(paste0(dir.path, sub.dir.path))  # returns path name
}



# Save Environment File ---------------------------------------------------
#   Description: Finds the latest R environment file
#   Input: path name to Robjects directory
#   Output: saves environment file with current date

save_env <- function(obj.path, ID = Sys.Date(), extra_info = NA){
    if(is.na(extra_info)){
        save.image(paste0(obj.path, "/environment", ID, "_ENV.RData"))
    } else{
  save.image(paste0(obj.path, "/environment_", extra_info, ID, "_ENV.RData"))
    }
}


# Find Latest Environment File --------------------------------------------
#   Description: Finds the latest R environment file
#   Input: directory
#   Output: Loads newest environment file or continues script


load_env <- function(obj.path) {
  
  out <- tryCatch(
    {
      files.robjects.dir <- file.info(list.files(obj.path, full.names = T))
      files.env <- files.robjects.dir %>%
        filter(grepl("ENV.RData", rownames(files.robjects.dir)))
      latest.ENV <- rownames(files.env)[which.max(files.env$mtime)]
      load(latest.ENV)
      print(paste0(latest.ENV, " was successfully loaded"))
    },
    
    error = function(e){          # Specifying error message
      message("No environment file found.
               Check that environment file is correctly located and named")
    },
    
    warning = function(w){        # Specifying warning message
      message("There was a warning message.")
    },
    
    finally = {                   # Specifying final message
      message("Proceed with script.")
    }
  )
}



#Ordination data table list ---------------------------------------------
#   Description: 
#   Input: 
#   Output: 

ord_dt_list <- function(model, physeq){
  return(
    lapply(model, function(model) {
      get.biplot.data(ps = physeq, ord = model)})
  )
  
}


# Sample Coordinate Datatable ---------------------------------------------
#   Description: 
#   Input: 
#   Output: 

samp_coord_dt <- function(ord.datatable.list, names){
  
  print(paste0("names: ", names))  # TEST
  
  sample.coord.dt <- lapply(names(ord.datatable.list), function(beta) {  # CHANGE
      ord.datatable.list[[beta]]$sample.coords[, Dist := beta]
      return(ord.datatable.list[[beta]]$sample.coords)
    }) %>% rbindlist()
    
  sample.coord.dt[
    , Dist := factor(Dist, levels = names(names)[c(1:length(names))]) ]
  
  return(sample.coord.dt)
}


# Axis Labels -------------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

axis_labels <- function(ord.datatable.list){
  return(axis.labels <- lapply(names(ord.datatable.list), function(beta) {
    data.table(
      Dist = beta,
      X.lab = ord.datatable.list[[beta]]$axes.labs[1],
      Y.lab = ord.datatable.list[[beta]]$axes.labs[2]
    )
  })  %>% rbindlist())
}







# Generate Multiple PCoA Plots (Beta) -----------------------------------------
#   Description: Generates alpha diversity plots of significant variables 
#   Input: data, significant variables, alpha score, facet conditions
#   Output: multiple plots displaying alpha diversity
#
#   To-Do: Handle interaction terms (eg, var1:var2, var2:var3, var1:var3)

gen_multi_beta_pcoa_plots <- function(data, axis, vars, methods, stats){
  
  # print(vars)
  var <- get_Variables(vars)
  tmp.stats <- stats
  # print(var)  # TEST
  # print(stats)  # TEST
  
  
  # Loop through sig vars -> loop beta metrics
  tmp.list <- lapply(var, function(v){ 
    
    # progress_Bar_Start(which(var == v), length(var))  # Progress bar
    
    lapply(names(methods), function(beta){
      
      progress_Bar_Start(which(methods == beta), length(methods))  # Progress bar
      
      statistic = round(tmp.stats[which(tmp.stats[2] == v &
                                tmp.stats["metric"] == methods[beta]),]$statistic, 3)
      p.val = tmp.stats[which(tmp.stats[2] == v &
                                tmp.stats["metric"] == methods[beta]),]$p.value
      
      
      if(length(p.val != 0)){
        p.val <- p.val
        statistic <- statistic
      } else {
        p.val <- ">> 0.05"
        statistic <- "..."
      }
      
      # Check for interaction terms
      if(grepl(":", v)){
        
        # Interaction Terms Plotting
        tmp.v <- unlist(strsplit(v, split = ":"))
        var.1 <- tmp.v[1]
        var.2 <- tmp.v[2]
        v.shape <- var.2
        v.label <- var.2
        
        # CENTROIDS
        ifelse(grepl(":", v),
               # If v == list, aka has interaction terms:
               centroid <- data[,.(CAP1 = mean(CAP1),
                                   CAP2 = mean(CAP2)),
                                by = c(var.1, var.2, "Dist")],
               # If v != list, aka single variable:
               centroid <- data[,.(CAP1 = mean(CAP1),
                                   CAP2 = mean(CAP2)),
                                by = c(var.1, "Dist")])
        
        # PLOT
        ggplot(data[Dist == methods[beta]], aes(x = CAP1, y = CAP2)) +
          geom_point(aes(fill = eval(parse(text = var.2)),
                         color = eval(parse(text = var.2)) ),  # v.shape set in if/else above
                     size = 2) +
          stat_ellipse(aes(color = eval(parse(text = var.1)) )) +
          geom_label_repel(data = centroid[Dist == methods[beta]],
                           aes(label = paste0(eval(parse(text = var.1)), ":",
                                              eval(parse(text = v.label)) ), 
                               size = 3)) +
          facet_wrap(~ Dist, scales = "free", ncol = 3) +
          labs(title = ifelse(grepl(":", v),
                              paste0(methods[beta], ": Beta diversity ~ ", var.1, ":", var.2),  # prints interaction terms if v == list
                              paste0(methods[beta], ": Beta diversity ~ ", var.1) ),
               caption = paste0("Statistic = ", statistic ,", P-value = ", p.val),
               x = axis[Dist == methods[beta]]$X.lab,
               y = axis[Dist == methods[beta]]$Y.lab)
        
      } else{
        
        # Plot for single variables
        var.1 <- v
        v.shape <- var.1  # Set shape to default (e.g., circle)
        v.label <- var.1
        
        # CENTROIDS
        centroid <- data[,.(CAP1 = mean(CAP1),
                            CAP2 = mean(CAP2)),
                         by = c(v, "Dist")]
        
        # PLOT
        ggplot(data[Dist == beta.methods[beta]], aes(x = CAP1, y = CAP2)) +
          geom_point(aes(fill = eval(parse(text = v)) ,
                         color = eval(parse(text = v)) ),
                     size = 2) +
          stat_ellipse(aes(color = eval(parse(text = v)) )) +
          geom_label(data = centroid[Dist == methods[beta]],
                     aes(label = eval(parse(text = v))),
                     size = 3) +
          facet_wrap(~ Dist, scales = "free", ncol = 3) +
          labs(title = paste0(methods[beta], ": Beta diversity ~ ", v),
               caption = paste0("Statistic = ", statistic ,", P-value = ", p.val),
               x = axis[Dist == methods[beta]]$X.lab,
               y = axis[Dist == methods[beta]]$Y.lab)
      }
      
    })
  })
  
  return(tmp.list)
}
  
  
#   function(data, axis, vars, methods){
#   
#   var <- get_Variables(vars)
#   
#   # Loop through sig vars -> loop beta metrics
#   tmp.list <- lapply(var, function(v){ 
#     
#     lapply(names(methods), function(beta){
#       
#       # Check for interaction terms
#       
#       if(grepl(":", v)){
#         tmp.v <- unlist(strsplit(v, split = ":"))
#         var.1 <- tmp.v[1]
#         var.2 <- tmp.v[2]
#         v.shape <- var.2
#         v.label <- var.2
#         
#       } else{
#         var.1 <- v
#         v.shape <- var.1  # Set shape to default (e.g., circle)
#         v.label <- var.1
#       }
#       
#       
#       # centroid labels
#       ifelse(grepl(":", v),
#              # If v == list, aka has interaction terms:
#              centroid <- data[,.(CAP1 = mean(CAP1),
#                                  CAP2 = mean(CAP2)),
#                               by = c(var.1, var.2, "Dist")],
#              # If v != list, aka single variable:
#              centroid <- data[,.(CAP1 = mean(CAP1),
#                                  CAP2 = mean(CAP2)),
#                               by = c(var.1, "Dist")])
#       
#       print(centroid)
#       
#       # Plot
#       ggplot(data[Dist == methods[beta]], aes(x = CAP1, y = CAP2)) +
#         geom_point(aes(fill = eval(parse(text = var.1)),
#                        color = eval(parse(text = var.1)),
#                        shape = eval(parse(text = v.shape)) ),  # v.shape set in if/else above
#                    size = 2) +
#         stat_ellipse(aes(color = eval(parse(text = var.1)) )) +
#         geom_label(data = centroid[Dist == methods[beta]],
#                    aes(label = eval(parse(text = v.label)), size = 3)) +
#         facet_wrap(~ Dist, scales = "free", ncol = 3) +
#         labs(title = ifelse(grepl(":", v),
#                             paste0(methods[beta], ": Beta diversity ~ ", var.1, ":", var.2),  # prints interaction terms if v == list
#                             paste0(methods[beta], ": Beta diversity ~ ", var.1) ),
#              # caption = "",
#              x = axis[Dist == methods[beta]]$X.lab,
#              y = axis[Dist == methods[beta]]$Y.lab)
#     })
#   })
#   
#   return(tmp.list)
# }


# gen_multi_beta_pcoa_plots <- function(data, axis, vars, methods){
#   
#   var <- get_Variables(vars)
#   
#   # Loop through sig vars -> loop beta metrics
#   tmp.list <- lapply(var, function(v){ 
#     lapply(names(methods), function(beta){
#       
#       # centroid labels
#       centroid <- data[,.(CAP1 = mean(CAP1),
#                           CAP2 = mean(CAP2)),
#                           by = c(v, "Dist")]
#       
#       # Plot
#       ggplot(data[Dist == beta.methods[beta]], aes(x = CAP1, y = CAP2)) + 
#         geom_point(aes(fill = eval(parse(text = v)) , 
#                        color = eval(parse(text = v)) ), 
#                        size = 2) + 
#         stat_ellipse(aes(color = eval(parse(text = v)) )) +
#         geom_label(data = centroid[Dist == methods[beta]], 
#                    aes(label = eval(parse(text = v))), 
#                    size = 3) +
#         facet_wrap(~ Dist, scales = "free", ncol = 3) + 
#         labs(title = paste0(methods[beta], ": Beta diversity ~ ", v),
#              # caption = "",
#              x = axis[Dist == methods[beta]]$X.lab, 
#              y = axis[Dist == methods[beta]]$Y.lab)
#     })
#   })
#   
#   return(tmp.list)
# }

# Generate Multiple PCoA Plots 2.0 (Beta) -----------------------------------------
#   Description: Same as other func, but checks for interaction terms
#   Input: 
#   Output: 
# 
# gen_multi_beta_pcoa_plots_2 <- function(data, axis, vars, methods, stats){
#   
#   var <- get_Variables(vars)
#   
#   
#   # Loop through sig vars -> loop beta metrics
#   tmp.list <- list()
#   
#   for (v in var) {
#     
#     
#   for (beta in methods) {
#     
#       p.val = stats[which(stats[2] == v &  # checking that v == term
#                             stats["metric"] == methods[beta]),]$p.value
#       
#       
#       print(paste0("P.val: ", p.val))
#       if(length(p.val) == 0){
#         
#         try(next
#         )
#         
#       } else{
#       
#       # Check for interaction terms
#       if(grepl(":", v)){
#         
#         # Interaction Terms Plotting
#         tmp.v <- unlist(strsplit(v, split = ":"))
#         var.1 <- tmp.v[1]
#         var.2 <- tmp.v[2]
#         v.shape <- var.2
#         v.label <- var.2
#         
#         # CENTROIDS
#         ifelse(grepl(":", v),
#                # If v == list, aka has interaction terms:
#                centroid <- data[,.(CAP1 = mean(CAP1),
#                                    CAP2 = mean(CAP2)),
#                                 by = c(var.1, var.2, "Dist")],
#                # If v != list, aka single variable:
#                centroid <- data[,.(CAP1 = mean(CAP1),
#                                    CAP2 = mean(CAP2)),
#                                 by = c(var.1, "Dist")])
#         
#         # PLOT
#         ggplot(data[Dist == methods[beta]], aes(x = CAP1, y = CAP2)) +
#           geom_point(aes(fill = eval(parse(text = var.2)),
#                          color = eval(parse(text = var.2)) ),  # v.shape set in if/else above
#                      size = 2) +
#           stat_ellipse(aes(color = eval(parse(text = var.1)) )) +
#           geom_label_repel(data = centroid[Dist == methods[beta]],
#                      aes(label = paste0(eval(parse(text = var.1)), ":",
#                                         eval(parse(text = v.label)) ), 
#                          size = 3)) +
#           facet_wrap(~ Dist, scales = "free", ncol = 3) +
#           labs(title = ifelse(grepl(":", v),
#                               paste0(methods[beta], ": Beta diversity ~ ", var.1, ":", var.2),  # prints interaction terms if v == list
#                               paste0(methods[beta], ": Beta diversity ~ ", var.1) ),
#                caption = paste0("P-value = ", p.val),
#                x = axis[Dist == methods[beta]]$X.lab,
#                y = axis[Dist == methods[beta]]$Y.lab)
#         
#       } else{
#         
#         # Plot for single variables
#         var.1 <- v
#         v.shape <- var.1  # Set shape to default (e.g., circle)
#         v.label <- var.1
#         
#         # CENTROIDS
#         centroid <- data[,.(CAP1 = mean(CAP1),
#                             CAP2 = mean(CAP2)),
#                             by = c(v, "Dist")]
# 
#         # PLOT
#         ggplot(data[Dist == beta.methods[beta]], aes(x = CAP1, y = CAP2)) +
#           geom_point(aes(fill = eval(parse(text = v)) ,
#                          color = eval(parse(text = v)) ),
#                          size = 2) +
#           stat_ellipse(aes(color = eval(parse(text = v)) )) +
#           geom_label(data = centroid[Dist == methods[beta]],
#                      aes(label = eval(parse(text = v))),
#                      size = 3) +
#           facet_wrap(~ Dist, scales = "free", ncol = 3) +
#           labs(title = paste0(methods[beta], ": Beta diversity ~ ", v),
#                caption = paste0("P-value = ", p.val),
#                x = axis[Dist == methods[beta]]$X.lab,
#                y = axis[Dist == methods[beta]]$Y.lab)
#       }
#       
#       }
#       tmp.list[v][] <-
#       }
#     
#   }
# 
#   return(tmp.list)
# }


# NA to Zero -----------------------------------------------------------
#   Description: Converts NAs to zeros
#   Input: data-table/frame and variables
#   Output: data-table/frame with select variabes containing NAs instead of zeros

na_to_zero <- function(data, vars){

  data <- as.data.frame(data)
  var <- get_Variables(vars)
  data[,var][is.na(data[,var])] <- 0
  return (data)
}


# Zero to NA -----------------------------------------------------------
#   Description: Converts zeros to NAs
#   Input: data-table/frame and variables
#   Output: data-table/frame with select variabes containing zeros instead of NA's

zero_to_na <- function(data, vars){
  data <- as.data.frame(data)
  var <- get_Variables(vars)
  data[,var] <- na_if(data[,var], 0)
  return(data)
}


# Differential Abundance -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 
#   To do: 

# diff_Abund_Maaslin2 <- function(f.eff, r.eff, ref, min.ab, min.prv, 
#                                 norm, trns, an.meth, corr, n.cores){
diff_Abund_Maaslin2 <- function(physeq, param, physeq.ID){
  
  print(param)
  
  maaslin.output <- paste0(output.path, "/MaAsLin2",
                           "_PSobj-", paste0(param$physeq.ID, collapse = "-"),
                           "_FixEf-", paste0(param$f.eff, collapse = "-"),
                           "_RanEf-", paste0(gsub("\\,", "-", param$r.eff), collapse = "-"),
                           "_Ref-", paste0(gsub("\\,", "", param$ref), collapse = "-"),
                           "_MinAb-", paste0(param$min.ab, collapse = "-"),
                           "_MinP-", paste0(param$min.prv, collapse = "-"),
                           "_N-", paste0(param$norm, collapse = "-"),
                           "_M-", paste0(param$an.meth, collapse = "-"),
                           "_C-", paste0(param$corr, collapse = "-"),
                           "_ID-", paste0(param$ID, collapse = "-"), 
                           "_", paste0(param$Add, collapse = "-")
                           )
  
  print(maaslin.output)
  
  fit_data = Maaslin2(input_data = otu.matrix(physeq),
                      input_metadata = sample.data.frame(physeq),
                      output = maaslin.output,
                      fixed_effects = param$f.eff,  # Explanatory Variables
                      random_effects = param$r.eff,  # Since same tanks were sampled multiple times
                      reference = param$ref,
                      min_abundance=param$min.ab,  # Conservative threshold
                      min_prevalence=param$min.prv,  # Default
                      normalization = param$norm,  # Data are rarefied/norm already
                      transform = param$trns,  # None, since using NegBin model
                      analysis_method = param$an.meth,
                      correction = param$corr,  # Multiple test correction
                      cores = param$n.cores)
  return(fit_data)
}



# Progress Bar -----------------------------------------------------------
#   Description: Displays a progress bar as your function proceeds
#   Input: 
#   Output: prints a progress bar to screen and "beeps" when finished if set T
#   Examples:
#     - progress_Bar_Start(which(<name_of_list> == <iterator>), length(<name_of_list>))  # for inside lapply loops
#     - progress_Bar_Start(iteration = i, n_iter = length(<name_of_list>)) 
#     - progress_Bar_Start(iteration = i, n_iter = <an_integar>) 

progress_Bar_Start <- function(iteration = 1, n_iter = 10, beep.sound = F){
  
  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = n_iter,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  setTxtProgressBar(pb, iteration)
  
  if(iteration == n_iter & beep.sound == TRUE){
    beep(10) # Random notification
    close(pb)
  } 
}


# Remove Outliers -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

check_Outliers_ps <- function(physeq, var, fill.var = NULL){
  # print(sample.data.table(physeq)[[var]])
  
  # Phyloseq Object, with BC containing no NA's
  ps.noNA <- subset_samples(physeq, !is.na(var))
  
  print(summary(sample.data.table(ps.noNA)[[var]]))
  
  lowerq = quantile(sample.data.table(ps.noNA)[[var]], na.rm = T)[2]
  upperq = quantile(sample.data.table(ps.noNA)[[var]], na.rm = T)[4]
  iqr = upperq - lowerq #Or use IQR(data)
  
  # Mild Outliers
  mild.threshold.upper = signif((iqr * 1.5) + upperq, digits = 4)
  mild.threshold.lower = signif(lowerq - (iqr * 1.5), digits = 4)
  
  print(paste0("Mild outliers. Lower: ", mild.threshold.lower, 
               ". Upper: ", mild.threshold.upper)) 

  # Extreme Outliers
  intermediate.threshold.upper = signif((iqr * 2) + upperq, digits = 4)
  intermediate.threshold.lower = signif(lowerq - (iqr * 2), digits = 4)
  
  print(paste0("Intermediate outliers. Lower: ", intermediate.threshold.lower, 
               ". Upper: ", intermediate.threshold.upper)) 
  
  # Extreme Outliers
  extreme.threshold.upper = signif((iqr * 3) + upperq, digits = 4)
  extreme.threshold.lower = signif(lowerq - (iqr * 3), digits = 4)
  
  print(paste0("Extreme outliers. Lower: ", extreme.threshold.lower, 
               ". Upper: ", extreme.threshold.upper)) 
  
  # Looks like there are outliers
  
  tmp.sam.data <- sample.data.frame(ps.noNA)
  
  outlier.plot <- {
    ggplot(tmp.sam.data, aes(x = eval(parse(text = var)), fill = eval(parse(text = fill.var)) )) + 
      geom_histogram() + 
      geom_vline(aes(xintercept = extreme.threshold.upper), color="red") + 
      geom_vline(aes(xintercept = extreme.threshold.lower), color="red") + 
      geom_vline(aes(xintercept = intermediate.threshold.upper), color="orange") + 
      geom_vline(aes(xintercept = intermediate.threshold.lower), color="orange") +
      geom_vline(aes(xintercept = mild.threshold.upper), color="green") + 
      geom_vline(aes(xintercept = mild.threshold.lower), color="green") + 
      labs(title = paste0(var, " outliers by ", fill.var),
           x = var,
           caption = str_wrap(paste0("Extreme outliers (3x IQR) below ", extreme.threshold.lower, 
                            " and above ", extreme.threshold.upper, " (red lines). ",
                            "Intermediate outliers (2x IQR) below ", intermediate.threshold.lower, 
                            " and above ", intermediate.threshold.upper, " (orange lines). ",
                            "Mild outliers (1.5x IQR) below ", mild.threshold.lower, 
                            " and above ", mild.threshold.upper, " (green lines).")), 80) + 
      theme(legend.position = "bottom") + 
      guides(size = F, alpha = F, fill=guide_legend(title=fill.var), color=guide_legend(title=fill.var)) + 
      scale_fill_brewer(palette = "Dark2") +
      scale_color_brewer(palette = "Dark2")
    }
    

  return(list(mild = c(mild.lower = mild.threshold.lower,
                       mild.upper = mild.threshold.upper),
              intermediate = c(intermediate.lower = intermediate.threshold.lower,
                               intermediate.upper = intermediate.threshold.upper),
              extreme = c(extreme.lower = extreme.threshold.lower,
                          extreme.upper = extreme.threshold.upper),
              plot = outlier.plot))
}


# KEATON symlink.fastqs -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

#' Symlink Fastqs
#'
#' A function to create symlinks with easier to parse names that the typical names that files have coming off the sequencing machine
#' @param seq.dir The path to the directory in which the raw sequence files are contained.
#' @param ids.tbl A data.frame or data.table containing sample names in one column and their corresponding file IDs in another (e.g. barcodes or sample IDs as used in the raw sequence files). Sample names and file IDs *may* be identical in some instances, please still provide a two column table.
#' @param smpl.id.col The name/number of the column in the `ids.tbl` that contains the name of the samples, which will be used in naming the symlinks. Default is 'Sample'.
#' @param file.id.col The name/number of the column in the `ids.tbl` that contains the corresponding file IDs for matching raw sequence files. Default is 2.
#' @param split.pattern A character string containg the pattern you want to use to separate the sample name and R1/R2 in the symlink names. Default is '--'.
#' @param quiet Logical, if TRUE, there will be no printing of progress. Default is FALSE
#' @export

symlink_fastqs_local <- function(
    seq.dir,
    ids.tbl,
    smpl.id.col = "Sample",
    file.id.col = 2,
    split.pattern = "--",
    quiet = FALSE
) {
    if ("data.table" %in% class(ids.tbl)) {
        ids.dt <- ids.tbl
    } else {
        ids.dt <- as.data.table(ids.dt)
    }
    setkeyv(ids.dt, smpl.id.col)
    for (smpl in ids.dt[[smpl.id.col]]) {
        file.id <- ids.dt[smpl][[file.id.col]]
        if (!quiet) {
            cat(paste0(file.id, " -> ", smpl, " ... "), sep = "")
        }
        seq.files <- list.files(
            path = seq.dir,
            pattern = file.id,
            full.names = T
        )
        if (length(seq.files) == 0) {
            stop("no files matched file ID")
        } else {
            read1.file <- seq.files[str_detect(seq.files, "[-_\\.]R1[-_\\.]")]
            read2.file <- seq.files[str_detect(seq.files, "[-_\\.]R2[-_\\.]")]
            lnName.r1 <- paste(smpl, "R1.fastq.gz", sep = split.pattern)
            lnName.r2 <- paste(smpl, "R2.fastq.gz", sep = split.pattern)
            cmd1 <- paste0("ln -s '", read1.file, "' ", lnName.r1)
            cmd2 <- paste0("ln -s '", read2.file,  "' ", lnName.r2)
            system(cmd1)
            system(cmd2)
            check.cmd1 <- paste("zcat", lnName.r1, "2>/dev/null | head -n 1")
            check1 <- length(system(check.cmd1, intern = TRUE)) > 0
            check.cmd2 <- paste("zcat", lnName.r2, "2>/dev/null | head -n 1")
            check2 <- length(system(check.cmd2, intern = TRUE)) > 0
            if (!quiet & check1 & check2) {
                cat("good", sep = "\n")
            } else if (!check1 & check2) {
                cat(paste0(lnName.r1, " is empty (links to", read1.file, ")"), sep = "\n")
            } else if (check1 & !check2) {
                cat(paste0(lnName.r2, " is empty (links to", read2.file, ")"), sep = "\n")
            } else {
                cat(
                    paste0(
                        lnName.r1, " and ", lnName.r2,
                        " are empty (link to", read1.file, " and ", read2.file, ", respectively)"
                    ),
                    sep = "\n"
                )
            }
        }
    }
}

# KEATON -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

#' Finish the dada2 pipeline
#'
#' This function completes the dada2 pipeline once the user has viewed the plots generated by \code{\link{dada2.upto.qualPlots}} and decided on the forward (and reverse) truncation lengths.
#'@param fastq.path Path to raw FASTQ files. This should be the same path supplied to `dada2.upto.qualPlots`.
#'@param truncFwd.len Truncate forward reads after this many bases. Reads shorter than this are discarded.
#'@param truncRev.len Truncate reverse reads after this many bases. Reads shorter than this are discarded.
#'@param taxa.db Path to the database file to be used to assign taxonomy to the ASVs.
#'@param metadata.file Path to file containing metadata for samples used in creating the output phyloseq object. Must be formatted as csv, tsv, or rds.
#'@param maxCores The maximum number of cores to utilize in steps that are parallelizabe. Default is 4.
#'@param build.tree Logical. Should a phylogenetic tree be created for the ASVs called by `dada2`? Defalt is TRUE.
#'#'@param guide.seqs.file (Optional) Path to file containing guide sequences for improved tree building. Only needed if above parameter `build.tree` is TRUE AND you *do not* want to use the provided files. Default is NULL.
#'@param alignment.template.file (Optional) Path to file containing aligment template used for NAST alignment of sequences prior to tree building. Only needed if above parameter `build.tree` is TRUE AND you *do not* want to use the provided files. Default is NULL.
#'@param fasttree.path (Optional) Path to FastTree executable. Only needed if above parameter `build.tree` is TRUE.
#'@param user.output.path Path to a directory if you want to use a specific output directory instead of the one programmatically generated by the pipeline. Primarily used if you are re-running a pipeline and want to save to previously generated output path. Default is NULL (which will autogenerate the output path).
#'@param paired Logical. Are these paired sequence data? Default is TRUE.
#'@param force Logical. By default the pipeline will skip steps if their corresponding output objects already exist. Set to TRUE to force the pipeline to run through all steps. Default is FALSE.
#'@seealso \code{\link{dada2}}, \code{\link{filterAndTrim}}, \code{\link{learnErrors}}, \code{\link{dada}}, \code{\link{mergePairs}}, \code{\link{removeBimeraDenovo}}, \code{\link{assignTaxonomy}}, \code{\link{phyloseq}}
#'@export

dada2.finish_local <- function(
    fastq.path,
    truncFwd.len,
    truncRev.len,
    taxa.db,
    metadata.file,
    maxCores = 4,
    build.tree = TRUE,
    guide.seqs.file = NULL,
    alignment.template.file = NULL,
    fasttree.path = NULL,
    user.output.path = NULL,
    paired = TRUE,
    force = FALSE
) {
    if (length(names(run.env)) == 0) {
        stop("Functions 'initiate.pipeline() and 'dada2.upto.qualPlots()' must be run first.")
    }
    if (is.null(user.output.path)) {
        output <- run.env$output.path
    } else {
        output <- user.output.path
    }
    if (!any(file.exists(list.files(path = output, pattern = "qualPlot.pdf", full.names = T)))) {
        stop("Function 'dada2.upto.qualPlots()' must be run first.")
    }
    if (build.tree) {
        if (length(suppressWarnings(system("which mothur", intern = T))) == 0) {
            stop(
                "It appears you are trying to build a phylogenetic tree, but mothur is not installed on your system. Please install mothur and try again."
            )
        }
        if (is.null(fasttree.path) | !file.exists(fasttree.path)) {
            stop(
                "It appears you are trying to build a phylogenetic tree, but you have not provided a viable path to FastTree."
            )
        }
        if (is.null(guide.seqs.file)) {
            writeLines(guide.seqs.lines, con = "guide_seqs.fasta")
            guide.seqs.file <- "guide_seqs.fasta"
        }
        if (is.null(alignment.template.file)) {
            writeLines(alignment.template.file.lines, con = "template.align")
            alignment.template.file <- "template.align"
        }
        if (!(
            guide.seqs.file %in% list.files() &
            alignment.template.file %in% list.files()
        )) {
            stop(
                paste(
                    "Files", guide.seqs.file, "and", alignment.template.file,
                    "must be in your current directory to build a tree."
                )
            )
        }
    }
    
    filtFs <- file.path(
        fastq.path,
        "Filtered",
        paste0(run.env$sample.names, "_F_filt.fastq.gz")
    )
    names(filtFs) <- run.env$sample.names
    if (paired) {
        filtRs <- file.path(
            fastq.path,
            "Filtered",
            paste0(run.env$sample.names, "_R_filt.fastq.gz")
        )
        names(filtRs) <- run.env$sample.names
    }
    
    out.file <- file.path(output, "filter_and_trim_numbers.rds")
    if (file.exists(out.file) & !force) {
        my.cat("Filtering and trimming completed previously, skipping...")
        out <- readRDS(out.file)
        my.cat("\tDONE:")
        print(head(out))
    } else {
        my.cat("Filtering and trimming...")
        if (paired) {
            out <- filterAndTrim(
                run.env$fnFs, filtFs,
                run.env$fnRs, filtRs,
                truncLen = c(truncFwd.len, truncRev.len),
                maxN = 0,
                maxEE = c(2, 2),
                truncQ = 2,
                rm.phix = TRUE,
                compress = TRUE,
                multithread = maxCores
            )
        } else {
            out <- filterAndTrim(
                run.env$fnFs, filtFs,
                truncLen = truncFwd.len,
                maxN = 0,
                maxEE = 2,
                truncQ = 2,
                rm.phix = TRUE,
                compress = TRUE,
                multithread = maxCores
            )
        }
        my.cat("\tDONE:")
        print(head(out))
        saveRDS(out, file = out.file)
    }
    
    
    filtFs <- filtFs[file.exists(filtFs)]
    if (paired) {
        filtRs <- filtRs[file.exists(filtRs)]
    }
    
    err.files <- file.path(output, "errF.rds")
    err.plot.files <- file.path(output, "errF_plot.pdf")
    if (paired) {
        err.files <- c(err.files, file.path(output, "errR.rds"))
        err.plot.files <- c(err.plot.files, file.path(output, "errR_plot.pdf"))
    }
    
    if (all(file.exists(err.files)) & !force) {
        my.cat("Errors learned previously, skipping...")
        errF <- readRDS(err.files[1])
        if (paired) {errR <- readRDS(err.files[2])}
    } else {
        my.cat("Learning errors and making error plots...")
        errF <- learnErrors(filtFs, multithread = maxCores)
        saveRDS(errF, file = err.files[1])
        errF.plot <- plotErrors(errF, nominalQ = TRUE)
        ggsave(errF.plot, file = err.plot.files[1])
        if (paired) {
            errR <- learnErrors(filtFs, multithread = maxCores)
            saveRDS(errR, file = err.files[2])
            errR.plot <- plotErrors(errR, nominalQ = TRUE)
            ggsave(errR.plot, file = err.plot.files[2])
        }
    }
    my.cat("\tDONE")
    
    dada.files <- file.path(output, "dadaFs.rds")
    if (paired) {
        dada.files <- c(dada.files, file.path(output, "dadaRs.rds"))
    }
    if (all(file.exists(dada.files)) & !force) {
        my.cat("dada-ing done previously, skipping...")
        dadaFs <- readRDS(dada.files[1])
        if (paired) {dadaRs <- readRDS(dada.files[2])}
    } else {
        my.cat("dada-ing...")
        dadaFs <- dada(filtFs, err = errF, multithread = maxCores)
        saveRDS(dadaFs, file = dada.files[1])
        if (paired) {
            dadaRs <- dada(filtRs, err = errR, multithread = maxCores)
            saveRDS(dadaRs, file = dada.files[2])
        }
    }
    my.cat("\tDONE")
    
    seqtab.file <- file.path(output, "seqtab.rds")
    if (paired) {
        mergers.file <- file.path(output, "mergers.rds")
    }
    
    if (file.exists(seqtab.file) & !force) {
        my.cat("Sequence table file exists, skipping...")
        seqtab <- readRDS(seqtab.file)
    } else {
        if (paired) {
            my.cat("Merging pairs...")
            mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
            merge.results <- sapply(mergers, getN)
            my.cat("\tDONE")
        } else {
            merge.results <- 1
        }
        
        if (sum(merge.results == 0) >= 3 | !paired) {
            if (paired) {
                my.cat("Multiple sequence merging fails, proceeding with FORWARD READS ONLY")
            }
            seqtab <- makeSequenceTable(dadaFs)
            saveRDS(seqtab, file = seqtab.file)
        } else {
            saveRDS(mergers, file = mergers.file)
            seqtab <- makeSequenceTable(mergers)
            saveRDS(seqtab, file = seqtab.file)
        }
    }
    
    seqtab.nochim.file <- file.path(output, "seqtab_nochim.rds")
    
    if (file.exists(seqtab.nochim.file) & !force) {
        my.cat("Bimeras previously removed, skipping...")
        seqtab.nochim <- readRDS(seqtab.nochim.file)
        my.cat(paste("\t\tPercent seqs kept:", sum(seqtab.nochim)/sum(seqtab)))
    } else {
        my.cat("Removing bimeras...")
        seqtab.nochim <- removeBimeraDenovo(
            seqtab,
            method = "consensus",
            multithread = maxCores,
            verbose = TRUE
        )
        my.cat("\tDONE")
        my.cat(paste("\t\tPercent seqs kept:", sum(seqtab.nochim)/sum(seqtab)))
        saveRDS(
            seqtab.nochim,
            file = seqtab.nochim.file
        )
    }
    track.file <- file.path(output, "track.rds")
    if (!file.exists(track.file)) {
        if (sum(merge.results == 0) >= 3) {
            track <- cbind(
                out,
                sapply(dadaFs, getN),
                sapply(dadaRs, getN),
                rowSums(seqtab.nochim)
            )
            colnames(track) <- c(
                "input",
                "filtered",
                "denoisedF",
                "denoisedR",
                "nonchimF"
            )
        } else {
            if (paired) {
                track <- cbind(
                    out,
                    sapply(dadaFs, getN),
                    sapply(dadaRs, getN),
                    merge.results,
                    rowSums(seqtab.nochim)
                )
                colnames(track) <- c(
                    "input",
                    "filtered",
                    "denoisedF",
                    "denoisedR",
                    "merged",
                    "nonchim"
                )
            } else {
                track <- cbind(
                    out,
                    sapply(dadaFs, getN),
                    rowSums(seqtab.nochim)
                )
                colnames(track) <- c(
                    "input",
                    "filtered",
                    "denoisedF",
                    "nonchimF"
                )
            }
        }
        rownames(track) <- run.env$sample.names
        saveRDS(track, file = track.file)
    } else {
        track <- readRDS(track.file)
    }
    my.cat(
        "See 'track.rds' in output directory for how many seqs made it through. First 6 samples look as follows:"
    )
    print(head(track))
    
    
    taxa.file <-  file.path(output, "taxa.rds")
    
    if (file.exists(taxa.file) & !force) {
        my.cat("Taxonomy previously assigned, skipping...")
        taxa <- readRDS(taxa.file)
    } else {
        my.cat("Assigning taxonomy...")
        taxa <- assignTaxonomy(
            seqtab.nochim,
            taxa.db,
            multithread = maxCores
        )
        saveRDS(taxa, file = taxa.file)
    }
    my.cat("\tDONE")
    
    smpl.tbl <- read.file(metadata.file)
    if ("data.table" %in% class(smpl.tbl)) {
        smpl.col <- smpl.tbl[
            ,
            sapply(
                .SD,
                function(col) {
                    sum(col %in% run.env$sample.names) == length(run.env$sample.names)
                }
            )
        ] %>%
            magrittr::extract(., .) %>%
            names()
        smpl.df <- as.data.frame(smpl.tbl)
        row.names(smpl.df) <- smpl.tbl$Sample
    } else {
        smpl.df <- smpl.tbl
    }
    
    smpl.df <- smpl.df[run.env$sample.names, , drop=F]
    ps0 <- phyloseq(
        otu_table(seqtab.nochim, taxa_are_rows = FALSE),
        sample_data(smpl.df),
        tax_table(taxa)
    )
    
    ps1 <- numbered.ASVs(
        ps = ps0,
        # prefix = paste0(proj.name, "_ASV"),
        save.dir = output,
        save.file = "asv_sequences"
    )
    asv.seqs <- readRDS(file.path(output, "asv_sequences.rds"))
    
    seqinr::write.fasta(
        sequences = as.list(asv.seqs),
        names = taxa_names(ps1),
        as.string = TRUE,
        file.out = file.path(output, "asv_sequences.fasta")
    )
    if (build.tree) {
        my.cat("Proceeding with phylogenetic tree:")
        asv.seqs.file <- file.path(output, "asv_sequences.fasta")
        asv.withguides.file  <- file.path(output, "asv_and_guide_seqs.fasta")
        asv.tree.rooted.file <- file.path(output, "asv_NASTaligned_seqs.nwk")
        
        cmd <- paste(
            "cat",
            asv.seqs.file,
            guide.seqs.file,
            ">",
            asv.withguides.file
        )
        system(cmd)
        
        my.cat("Aligning sequences...")
        cmd <- paste0(
            "mothur \"#align.seqs( fasta=",
            asv.withguides.file,
            ", reference=",
            alignment.template.file,
            ", flip=t",
            ", keepdots=t",
            ", processors=", maxCores,
            ", outputdir=",
            output,
            "/ )\""
        )
        system(cmd)
        my.cat("\tDONE")
        
        mothur.output.file <- file.path(output, "asv_and_guide_seqs.align")
        fasttree.log.file <- file.path(output, "fasttree.log")
        fasttree.output.file <- file.path(output, "asv_and_guide_seqs.nwk")
        
        my.cat("Building phylogenetic tree...")
        cmd <- paste0(
            "export OMP_NUM_THREADS=",
            maxCores, "; '",
            fasttree.path, "' -nt -nosupport -quote -gtr -gamma -log ",
            fasttree.log.file,
            " ",
            mothur.output.file,
            " > ",
            fasttree.output.file
            
        )
        system(cmd)
        my.cat("\tDONE")
        asvs.and.guides.tree <- read_tree(fasttree.output.file)
        asvs.and.guides.tree.rooted <- phangorn::midpoint(asvs.and.guides.tree)
        
        guides <- scan(guide.seqs.file, what = "character" )
        guide.ids <- guides[stringr::str_detect(guides, ">" )]
        guide.ids <- stringr::str_remove(guide.ids, ">")
        
        asvs.tree.rooted <- ape::drop.tip(asvs.and.guides.tree.rooted, guide.ids)
        write.tree(asvs.tree.rooted, file = asv.tree.rooted.file)
        
        phy_tree(ps1) <- phy_tree(asvs.tree.rooted)
        system(paste("mv mothur*", output))
    }
    saveRDS(ps1, file = file.path(output, "phyloseq.rds"))
    for (file in c("tmp.txt", guide.seqs.file, alignment.template.file)) {
        if (file.exists(file)) { file.remove(file) }
    }
    my.cat("\tDONE and DONE")
}

# Load RData -----------------------------------------------------------
#   Description: 
#   Input: filename for .RData file
#   Output: Object stored from .RData file


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


# GGSave -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

save_plots_2 <- function(out.dir = output.path, sub.dir = "Figures", obj, name = NULL, wdt = 8, ht = 5, ID = Sys.Date(), file.type = ".pdf"){
  ggsave(filename = paste0(out.dir, "/", sub.dir, "/", ifelse(is.null(name), paste0(wdt,"w",ht,"h"), name), "_", ID, file.type),
         plot = marrangeGrob(obj, nrow = 1, ncol = 1, top=NULL),
         width = wdt, height = ht)
}

# Set Names Beta ----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

set_names_beta <- function(x) {
  setNames(x, methods.beta)
  }

# Flex Table fit to page -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

FitFlextableToPage <- function(ft, pgwidth = 6){
  
  ft_out <- ft %>% autofit()
  
  ft_out <- width(ft_out, width = dim(ft_out)$widths*pgwidth /(flextable_dim(ft_out)$widths))
  return(ft_out)
}

# Extract Stats  -----------------------------------------------------------
#   Description: helper function for ggstatsplot::extract_stats
#   Input: takes in a ggstatsplot object
#   Output: returns a list of tibbles or one dataframe with statistical information
#   Improve: More granularity of what to return (ie, subtitle data or other data, or both?)

extract_stats2 <- function(stats.obj, return.all = F){
  
  # Extract stats information
  tmp.table <- stats.obj %>%
    extract_stats()
  
  # Removes empty list elements
  tmp.table[sapply(tmp.table, is.null)] <- NULL
  
  # Check if return all or not
  if(isTRUE(return.all)){
    
    # Returns a list of tibbles
    return(tmp.table)
    
    # If return.all == F, then returns a dataframe of stats info
  } else {
    
    # Take last list element, convert to dataframe
    tmp.df <- as.data.frame(head(tmp.table, 1)[[1]]) 
    
    # Remove "expression" column
    tmp.df <- tmp.df[1:(length(tmp.df)-1)]
    
    # Returns one dataframe
    return(tmp.df)
    }
  
}


# Display Stats -----------------------------------------------------------
#   Description: Helper function to display stats from extract_stats2()
#   Input: 
#   Output: 

disp_stats <- function (stats.table){
  
  return(
    extract_stats2(stats.table) %>%
      kbl() %>%
      kable_styling(bootstrap_options = c("striped", "hover"))
  )
}

# ggcoefstats_2 -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

ggcoefstats_2 <- function(tmp.mod){
  
  tmp.plot <-combine_plots(
                tmp.mod,  # list(p1, p2, p3, p4),
                plotgrid.args = list(nrow = 2)#,
                # annotation.args = list(
                #   title = "Comparison of life expectancy between 1957 and 2007",
                #   caption = "Source: Gapminder Foundation"
                # )
                )
  
  return(tmp.plot)
}

# Extract Model Name and Family -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

extract_modForm <- function(tmp.mod){
  tmp.form = deparse(as.list(tmp.mod$call)$formula)
  
  return(tmp.form)
}






# Maaslin2 Taxa Names -----------------------------------------------------
#   Description: add taxa names to maaslin2 results output dataframe
#   Input: 
#   Output: 

maaslin2_taxaNames <- function(data, physeq, filt.sig = T){
  tmp.data <- data$results %>%
    arrange(feature) 
  
  # Check for identical "Taxon" column and delete second one
  taxa.names <- taxa.data.table(physeq) %>%
    subset(., select = which(!duplicated(names(.))))
  
  # Assign taxonomy to features and check that they equal each other
  tmp.data.nm <- tmp.data %>%
    arrange(feature) %>%  # Arrange numerically
    mutate(temp <- merge(taxa.names, tmp.data, by.x='Taxon', by.y='feature')[, c(3,6,7,1)]) %>% # add 1 to array to check ASVs match
    arrange(qval) %>% # Rearrange by most significant by adjusted pval (aka qval)
    mutate(feat.equals.taxon = sapply(1:nrow(tmp.data), function(i){
      if(feature[i] == Taxon[i]){TRUE} else {FALSE}
    }))

  # Return dataframe with Maaslin2 results and taxa names
  tmp.data.2 <- tmp.data.nm[1:13]

  # Removes any taxa below 0.05 p-value threshold
  if(isTRUE(filt.sig)){
    tmp.data.2 <- tmp.data.2 %>%
      filter(qval < 0.05)
  }


  return(tmp.data.2)
}



# Maaslin Heatmap ---------------------------------------------------------
#   Description: creates a heat map of the most n sig features
#   Input: 
#   Output: 
#   To do: simplify, erase the matrix building stuff and make that a sep func

maaslin2_heatmap <-
  function(
    # output_results,
    df,
    sig.results.df,
    first_n,
    title = NA,
    cell_value = 'qval',
    data_label = 'data',
    metadata_label = 'metadata',
    border_color = 'grey93',
    color = colorRampPalette(c(pal.dark2[3], "grey90", pal.dark2[2])),
    col_rotate = 90
    ) {
    
    # print(paste0("2) first_n: ", first_n)) # TEST
    # first_n = first_n
    # read MaAsLin output
    # df <- read.table(
    #   output_results,
    #   header = TRUE,
    #   sep = ",",
    #   fill = TRUE,
    #   comment.char = "" ,
    #   check.names = FALSE
    # )
    df <- df$results
    
    # test.count <- 0  # TEST
    # print(paste0("Test.count: ", test.count))  # TEST 0
    # test.count <- test.count + 1  # TEST
    
    title_additional <- ""
    
    title_additional <- ""
    if (!is.na(first_n) & first_n > 0 & first_n < dim(df)[1]) {
      if (cell_value == 'coef') {
        df <- df[order(-abs(df[cell_value])) , ]
      } else{
        df <- df[order(df[cell_value]), ]
      }
      # get the top n features with significant associations
      df_sub <- df[1:first_n,]
      for (first_n_index in seq(first_n, dim(df)[1]))
      {
        if (length(unique(df_sub$feature)) == first_n)
        {
          break
        }
        df_sub <- df[1:first_n_index,]
      }
      # get all rows that have the top N features
      df <- df[which(df$feature %in% df_sub$feature),]
      title_additional <- paste("Top", first_n, sep=" ")
    }
    
    # print(paste0("Test.count: ", test.count))  # TEST 1
    # test.count <- test.count + 1  # TEST
    
    if (dim(df)[1] < 2) {
      print('There are no associations to plot!')
      return(NULL)
    }
    
    metadata <- df$metadata
    data <- df$feature
    dfvalue <- df$value
    value <- NA
    
    # values to use for coloring the heatmap
    # and set the colorbar boundaries
    if (cell_value == "pval") {
      value <- -log(df$pval) * sign(df$coef)
      value <- pmax(-20, pmin(20, value))
      if (is.null(title))
        title <- "(-log(pval)*sign(coeff))"
    } else if (cell_value == "qval") {
      value <- -log(df$qval) * sign(df$coef)
      value <- pmax(-20, pmin(20, value))
      if (is.null(title))
        title <- "(-log(qval)*sign(coeff))"
    } else if (cell_value == "coef") {
      value <- df$coef
      if (is.null(title))
        title <- "(coeff)"
    }
    
    # print(paste0("Test.count: ", test.count))  # TEST 2
    # test.count <- test.count + 1  # TEST
    
    if (title_additional!="") {
      title <- paste(title_additional, "features with significant associations", title, sep=" ")
    } else {
      title <- paste("Significant associations", title, sep=" ")
    }
    
    # identify variables with more than one level present
    verbose_metadata <- c()
    metadata_multi_level <- c()
    for (i in unique(metadata)) {
      levels <- unique(df$value[df$metadata == i])
      if (length(levels) > 1) {
        metadata_multi_level <- c(metadata_multi_level, i)
        for (j in levels) {
          verbose_metadata <- c(verbose_metadata, paste(i, j))
        }
      } else {
        verbose_metadata <- c(verbose_metadata, i)
      }
    }
    
    # print(paste0("Test.count: ", test.count))  # TEST 3
    # test.count <- test.count + 1  # TEST
    
    n <- length(unique(data))
    m <- length(unique(verbose_metadata))
    
    if (n < 2) {
      print(
        paste(
          "There is not enough features in the associations",
          "to create a heatmap plot.",
          "Please review the associations in text output file.")
      )
      return(NULL)
    }
    
    if (m < 2) {
      print(
        paste(
          "There is not enough metadata in the associations",
          "to create a heatmap plot.",
          "Please review the associations in text output file.")
      )
      return(NULL)
    }
    
    a = matrix(0, nrow = n, ncol = m)
    a <- as.data.frame(a)
    
    # print(paste0("Test.count: ", test.count))  # TEST 4
    # test.count <- test.count + 1  # TEST
    
    rownames(a) <- unique(data)
    colnames(a) <- unique(verbose_metadata)
    
    # print(paste0("Test.count: ", test.count))  # TEST 5
    # test.count <- test.count + 1  # TEST
    
    tryCatch({
      for (i in seq_len(dim(df)[1])) {
        current_metadata <- metadata[i]
        if (current_metadata %in% metadata_multi_level) {
          current_metadata <- paste(metadata[i], dfvalue[i])
        }
        if (abs(a[as.character(data[i]), 
                  as.character(current_metadata)]) > abs(value[i]))
          next
        a[as.character(data[i]), as.character(current_metadata)] <- value[i]
      }
      }, warning = function(w) {
        message("Warning: Heatmap > setting up matrix")
        return(NULL)
      })
    
    # for (i in seq_len(dim(df)[1])) {
    #   current_metadata <- metadata[i]
    #   if (current_metadata %in% metadata_multi_level) {
    #     current_metadata <- paste(metadata[i], dfvalue[i])
    #   }
    #   if (abs(a[as.character(data[i]), 
    #             as.character(current_metadata)]) > abs(value[i]))
    #     next
    #   a[as.character(data[i]), as.character(current_metadata)] <- value[i]
    # }
    
    # print(paste0("Test.count: ", test.count))  # TEST 6
    # test.count <- test.count + 1  # TEST
    
    
    # Appends the Phylum and Genus next to the ASV
    tmp.row.names <- paste0(sig.results.df$feature, 
                        " | ",
                        sig.results.df$Phylum, " - ",
                        sig.results.df$Genus)[1:first_n]
    

    
    
    # print(paste0("Test.count: ", test.count))  # TEST 7
    # test.count <- test.count + 1  # TEST
    
    # get the range for the colorbar
    max_value <- ceiling(max(a[,1:(ncol(a)-1)]))
    min_value <- ceiling(min(a[,1:(ncol(a)-1)]))
    range_value <- max(c(abs(max_value),abs(min_value)))
    breaks <- seq(-1*range_value, range_value, by = 1)
    
    # print(paste0("Test.count: ", test.count))  # TEST 8
    # test.count <- test.count + 1  # TEST
    
    # TEST
    # saveRDS(a, file = file.path(paste0(objects.path, "/diffAbundHeatMap", Sys.Date(),".rds")))
    
    p <- NULL
    tryCatch({
      p <-
        pheatmap::pheatmap(
          a,
          cellwidth = 10,
          cellheight = 10,
          main = title,
          fontsize = 6,
          kmeans_k = NA,
          border = TRUE,
          show_rownames = TRUE,
          show_colnames = TRUE,
          scale = "none",
          cluster_rows = FALSE,
          labels_row = tmp.row.names,  # TEST
          cluster_cols = TRUE,
          clustering_distance_rows = "euclidean",
          clustering_distance_cols = "euclidean",
          legend = TRUE,
          border_color = border_color,
          color = color(range_value*2),
          breaks = breaks,
          treeheight_row = 0,
          treeheight_col = 0,
          display_numbers = matrix(ifelse(
            a > 0.0, "+", ifelse(a < 0.0, "-", "")), nrow(a)),
          fontsize_number = 6,
          angle_col = 315
        )
    }, warning = function(w) {
      message("Error in heatmap")
    })
    return(p)
  }



# Maaslin melt df heatmap -------------------------------------------------
#   Description: returns a matrix of the top correlation values for n asvs
#   Input: 
#   Output: 


# # Example Code
# test <- maaslin2_topASVmatrix(
#   df = diffAbund,
#   sig.results.df = diffAbund.sig,
#   title = "",
#   first_n = 20)

maaslin2_topASVmatrix <-
  function(
    output_results,
    df,
    sig.results.df,
    first_n,
    title = NA,
    cell_value = 'qval',
    data_label = 'data',
    metadata_label = 'metadata'
  ) {
    
    df <- df$results
    
    
    title_additional <- ""
    if (!is.na(first_n) & first_n > 0 & first_n < dim(df)[1]) {
      if (cell_value == 'coef') {
        df <- df[order(-abs(df[cell_value])) , ]
      } else{
        df <- df[order(df[cell_value]), ]
      }
      # get the top n features with significant associations
      df_sub <- df[1:first_n,]
      for (first_n_index in seq(first_n, dim(df)[1]))
      {
        if (length(unique(df_sub$feature)) == first_n)
        {
          break
        }
        df_sub <- df[1:first_n_index,]
      }
      # get all rows that have the top N features
      df <- df[which(df$feature %in% df_sub$feature),]
      title_additional <- paste("Top", first_n, sep=" ")
    }
    
    if (dim(df)[1] < 2) {
      print('There are no associations to plot!')
      return(NULL)
    }
    
    metadata <- df$metadata
    data <- df$feature
    dfvalue <- df$value
    value <- NA
    
    # values to use for coloring the heatmap
    # and set the colorbar boundaries
    if (cell_value == "pval") {
      value <- -log(df$pval) * sign(df$coef)
      value <- pmax(-20, pmin(20, value))
      if (is.null(title))
        title <- "(-log(pval)*sign(coeff))"
    } else if (cell_value == "qval") {
      value <- -log(df$qval) * sign(df$coef)
      value <- pmax(-20, pmin(20, value))
      if (is.null(title))
        title <- "(-log(qval)*sign(coeff))"
    } else if (cell_value == "coef") {
      value <- df$coef
      if (is.null(title))
        title <- "(coeff)"
    }
    
    if (title_additional!="") {
      title <- paste(title_additional, "features with significant associations", title, sep=" ")
    } else {
      title <- paste("Significant associations", title, sep=" ")
    }
    
    # identify variables with more than one level present
    verbose_metadata <- c()
    metadata_multi_level <- c()
    for (i in unique(metadata)) {
      levels <- unique(df$value[df$metadata == i])
      if (length(levels) > 1) {
        metadata_multi_level <- c(metadata_multi_level, i)
        for (j in levels) {
          verbose_metadata <- c(verbose_metadata, paste(i, j))
        }
      } else {
        verbose_metadata <- c(verbose_metadata, i)
      }
    }
    
    n <- length(unique(data))
    m <- length(unique(verbose_metadata))
    
    if (n < 2) {
      print(
        paste(
          "There is not enough features in the associations",
          "to create a heatmap plot.",
          "Please review the associations in text output file.")
      )
      return(NULL)
    }
    
    if (m < 2) {
      print(
        paste(
          "There is not enough metadata in the associations",
          "to create a heatmap plot.",
          "Please review the associations in text output file.")
      )
      return(NULL)
    }
    
    a = matrix(0, nrow = n, ncol = m)
    a <- as.data.frame(a)
    
    rownames(a) <- unique(data)
    colnames(a) <- unique(verbose_metadata)
    
    for (i in seq_len(dim(df)[1])) {
      current_metadata <- metadata[i]
      if (current_metadata %in% metadata_multi_level) {
        current_metadata <- paste(metadata[i], dfvalue[i])
      }
      if (abs(a[as.character(data[i]), 
                as.character(current_metadata)]) > abs(value[i]))
        next
      a[as.character(data[i]), as.character(current_metadata)] <- value[i]
    }
    
    # Appends the Phylum and Genus next to the ASV
    rownames(a) <- paste0(unique(sig.results.df$feature), 
                          " | ", 
                          sig.results.df$Phylum, " - ", 
                          sig.results.df$Genus)[1:first_n]
    
    # get the range for the colorbar
    max_value <- ceiling(max(a))
    min_value <- ceiling(min(a))
    range_value <- max(c(abs(max_value),abs(min_value)))
    breaks <- seq(-1*range_value, range_value, by = 1)
    
    return(a)
  }


# Save Maaslin2 Heatmap -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

# save_heatmap <-
#   function(
#     results_file,
#     sig.results.df,
#     heatmap_file,
#     figures_folder,
#     first_n,
#     title = NULL,
#     cell_value = "qval",
#     data_label = 'data',
#     metadata_label = 'metadata',
#     border_color = "grey93",
#     color = colorRampPalette(c(pal.dark2[3], "grey90", pal.dark2[2])),
#     ) {
#     
#     print(paste0("1) first_n: ", first_n)) # TEST
#     
#     # generate a heatmap and save it to a pdf and as a png
#     heatmap <-
#       maaslin2_heatmap(
#         results_file,
#         sig.results.df,
#         title,
#         cell_value,
#         data_label,
#         metadata_label,
#         border_color,
#         color,
#         first_n = first_n)
#     
#     if (!is.null(heatmap)) {
#       pdf(heatmap_file, width = 5, height = 8)
#       print(heatmap)
#       dev.off()
#       
#       png_file <- file.path(paste0(output.path, "/Figures/diffAbund-Plot_", 
#                                    paste0(c(diff.abund.param[c(1,11)]), collapse="-") , ".png"))
#       png(png_file, res = 150, height = 1500, width = 1100)
#       print(heatmap)
#       dev.off()
#     }
#     
#     return(heatmap)
#     
#     
#     
#   }



# Differential Abundance Plot ---------------------------------------------
#   Description: 
#   Input: 
#   Output: 

plot_diffAbund <- function(data, x.var = "Phylum", fill.var = "Genus", title = ""){
  data %>%
    ggplot(aes(y = eval(parse(text = x.var)) )) +
    geom_bar(aes(fill = eval(parse(text = fill.var)) ), 
             position = position_dodge2(width = 0.9, preserve = "single")) + 
    facet_wrap(.~value, scale="fixed", ncol = 1) + 
    scale_fill_manual(values = c(pal.dark2, pal.paired)) + 

    labs(title = ifelse(title != "", title, "Counts of significantly abundant taxa relative to their reference variables"),
         # caption = "",
         y = x.var,
         x = "Counts"
         ) + scale_x_continuous(breaks = scales::breaks_pretty(10)) +
    guides(fill=guide_legend(title=fill.var))
}

plot_diffAbund_2 <- function(data, x.var="Phylum", fill.var="coef"){
  subset(data, !is.na(Genus)) %>%
  ggplot(aes(x=fct_rev(eval(parse(text = x.var)) ), 
             fill= Genus, # factor(sign(eval(parse(text = fill.var)) )),
             color = factor(sign(eval(parse(text = fill.var)) )))) + # Genus
    #geom_bar(position="fill", size = 1.5) + 
    geom_bar(aes(fill = Genus ), 
             position = position_dodge2(width = 0.9, preserve = "single"), 
             size = ifelse(length(unique(data$Genus)) < 8,1,.5) ) + 
    
    facet_wrap(.~value, ncol = 1) +
    scale_fill_manual(values = c(pal.dark2, pal.paired, pal.set1)) +
    scale_color_manual(name="Correlation", 
                      labels = c("Negative","Positive"), 
                      values=c(pal.accent[3], pal.accent[2])) +
    theme(axis.text.x = element_text(angle = 0))+ #,
    # axis.text.y = element_text(color = rev(relAbund_colors$Color))) +  
    coord_flip() + labs(y="Relative Abundance Significant Taxa", 
                        x="Genus (color coded by Phylum)") 
}




# Differential Abundance - Maaslin -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

diffAbund_Maaslin2 <- function(physeq,
                               physeq.ID,
                               fixed.effects,  # Explanatory Variables
                               random.effects = c(""),  # Since same tanks were sampled multiple times
                               reference.vars = c(""), #c("Diet,Gemma", "Timepoint,Initial", "Sex,M", "Exposure,Unexposed"),
                               min.abundance = 0.01, #0.0001,  # Conservative threshold 
                               min.preval = 0.1,  # Default 
                               normalization = "NONE",  # Data are rarefied/norm already
                               transformation = "NONE",   # None, since using NegBin model
                               analysis.method = "NEGBIN",  # negative binomial glm
                               correction = "BH",  # Multiple test correction
                               n.cores = num.cores, 
                               ID = analysis.ID,
                               add.info = "",
                               num_taxa = NULL) {
  
  ps.obj <- physeq
  
  # Instantiates a list with empty elements, for consistency and in case of errors
  tmp.obj <- list("Parameters" = vector(mode='list', length=16),
                  "Results" = data.frame(), 
                  "Sig.Res" = NULL, 
                  "Heatmap" = NULL)
  
  tmp.obj[["Parameters"]] <- list(physeq.ID = physeq.ID,  # Phyloseq ID ("ps.T0")
                                f.eff = fixed.effects,  # Explanatory Variables
                                r.eff = random.effects,  # Since same tanks were sampled multiple times
                                ref = reference.vars, #c("Diet,Gemma", "Timepoint,Initial", "Sex,M", "Exposure,Unexposed"),
                                min.ab = min.abundance,  # Conservative threshold 
                                min.prv = min.preval,  # Default 
                                norm = normalization,  # Data are rarefied/norm already
                                trns = transformation,   # None, since using NegBin model
                                an.meth = analysis.method,  # negative binomial glm
                                corr = correction,  # Multiple test correction
                                n.cores = num.cores, 
                                ID = analysis.ID,
                                Add = add.info)
  
  
    tmp.obj[["Results"]] <- #tryCatch({
    # Insert Code
    diff_Abund_Maaslin2(physeq = ps.obj, param = tmp.obj[["Parameters"]], physeq.ID)
    
  # }, warning = function(w) {
  #   message("Warning: Results")
  #   return(NULL)
  # })
  

  
    tmp.obj[["Sig.Res"]] <- #tryCatch({
    # Insert Code
    maaslin2_taxaNames(data = tmp.obj[["Results"]], 
                                               physeq = ps.obj,
                                               filt.sig = T)
    
  # }, warning = function(w) {
  #   message("Warning: Sig.Res")
  #   return(NULL)
  # })
  
  
  n_taxa <- ifelse(is.null(num_taxa), nrow(tmp.obj[["Sig.Res"]]), num_taxa)
  
  # tmp.obj[["Heatmap"]] <-  #tryCatch({
  #   # Insert Code
  #   maaslin2_heatmap(
  #     df = tmp.obj[["Results"]],
  #     sig.results.df = tmp.obj[["Sig.Res"]],
  #     title = "",
  #     first_n = n_taxa)
    
  # }, warning = function(w) {
  #   message("Warning: Heatmap")
  #   return(NULL)
  # })
  
  
  tmp.obj[["Stats.Table"]] <- #tryCatch({
    # Insert Code
    tmp.obj[["Sig.Res"]] %>%
      head(n = n_taxa) %>%
      arrange(coef) %>%
      flextable() %>%
      colformat_double(j = c(4,5), digits = 4, ) %>%
      set_formatter(values = list("pval" = p_val_format) )%>%
      set_formatter(values = list("qval" = p_val_format) )%>%
      autofit()
    
  # }, warning = function(w) {
  #   message("Warning: Stats.Table")
  #   return(NULL)
  # })

  
  return(tmp.obj)
}







# Try Catch Errors -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

# tryCatch({
#   # Insert Code
#    
#    
# }, warning = function(w) {
#   message("Warning: ")
# })



# Interaction Matrix -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

create_intMatrix <- function(data, vars, terms = 2, set.samp.name = T){
  
  # Create interaction matrix with (0's and 1's)
  tmp.data <- as.data.table(
    model.matrix( 
    as.formula( paste0("~ (", paste0(vars, collapse = "+"), ")^", terms, "- 1") ), data) ) 
  
  tmp.data <- tmp.data %>%
    dplyr::arrange(.by_group = T)

  if(isTRUE(set.samp.name)){
    tmp.data$Sample <- data$Sample
  }

  return(tmp.data)
}




# Merge Datatable and IntMatrix -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

merge_DTandIntMatrix <- function(datatable, int.matrix, vars){
  
  # Arranges datatable by variables of interest
  tmp.dt <- datatable %>%
    dplyr::arrange(!!!vars, .by_group = T)
  
  # Merges 
  tmp.merge.dt <- tmp.dt %>%
    dplyr::left_join(int.matrix, by = "Sample")
  
}




# Interaction Matrix Datatable -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

get_intMatrixDT <- function(tmp.int.dt, tmp.int.vars, tmp.terms = 2){
  # Create an interaction matrix 
  tmp.int.mat <- create_intMatrix(data = tmp.int.dt,
                                  vars = tmp.int.vars,
                                  terms = tmp.terms,
                                  set.samp.name = T)
  
  # Datatable including interaction matrix
  tmp.dt.int.mat <- merge_DTandIntMatrix(datatable = tmp.int.dt, 
                                         int.matrix = tmp.int.mat, 
                                         vars = tmp.int.vars
                                          )
  
  return(tmp.dt.int.mat)
}




# Merge PS and Int Matrix DT -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

mergePSandIntMat <- function(physeq, int.matrix){
  
  tmp.ps <- physeq
  
  
  for(col in colnames(int.matrix)){
    sample_data(tmp.ps)[[col]] <- (int.matrix[[col]])
  }
  
  return(tmp.ps)
}




# Function Name -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

get_intMatPSobj <- function(physeq, data, vars){
  
  tmp.ps <- physeq
  
  tmp.int.matrix <- get_intMatrixDT(tmp.int.dt = data,
                          tmp.int.vars = vars,
                          tmp.terms = 2
                          )
  
  tmp.ps.merged <- mergePSandIntMat(physeq = tmp.ps, 
                              int.matrix = tmp.int.matrix)
  
  
  return(tmp.ps.merged)
}




# ANCOM -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

# OTU table should be a matrix/data.frame with each feature in rows and sample in columns. 
# Metadata should be a matrix/data.frame containing the sample identifier. 

# Data Pre-Processing
feature_table_pre_process = function(feature_table, meta_data, sample_var, group_var = NULL, 
                                     out_cut = 0.05, zero_cut = 0.90, lib_cut, neg_lb){
  feature_table = data.frame(feature_table, check.names = FALSE)
  meta_data = data.frame(meta_data, check.names = FALSE)
  # Drop unused levels
  meta_data[] = lapply(meta_data, function(x) if(is.factor(x)) factor(x) else x)
  # Match sample IDs between metadata and feature table
  sample_ID = data.table::intersect(meta_data[, sample_var], colnames(feature_table))
  feature_table = feature_table[, sample_ID]
  meta_data = meta_data[match(sample_ID, meta_data[, sample_var]), ]
  
  # 1. Identify outliers within each taxon
  if (!is.null(group_var)) {
    group = meta_data[, group_var]
    z = feature_table + 1 # Add pseudo-count (1) 
    f = log(z); f[f == 0] = NA; f = colMeans(f, na.rm = T)
    f_fit = lm(f ~ group)
    e = rep(0, length(f)); e[!is.na(group)] = residuals(f_fit)
    y = t(t(z) - e)
    
    outlier_check = function(x){
      # Fitting the mixture model using the algorithm of Peddada, S. Das, and JT Gene Hwang (2002)
      mu1 = quantile(x, 0.25, na.rm = T); mu2 = quantile(x, 0.75, na.rm = T)
      sigma1 = quantile(x, 0.75, na.rm = T) - quantile(x, 0.25, na.rm = T); sigma2 = sigma1
      pi = 0.75
      n = length(x)
      epsilon = 100
      tol = 1e-5
      score = pi*dnorm(x, mean = mu1, sd = sigma1)/((1 - pi)*dnorm(x, mean = mu2, sd = sigma2))
      while (epsilon > tol) {
        grp1_ind = (score >= 1)
        mu1_new = mean(x[grp1_ind]); mu2_new = mean(x[!grp1_ind])
        sigma1_new = sd(x[grp1_ind]); if(is.na(sigma1_new)) sigma1_new = 0
        sigma2_new = sd(x[!grp1_ind]); if(is.na(sigma2_new)) sigma2_new = 0
        pi_new = sum(grp1_ind)/n
        
        para = c(mu1_new, mu2_new, sigma1_new, sigma2_new, pi_new)
        if(any(is.na(para))) break
        
        score = pi_new * dnorm(x, mean = mu1_new, sd = sigma1_new)/
          ((1-pi_new) * dnorm(x, mean = mu2_new, sd = sigma2_new))
        
        epsilon = sqrt((mu1 - mu1_new)^2 + (mu2 - mu2_new)^2 + 
                         (sigma1 - sigma1_new)^2 + (sigma2 - sigma2_new)^2 + (pi - pi_new)^2)
        mu1 = mu1_new; mu2 = mu2_new; sigma1 = sigma1_new; sigma2 = sigma2_new; pi = pi_new
      }
      
      if(mu1 + 1.96 * sigma1 < mu2 - 1.96 * sigma2){
        if(pi < out_cut){
          out_ind = grp1_ind
        }else if(pi > 1 - out_cut){
          out_ind = (!grp1_ind)
        }else{
          out_ind = rep(FALSE, n)
        }
      }else{
        out_ind = rep(FALSE, n)
      }
      return(out_ind)
    }
    out_ind = matrix(FALSE, nrow = nrow(feature_table), ncol = ncol(feature_table))
    out_ind[, !is.na(group)] = t(apply(y, 1, function(i) 
      unlist(tapply(i, group, function(j) outlier_check(j)))))
    
    feature_table[out_ind] = NA
  }
  
  # 2. Discard taxa with zeros  >=  zero_cut
  zero_prop = apply(feature_table, 1, function(x) sum(x == 0, na.rm = T)/length(x[!is.na(x)]))
  taxa_del = which(zero_prop >= zero_cut)
  if(length(taxa_del) > 0){
    feature_table = feature_table[- taxa_del, ]
  }
  
  # 3. Discard samples with library size < lib_cut
  lib_size = colSums(feature_table, na.rm = T)
  if(any(lib_size < lib_cut)){
    subj_del = which(lib_size < lib_cut)
    feature_table = feature_table[, - subj_del]
    meta_data = meta_data[- subj_del, ]
  }
  
  # 4. Identify taxa with structure zeros
  if (!is.null(group_var)) {
    group = factor(meta_data[, group_var])
    present_table = as.matrix(feature_table)
    present_table[is.na(present_table)] = 0
    present_table[present_table != 0] = 1
    
    p_hat = t(apply(present_table, 1, function(x)
      unlist(tapply(x, group, function(y) mean(y, na.rm = T)))))
    samp_size = t(apply(feature_table, 1, function(x)
      unlist(tapply(x, group, function(y) length(y[!is.na(y)])))))
    p_hat_lo = p_hat - 1.96 * sqrt(p_hat * (1 - p_hat)/samp_size)
    
    struc_zero = (p_hat == 0) * 1
    # Whether we need to classify a taxon into structural zero by its negative lower bound?
    if(neg_lb) struc_zero[p_hat_lo <= 0] = 1
    
    # Entries considered to be structural zeros are set to be 0s
    struc_ind = struc_zero[, group]
    feature_table = feature_table * (1 - struc_ind)
    
    colnames(struc_zero) = paste0("structural_zero (", colnames(struc_zero), ")")
  }else{
    struc_zero = NULL
  }
  
  # 5. Return results
  res = list(feature_table = feature_table, meta_data = meta_data, structure_zeros = struc_zero)
  return(res)
}

# ANCOM main function
ANCOM = function(feature_table, meta_data, struc_zero = NULL, main_var, p_adj_method = "BH", 
                 alpha = 0.05, adj_formula = NULL, rand_formula = NULL, ...){
  # OTU table transformation: 
  # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
  if (!is.null(struc_zero)) {
    num_struc_zero = apply(struc_zero, 1, sum)
    comp_table = feature_table[num_struc_zero == 0, ]
  }else{
    comp_table = feature_table
  }
  comp_table = log(as.matrix(comp_table) + 1)
  n_taxa = dim(comp_table)[1]
  taxa_id = rownames(comp_table)
  n_samp = dim(comp_table)[2]
  
  # Determine the type of statistical test and its formula.
  if (is.null(rand_formula) & is.null(adj_formula)) {
    # Basic model
    # Whether the main variable of interest has two levels or more?
    if (length(unique(meta_data%>%pull(main_var))) == 2) {
      # Two levels: Wilcoxon rank-sum test
      tfun = stats::wilcox.test
    } else{
      # More than two levels: Kruskal-Wallis test
      tfun = stats::kruskal.test
    }
    # Formula
    tformula = formula(paste("x ~", main_var, sep = " "))
  }else if (is.null(rand_formula) & !is.null(adj_formula)) {
    # Model: ANOVA
    tfun = stats::aov
    # Formula
    tformula = formula(paste("x ~", main_var, "+", adj_formula, sep = " "))
  }else if (!is.null(rand_formula)) {
    # Model: Mixed-effects model
    tfun = nlme::lme
    # Formula
    if (is.null(adj_formula)) {
      # Random intercept model
      tformula = formula(paste("x ~", main_var))
    }else {
      # Random coefficients/slope model
      tformula = formula(paste("x ~", main_var, "+", adj_formula))
    }
  }
  
  # Calculate the p-value for each pairwise comparison of taxa.
  p_data = matrix(NA, nrow = n_taxa, ncol = n_taxa)
  colnames(p_data) = taxa_id
  rownames(p_data) = taxa_id
  for (i in 1:(n_taxa - 1)) {
    # Loop through each taxon.
    # For each taxon i, additive log ratio (alr) transform the OTU table using taxon i as the reference.
    # e.g. the first alr matrix will be the log abundance data (comp_table) recursively subtracted 
    # by the log abundance of 1st taxon (1st column) column-wisely, and remove the first i columns since:
    # the first (i - 1) columns were calculated by previous iterations, and
    # the i^th column contains all zeros.
    alr_data = apply(comp_table, 1, function(x) x - comp_table[i, ]) 
    # apply(...) allows crossing the data in a number of ways and avoid explicit use of loop constructs.
    # Here, we basically want to iteratively subtract each column of the comp_table by its i^th column.
    alr_data = alr_data[, - (1:i), drop = FALSE]
    n_lr = dim(alr_data)[2] # number of log-ratios (lr)
    alr_data = cbind(alr_data, meta_data) # merge with the metadata
    
    # P-values
    if (is.null(rand_formula) & is.null(adj_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        suppressWarnings(tfun(tformula, 
                              data = data.frame(x, alr_data, 
                                                check.names = FALSE))$p.value)
      }
      ) 
    }else if (is.null(rand_formula) & !is.null(adj_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        fit = tfun(tformula, 
                   data = data.frame(x, alr_data, check.names = FALSE), 
                   na.action = na.omit)
        summary(fit)[[1]][main_var, "Pr(>F)"]
      }
      )
    }else if (!is.null(rand_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        fit = tfun(fixed = tformula, 
                   data = data.frame(x, alr_data, check.names = FALSE),
                   random = formula(rand_formula),
                   na.action = na.omit, ...)
        anova(fit)[main_var, "p-value"]
      }
      ) 
    }
  }
  # Complete the p-value matrix.
  # What we got from above iterations is a lower triangle matrix of p-values.
  p_data[upper.tri(p_data)] = t(p_data)[upper.tri(p_data)]
  diag(p_data) = 1 # let p-values on diagonal equal to 1
  p_data[is.na(p_data)] = 1 # let p-values of NA equal to 1
  
  # Multiple comparisons correction.
  q_data = apply(p_data, 2, function(x) p.adjust(x, method = p_adj_method))
  
  # Calculate the W statistic of ANCOM.
  # For each taxon, count the number of q-values < alpha.
  W = apply(q_data, 2, function(x) sum(x < alpha))
  
  # Organize outputs
  out_comp = data.frame(taxa_id, W, row.names = NULL, check.names = FALSE)
  # Declare a taxon to be differentially abundant based on the quantile of W statistic.
  # We perform (n_taxa - 1) hypothesis testings on each taxon, so the maximum number of rejections is (n_taxa - 1).
  out_comp = out_comp%>%mutate(detected_0.9 = ifelse(W > 0.9 * (n_taxa -1), TRUE, FALSE),
                               detected_0.8 = ifelse(W > 0.8 * (n_taxa -1), TRUE, FALSE),
                               detected_0.7 = ifelse(W > 0.7 * (n_taxa -1), TRUE, FALSE),
                               detected_0.6 = ifelse(W > 0.6 * (n_taxa -1), TRUE, FALSE))
  
  # Taxa with structural zeros are automatically declared to be differentially abundant
  if (!is.null(struc_zero)){
    out = data.frame(taxa_id = rownames(struc_zero), W = Inf, detected_0.9 = TRUE, 
                     detected_0.8 = TRUE, detected_0.7 = TRUE, detected_0.6 = TRUE, 
                     row.names = NULL, check.names = FALSE)
    out[match(taxa_id, out$taxa_id), ] = out_comp
  }else{
    out = out_comp
  }
  
  # Draw volcano plot
  # Calculate clr
  clr_table = apply(feature_table, 2, clr)
  # Calculate clr mean difference
  eff_size = apply(clr_table, 1, function(y) 
    lm(y ~ x, data = data.frame(y = y, 
                                x = meta_data %>% pull(main_var),
                                check.names = FALSE))$coef[-1])
  
  if (is.matrix(eff_size)){
    # Data frame for the figure
    dat_fig = data.frame(taxa_id = out$taxa_id, t(eff_size), y = out$W, check.names = FALSE) %>% 
      mutate(zero_ind = factor(ifelse(is.infinite(y), "Yes", "No"), levels = c("Yes", "No"))) %>%
      gather(key = group, value = x, rownames(eff_size))
    # Replcace "x" to the name of covariate
    dat_fig$group = sapply(dat_fig$group, function(x) gsub("x", paste0(main_var, " = "), x))
    # Replace Inf by (n_taxa - 1) for structural zeros
    dat_fig$y = replace(dat_fig$y, is.infinite(dat_fig$y), n_taxa - 1)
    
    fig = ggplot(data = dat_fig) + aes(x = x, y = y) + 
      geom_point(aes(color = zero_ind)) + 
      facet_wrap(~ group) +
      labs(x = "CLR mean difference", y = "W statistic") +
      scale_color_discrete(name = "Structural zero", drop = FALSE) + 
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "top",
            strip.background = element_rect(fill = "white"))
    fig  
  } else{
    # Data frame for the figure
    dat_fig = data.frame(taxa_id = out$taxa_id, x = eff_size, y = out$W) %>% 
      mutate(zero_ind = factor(ifelse(is.infinite(y), "Yes", "No"), levels = c("Yes", "No")))
    # Replace Inf by (n_taxa - 1) for structural zeros
    dat_fig$y = replace(dat_fig$y, is.infinite(dat_fig$y), n_taxa - 1)
    
    fig = ggplot(data = dat_fig) + aes(x = x, y = y) + 
      geom_point(aes(color = zero_ind)) + 
      labs(x = "CLR mean difference", y = "W statistic") +
      scale_color_discrete(name = "Structural zero", drop = FALSE) + 
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "top")
    fig
  }
  
  res = list(out = out, fig = fig)
  return(res)
}





# Function Name -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

agglom_taxa <- function(physeq,
                        rank.name = "Genus"){
  
  tmp.ps <- physeq
  tmp.ps <- tax_glom(tmp.ps, taxrank = rank.name )  # "compress" or agglom taxa of similar rank
  taxa_names(tmp.ps) <- taxa.data.table(tmp.ps)[[rank.name]]  # rename asvs to taxa rank level
  
  warning(paste0(ntaxa(tmp.ps), " were aggolomerated."))
  
  return(tmp.ps)
}




# Wilcoxon tests (one and two vars, custom) -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

pair.wilcox.test_1 <- function(
  data,
  respVar,
  expVar1,
  tmp.list = c(""),
  subset = NULL
){
  test.res <- lapply(tmp.list, function(x){
    data <- ifelse(is.null(subset), subset(data, eval(parse(text = subset)) == x), data )
      pairwise.wilcox.test(x=eval(parse(text = paste0("data$", respVar))),
                           g=eval(parse(text = paste0("data$", expVar1))),
                           p.adjust.method = "BH") %>% tidy() %>% 
        mutate(metric = x, .before = 1) %>%
        mutate(sig = ifelse(p.value <= 0.05, "*", ""))
  }) %>% bind_rows() %>% arrange(metric) %>%
    flextable() %>%
    # merge_v(j = c(1)) %>% 
    set_caption(paste0("Pairwise wilcoxon, p.adj: BH. ", respVar, " ~ ", expVar1)) %>%
    set_formatter(values = list("p.value" = p_val_format) )
  
  return(test.res)
}

pair.wilcox.test_2 <- function(
  data,
  respVar,
  expVar1,
  expVar2,
  tmp.list,
  subset
){
  test.res <- lapply(unique(data[[expVar1]]), function(y){
                data <- subset(data, eval(parse(text = expVar1)) == y)
                lapply(tmp.list, function(x){
                  data <- subset(data, eval(parse(text = subset)) == x)
                  pairwise.wilcox.test(x=eval(parse(text = paste0("data$", respVar))),
                                       g=eval(parse(text = paste0("data$", expVar2))),
                                       p.adjust.method = "BH") %>% tidy() %>% 
                    mutate(metric = x, .before = 1) %>%
                    mutate(sig = ifelse(p.value <= 0.05, "*", ""))
                }) %>% bind_rows() %>% 
                  mutate(subset = y, .after = 1)
              }) %>% bind_rows() %>% arrange(metric) %>%
                flextable() %>%
                # merge_v(j = c(1)) %>% 
                set_caption(paste0("Pairwise wilcoxon, p.adj: BH. ", respVar, " ~ ", expVar1, ":", expVar2)) %>%
                set_formatter(values = list("p.value" = p_val_format) )
  
  return(test.res)
}




# Get Prevalance -----------------------------------------------------------
#   Description: Chris Gaulke 
#   Source: https://github.com/chrisgaulke/pcap_2018/blob/master/gaulke_et_al_pcap_microbiome_models.R
#   Input: 
#   Output: 


get_prevalence <- function(df, id.col, value.col){
  #this function will take an long form data frame and
  #calculate the prevalence for each taxa in id.col
  #id.col: the name of the column that contains taxa
  #value.col: the column that contains the taxa abd information
  
  prev.vec <- NULL
  seq_id <- unique(df[[id.col]])
  for(i in 1:length(seq_id)){
    sub.vec <- unlist(df[which(df[,id.col] == seq_id[i]),value.col])
    #turn into a logical (i.e., T/F). Since T=1 and F=0 summing gives us an
    #count of the individuals with non-zero abd for this taxa
    seq.prev <- sum(as.logical(sub.vec)) / length(sub.vec)
    prev.vec <- c(prev.vec, seq.prev)
  }
  names(prev.vec) <- seq_id
  return(prev.vec)
}


# Jittered Dist Plot -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

dist_plot_jittered <- function(
  
){
  
  
  
}



# Function Name -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 

# Function Name -----------------------------------------------------------
#   Description: 
#   Input: 
#   Output: 











