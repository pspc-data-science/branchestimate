
#' Extractor of COVID-19 data from the NYT.
#'
#' This function extracts county/state level COVID-19 cumulative count data from the NYT: https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv
#' The function filters only those regions that exhibit exponential growth. The filtering occurs through a sign check
#' on a quadratic regression over the selected time window.
#'
#' @param date_start The start date of the time window. Default 2020-05-01
#' @param date_end The end date of the time window. Default 2020-08-20
#' @param county A logical TRUE = county level data; FALSE = state level data
#'
#' @return A tibble containing the regional data including fips number for county level data.
#'
#' @export
US_covid_data_creator<- function(date_start = "2020-05-01" , date_end = "2020-08-20", county = TRUE){

  covid_tibble_US<- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv") %>%
    unite(county, county, state, sep = "; ") %>%
    select(date, region = county, n_infected = cases, deaths, fips) %>% arrange(region, date) %>%
    filter(!str_detect(region, "Unknown")) %>% group_by(region) %>%
    filter(date> as_date(date_start)) %>%
    filter(date< as_date(date_end)) %>%
    nest() %>% ungroup()

  day_fun<- function(dF){

    x<- dF$date - dF$date[1]
    x<- as.numeric(x)

    dF$time<- x
    return(dF)

  }

  make_regressions_old<- function(dF){

    mod_1<- lm(dF$n_infected ~ dF$time) %>% summary()
    mod_2<- lm(log(dF$n_infected) ~ dF$time) %>% summary()

    if(is.nan(mod_1$adj.r.squared) | is.nan(mod_2$adj.r.squared) ){
      return(FALSE)
    }else{
      if(mod_2$adj.r.squared/mod_1$adj.r.squared > 0.98){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }
  }


  make_regressions<- function(dF){

    mod_1<- lm(n_infected~ time + I(time^2), data = dF)

    if(coef(mod_1)[3]>0){
      return(TRUE)}else{
        return(FALSE)}

  }




  if(county == TRUE){
    covid_tibble_US_counties<- covid_tibble_US %>% mutate(data = map(data, day_fun)) %>%
      mutate(num_records = map_int(data, nrow)) %>%
      filter(num_records>5) %>%
      mutate(exp_true = map_lgl(data, make_regressions)) %>%
      filter(exp_true == TRUE) %>%
      unnest(data) %>% ungroup() %>%
      select(region, time, n_infected, fips) %>% mutate(n_infected = log(n_infected))

    region_names<- unique(covid_tibble_US_counties$region)
    region_index<- tibble(region = region_names, id_sim = 1:length(region_names))

    covid_tibble_US_counties<- covid_tibble_US_counties %>% left_join(region_index, by = "region" )
  }

  if(county == FALSE){
    covid_tibble_US_state<- covid_tibble_US %>% unnest(data) %>%
      separate(region, into = c("county_name", "region"), sep = "; ") %>%
      group_by(region, date) %>%
      summarise(n_infected = sum(n_infected), .groups = "drop") %>%
      group_by(region) %>%
      nest() %>% mutate(data= map(data, day_fun)) %>%
      mutate(num_records = map_int(data, nrow)) %>%
      filter(num_records>5) %>%
      mutate(exp_true = map_lgl(data, make_regressions)) %>%
      filter(exp_true == TRUE) %>%
      unnest(data) %>% ungroup() %>%
      mutate(n_infected = log(n_infected))

    region_names<- unique(covid_tibble_US_state$region)
    region_index<- tibble(region = region_names, id_sim = 1:length(region_names))

    covid_tibble_US_state<- covid_tibble_US_state %>% left_join(region_index, by = "region" )

  }

  if(county == TRUE){return(covid_tibble_US_counties)}else{
    return(covid_tibble_US_state)
  }

}



#' A JAGS wrapper for Malthusian parameter estimation based on a reginal random effect.
#'
#' This function creates the MCMC posterior output of JAGS running on the US regional COVID-19 data.
#' The output is based on a linear mixed model of the log transformed count data providing both the fixed effects
#' and the random effects.
#'
#' @param master_data_2 The result of US_covid_data_creator(): county or state regions.
#'
#' @return A list containing a BUGS mcmc object and the original data.
#'
#' @export
run_jag_estimator = function(master_data_2){

  X<- model.matrix(~time, data = master_data_2)
  K<- ncol(X)
  N<- nrow(master_data_2)
  re<- master_data_2$id_sim
  Nre<- length(unique(re))

  meta_data<- list(infected = master_data_2$n_infected,
                   X = X,
                   N = N,
                   K = K,
                   re = re,
                   Nre = Nre)

  sink("lmm.txt")
  cat("
      model{
      #1A Diffuse priors for regression coefs
      for (i in 1:K) {beta[i] ~ dnorm(0,0.0001)}
      #1B Priors for random effects
      for(i in 1:Nre){
      a[i] ~ dnorm(0, tau.ri)
      b[i] ~ dnorm(0, tau.rs)
      }

      #Priors for variances
      num.ri ~ dnorm(0,0.0016)
      num.rs ~ dnorm(0,0.0016)
      num.eps ~ dnorm(0,0.0016)
      denom.ri ~ dnorm(0,1)
      denom.rs ~ dnorm(0,1)
      denom.eps ~ dnorm(0,1)


      sigma.ri<- abs(num.ri/denom.ri)
      sigma.rs<- abs(num.rs/denom.rs)
      sigma.eps<- abs(num.eps/denom.eps)
      #sigma.eps ~ dunif(0.001,10)
      tau.ri<- 1/(sigma.ri*sigma.ri)
      tau.rs<- 1/(sigma.rs*sigma.rs)
      tau.eps<- 1/(sigma.eps*sigma.eps)


      #Likelihood
      for(i in 1:N){
      infected[i] ~ dnorm(mu[i], tau.eps)
      mu[i]<- eta[i]
      eta[i]<- inprod(beta[], X[i,]) + a[re[i]] + b[re[i]]*X[i,2]
      }
      }
      ")
  sink()
  inits<- function(){
    list(beta = rnorm(K, 0, 0.01),
         a =rnorm(Nre,0,.1),
         b = rnorm(Nre,0,.1),
         num.ri = rnorm(1,0,25),
         num.rs = rnorm(1,0,25),
         num.eps = rnorm(1,0,25),
         denom.ri = rnorm(1,0,1),
         denom.rs = rnorm(1,0,1),
         denom.eps = rnorm(1,0,1)
    )
  }


  params<- c("beta", "a", "b", "sigma.ri", "sigma.rs", "sigma.eps")


  malthusian_mcmc<- jags.parallel(data = meta_data,
                                  inits = inits,
                                  parameters.to.save = params,
                                  model.file = "lmm.txt",
                                  n.chains = 10,
                                  n.iter = 5000,
                                  n.thin = 3,
                                  n.burnin = 100
  )


  mcmc_list<- list(malthusian_mcmc, master_data_2)
  return(mcmc_list)

}


#' Creates a tibble for choroplethR function county_choropleth.
#'
#' This function creates the tibble that county_choropleth() requires in the correct format. The value
#' is the R effective for the region.
#'
#' @param mcmc_object The list provided by run_jags_estimator() for county regions.
#'
#' @return A tibble for county_choropleth
#'
#' @export
make_county_choro_data<- function(mcmc_object){

  mcmc_malthusian<- mcmc_object[[1]]
  covid_data<- mcmc_object[[2]]

  gamma_pars<- branchsim::find_gamma_parameters()
  a<- gamma_pars[1]
  b<- gamma_pars[2]

  Alpha<- mcmc_malthusian$BUGSoutput$sims.list$beta[,2]

  #m<- Alpha/(1- (b/(Alpha+ b))^a)


  posterior_results<- mcmc_malthusian$BUGSoutput$sims.matrix %>% as_tibble()

  new_posterior_col_names<- names(posterior_results) %>% str_replace_all("(\\[|\\])","_")

  posterior_results<- setNames(posterior_results,new_posterior_col_names)%>%
    select(starts_with("b")) %>% select(ends_with("_")) %>% select(-beta_1_, -beta_2_)

  posterior_id_sim<- names(posterior_results) %>% str_sub(3,100) %>% str_replace("_","") %>% as.integer()

  posterior_means<- apply(posterior_results,2,mean)

  posterior_lower<- apply(posterior_results + Alpha,2,quantile,probs=0.025)
  posterior_upper<- apply(posterior_results + Alpha,2,quantile,probs=0.975)


  posterior_id_sim_means<- tibble(id_sim = posterior_id_sim, slope = posterior_means + mean(Alpha),
                                  lower = posterior_lower, upper = posterior_upper)

  county_id_sim<-covid_data %>% group_by(region) %>% nest() %>% select(region) %>% ungroup() %>%
    mutate(id_sim = 1:nrow(.))
  county_fips<- covid_data %>% select(region, fips) %>% distinct()
  county_fips_id_sim<-left_join(county_id_sim,county_fips, by = "region")

  #fix NYC
  nyc_id_sim<- county_fips_id_sim %>% filter(str_detect(region, "New York City")) %>% pull(id_sim)

  if(length(nyc_id_sim)>0){
    nyc_fips<- tibble(region = "New York City; New York", id_sim = nyc_id_sim,
                      fips = as.character(c(36061, 36081, 36005, 36047, 36085)))
    county_fips_id_sim<- county_fips_id_sim %>% bind_rows(nyc_fips)
  }
  #data("df_pop_county")

  posterior_choro_data<- left_join(posterior_id_sim_means,county_fips_id_sim, by = "id_sim") %>%
    setNames(c("id_sim", "value", "lower", "upper", "county_name", "region"))

  county_population<- branchestimate::df_pop_county %>% as_tibble() %>% setNames(c("region", "population"))
  posterior_choro_data$region<- as.numeric(posterior_choro_data$region)

  posterior_choro_data<- left_join(county_population,posterior_choro_data, by = "region")


  compute_R<- function(x,A=a,B=b){

    result<- x/(1- (B/(x+ B))^A)*A/B
    return(result)

  }

  posterior_choro_data<- posterior_choro_data %>% mutate(R = compute_R(value), R_lower = compute_R(lower),
                                                         R_upper = compute_R(upper)) %>% na.omit() %>% filter(value>0)
  posterior_choro_data$alpha<- posterior_choro_data$value
  posterior_choro_data$value<- posterior_choro_data$R

  return(posterior_choro_data)



}


#' Creates a tibble for choroplethR function state_choropleth.
#'
#' This function creates the tibble that state_choropleth() requires in the correct format. The value
#' is the R effective for the region
#'
#' @param mcmc_object The list provided by run_jags_estimator() for state regions.
#'
#' @return A tibble for state_choropleth
#'
#' @export
make_state_choro_data<- function(mcmc_object){

  mcmc_malthusian<- mcmc_object[[1]]
  covid_data<- mcmc_object[[2]]

  gamma_pars<- branchsim::find_gamma_parameters()
  a<- gamma_pars[1]
  b<- gamma_pars[2]

  Alpha<- mcmc_malthusian$BUGSoutput$sims.list$beta[,2]

  #m<- Alpha/(1- (b/(Alpha+ b))^a)


  posterior_results<- mcmc_malthusian$BUGSoutput$sims.matrix %>% as_tibble()

  new_posterior_col_names<- names(posterior_results) %>% str_replace_all("(\\[|\\])","_")

  posterior_results<- setNames(posterior_results,new_posterior_col_names)%>%
    select(starts_with("b")) %>% select(ends_with("_")) %>% select(-beta_1_, -beta_2_)

  posterior_id_sim<- names(posterior_results) %>% str_sub(3,100) %>% str_replace("_","") %>% as.integer()

  posterior_means<- apply(posterior_results,2,mean)

  posterior_lower<- apply(posterior_results + Alpha,2,quantile,probs=0.025)
  posterior_upper<- apply(posterior_results + Alpha,2,quantile,probs=0.975)

  posterior_id_sim_means<- tibble(id_sim = posterior_id_sim, slope = posterior_means + mean(Alpha),
                                  lower = posterior_lower, upper = posterior_upper)

  state_id_sim<-covid_data %>% group_by(region) %>% nest() %>% select(region) %>% ungroup() %>%
    mutate(id_sim = 1:nrow(.))
  states_id_sim<- covid_data %>% select(region, id_sim) %>% distinct()


  posterior_choro_data<- left_join(posterior_id_sim_means,states_id_sim, by = "id_sim") %>%
    setNames(c("id_sim", "value", "lower", "upper", "region")) %>% mutate(region = tolower(region))

  #data("df_pop_state")



  state_population<- branchestimate::df_pop_state %>% as_tibble() %>% setNames(c("region", "population"))


  posterior_choro_data<- left_join(state_population,posterior_choro_data, by = "region")


  compute_R<- function(x,A=a,B=b){

    result<- x/(1- (B/(x+ B))^A)*A/B
    return(result)

  }

  posterior_choro_data<- posterior_choro_data %>% mutate(R = compute_R(value), R_lower = compute_R(lower),
                                                         R_upper = compute_R(upper)) %>% na.omit() %>% filter(value>0)

  posterior_choro_data$alpha<- posterior_choro_data$value
  posterior_choro_data$value<- posterior_choro_data$R

  return(posterior_choro_data)

}



#' Renders a state choropleth for R effective.
#'
#' This function creates a state choropleth depending on the output of make_state_choro_data().
#'
#' @param choro_data The tibble provided by make_state_choro_data().
#'
#' @return A state level choropleth for R effective in those states seeing an epidemic outbreak (R_eff>1).
#'
#' @export
make_state_map<- function(choro_data){

  #col.pal<- (brewer.pal(7,"Set2"))
  col.pal<- (brewer.pal(7,"YlOrRd"))

  choro1<- StateChoropleth$new(choro_data)
  choro1$ggplot_scale <- scale_fill_manual(name="R effective\n(NA: sub-exponential growth)",values=col.pal, drop=FALSE)
  choro1$add_state_outline<- TRUE
  choro1$render() + ggtitle("US counties with R effective > 1") +
    theme(plot.title = element_text(hjust = 0.5)) + theme(text = element_text(size=15))

}


#' Renders a county choropleth for R effective.
#'
#' This function creates a county choropleth depending on the output of make_county_choro_data().
#'
#' @param choro_data The tibble provided by make_county_choro_data().
#' @param zoom The state to zoom on. If NA, show the entire US.
#'
#' @return A state level choropleth for R effective in those counties seeing an epidemic outbreak (R_eff>1).
#'
#' @export
make_county_map<- function(choro_data, zoom_state = NA){

    #col.pal<- (brewer.pal(7,"Set2"))
    col.pal<- (brewer.pal(7,"YlOrRd"))

    choro1<- CountyChoropleth$new(choro_data)
    if (!is.na(zoom_state) & zoom_state != "NA") {
        region_label <- str_to_title(zoom_state)
        choro1$set_zoom(zoom_state)
    } else {
        region_label <- "US"
    }
  choro1$ggplot_scale <- scale_fill_manual(name="R effective\n(NA: sub-exponential growth)",values=col.pal, drop=FALSE)
  choro1$add_state_outline<- TRUE
  choro1$render() + ggtitle(str_c(region_label, " counties with R effective > 1")) +
    theme(plot.title = element_text(hjust = 0.5)) + theme(text = element_text(size=15))

}
