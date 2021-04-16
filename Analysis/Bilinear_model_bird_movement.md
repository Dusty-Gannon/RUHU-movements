Sender-receiver Poisson model fit
================
D. G. Gannon
October 2020

Let ![y\_{i,j,k} \\in
\\mathbb{N}^+](https://latex.codecogs.com/png.latex?y_%7Bi%2Cj%2Ck%7D%20%5Cin%20%5Cmathbb%7BN%7D%5E%2B
"y_{i,j,k} \\in \\mathbb{N}^+") be the number of movements detected
between readers ![i](https://latex.codecogs.com/png.latex?i "i") and
![j](https://latex.codecogs.com/png.latex?j "j") by hummingbirds in year
![k](https://latex.codecogs.com/png.latex?k "k"), where ![i \\ne j
= 1,2,...,R](https://latex.codecogs.com/png.latex?i%20%5Cne%20j%20%3D%201%2C2%2C...%2CR
"i \\ne j = 1,2,...,R"). Following Hoff (2005), we model the number of
movements in a bilinear model:

  
![&#10; \\log (\\lambda\_{i,j,k}) = {\\bf x}\_{ij}\\boldsymbol \\beta +
u\_i + w\_j + \\gamma\_{ij} + \\eta\_k +\\log(b\_k) +
\\log(d\_k)&#10;](https://latex.codecogs.com/png.latex?%0A%20%5Clog%20%28%5Clambda_%7Bi%2Cj%2Ck%7D%29%20%3D%20%20%7B%5Cbf%20x%7D_%7Bij%7D%5Cboldsymbol%20%5Cbeta%20%2B%20u_i%20%2B%20w_j%20%2B%20%5Cgamma_%7Bij%7D%20%2B%20%5Ceta_k%20%2B%5Clog%28b_k%29%20%2B%20%5Clog%28d_k%29%0A
"
 \\log (\\lambda_{i,j,k}) =  {\\bf x}_{ij}\\boldsymbol \\beta + u_i + w_j + \\gamma_{ij} + \\eta_k +\\log(b_k) + \\log(d_k)
")  
where ![{\\boldsymbol
\\beta}](https://latex.codecogs.com/png.latex?%7B%5Cboldsymbol%20%5Cbeta%7D
"{\\boldsymbol \\beta}") is a vector of dyad specific effects with
![{\\bf
x}\_{i,j}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20x%7D_%7Bi%2Cj%7D
"{\\bf x}_{i,j}") the vector of dyad-specific regressors for dyad
![\\{i,j\\}](https://latex.codecogs.com/png.latex?%5C%7Bi%2Cj%5C%7D
"\\{i,j\\}"). The effects
![u\_i](https://latex.codecogs.com/png.latex?u_i "u_i") and
![w\_i](https://latex.codecogs.com/png.latex?w_i "w_i") are the average
effects of reader ![i](https://latex.codecogs.com/png.latex?i "i") as a
sender and receiver (respectively), and the term
![\\gamma\_{i,j}](https://latex.codecogs.com/png.latex?%5Cgamma_%7Bi%2Cj%7D
"\\gamma_{i,j}") is an average effect on movement for the pair of
readers
![\\{i,j\\}](https://latex.codecogs.com/png.latex?%5C%7Bi%2Cj%5C%7D
"\\{i,j\\}"). Finally, ![b\_k](https://latex.codecogs.com/png.latex?b_k
"b_k") and ![d\_k](https://latex.codecogs.com/png.latex?d_k "d_k") are
offsets for the cumulative number of birds that had been implanted with
RFID tags in year ![k](https://latex.codecogs.com/png.latex?k "k") and
the number of days the readers were maintained in year
![k](https://latex.codecogs.com/png.latex?k "k") (respectively).

To induce dependence among observations that involve reader
![i](https://latex.codecogs.com/png.latex?i "i"), we assume

  
![&#10;\\begin{aligned}&#10;\\begin{bmatrix}&#10;u\_i
\\\\&#10;w\_i&#10;\\end{bmatrix} &\\sim \\mathcal{N}\\left({\\bf 0},\\
{\\boldsymbol \\Sigma}\_{uw} \\right),\\\\&#10;{\\boldsymbol
\\Sigma}\_{uw} &= \\begin{bmatrix}&#10;\\sigma\_{u}^2 &
\\sigma\_{uw}\\\\&#10;\\sigma\_{uw} &
\\sigma\_w^2&#10;\\end{bmatrix},&#10;\\end{split}&#10;](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0A%5Cbegin%7Bbmatrix%7D%0Au_i%20%5C%5C%0Aw_i%0A%5Cend%7Bbmatrix%7D%20%26%5Csim%20%5Cmathcal%7BN%7D%5Cleft%28%7B%5Cbf%200%7D%2C%5C%20%7B%5Cboldsymbol%20%5CSigma%7D_%7Buw%7D%20%5Cright%29%2C%5C%5C%0A%7B%5Cboldsymbol%20%5CSigma%7D_%7Buw%7D%20%26%3D%20%5Cbegin%7Bbmatrix%7D%0A%5Csigma_%7Bu%7D%5E2%20%26%20%5Csigma_%7Buw%7D%5C%5C%0A%5Csigma_%7Buw%7D%20%26%20%5Csigma_w%5E2%0A%5Cend%7Bbmatrix%7D%2C%0A%5Cend%7Bsplit%7D%0A
"
\\begin{aligned}
\\begin{bmatrix}
u_i \\\\
w_i
\\end{bmatrix} &\\sim \\mathcal{N}\\left({\\bf 0},\\ {\\boldsymbol \\Sigma}_{uw} \\right),\\\\
{\\boldsymbol \\Sigma}_{uw} &= \\begin{bmatrix}
\\sigma_{u}^2 & \\sigma_{uw}\\\\
\\sigma_{uw} & \\sigma_w^2
\\end{bmatrix},
\\end{split}
")  

such that observations that have a common sender or receiver may be
correlated. Importantly, this allows for the potential that good senders
may be poor receivers, so negative correlation is a possibility.
Finally, let

  
![&#10;\\begin{aligned}&#10;\\begin{bmatrix}&#10;\\gamma\_{i,j}\\\\&#10;\\gamma\_{j,i}&#10;\\end{bmatrix}
&\\sim \\mathcal{N}({\\bf 0},\\ {\\boldsymbol
\\Sigma}\_\\gamma),\\\\&#10;{\\boldsymbol \\Sigma}\_\\gamma &=
\\begin{bmatrix}&#10;\\sigma\_\\gamma^2 &
\\rho\\sigma\_\\gamma^2\\\\&#10;\\rho\\sigma\_\\gamma^2 &
\\sigma\_\\gamma^2&#10;\\end{bmatrix},&#10;\\end{aligned}&#10;](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0A%5Cbegin%7Bbmatrix%7D%0A%5Cgamma_%7Bi%2Cj%7D%5C%5C%0A%5Cgamma_%7Bj%2Ci%7D%0A%5Cend%7Bbmatrix%7D%20%26%5Csim%20%5Cmathcal%7BN%7D%28%7B%5Cbf%200%7D%2C%5C%20%7B%5Cboldsymbol%20%5CSigma%7D_%5Cgamma%29%2C%5C%5C%0A%7B%5Cboldsymbol%20%5CSigma%7D_%5Cgamma%20%26%3D%20%5Cbegin%7Bbmatrix%7D%0A%5Csigma_%5Cgamma%5E2%20%26%20%5Crho%5Csigma_%5Cgamma%5E2%5C%5C%0A%5Crho%5Csigma_%5Cgamma%5E2%20%26%20%5Csigma_%5Cgamma%5E2%0A%5Cend%7Bbmatrix%7D%2C%0A%5Cend%7Baligned%7D%0A
"
\\begin{aligned}
\\begin{bmatrix}
\\gamma_{i,j}\\\\
\\gamma_{j,i}
\\end{bmatrix} &\\sim \\mathcal{N}({\\bf 0},\\ {\\boldsymbol \\Sigma}_\\gamma),\\\\
{\\boldsymbol \\Sigma}_\\gamma &= \\begin{bmatrix}
\\sigma_\\gamma^2 & \\rho\\sigma_\\gamma^2\\\\
\\rho\\sigma_\\gamma^2 & \\sigma_\\gamma^2
\\end{bmatrix},
\\end{aligned}
")  

where ![\\rho](https://latex.codecogs.com/png.latex?%5Crho "\\rho") is
the correlation between the rate of movement from ![i\\to
j](https://latex.codecogs.com/png.latex?i%5Cto%20j "i\\to j") and ![j
\\to i](https://latex.codecogs.com/png.latex?j%20%5Cto%20i "j \\to i").
These differ from standard random effects models because they allow for
negative correlation of the observations within a dyad.

## Stan Model Code

``` stan

data{

  int<lower=0> N;             //number of observations (edges*years)
  int<lower=0> R;             //number of nodes (readers)
  int<lower=0> D;             //number of dyads
  int<lower=0> P;             //number of fixed effects
  
  int<lower=0,upper=R> r1[N];     //reader 1 index
  int<lower=0,upper=R> r2[N];     //reader 2 index
  int<lower=0,upper=D> dyad[N];   //index for the dyad
  int<lower=0> order[N];          //order for the dyad
  
  matrix[N,P] X;                  //fixed effects model matrix
  real<lower=0> os_d[N];          //number of days connections were active
  real<lower=0> os_b[N];          //cumulative number of birds with RFID tags
  int<lower=0> y[N];              //number of movements

}

parameters{

  vector[P] beta;                //regression parameters
  vector[2] u_raw[R];            //random effects for sender/receiver
  vector[2] gamma_raw[D];        //random effects for the dyad

  vector<lower=0>[2] sigma_u;            //scale for sender/receiver effects
  cholesky_factor_corr[2] L_omega_u;
  
  real<lower=0> sigma_gamma;        //scale for dyad effects
  cholesky_factor_corr[2] L_omega_gamma;
  

}

transformed parameters{

  vector[2] u[R];                
  vector[2] gamma[D];
  real<lower=0> lambda[N];
  
  for(r in 1:R){
    u[r] = diag_pre_multiply(sigma_u, L_omega_u)*u_raw[r];
  }
  
  for(d in 1:D){
    gamma[d] = diag_pre_multiply(rep_vector(sigma_gamma,2),
                  L_omega_gamma)*gamma_raw[d];
  }
  
  //create linear predictor
  for(i in 1:N){
    lambda[i] = exp(X[i,]*beta + u[r1[i],1] + u[r2[i],2] +
                    gamma[dyad[i],order[i]] 
                    + log(os_d[i]) + log(os_b[i]));
  }
  
}

model{

  //priors
  beta ~ normal(0,1);
  sigma_u ~ normal(0,2);
  sigma_gamma ~ normal(0,2);
  L_omega_u ~ lkj_corr_cholesky(5);
  L_omega_gamma ~ lkj_corr_cholesky(5);
  
  
  //model
  for(r in 1:R){
    u_raw[r] ~ std_normal();
  }

  for(d in 1:D){
    gamma_raw[d] ~ std_normal(); 
  }
  
  
  for(i in 1:N){
    y[i] ~ poisson(lambda[i]);
  }


}

generated quantities{

  vector[N] log_lik;
  vector[N] y_hat;
  matrix[2,2] omega_u;
  matrix[2,2] omega_gamma;
  
  omega_u = L_omega_u*L_omega_u';
  omega_gamma = L_omega_gamma*L_omega_gamma';
  
  for(i in 1:N){
    log_lik[i] = poisson_lpmf(y[i] | lambda[i]);
    y_hat[i] = poisson_rng(lambda[i]);
  }

}
```

### Load data

``` r
# number of movements for each year
  mvmnts <- read.table(file = here("Data","total_movements_by_year.tsv"), 
                       header = T, as.is = T, sep = "\t")
  
# reader locations
  rlocs <- read.table(file = here("Data", "hja_reader_locations_centercomb.txt"),
                      header = T, sep = "\t")
  rdrs <- rlocs$reader
  
# bird capture data
  caps <- read_csv(file = here("Data", "hja_hummingbird_captures_rfid.csv"))
```

    ## Parsed with column specification:
    ## cols(
    ##   bdr = col_character(),
    ##   loc = col_character(),
    ##   yr = col_double(),
    ##   mo = col_double(),
    ##   day = col_double(),
    ##   time = col_time(format = ""),
    ##   status = col_character(),
    ##   band = col_character(),
    ##   species = col_character(),
    ##   sex = col_character(),
    ##   age = col_character(),
    ##   rfid = col_character()
    ## )

``` r
  birds_yr <- group_by(caps, yr) %>%
    summarise(n=n())
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
  birds_yr$cum_birds <- cumsum(birds_yr$n)
  
# merge back with movement data
  mvmnts <- merge(mvmnts, birds_yr, by.x="year", by.y="yr")

# create a column for the order
  #if the position for reader 1 > reader 2, order=1
 mvmnts$order <- NA
 for(i in 1:nrow(mvmnts)){
   r1 <- which(rdrs == mvmnts$reader_i[i])
   r2 <- which(rdrs == mvmnts$reader_j[i])
   if(r1 < r2){
     mvmnts$order[i] <- 1
   } else{
     mvmnts$order[i] <- 2
   }
 }
 
# indicator variables for whether a reader is inside the forest
 mvmnts$forest_i <- 0
 mvmnts$forest_i[which(mvmnts$hab_i == "forest")] <- 1
 
 mvmnts$forest_j <- 0
 mvmnts$forest_j[which(mvmnts$hab_j == "forest")] <- 1
 
 mvmnts$forest_ij <- with(mvmnts,
                          (forest_i + forest_j)/2)
```

# Compile data

``` r
  N <- nrow(mvmnts)
  y <- mvmnts$count
  X <- model.matrix(~ dist_ij + pfor_ij + forest_ij, data = mvmnts)
  
  P <- dim(X)[2]
  ord <- mvmnts$order
  dyad <- as.integer(as.factor(mvmnts$dyad))
  D <- length(unique(dyad))
  
  r1 <- as.integer(factor(mvmnts$reader_i, levels = rdrs))
  r2 <- as.integer(factor(mvmnts$reader_j, levels = rdrs))
  R <- length(unique(r1))
  os_d <- (mvmnts$active_days)/7
  os_b <- mvmnts$cum_birds
  
# compile data into list
  mod.data <- list(N=N, P=P, R=R, D=D, r1=r1,
                   r2=r2, dyad=dyad, order=ord,
                   X=X, os_d=os_d, os_b=os_b,
                   y=y)
```

# Fitting the model

``` r
# fit the model
  fit <- sampling(soc_net_poisson, data=mod.data, chains=2, 
                  control=list(adapt_delta=0.9), init=0)
```

<div id="refs" class="references">

<div id="ref-hoff2005b">

Hoff, Peter D. 2005. “Bilinear Mixed-Effects Models for Dyadic Data.”
*Journal of the American Statistical Association* 100 (469): 286–95.

</div>

</div>
