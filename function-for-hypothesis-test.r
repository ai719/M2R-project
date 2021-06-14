library(spatstat)

### CSR function ###
csr <- function(xMin, xMax, yMin, yMax, lambda) {
  xDelta=xMax-xMin;yDelta=yMax-yMin; # rectangle dimensions
  areaTotal=xDelta*yDelta; # study area
  
  numbPoints=rpois(1,areaTotal*lambda);#Poisson number of points
  xx=xDelta*runif(numbPoints)+xMin; # x coordinates of Poisson points
  yy=yDelta*runif(numbPoints)+yMin; # y coordinates of Poisson points
  df = data.frame(x=xx, y=yy) # dataframe of x and y points
  return(df)
}


### CSR test ###
testing_csr <- function(x,y,nsim,alpha) 
  {
  # min and max x and y valuees
  xMin = min(x); xMax = max(x);
  yMin = min(y); yMax = max(y);
  # area of the study area
  A = (xMax - xMin) * (yMax - yMin);
  # number of events
  N = length(x) ;
  # estimated intensity for for csr
  lambda = N/A ;
  # pp of our data we want to test
  data_pp = ppp(x, y, c(xMin, xMax), c(yMin, yMax)) ;
  # l function on data
  data_l_func <- Lest(data_pp, correction = "Ripley") ;
  # test statistic for data
  t = max(abs(data_l_func[["iso"]]-data_l_func[["r"]])) ;
  # vector of t_i for the nsim csr s
  t_sim_vals = numeric(nsim) ;
  for (i in 1:nsim) {
    # calling csr function
    csr_i = csr(xMin, xMax, yMin, yMax, lambda)
    # pp for simulation
    i_hpp = ppp(csr_i$x, csr_i$y, c(xMin, xMax), c(yMin, yMax))
    # l function of simulation
    l_i = Lest(i_hpp, correction = "Ripley")
    # t_i statistic assigned to i th value of vector
    t_sim_vals[i] = max(abs(l_i[["iso"]]-l_i[["r"]]))
  } ;
  # sort so in order
  sort(t_sim_vals);
  # as 2 tailed test
  percentile = floor(nsim*(1 - alpha/2)) + 1 ;
  c_alpha_over_2 = t_sim_vals[percentile] ;
  # output if in rekection region or not
  if(c_alpha_over_2<t){
    print("Reject H_0")
  } else {
    print("Failed to Reject H_0")
  }
  # a data from of thetest statistic for data & c_a/2 
  output = data.frame("Test Statistic" = t,
                      c_alpha_over_2 = c_alpha_over_2)
  return(output)
}
