
# f = (x, x, x, 1-3x)

f13.def2 = function(alpha=0.025, beta=0.01, delta=0.005, sigma=0.013, u=c(1,1,1,1), b=0.005/4, set.pow=0.8)   
  
{
  z.alpha = qnorm(1-alpha)
  
  z.beta = qnorm(1-beta)
  
  ## Total number of patients for placebo group to detect overall treatment effect delta 
  ## with type I error=alpha, power=1-beta with 2 (trxt):1(placebo) allocation
  
  N = 3*sigma^2*(z.alpha+z.beta)^2/2/delta^2	
  
  N = round(N)
  
  s = 4
  
  un.pow = function(x, a)   # x=f[1] a=c(s, u, b, delta, sigma, N, set.pow)
  {
    s = a[1]
    u = a[2:(s+1)]
    b = a[s+2]
    delta = a[s+3]
    sigma = a[s+4]
    N = a[s+5]
    set.pow = a[s+6]
    
    mean = u*delta
    
    f = c(x, x, x, 1-3*x)
    
    sgm = 3*sigma^2/N/f/2
    
    # unconditional power - set.pow
    
    prod(pnorm(b, mean = mean, sd = sqrt(sgm), lower.tail = FALSE))-set.pow		
  }
  
  ## minimum required f_1 for uncondition power=set.pow
  
  f1.un = uniroot(un.pow, c(0.05, 1/s), tol=0.0001, a=c(s, u, b, delta, sigma, N, set.pow))[1]
  
  f1.un = round(as.numeric(f1.un), digits=2)
  
  con.pow=function(x, a)	# x=f[1] a=c(s, u, b, delta, sigma, N, z.alpha, set.pow)
  {
    s = a[1]
    u = a[2:(s+1)]
    b = a[s+2]
    delta = a[s+3]
    sigma = a[s+4]
    N = a[s+5]
    z.alpha = a[s+6]
    set.pow = a[s+7]
    
    mean = c(u*delta, delta)
    
    f = c(x, x, x, 1-3*x)
    
    sgm = matrix(rep(0, (s+1)*(s+1)), (s+1), (s+1))
    
    for (i in 1:s) sgm[i,i] = 3*sigma^2/N/f[i]/2
    
    sgm[(s+1),] = 3*sigma^2/N/2
    
    sgm[,(s+1)] = 3*sigma^2/N/2
    
    lower = c(rep(b, s), z.alpha*sqrt(3*sigma^2/N/2))
    
    upper = rep(Inf, (s+1))
    
    # conditional power-set.pow
    
    as.numeric(pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sgm)[1])/(1-pnorm(z.alpha-delta/sqrt(3*sigma^2/N/2)))-set.pow		
  }
  
  ## minimum required f_1 for condition power=set.pow
  
  f1.con = uniroot(con.pow, c(0.05, 1/s), tol=0.0001, a=c(s, u, b, delta, sigma, N, z.alpha, set.pow))[1]
  
  f1.con = round(as.numeric(f1.con), digits=2)
  
  
  return(list(uncon.f1 = f1.un, con.f1 = f1.con))
  
}
