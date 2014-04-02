#Function to calculate Proportion of Japanese Patients according to 

#Sample size and proportion of Japanese patients in multi-regional trials
#Kimitoshi Ikeda and Frank Bretz
#Pharmaceutical Statistics 9: 207–216 (2010)

#to ensure delta_Japan / delta_All > omega
#and delta_All > 0

prop.jap1 <- function(omega=0.5,beta=0.8,alpha=0.025,gamma=0.8){
  Zalpha <- qnorm(p=1-alpha)
  Zbeta  <- qnorm(p=1-beta)
  Zgamma <- qnorm(p=1-gamma)
  p <- Zgamma^2/((1-omega)^2 * (Zalpha+Zbeta)^2 - (Zgamma^2*omega*(omega-2)))
  return(list(Proportion=p))
}



