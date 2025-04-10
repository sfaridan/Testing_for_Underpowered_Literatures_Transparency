sigma_Y <- 1
c <- sqrt(2)

q <- log( (1+sigma_Y^2)/(1+sigma_Y^2-c^(-2)) )/log((1+sigma_Y^2)/sigma_Y^2)
eta <- (1+sigma_Y^2-c^(-2))/(1+sigma_Y^2)
lambda <- (sigma_Y^2)/(1+sigma_Y^2-c^(-2))

#lambda^{2\beta} =eta^2
#lambda^{beta} = eta
#beta =log_{lambda}(eta)
#beta = log(eta)/log(lambda)
beta <- log(eta)/log(lambda)
wedge <- min(c(beta,2))

q # my MSE
wedge/(wedge+2)  #carrasco MISE
