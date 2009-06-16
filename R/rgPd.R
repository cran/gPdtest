# This function generates Generalized Pareto random numbers
rgPd <- function(n,shape,scale=1)
{
	if(shape!=0)
		return((scale/shape)*(runif(n)**(-shape)-1))
	else return(rexp(n,scale))
}
