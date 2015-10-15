# Levenberg Marquardt curve-fitting

minimize sum of weighted squared residuals. Javascript version of matlab library from Henri Gavin.

See example for usage

### ---------  INPUT  VARIABLES  -----------
 ``
 func   = function of n independent variables, 't', and m parameters, 'p',
 
 returning the simulated model: y_hat = func(t,p,c)
 
 p      = n-vector of initial guess of parameter values
 
 t      = m-vectors or matrix of independent variables (used as arg to func)
 
 y_dat  = m-vectors or matrix of data to be fit by func(t,p)
 
 weight = weighting vector for least squares fit ( weight >= 0 ) ...
 
 inverse of the standard measurement errors
 
 Default:  sqrt(d.o.f. / ( y_dat' * y_dat ))
 
 dp     = fractional increment of 'p' for numerical derivatives
 
 dp(j)>0 central differences calculated
 
 dp(j)<0 one sided 'backwards' differences calculated
 
 dp(j)=0 sets corresponding partials to zero; i.e. holds p(j) fixed
 
 Default:  0.001;
 
 p_min  = n-vector of lower bounds for parameter values
 
 p_max  = n-vector of upper bounds for parameter values
 
 c      = an optional matrix of values passed to func(t,p,c)
 
 opts   = vector of algorithmic parameters
 
 parameter    defaults    meaning
 
 opts(1)  =  prnt            3        >1 intermediate results; >2 plots
 
 opts(2)  =  MaxIter      10*Npar     maximum number of iterations
 
 opts(3)  =  epsilon_1       1e-3     convergence tolerance for gradient
 
 opts(4)  =  epsilon_2       1e-3     convergence tolerance for parameters
 
 opts(5)  =  epsilon_3       1e-3     convergence tolerance for Chi-square
 
 opts(6)  =  epsilon_4       1e-2     determines acceptance of a L-M step
 
 opts(7)  =  lambda_0        1e-2     initial value of L-M paramter
 
 opts(8)  =  lambda_UP_fac   11       factor for increasing lambda
 
 opts(9)  =  lambda_DN_fac    9       factor for decreasing lambda
 
 opts(10) =  Update_Type      1       1: Levenberg-Marquardt lambda update
 
 2: Quadratic update
 
 3: Nielsen's lambda update equations
 
```
### ----------  OUTPUT  VARIABLES  -----------

 p       = least-squares optimal estimate of the parameter values
 
 X2      = Chi squared criteria
 

 Henri Gavin, Dept. Civil & Environ. Engineering, Duke Univ. 22 Sep 2013
 modified from: [http://octave.sourceforge.net/optim/function/leasqr.html]
 using references by
 
 Press, et al., Numerical Recipes, Cambridge Univ. Press, 1992, Chapter 15.
 
 Sam Roweis      [http://www.cs.toronto.edu/~roweis/notes/lm.pdf]
 
 Manolis Lourakis [http://www.ics.forth.gr/~lourakis/levmar/levmar.pdf]
 
 Hans Nielson     [http://www2.imm.dtu.dk/~hbn/publ/TR9905.ps]
 
 Mathworks        optimization toolbox reference manual
 
 K. Madsen, H.B., Nielsen, and O. Tingleff
 
 [http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf]
