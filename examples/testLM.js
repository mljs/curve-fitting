/**
 * Created by acastillo on 8/6/15.
 */
var Matrix = require("ml-matrix");
var math = require("../src/algebra");
var LM = require("../src/LM");

//p is a column
var lm_func = function(t,p,c){
    //console.log(t);
    //y_hat = p(1)*exp(-t/p(2)) + p(3)*t.*exp(-t/p(4));
    var value = math.add(math.multiply(math.exp(math.multiply(t,-1/p[1][0])),p[0][0]), math.multiply(p[2][0],math.dotMultiply(t,math.exp(math.multiply(t,-1/p[3][0])))));
    return value;
};
// Levenberg-Marquardt example/test

//   Henri Gavin, Dept. Civil & Environ. Engineering, Duke Univ. November 2005
//   modified from: ftp://fly.cnuce.cnr.it/pub/software/octave/leasqr/
    //   Press, et al., Numerical Recipes, Cambridge Univ. Press, 1992, Chapter 15.

//randn('seed',0);	// specify a particular random sequence for msmnt error

//epsPlots = 0;  formatPlot(epsPlots);		// 1: make .eps plots, 0: don't

// *** For this demonstration example, simulate some artificial measurements by
// *** adding random errors to the curve-fit equation.
var example_number = 1;			  // which example to run.

var consts = [ ];   // optional vector of constants

var Npnt = 100;				  // number of data points

var t = math.matrix(Npnt,1);//[1:Npnt]';				  // independent variable

for(var i=0;i<Npnt;i++)
    t[i][0]=i+1;
//t = math.transpose(math.matrix(t));
// true value of parameters ...
var p_true = [];
if(example_number == 1) p_true  = math.matrix([ [20],   [10],   [1],  [50] ]);
if(example_number == 2) p_true  = math.matrix([ [20],  [-24],  [30], [-40] ]);
if(example_number == 3) p_true  = math.matrix([  [6],   [20],   [1],   [5] ]);

var y_dat = lm_func(t, p_true, consts);

//console.log("here "+math.size(y_dat)+" "+math.size(t)+" "+math.size(math.random([Npnt,1], 0, 0.1)));
var y_dat = math.add(y_dat, math.multiply(math.random(Npnt,1),0.1));//0.1*Math.randn(Npnt,1);	  // add random noise

// range of values for basic paramter search
/*var p1 = 0.1*p_true(1):0.2*p_true(1):2*p_true(1);
var p2 = 0.1*p_true(2):0.2*p_true(2):2*p_true(2);
var p3 = 0.1*p_true(3):0.2*p_true(3):2*p_true(3);
var p4 = 0.1*p_true(4):0.2*p_true(4):2*p_true(4);*/
var p1 = new Array();
var p2 = new Array();
var p3 = new Array();
var p4 = new Array();
for(var i=0.1;i<=2;i+=0.2){
    p1.push([p_true[0][0]*i]);
    p2.push([p_true[1][0]*i]);
    p3.push([p_true[2][0]*i]);
    p4.push([p_true[3][0]*i]);
}

// parameter search
var X2 = new Array(p2.length)
for(var ip2 = 0;ip2<p2.length;ip2++){
    X2[ip2]= new Array(p4.length);
    for(var ip4 = 0;ip4<p4.length;ip4++){
        pt = [ p_true[1],  p2[ip2], p_true[3], p4[ip4] ];
        delta_y = math.subtract(y_dat,lm_func(t,pt,consts));
        X2[ip2][ip4] = math.multiply(math.multiply(math.transpose(delta_y), delta_y),1/2);
    }
}
/*figure(1); // ------------ plot shape of Chi-squared objective function
clf
mesh(p2,p4,log10(X2))
xlabel('p_2')
ylabel('p_4')
zlabel('log_{10}(\chi^2)')
plotfile = ['lm_exampA',int2str(example_number),'.eps'];
//if epsPlots, print(plotfile,'-solid','-color','-deps','-F:28'); end
*/
// *** Replace the lines above with a read-in of some
// *** experimentally-measured data.

// initial guess parameters  ...
var p_init = [];
if(example_number == 1) p_init  = math.matrix([[5], [2], [0.2], [10] ]);
if(example_number == 2) p_init  = math.matrix([[4], [-5], [6], [10] ]);
if(example_number == 3) p_init  =math.matrix([ [10], [50], [5], [5.6] ]);

//weight = Npnt/sqrt(y_dat'*y_dat);	  // sqrt of sum of data squared
var weight = [Npnt / math.sqrt(math.multiply(math.transpose(y_dat),y_dat))];


var p_min = math.multiply(math.abs(p_init),-10);
var p_max = math.multiply(math.abs(p_init),10);


// Algorithmic Parameters
//         prnt MaxIter  eps1  eps2  epx3  eps4  lam0  lamUP lamDN UpdateType
var opts = [  3,    100, 1e-3, 1e-3, 1e-3, 1e-2, 1e-2,    11,    9,        1 ];

//[p_fit,Chi_sq,sigma_p,sigma_y,corr,R2,cvg_hst] =  lm('lm_func',p_init,t,y_dat,weight,-0.01,p_min,p_max,consts,opts);
var iterations = 0;
console.log(LM);
var p_fit = LM.optimize(lm_func,p_init,t,y_dat,weight,-0.01,p_min,p_max,consts,opts);
p_fit = p_fit.p;
var y_fit = lm_func(t,p_fit,consts);

console.log('    initial    true       fit');
console.log(' -----------------------------------');
console.log([ p_init,  p_true,  p_fit])//, sigma_p, 100*Math.abs(sigma_p./p_fit) ]);

