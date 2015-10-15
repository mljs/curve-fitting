/**
 * Created by acastillo on 15/10/15.
 */
/**
 * Trilateration algorithm. https://en.wikipedia.org/wiki/Trilateration
 * @type {*|exports|module.exports}
 */
var Matrix = require("ml-matrix");
var math = require("../src/algebra");
var LM = require("../src/LM");

//Distance function. p is the guessed point.
var euclidean = function(t,p,c){
    var rows = t.rows;
    var result = new Matrix(t.rows, 1);
    for(var i=0;i<rows;i++){
        result[i][0] = Math.sqrt(Math.pow(t[i][0]-p[0][0],2)+Math.pow(t[i][1]-p[1][0],2));
    }

    return result;
};

var data =  [[0.0, 0.0, 10.0], [10.0, 10.0, 10], [10.0, 0.0, 14.142135]];
var nbPoints = data.length;
var t = math.matrix(nbPoints,2);//[1:Npnt]';                              // independent variable

var y_data = math.matrix(nbPoints, 1);
var sum = 0;
var maxY=0;

for(var i=0;i<nbPoints;i++){
    t[i][0] = data[i][0];
    t[i][1] = data[i][1];
    y_data[i][0]=data[i][2];
}
//console.log(nbPoints);
var weight = [1];
//console.log("weight: "+weight);
var opts = [  2,    100, 1e-3, 1e-3, 1e-3, 1e-2, 1e-2, 11, 9, 1 ];
var consts = [];
var p_init = math.matrix([[5],[5]]);
var p_min = math.matrix([[-20],[-20]]);
var p_max = math.matrix([[20],[25]]);
var p_fit = LM.optimize(euclidean,p_init,t,y_data,weight,-0.01,p_min,p_max,consts,opts);
p_fit = p_fit.p;
console.log("Optimus: ");
console.log(p_fit);
console.log("Distance to optimus: ")
console.log(euclidean(t,p_fit,consts));
