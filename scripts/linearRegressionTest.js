// simple linear regression
var line = require("./line.js");

var LM = require("../src/index");
var Matrix = LM.Matrix;
var math = Matrix.algebra;

var t = math.matrix([[0], [1], [2], [3], [4], [5]]);
var y = math.matrix([[1], [3], [5], [7], [9], [11]]);
var y2 = math.matrix([[1], [4], [7], [10], [13], [16]]);

var p_init = math.matrix([[2], [0]]);

var p_min = math.matrix([[-5], [-5]]);
var p_max = math.matrix([[5], [5]]);

var weight = [1];
var dp = -0.01;

console.log({ t, y, p_init });
var fit = LM.optimize(line, p_init, t, y, weight, dp, p_min, p_max);

// p should be 2, 1
console.log("LINEAR REGRESSION");
console.log(fit.p);
