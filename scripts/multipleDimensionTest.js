// curve fitting with 2 dimensions data
// we want to find the transformation matrix from source to destination points
var { affineTransform, getTransformMatrix } = require("./affineTransform.js");

var LM = require("../src/index");
var Matrix = LM.Matrix;
var math = Matrix.algebra;

var t = math.matrix([
  [1, 5, 1],
  [2, 6, 1],
  [3, 7, 1],
  [4, 8, 1],
]);

var y = math.matrix([
  [1, 5, 1],
  [2, 6, 1],
  [3, 7, 1],
  [4, 8, 1],
]);

var p_init = math.matrix([[0], [0], [0], [2]]);

var p_min = math.matrix([[0], [0], [0], [0]]);
var p_max = math.matrix([[0], [0], [0], [0]]);

// weight of each data point
var weight = [1];

var dp = -0.01;

console.log({ t, y, p_init });
var fit = LM.optimize(affineTransform, p_init, t, y, weight, dp, p_min, p_max);

console.log(fit.p);
