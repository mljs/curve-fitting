// curve fitting with 2 dimensions data
// we want to find the transformation matrix from source to destination points
var { affineTransform, getTransformMatrix } = require("./affineTransform.js");

var LM = require("../src/index");
var Matrix = LM.Matrix;
var math = Matrix.algebra;

var t = math.matrix([
  [0, 0, 0],
  [1, 0, 0],
  [2, 2, 0],
  [3, 1, 0],
]);

var y = math.matrix([
  [0, 0, 0],
  [1, 0, 0],
  [2, 2, 0],
  [3, 1, 0],
]);

var p_init = getTransformMatrix([0, 0, 0, 0]);
console.log({ p_init });
console.log({ t, y, p_init });
var fit = LM.optimize(affineTransform, p_init, t, y);

console.log(fit.p);
