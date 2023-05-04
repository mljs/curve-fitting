var LM = require("../../src/index");
var Matrix = LM.Matrix;
var math = Matrix.algebra;

var { affineTransform } = require("../affineTransform.js");

// padding third dimension with zeros
var t = math.matrix([
  [0, 0, 0],
  [1, 0, 0],
  [2, 2, 0],
  [3, 1, 0],
]);

var p = math.matrix([[0], [0], [0], [2]]);

var y = affineTransform(t, p);

console.log(y);
// should be
// [
//     [0, 0, 0],
//     [2, 0, 0],
//     [4, 4, 0],
//     [6, 2, 0],
//   ]
