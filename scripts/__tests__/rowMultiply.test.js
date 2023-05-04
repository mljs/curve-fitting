var LM = require("../../src/index");
var Matrix = LM.Matrix;
var math = Matrix.algebra;

var A = math.matrix([
  [1, 0, 0],
  [2, 2, 0],
  [3, 1, 0],
]);

var B = math.matrix([[1], [2], [3]]);

var result = math.rowMultiply(A, B);
console.log(result);
