var LM = require("../../src/index");
var Matrix = LM.Matrix;
var math = Matrix.algebra;

var A = math.matrix([
  [0, 0, 0],
  [1, 0, 0],
  [2, 2, 0],
  [3, 1, 0],
]);

var B = math.matrix([
  [0, 0, 0],
  [1, 0, 0],
  [2, 2, 0],
  [3, 1, 0],
]);

// equivalent to .* in matlab
var result = math.dotMultiply(A, B);

console.log(result);
