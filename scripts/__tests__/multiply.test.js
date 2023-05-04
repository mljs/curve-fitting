var LM = require("../../src/index");
var Matrix = LM.Matrix;
var math = Matrix.algebra;

var A = math.matrix([[1], [2], [3]]);

var B = math.matrix([[2, 1]]);

// equivalent to .* in matlab
var result = math.multiply(A, B);

console.log(result);
