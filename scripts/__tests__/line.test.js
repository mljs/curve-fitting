var LM = require("../../src/index");
var Matrix = LM.Matrix;
var math = Matrix.algebra;

var line = require("../line.js");

var t = math.matrix([[0], [1], [2], [3], [4], [5]]);
var p = math.matrix([[2], [1]]);

var y = line(t, p);

console.log(y);
// should be [1,3,5,7,9,11]
