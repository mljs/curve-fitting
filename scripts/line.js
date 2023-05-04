var LM = require("../src/index");
var Matrix = LM.Matrix;
var math = Matrix.algebra;

function line(t, p) {
  var result = math.matrix(t.rows, 1);
  for (var i = 0; i < t.rows; i++) {
    result[i][0] = p[0][0] * t[i][0] + p[1][0];
  }
  return result;
}

module.exports = line;
