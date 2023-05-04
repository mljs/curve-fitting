var LM = require("../src/index");
var Matrix = LM.Matrix;
var math = Matrix.algebra;
// affine transformation consisting of rotation, scaling and translation for 2D points with 3rd dimension padded with zeros
function getTransformMatrix([xOffset, yOffset, angle, scale]) {
  return math.matrix([
    [Math.cos(angle) * scale, -Math.sin(angle) * scale, xOffset],
    [Math.sin(angle) * scale, Math.cos(angle) * scale, yOffset],
    [0, 0, 1],
  ]);
}

function affineTransform(t, matrix) {
  return math.multiply(t, matrix);
}

module.exports = { affineTransform, getTransformMatrix };
