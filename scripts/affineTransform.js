var LM = require("../src/index");
var Matrix = LM.Matrix;
var math = Matrix.algebra;
/**
 * Affine transformation consisting of rotation, scaling and translation for 2D points with 3rd dimension padded with zeros
 * p - Parameters: [xOffset, yOffset, angle, scale]
 */
function affineTransform(t, p) {
  var xOffset = p[0][0];
  var yOffset = p[1][0];
  var angle = p[2][0];
  var scale = p[3][0];

  var transformMatrix = math.matrix([
    [Math.cos(angle) * scale, -Math.sin(angle) * scale, xOffset],
    [Math.sin(angle) * scale, Math.cos(angle) * scale, yOffset],
    [0, 0, 1],
  ]);
  return math.multiply(t, transformMatrix);
}

module.exports = { affineTransform };
