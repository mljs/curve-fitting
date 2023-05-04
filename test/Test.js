"use strict";

var LM = require("../src/index");
var Matrix = LM.Matrix;
var math = Matrix.algebra;

describe("Levenberg-Marquardt", function () {
  // Test case obtained from Pag 443, Chap 8.
  it("Convergence with good initial aproximation", function () {
    //Lorentzian function
    //p is a column where p[0][0] is the center of the distribution, and p[1][0] is the parameter specifying the width.
    var lm_func = function (t, p, c) {
      var factor = p[2][0] * Math.pow(p[1][0], 2);
      var rows = t.rows;
      var result = math.matrix(t.rows, t.columns);
      // var tmp = math.add(math.dotPow(math.subtract(t,p[0][0]),2),Math.pow(p[1][0],2));
      for (var i = 0; i < rows; i++) {
        result[i][0] =
          p[3][0] +
          factor / (Math.pow(t[i][0] - p[0][0], 2) + Math.pow(p[1][0], 2));
      }
      return result;
    };

    var data = [
      [3.8436861476749726, 6631317.628571428, 0.00170772075596215],
      [3.847101589186897, 6746573.571428571, 0.002561581133943225],
      [3.8547863325887266, 7012289.457142857, 0.00170772075596215],
      [3.858201774100651, 7137891.571428571, 0.002561581133943225],
      [3.8633249363685374, 7331528.342857143, 0.0034154415119243],
      [3.8675942382584427, 7491647.971428571, 0.00170772075596215],
      [3.871863540148348, 7656913.371428572, 0.00170772075596215],
      [3.8752789816602724, 7794357.685714286, 0.00170772075596215],
      [3.882963725062102, 8115994, 0.0034154415119243],
      [3.8906484684639318, 8444931.057142857, 0.00170772075596215],
      [3.894917770353837, 8636644.514285713, 0.00170772075596215],
      [3.9000409326217236, 8876069.142857142, 0.00170772075596215],
      [3.904310234511629, 9079192.657142857, 0.002561581133943225],
      [3.9077256760235533, 9244712.971428571, 0.00170772075596215],
      [3.9137026986694208, 9544687.028571429, 0.00170772075596215],
      [3.917118140181345, 9719034.485714285, 0.002561581133943225],
      [3.9239490232051937, 10075249.314285714, 0.002561581133943225],
      [3.929926045851061, 10397994.685714286, 0.00170772075596215],
      [3.939318510008853, 10920887.657142857, 0.00170772075596215],
      [3.9452955326547205, 11267839.2, 0.00170772075596215],
      [3.9555418571904934, 11875819.228571428, 0.00170772075596215],
      [3.962372740214342, 12297532.457142856, 0.00170772075596215],
      [3.9666420421042474, 12563071.57142857, 0.0034154415119243],
      [3.9700574836161717, 12773725.2, 0.00170772075596215],
      [3.974326785506077, 13040415.142857142, 0.00170772075596215],
      [3.9768883666400203, 13201829.514285713, 0.00170772075596215],
      [3.98457311004185, 13687291.028571429, 0.002561581133943225],
      [3.987134691175793, 13849282.342857143, 0.00170772075596215],
      [3.994819434577623, 14334897.971428571, 0.002561581133943225],
      [3.9990887364675283, 14602626.885714285, 0.002561581133943225],
      [4.0042118987354165, 14917291.371428572, 0.003415441511926076],
      [4.009335061003306, 15227003.457142856, 0.001707720755963038],
      [4.014458223271195, 15531917.885714285, 0.004269301889907595],
      [4.018727525161102, 15776362.028571429, 0.002561581133944557],
      [4.02299682705101, 16014102.17142857, 0.005977022645870633],
      [4.0255584081849545, 16151113.4, 0.005977022645870633],
      [4.030681570452844, 16413148.542857142, 0.002561581133944557],
      [4.036658593098714, 16703256.257142857, 0.005123162267889114],
      [4.041781755366603, 16931194.057142857, 0.011100184913759747],
      [4.0469049176344924, 17138768.22857143, 0.011100184913759747],
      [4.049466498768437, 17233456.542857144, 0.011100184913759747],
      [4.052881940280363, 17351988.914285712, 0.007684743401833671],
      [4.058005102548252, 17506438.514285713, 0.007684743401833671],
      [4.061420544060178, 17594415.285714284, 0.019638788693574938],
      [4.063982125194123, 17655829.42857143, 0.019638788693574938],
      [4.066543706328067, 17709095.4, 0.019638788693574938],
      [4.070813008217975, 17781712.885714285, 0.019638788693574938],
      [4.0733745893519195, 17815329.314285714, 0.019638788693574938],
      [4.077643891241827, 17852089.171428572, 0.019638788693574938],
      [4.081913193131735, 17869869.028571427, 0.04440073965503899],
      [4.087036355399624, 17862016.942857143, 0.04440073965503899],
      [4.092159517667513, 17825243.42857143, 0.04440073965503899],
      [4.095574959179439, 17782817.4, 0.04440073965503899],
      [4.100698121447328, 17690851.8, 0.04440073965503899],
      [4.104967423337236, 17592028.085714284, 0.04440073965503899],
      [4.109236725227143, 17473718.114285715, 0.04440073965503899],
      [4.112652166739069, 17367166.62857143, 0.04440073965503899],
      [4.1160676082509955, 17246314.77142857, 0.04440073965503899],
      [4.11862918938494, 17148483.657142855, 0.04440073965503899],
      [4.122898491274848, 16972286.314285714, 0.04440073965503899],
      [4.128021653542737, 16735571.285714285, 0.00853860377981519],
      [4.130583234676681, 16609499.799999999, 0.00853860377981519],
      [4.133998676188607, 16432863.142857142, 0.00853860377981519],
      [4.1374141177005335, 16247389.2, 0.005977022645870633],
      [4.1408295592124595, 16057635.657142857, 0.005977022645870633],
      [4.145952721480349, 15756797.714285715, 0.004269301889907595],
      [4.151929744126219, 15388700.17142857, 0.003415441511926076],
      [4.156199046016127, 15112549.228571428, 0.003415441511926076],
      [4.161322208284016, 14775590.828571428, 0.002561581133944557],
      [4.164737649795942, 14548151.542857142, 0.002561581133944557],
      [4.168153091307868, 14314625.914285714, 0.002561581133944557],
      [4.173276253575757, 13959890.971428571, 0.001707720755963038],
      [4.179253276221628, 13544738, 0.001707720755963038],
      [4.1835225781115355, 13246886.799999999, 0.001707720755963038],
      [4.187791880001443, 12949848.6, 0.002561581133944557],
      [4.194622763025295, 12475132, 0.001707720755963038],
      [4.198892064915203, 12180167.342857143, 0.001707720755963038],
      [4.204015227183092, 11830441.799999999, 0.001707720755963038],
      [4.21426155171887, 11146941.77142857, 0.002561581133944557],
      [4.221946295120704, 10647255, 0.001707720755963038],
      [4.2296310385225375, 10164768.42857143, 0.001707720755963038],
      [4.23646192154639, 9750976.085714286, 0.001707720755963038],
      [4.240731223436297, 9501366.628571428, 0.002561581133944557],
      [4.2509775479720755, 8918881.285714285, 0.001707720755963038],
      [4.254392989484002, 8731254.799999999, 0.001707720755963038],
      [4.258662291373909, 8503638.057142857, 0.001707720755963038],
      [4.273177917799595, 7770267.342857143, 0.001707720755963038],
      [4.277447219689503, 7568759.085714285, 0.001707720755963038],
      [4.285985823469318, 7180737.257142857, 0.003415441511926076],
      [4.2936705668711515, 6850695.657142857, 0.001707720755963038],
      [4.297939868761059, 6675134.171428571, 0.002561581133944557],
    ];

    var nbPoints = data.length;
    var t = math.matrix(nbPoints, 1); //[1:Npnt]';				  // independent variable
    var y_data = math.matrix(nbPoints, 1);
    var sum = 0;
    t[0][0] = data[0][0];
    y_data[0][0] = data[0][1];
    var maxY = 0;
    for (var i = 0; i < nbPoints; i++) {
      t[i][0] = data[i][0];
      y_data[i][0] = data[i][1];
      if (data[i][1] > maxY) maxY = data[i][1];
    }
    for (var i = 0; i < nbPoints; i++) {
      y_data[i][0] /= maxY;
    }
    var weight = [nbPoints / math.sqrt(y_data.dot(y_data))];
    var opts = [3, 100, 1e-3, 1e-3, 1e-3, 1e-2, 1e-2, 11, 9, 1];
    var consts = []; // optional vector of constants
    var p_init = math.matrix([
      [(t[0][0] + t[nbPoints - 1][0]) / 2],
      [Math.abs(t[0][0] - t[nbPoints - 1][0]) / 2],
      [1],
      [0],
    ]);
    var p_min = math.matrix([[t[0][0]], [0.0], [0], [0]]);
    var p_max = math.matrix([
      [t[nbPoints - 1][0]],
      [Math.abs(t[0][0] - t[nbPoints - 1][0])],
      [1.5],
      [0.5],
    ]);
    // console.log({ t, y_data, p_init, weight });

    var p_fit = LM.optimize(
      lm_func,
      p_init,
      t,
      y_data,
      weight,
      -0.01,
      p_min,
      p_max,
      consts,
      opts
    );
    p_fit = p_fit.p;
    p_fit[0][0].should.approximatelyDeep(4.08, 0.01); //Center
    p_fit[1][0].should.approximatelyDeep(0.176, 0.01); //hwhh
    p_fit[2][0].should.approximatelyDeep(1, 0.01); //height
    p_fit[3][0].should.approximatelyDeep(0, 0.01); //constant
  });

  it("linear regression", function () {
    function line(t, p) {
      var result = math.matrix(t.rows, 1);

      for (var i = 0; i < t.rows; t++) {
        result[i][0] = p[0] * t[i][0] + p[1];
      }
      return result;
    }

    var t = math.matrix([[0], [1], [2], [3], [4], [5]]);
    var y = math.matrix([[1], [3], [5], [7], [9], [11]]);

    var p_init = math.matrix([[0], [0]]);

    var p_min = math.matrix([[-5], [-5]]);
    var p_max = math.matrix([[5], [5]]);
    var weight = [1];

    console.log({ t, y, p_init });
    var fit = LM.optimize(line, p_init, t, y, weight, 0.001, p_min, p_max);

    console.log(fit);
  });
});
