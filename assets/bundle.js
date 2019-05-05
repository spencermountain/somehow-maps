(function(){function r(e,n,t){function o(i,f){if(!n[i]){if(!e[i]){var c="function"==typeof require&&require;if(!f&&c)return c(i,!0);if(u)return u(i,!0);var a=new Error("Cannot find module '"+i+"'");throw a.code="MODULE_NOT_FOUND",a}var p=n[i]={exports:{}};e[i][0].call(p.exports,function(r){var n=e[i][1][r];return o(n||r)},p,p.exports,r,e,n,t)}return n[i].exports}for(var u="function"==typeof require&&require,i=0;i<t.length;i++)o(t[i]);return o}return r})()({1:[function(_dereq_,module,exports){
"use strict";

var somehowMaps = _dereq_('./src');

var w = somehowMaps({
  height: 300,
  aspect: 'widescreen'
});
w.shape({
  shape: 'great-lakes'
}).fill('lightblue');
w.shape({
  shape: 'provinces'
}).color('lightgrey');
w.shape({
  shape: 'north-america'
}); // w.shape({
//   shape: 'rivers'
// }).stroke('blue')

w.line().from('toronto').to('winnipeg').color('red');
w.dot().at('barrie').color('blue');
w.text('Toronto').at('toronto').color('red');
w.clip(false);
document.querySelector('#stage').innerHTML = w.build();

},{"./src":22}],2:[function(_dereq_,module,exports){
// https://d3js.org/d3-array/ v1.2.4 Copyright 2018 Mike Bostock
(function (global, factory) {
typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
typeof define === 'function' && define.amd ? define(['exports'], factory) :
(factory((global.d3 = global.d3 || {})));
}(this, (function (exports) { 'use strict';

function ascending(a, b) {
  return a < b ? -1 : a > b ? 1 : a >= b ? 0 : NaN;
}

function bisector(compare) {
  if (compare.length === 1) compare = ascendingComparator(compare);
  return {
    left: function(a, x, lo, hi) {
      if (lo == null) lo = 0;
      if (hi == null) hi = a.length;
      while (lo < hi) {
        var mid = lo + hi >>> 1;
        if (compare(a[mid], x) < 0) lo = mid + 1;
        else hi = mid;
      }
      return lo;
    },
    right: function(a, x, lo, hi) {
      if (lo == null) lo = 0;
      if (hi == null) hi = a.length;
      while (lo < hi) {
        var mid = lo + hi >>> 1;
        if (compare(a[mid], x) > 0) hi = mid;
        else lo = mid + 1;
      }
      return lo;
    }
  };
}

function ascendingComparator(f) {
  return function(d, x) {
    return ascending(f(d), x);
  };
}

var ascendingBisect = bisector(ascending);
var bisectRight = ascendingBisect.right;
var bisectLeft = ascendingBisect.left;

function pairs(array, f) {
  if (f == null) f = pair;
  var i = 0, n = array.length - 1, p = array[0], pairs = new Array(n < 0 ? 0 : n);
  while (i < n) pairs[i] = f(p, p = array[++i]);
  return pairs;
}

function pair(a, b) {
  return [a, b];
}

function cross(values0, values1, reduce) {
  var n0 = values0.length,
      n1 = values1.length,
      values = new Array(n0 * n1),
      i0,
      i1,
      i,
      value0;

  if (reduce == null) reduce = pair;

  for (i0 = i = 0; i0 < n0; ++i0) {
    for (value0 = values0[i0], i1 = 0; i1 < n1; ++i1, ++i) {
      values[i] = reduce(value0, values1[i1]);
    }
  }

  return values;
}

function descending(a, b) {
  return b < a ? -1 : b > a ? 1 : b >= a ? 0 : NaN;
}

function number(x) {
  return x === null ? NaN : +x;
}

function variance(values, valueof) {
  var n = values.length,
      m = 0,
      i = -1,
      mean = 0,
      value,
      delta,
      sum = 0;

  if (valueof == null) {
    while (++i < n) {
      if (!isNaN(value = number(values[i]))) {
        delta = value - mean;
        mean += delta / ++m;
        sum += delta * (value - mean);
      }
    }
  }

  else {
    while (++i < n) {
      if (!isNaN(value = number(valueof(values[i], i, values)))) {
        delta = value - mean;
        mean += delta / ++m;
        sum += delta * (value - mean);
      }
    }
  }

  if (m > 1) return sum / (m - 1);
}

function deviation(array, f) {
  var v = variance(array, f);
  return v ? Math.sqrt(v) : v;
}

function extent(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      min,
      max;

  if (valueof == null) {
    while (++i < n) { // Find the first comparable value.
      if ((value = values[i]) != null && value >= value) {
        min = max = value;
        while (++i < n) { // Compare the remaining values.
          if ((value = values[i]) != null) {
            if (min > value) min = value;
            if (max < value) max = value;
          }
        }
      }
    }
  }

  else {
    while (++i < n) { // Find the first comparable value.
      if ((value = valueof(values[i], i, values)) != null && value >= value) {
        min = max = value;
        while (++i < n) { // Compare the remaining values.
          if ((value = valueof(values[i], i, values)) != null) {
            if (min > value) min = value;
            if (max < value) max = value;
          }
        }
      }
    }
  }

  return [min, max];
}

var array = Array.prototype;

var slice = array.slice;
var map = array.map;

function constant(x) {
  return function() {
    return x;
  };
}

function identity(x) {
  return x;
}

function range(start, stop, step) {
  start = +start, stop = +stop, step = (n = arguments.length) < 2 ? (stop = start, start = 0, 1) : n < 3 ? 1 : +step;

  var i = -1,
      n = Math.max(0, Math.ceil((stop - start) / step)) | 0,
      range = new Array(n);

  while (++i < n) {
    range[i] = start + i * step;
  }

  return range;
}

var e10 = Math.sqrt(50),
    e5 = Math.sqrt(10),
    e2 = Math.sqrt(2);

function ticks(start, stop, count) {
  var reverse,
      i = -1,
      n,
      ticks,
      step;

  stop = +stop, start = +start, count = +count;
  if (start === stop && count > 0) return [start];
  if (reverse = stop < start) n = start, start = stop, stop = n;
  if ((step = tickIncrement(start, stop, count)) === 0 || !isFinite(step)) return [];

  if (step > 0) {
    start = Math.ceil(start / step);
    stop = Math.floor(stop / step);
    ticks = new Array(n = Math.ceil(stop - start + 1));
    while (++i < n) ticks[i] = (start + i) * step;
  } else {
    start = Math.floor(start * step);
    stop = Math.ceil(stop * step);
    ticks = new Array(n = Math.ceil(start - stop + 1));
    while (++i < n) ticks[i] = (start - i) / step;
  }

  if (reverse) ticks.reverse();

  return ticks;
}

function tickIncrement(start, stop, count) {
  var step = (stop - start) / Math.max(0, count),
      power = Math.floor(Math.log(step) / Math.LN10),
      error = step / Math.pow(10, power);
  return power >= 0
      ? (error >= e10 ? 10 : error >= e5 ? 5 : error >= e2 ? 2 : 1) * Math.pow(10, power)
      : -Math.pow(10, -power) / (error >= e10 ? 10 : error >= e5 ? 5 : error >= e2 ? 2 : 1);
}

function tickStep(start, stop, count) {
  var step0 = Math.abs(stop - start) / Math.max(0, count),
      step1 = Math.pow(10, Math.floor(Math.log(step0) / Math.LN10)),
      error = step0 / step1;
  if (error >= e10) step1 *= 10;
  else if (error >= e5) step1 *= 5;
  else if (error >= e2) step1 *= 2;
  return stop < start ? -step1 : step1;
}

function sturges(values) {
  return Math.ceil(Math.log(values.length) / Math.LN2) + 1;
}

function histogram() {
  var value = identity,
      domain = extent,
      threshold = sturges;

  function histogram(data) {
    var i,
        n = data.length,
        x,
        values = new Array(n);

    for (i = 0; i < n; ++i) {
      values[i] = value(data[i], i, data);
    }

    var xz = domain(values),
        x0 = xz[0],
        x1 = xz[1],
        tz = threshold(values, x0, x1);

    // Convert number of thresholds into uniform thresholds.
    if (!Array.isArray(tz)) {
      tz = tickStep(x0, x1, tz);
      tz = range(Math.ceil(x0 / tz) * tz, x1, tz); // exclusive
    }

    // Remove any thresholds outside the domain.
    var m = tz.length;
    while (tz[0] <= x0) tz.shift(), --m;
    while (tz[m - 1] > x1) tz.pop(), --m;

    var bins = new Array(m + 1),
        bin;

    // Initialize bins.
    for (i = 0; i <= m; ++i) {
      bin = bins[i] = [];
      bin.x0 = i > 0 ? tz[i - 1] : x0;
      bin.x1 = i < m ? tz[i] : x1;
    }

    // Assign data to bins by value, ignoring any outside the domain.
    for (i = 0; i < n; ++i) {
      x = values[i];
      if (x0 <= x && x <= x1) {
        bins[bisectRight(tz, x, 0, m)].push(data[i]);
      }
    }

    return bins;
  }

  histogram.value = function(_) {
    return arguments.length ? (value = typeof _ === "function" ? _ : constant(_), histogram) : value;
  };

  histogram.domain = function(_) {
    return arguments.length ? (domain = typeof _ === "function" ? _ : constant([_[0], _[1]]), histogram) : domain;
  };

  histogram.thresholds = function(_) {
    return arguments.length ? (threshold = typeof _ === "function" ? _ : Array.isArray(_) ? constant(slice.call(_)) : constant(_), histogram) : threshold;
  };

  return histogram;
}

function quantile(values, p, valueof) {
  if (valueof == null) valueof = number;
  if (!(n = values.length)) return;
  if ((p = +p) <= 0 || n < 2) return +valueof(values[0], 0, values);
  if (p >= 1) return +valueof(values[n - 1], n - 1, values);
  var n,
      i = (n - 1) * p,
      i0 = Math.floor(i),
      value0 = +valueof(values[i0], i0, values),
      value1 = +valueof(values[i0 + 1], i0 + 1, values);
  return value0 + (value1 - value0) * (i - i0);
}

function freedmanDiaconis(values, min, max) {
  values = map.call(values, number).sort(ascending);
  return Math.ceil((max - min) / (2 * (quantile(values, 0.75) - quantile(values, 0.25)) * Math.pow(values.length, -1 / 3)));
}

function scott(values, min, max) {
  return Math.ceil((max - min) / (3.5 * deviation(values) * Math.pow(values.length, -1 / 3)));
}

function max(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      max;

  if (valueof == null) {
    while (++i < n) { // Find the first comparable value.
      if ((value = values[i]) != null && value >= value) {
        max = value;
        while (++i < n) { // Compare the remaining values.
          if ((value = values[i]) != null && value > max) {
            max = value;
          }
        }
      }
    }
  }

  else {
    while (++i < n) { // Find the first comparable value.
      if ((value = valueof(values[i], i, values)) != null && value >= value) {
        max = value;
        while (++i < n) { // Compare the remaining values.
          if ((value = valueof(values[i], i, values)) != null && value > max) {
            max = value;
          }
        }
      }
    }
  }

  return max;
}

function mean(values, valueof) {
  var n = values.length,
      m = n,
      i = -1,
      value,
      sum = 0;

  if (valueof == null) {
    while (++i < n) {
      if (!isNaN(value = number(values[i]))) sum += value;
      else --m;
    }
  }

  else {
    while (++i < n) {
      if (!isNaN(value = number(valueof(values[i], i, values)))) sum += value;
      else --m;
    }
  }

  if (m) return sum / m;
}

function median(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      numbers = [];

  if (valueof == null) {
    while (++i < n) {
      if (!isNaN(value = number(values[i]))) {
        numbers.push(value);
      }
    }
  }

  else {
    while (++i < n) {
      if (!isNaN(value = number(valueof(values[i], i, values)))) {
        numbers.push(value);
      }
    }
  }

  return quantile(numbers.sort(ascending), 0.5);
}

function merge(arrays) {
  var n = arrays.length,
      m,
      i = -1,
      j = 0,
      merged,
      array;

  while (++i < n) j += arrays[i].length;
  merged = new Array(j);

  while (--n >= 0) {
    array = arrays[n];
    m = array.length;
    while (--m >= 0) {
      merged[--j] = array[m];
    }
  }

  return merged;
}

function min(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      min;

  if (valueof == null) {
    while (++i < n) { // Find the first comparable value.
      if ((value = values[i]) != null && value >= value) {
        min = value;
        while (++i < n) { // Compare the remaining values.
          if ((value = values[i]) != null && min > value) {
            min = value;
          }
        }
      }
    }
  }

  else {
    while (++i < n) { // Find the first comparable value.
      if ((value = valueof(values[i], i, values)) != null && value >= value) {
        min = value;
        while (++i < n) { // Compare the remaining values.
          if ((value = valueof(values[i], i, values)) != null && min > value) {
            min = value;
          }
        }
      }
    }
  }

  return min;
}

function permute(array, indexes) {
  var i = indexes.length, permutes = new Array(i);
  while (i--) permutes[i] = array[indexes[i]];
  return permutes;
}

function scan(values, compare) {
  if (!(n = values.length)) return;
  var n,
      i = 0,
      j = 0,
      xi,
      xj = values[j];

  if (compare == null) compare = ascending;

  while (++i < n) {
    if (compare(xi = values[i], xj) < 0 || compare(xj, xj) !== 0) {
      xj = xi, j = i;
    }
  }

  if (compare(xj, xj) === 0) return j;
}

function shuffle(array, i0, i1) {
  var m = (i1 == null ? array.length : i1) - (i0 = i0 == null ? 0 : +i0),
      t,
      i;

  while (m) {
    i = Math.random() * m-- | 0;
    t = array[m + i0];
    array[m + i0] = array[i + i0];
    array[i + i0] = t;
  }

  return array;
}

function sum(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      sum = 0;

  if (valueof == null) {
    while (++i < n) {
      if (value = +values[i]) sum += value; // Note: zero and null are equivalent.
    }
  }

  else {
    while (++i < n) {
      if (value = +valueof(values[i], i, values)) sum += value;
    }
  }

  return sum;
}

function transpose(matrix) {
  if (!(n = matrix.length)) return [];
  for (var i = -1, m = min(matrix, length), transpose = new Array(m); ++i < m;) {
    for (var j = -1, n, row = transpose[i] = new Array(n); ++j < n;) {
      row[j] = matrix[j][i];
    }
  }
  return transpose;
}

function length(d) {
  return d.length;
}

function zip() {
  return transpose(arguments);
}

exports.bisect = bisectRight;
exports.bisectRight = bisectRight;
exports.bisectLeft = bisectLeft;
exports.ascending = ascending;
exports.bisector = bisector;
exports.cross = cross;
exports.descending = descending;
exports.deviation = deviation;
exports.extent = extent;
exports.histogram = histogram;
exports.thresholdFreedmanDiaconis = freedmanDiaconis;
exports.thresholdScott = scott;
exports.thresholdSturges = sturges;
exports.max = max;
exports.mean = mean;
exports.median = median;
exports.merge = merge;
exports.min = min;
exports.pairs = pairs;
exports.permute = permute;
exports.quantile = quantile;
exports.range = range;
exports.scan = scan;
exports.shuffle = shuffle;
exports.sum = sum;
exports.ticks = ticks;
exports.tickIncrement = tickIncrement;
exports.tickStep = tickStep;
exports.transpose = transpose;
exports.variance = variance;
exports.zip = zip;

Object.defineProperty(exports, '__esModule', { value: true });

})));

},{}],3:[function(_dereq_,module,exports){
// https://d3js.org/d3-geo/ v1.11.3 Copyright 2018 Mike Bostock
(function (global, factory) {
typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports, _dereq_('d3-array')) :
typeof define === 'function' && define.amd ? define(['exports', 'd3-array'], factory) :
(factory((global.d3 = global.d3 || {}),global.d3));
}(this, (function (exports,d3Array) { 'use strict';

// Adds floating point numbers with twice the normal precision.
// Reference: J. R. Shewchuk, Adaptive Precision Floating-Point Arithmetic and
// Fast Robust Geometric Predicates, Discrete & Computational Geometry 18(3)
// 305–363 (1997).
// Code adapted from GeographicLib by Charles F. F. Karney,
// http://geographiclib.sourceforge.net/

function adder() {
  return new Adder;
}

function Adder() {
  this.reset();
}

Adder.prototype = {
  constructor: Adder,
  reset: function() {
    this.s = // rounded value
    this.t = 0; // exact error
  },
  add: function(y) {
    add(temp, y, this.t);
    add(this, temp.s, this.s);
    if (this.s) this.t += temp.t;
    else this.s = temp.t;
  },
  valueOf: function() {
    return this.s;
  }
};

var temp = new Adder;

function add(adder, a, b) {
  var x = adder.s = a + b,
      bv = x - a,
      av = x - bv;
  adder.t = (a - av) + (b - bv);
}

var epsilon = 1e-6;
var epsilon2 = 1e-12;
var pi = Math.PI;
var halfPi = pi / 2;
var quarterPi = pi / 4;
var tau = pi * 2;

var degrees = 180 / pi;
var radians = pi / 180;

var abs = Math.abs;
var atan = Math.atan;
var atan2 = Math.atan2;
var cos = Math.cos;
var ceil = Math.ceil;
var exp = Math.exp;
var log = Math.log;
var pow = Math.pow;
var sin = Math.sin;
var sign = Math.sign || function(x) { return x > 0 ? 1 : x < 0 ? -1 : 0; };
var sqrt = Math.sqrt;
var tan = Math.tan;

function acos(x) {
  return x > 1 ? 0 : x < -1 ? pi : Math.acos(x);
}

function asin(x) {
  return x > 1 ? halfPi : x < -1 ? -halfPi : Math.asin(x);
}

function haversin(x) {
  return (x = sin(x / 2)) * x;
}

function noop() {}

function streamGeometry(geometry, stream) {
  if (geometry && streamGeometryType.hasOwnProperty(geometry.type)) {
    streamGeometryType[geometry.type](geometry, stream);
  }
}

var streamObjectType = {
  Feature: function(object, stream) {
    streamGeometry(object.geometry, stream);
  },
  FeatureCollection: function(object, stream) {
    var features = object.features, i = -1, n = features.length;
    while (++i < n) streamGeometry(features[i].geometry, stream);
  }
};

var streamGeometryType = {
  Sphere: function(object, stream) {
    stream.sphere();
  },
  Point: function(object, stream) {
    object = object.coordinates;
    stream.point(object[0], object[1], object[2]);
  },
  MultiPoint: function(object, stream) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n) object = coordinates[i], stream.point(object[0], object[1], object[2]);
  },
  LineString: function(object, stream) {
    streamLine(object.coordinates, stream, 0);
  },
  MultiLineString: function(object, stream) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n) streamLine(coordinates[i], stream, 0);
  },
  Polygon: function(object, stream) {
    streamPolygon(object.coordinates, stream);
  },
  MultiPolygon: function(object, stream) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n) streamPolygon(coordinates[i], stream);
  },
  GeometryCollection: function(object, stream) {
    var geometries = object.geometries, i = -1, n = geometries.length;
    while (++i < n) streamGeometry(geometries[i], stream);
  }
};

function streamLine(coordinates, stream, closed) {
  var i = -1, n = coordinates.length - closed, coordinate;
  stream.lineStart();
  while (++i < n) coordinate = coordinates[i], stream.point(coordinate[0], coordinate[1], coordinate[2]);
  stream.lineEnd();
}

function streamPolygon(coordinates, stream) {
  var i = -1, n = coordinates.length;
  stream.polygonStart();
  while (++i < n) streamLine(coordinates[i], stream, 1);
  stream.polygonEnd();
}

function geoStream(object, stream) {
  if (object && streamObjectType.hasOwnProperty(object.type)) {
    streamObjectType[object.type](object, stream);
  } else {
    streamGeometry(object, stream);
  }
}

var areaRingSum = adder();

var areaSum = adder(),
    lambda00,
    phi00,
    lambda0,
    cosPhi0,
    sinPhi0;

var areaStream = {
  point: noop,
  lineStart: noop,
  lineEnd: noop,
  polygonStart: function() {
    areaRingSum.reset();
    areaStream.lineStart = areaRingStart;
    areaStream.lineEnd = areaRingEnd;
  },
  polygonEnd: function() {
    var areaRing = +areaRingSum;
    areaSum.add(areaRing < 0 ? tau + areaRing : areaRing);
    this.lineStart = this.lineEnd = this.point = noop;
  },
  sphere: function() {
    areaSum.add(tau);
  }
};

function areaRingStart() {
  areaStream.point = areaPointFirst;
}

function areaRingEnd() {
  areaPoint(lambda00, phi00);
}

function areaPointFirst(lambda, phi) {
  areaStream.point = areaPoint;
  lambda00 = lambda, phi00 = phi;
  lambda *= radians, phi *= radians;
  lambda0 = lambda, cosPhi0 = cos(phi = phi / 2 + quarterPi), sinPhi0 = sin(phi);
}

function areaPoint(lambda, phi) {
  lambda *= radians, phi *= radians;
  phi = phi / 2 + quarterPi; // half the angular distance from south pole

  // Spherical excess E for a spherical triangle with vertices: south pole,
  // previous point, current point.  Uses a formula derived from Cagnoli’s
  // theorem.  See Todhunter, Spherical Trig. (1871), Sec. 103, Eq. (2).
  var dLambda = lambda - lambda0,
      sdLambda = dLambda >= 0 ? 1 : -1,
      adLambda = sdLambda * dLambda,
      cosPhi = cos(phi),
      sinPhi = sin(phi),
      k = sinPhi0 * sinPhi,
      u = cosPhi0 * cosPhi + k * cos(adLambda),
      v = k * sdLambda * sin(adLambda);
  areaRingSum.add(atan2(v, u));

  // Advance the previous points.
  lambda0 = lambda, cosPhi0 = cosPhi, sinPhi0 = sinPhi;
}

function area(object) {
  areaSum.reset();
  geoStream(object, areaStream);
  return areaSum * 2;
}

function spherical(cartesian) {
  return [atan2(cartesian[1], cartesian[0]), asin(cartesian[2])];
}

function cartesian(spherical) {
  var lambda = spherical[0], phi = spherical[1], cosPhi = cos(phi);
  return [cosPhi * cos(lambda), cosPhi * sin(lambda), sin(phi)];
}

function cartesianDot(a, b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function cartesianCross(a, b) {
  return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}

// TODO return a
function cartesianAddInPlace(a, b) {
  a[0] += b[0], a[1] += b[1], a[2] += b[2];
}

function cartesianScale(vector, k) {
  return [vector[0] * k, vector[1] * k, vector[2] * k];
}

// TODO return d
function cartesianNormalizeInPlace(d) {
  var l = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  d[0] /= l, d[1] /= l, d[2] /= l;
}

var lambda0$1, phi0, lambda1, phi1, // bounds
    lambda2, // previous lambda-coordinate
    lambda00$1, phi00$1, // first point
    p0, // previous 3D point
    deltaSum = adder(),
    ranges,
    range;

var boundsStream = {
  point: boundsPoint,
  lineStart: boundsLineStart,
  lineEnd: boundsLineEnd,
  polygonStart: function() {
    boundsStream.point = boundsRingPoint;
    boundsStream.lineStart = boundsRingStart;
    boundsStream.lineEnd = boundsRingEnd;
    deltaSum.reset();
    areaStream.polygonStart();
  },
  polygonEnd: function() {
    areaStream.polygonEnd();
    boundsStream.point = boundsPoint;
    boundsStream.lineStart = boundsLineStart;
    boundsStream.lineEnd = boundsLineEnd;
    if (areaRingSum < 0) lambda0$1 = -(lambda1 = 180), phi0 = -(phi1 = 90);
    else if (deltaSum > epsilon) phi1 = 90;
    else if (deltaSum < -epsilon) phi0 = -90;
    range[0] = lambda0$1, range[1] = lambda1;
  }
};

function boundsPoint(lambda, phi) {
  ranges.push(range = [lambda0$1 = lambda, lambda1 = lambda]);
  if (phi < phi0) phi0 = phi;
  if (phi > phi1) phi1 = phi;
}

function linePoint(lambda, phi) {
  var p = cartesian([lambda * radians, phi * radians]);
  if (p0) {
    var normal = cartesianCross(p0, p),
        equatorial = [normal[1], -normal[0], 0],
        inflection = cartesianCross(equatorial, normal);
    cartesianNormalizeInPlace(inflection);
    inflection = spherical(inflection);
    var delta = lambda - lambda2,
        sign$$1 = delta > 0 ? 1 : -1,
        lambdai = inflection[0] * degrees * sign$$1,
        phii,
        antimeridian = abs(delta) > 180;
    if (antimeridian ^ (sign$$1 * lambda2 < lambdai && lambdai < sign$$1 * lambda)) {
      phii = inflection[1] * degrees;
      if (phii > phi1) phi1 = phii;
    } else if (lambdai = (lambdai + 360) % 360 - 180, antimeridian ^ (sign$$1 * lambda2 < lambdai && lambdai < sign$$1 * lambda)) {
      phii = -inflection[1] * degrees;
      if (phii < phi0) phi0 = phii;
    } else {
      if (phi < phi0) phi0 = phi;
      if (phi > phi1) phi1 = phi;
    }
    if (antimeridian) {
      if (lambda < lambda2) {
        if (angle(lambda0$1, lambda) > angle(lambda0$1, lambda1)) lambda1 = lambda;
      } else {
        if (angle(lambda, lambda1) > angle(lambda0$1, lambda1)) lambda0$1 = lambda;
      }
    } else {
      if (lambda1 >= lambda0$1) {
        if (lambda < lambda0$1) lambda0$1 = lambda;
        if (lambda > lambda1) lambda1 = lambda;
      } else {
        if (lambda > lambda2) {
          if (angle(lambda0$1, lambda) > angle(lambda0$1, lambda1)) lambda1 = lambda;
        } else {
          if (angle(lambda, lambda1) > angle(lambda0$1, lambda1)) lambda0$1 = lambda;
        }
      }
    }
  } else {
    ranges.push(range = [lambda0$1 = lambda, lambda1 = lambda]);
  }
  if (phi < phi0) phi0 = phi;
  if (phi > phi1) phi1 = phi;
  p0 = p, lambda2 = lambda;
}

function boundsLineStart() {
  boundsStream.point = linePoint;
}

function boundsLineEnd() {
  range[0] = lambda0$1, range[1] = lambda1;
  boundsStream.point = boundsPoint;
  p0 = null;
}

function boundsRingPoint(lambda, phi) {
  if (p0) {
    var delta = lambda - lambda2;
    deltaSum.add(abs(delta) > 180 ? delta + (delta > 0 ? 360 : -360) : delta);
  } else {
    lambda00$1 = lambda, phi00$1 = phi;
  }
  areaStream.point(lambda, phi);
  linePoint(lambda, phi);
}

function boundsRingStart() {
  areaStream.lineStart();
}

function boundsRingEnd() {
  boundsRingPoint(lambda00$1, phi00$1);
  areaStream.lineEnd();
  if (abs(deltaSum) > epsilon) lambda0$1 = -(lambda1 = 180);
  range[0] = lambda0$1, range[1] = lambda1;
  p0 = null;
}

// Finds the left-right distance between two longitudes.
// This is almost the same as (lambda1 - lambda0 + 360°) % 360°, except that we want
// the distance between ±180° to be 360°.
function angle(lambda0, lambda1) {
  return (lambda1 -= lambda0) < 0 ? lambda1 + 360 : lambda1;
}

function rangeCompare(a, b) {
  return a[0] - b[0];
}

function rangeContains(range, x) {
  return range[0] <= range[1] ? range[0] <= x && x <= range[1] : x < range[0] || range[1] < x;
}

function bounds(feature) {
  var i, n, a, b, merged, deltaMax, delta;

  phi1 = lambda1 = -(lambda0$1 = phi0 = Infinity);
  ranges = [];
  geoStream(feature, boundsStream);

  // First, sort ranges by their minimum longitudes.
  if (n = ranges.length) {
    ranges.sort(rangeCompare);

    // Then, merge any ranges that overlap.
    for (i = 1, a = ranges[0], merged = [a]; i < n; ++i) {
      b = ranges[i];
      if (rangeContains(a, b[0]) || rangeContains(a, b[1])) {
        if (angle(a[0], b[1]) > angle(a[0], a[1])) a[1] = b[1];
        if (angle(b[0], a[1]) > angle(a[0], a[1])) a[0] = b[0];
      } else {
        merged.push(a = b);
      }
    }

    // Finally, find the largest gap between the merged ranges.
    // The final bounding box will be the inverse of this gap.
    for (deltaMax = -Infinity, n = merged.length - 1, i = 0, a = merged[n]; i <= n; a = b, ++i) {
      b = merged[i];
      if ((delta = angle(a[1], b[0])) > deltaMax) deltaMax = delta, lambda0$1 = b[0], lambda1 = a[1];
    }
  }

  ranges = range = null;

  return lambda0$1 === Infinity || phi0 === Infinity
      ? [[NaN, NaN], [NaN, NaN]]
      : [[lambda0$1, phi0], [lambda1, phi1]];
}

var W0, W1,
    X0, Y0, Z0,
    X1, Y1, Z1,
    X2, Y2, Z2,
    lambda00$2, phi00$2, // first point
    x0, y0, z0; // previous point

var centroidStream = {
  sphere: noop,
  point: centroidPoint,
  lineStart: centroidLineStart,
  lineEnd: centroidLineEnd,
  polygonStart: function() {
    centroidStream.lineStart = centroidRingStart;
    centroidStream.lineEnd = centroidRingEnd;
  },
  polygonEnd: function() {
    centroidStream.lineStart = centroidLineStart;
    centroidStream.lineEnd = centroidLineEnd;
  }
};

// Arithmetic mean of Cartesian vectors.
function centroidPoint(lambda, phi) {
  lambda *= radians, phi *= radians;
  var cosPhi = cos(phi);
  centroidPointCartesian(cosPhi * cos(lambda), cosPhi * sin(lambda), sin(phi));
}

function centroidPointCartesian(x, y, z) {
  ++W0;
  X0 += (x - X0) / W0;
  Y0 += (y - Y0) / W0;
  Z0 += (z - Z0) / W0;
}

function centroidLineStart() {
  centroidStream.point = centroidLinePointFirst;
}

function centroidLinePointFirst(lambda, phi) {
  lambda *= radians, phi *= radians;
  var cosPhi = cos(phi);
  x0 = cosPhi * cos(lambda);
  y0 = cosPhi * sin(lambda);
  z0 = sin(phi);
  centroidStream.point = centroidLinePoint;
  centroidPointCartesian(x0, y0, z0);
}

function centroidLinePoint(lambda, phi) {
  lambda *= radians, phi *= radians;
  var cosPhi = cos(phi),
      x = cosPhi * cos(lambda),
      y = cosPhi * sin(lambda),
      z = sin(phi),
      w = atan2(sqrt((w = y0 * z - z0 * y) * w + (w = z0 * x - x0 * z) * w + (w = x0 * y - y0 * x) * w), x0 * x + y0 * y + z0 * z);
  W1 += w;
  X1 += w * (x0 + (x0 = x));
  Y1 += w * (y0 + (y0 = y));
  Z1 += w * (z0 + (z0 = z));
  centroidPointCartesian(x0, y0, z0);
}

function centroidLineEnd() {
  centroidStream.point = centroidPoint;
}

// See J. E. Brock, The Inertia Tensor for a Spherical Triangle,
// J. Applied Mechanics 42, 239 (1975).
function centroidRingStart() {
  centroidStream.point = centroidRingPointFirst;
}

function centroidRingEnd() {
  centroidRingPoint(lambda00$2, phi00$2);
  centroidStream.point = centroidPoint;
}

function centroidRingPointFirst(lambda, phi) {
  lambda00$2 = lambda, phi00$2 = phi;
  lambda *= radians, phi *= radians;
  centroidStream.point = centroidRingPoint;
  var cosPhi = cos(phi);
  x0 = cosPhi * cos(lambda);
  y0 = cosPhi * sin(lambda);
  z0 = sin(phi);
  centroidPointCartesian(x0, y0, z0);
}

function centroidRingPoint(lambda, phi) {
  lambda *= radians, phi *= radians;
  var cosPhi = cos(phi),
      x = cosPhi * cos(lambda),
      y = cosPhi * sin(lambda),
      z = sin(phi),
      cx = y0 * z - z0 * y,
      cy = z0 * x - x0 * z,
      cz = x0 * y - y0 * x,
      m = sqrt(cx * cx + cy * cy + cz * cz),
      w = asin(m), // line weight = angle
      v = m && -w / m; // area weight multiplier
  X2 += v * cx;
  Y2 += v * cy;
  Z2 += v * cz;
  W1 += w;
  X1 += w * (x0 + (x0 = x));
  Y1 += w * (y0 + (y0 = y));
  Z1 += w * (z0 + (z0 = z));
  centroidPointCartesian(x0, y0, z0);
}

function centroid(object) {
  W0 = W1 =
  X0 = Y0 = Z0 =
  X1 = Y1 = Z1 =
  X2 = Y2 = Z2 = 0;
  geoStream(object, centroidStream);

  var x = X2,
      y = Y2,
      z = Z2,
      m = x * x + y * y + z * z;

  // If the area-weighted ccentroid is undefined, fall back to length-weighted ccentroid.
  if (m < epsilon2) {
    x = X1, y = Y1, z = Z1;
    // If the feature has zero length, fall back to arithmetic mean of point vectors.
    if (W1 < epsilon) x = X0, y = Y0, z = Z0;
    m = x * x + y * y + z * z;
    // If the feature still has an undefined ccentroid, then return.
    if (m < epsilon2) return [NaN, NaN];
  }

  return [atan2(y, x) * degrees, asin(z / sqrt(m)) * degrees];
}

function constant(x) {
  return function() {
    return x;
  };
}

function compose(a, b) {

  function compose(x, y) {
    return x = a(x, y), b(x[0], x[1]);
  }

  if (a.invert && b.invert) compose.invert = function(x, y) {
    return x = b.invert(x, y), x && a.invert(x[0], x[1]);
  };

  return compose;
}

function rotationIdentity(lambda, phi) {
  return [abs(lambda) > pi ? lambda + Math.round(-lambda / tau) * tau : lambda, phi];
}

rotationIdentity.invert = rotationIdentity;

function rotateRadians(deltaLambda, deltaPhi, deltaGamma) {
  return (deltaLambda %= tau) ? (deltaPhi || deltaGamma ? compose(rotationLambda(deltaLambda), rotationPhiGamma(deltaPhi, deltaGamma))
    : rotationLambda(deltaLambda))
    : (deltaPhi || deltaGamma ? rotationPhiGamma(deltaPhi, deltaGamma)
    : rotationIdentity);
}

function forwardRotationLambda(deltaLambda) {
  return function(lambda, phi) {
    return lambda += deltaLambda, [lambda > pi ? lambda - tau : lambda < -pi ? lambda + tau : lambda, phi];
  };
}

function rotationLambda(deltaLambda) {
  var rotation = forwardRotationLambda(deltaLambda);
  rotation.invert = forwardRotationLambda(-deltaLambda);
  return rotation;
}

function rotationPhiGamma(deltaPhi, deltaGamma) {
  var cosDeltaPhi = cos(deltaPhi),
      sinDeltaPhi = sin(deltaPhi),
      cosDeltaGamma = cos(deltaGamma),
      sinDeltaGamma = sin(deltaGamma);

  function rotation(lambda, phi) {
    var cosPhi = cos(phi),
        x = cos(lambda) * cosPhi,
        y = sin(lambda) * cosPhi,
        z = sin(phi),
        k = z * cosDeltaPhi + x * sinDeltaPhi;
    return [
      atan2(y * cosDeltaGamma - k * sinDeltaGamma, x * cosDeltaPhi - z * sinDeltaPhi),
      asin(k * cosDeltaGamma + y * sinDeltaGamma)
    ];
  }

  rotation.invert = function(lambda, phi) {
    var cosPhi = cos(phi),
        x = cos(lambda) * cosPhi,
        y = sin(lambda) * cosPhi,
        z = sin(phi),
        k = z * cosDeltaGamma - y * sinDeltaGamma;
    return [
      atan2(y * cosDeltaGamma + z * sinDeltaGamma, x * cosDeltaPhi + k * sinDeltaPhi),
      asin(k * cosDeltaPhi - x * sinDeltaPhi)
    ];
  };

  return rotation;
}

function rotation(rotate) {
  rotate = rotateRadians(rotate[0] * radians, rotate[1] * radians, rotate.length > 2 ? rotate[2] * radians : 0);

  function forward(coordinates) {
    coordinates = rotate(coordinates[0] * radians, coordinates[1] * radians);
    return coordinates[0] *= degrees, coordinates[1] *= degrees, coordinates;
  }

  forward.invert = function(coordinates) {
    coordinates = rotate.invert(coordinates[0] * radians, coordinates[1] * radians);
    return coordinates[0] *= degrees, coordinates[1] *= degrees, coordinates;
  };

  return forward;
}

// Generates a circle centered at [0°, 0°], with a given radius and precision.
function circleStream(stream, radius, delta, direction, t0, t1) {
  if (!delta) return;
  var cosRadius = cos(radius),
      sinRadius = sin(radius),
      step = direction * delta;
  if (t0 == null) {
    t0 = radius + direction * tau;
    t1 = radius - step / 2;
  } else {
    t0 = circleRadius(cosRadius, t0);
    t1 = circleRadius(cosRadius, t1);
    if (direction > 0 ? t0 < t1 : t0 > t1) t0 += direction * tau;
  }
  for (var point, t = t0; direction > 0 ? t > t1 : t < t1; t -= step) {
    point = spherical([cosRadius, -sinRadius * cos(t), -sinRadius * sin(t)]);
    stream.point(point[0], point[1]);
  }
}

// Returns the signed angle of a cartesian point relative to [cosRadius, 0, 0].
function circleRadius(cosRadius, point) {
  point = cartesian(point), point[0] -= cosRadius;
  cartesianNormalizeInPlace(point);
  var radius = acos(-point[1]);
  return ((-point[2] < 0 ? -radius : radius) + tau - epsilon) % tau;
}

function circle() {
  var center = constant([0, 0]),
      radius = constant(90),
      precision = constant(6),
      ring,
      rotate,
      stream = {point: point};

  function point(x, y) {
    ring.push(x = rotate(x, y));
    x[0] *= degrees, x[1] *= degrees;
  }

  function circle() {
    var c = center.apply(this, arguments),
        r = radius.apply(this, arguments) * radians,
        p = precision.apply(this, arguments) * radians;
    ring = [];
    rotate = rotateRadians(-c[0] * radians, -c[1] * radians, 0).invert;
    circleStream(stream, r, p, 1);
    c = {type: "Polygon", coordinates: [ring]};
    ring = rotate = null;
    return c;
  }

  circle.center = function(_) {
    return arguments.length ? (center = typeof _ === "function" ? _ : constant([+_[0], +_[1]]), circle) : center;
  };

  circle.radius = function(_) {
    return arguments.length ? (radius = typeof _ === "function" ? _ : constant(+_), circle) : radius;
  };

  circle.precision = function(_) {
    return arguments.length ? (precision = typeof _ === "function" ? _ : constant(+_), circle) : precision;
  };

  return circle;
}

function clipBuffer() {
  var lines = [],
      line;
  return {
    point: function(x, y) {
      line.push([x, y]);
    },
    lineStart: function() {
      lines.push(line = []);
    },
    lineEnd: noop,
    rejoin: function() {
      if (lines.length > 1) lines.push(lines.pop().concat(lines.shift()));
    },
    result: function() {
      var result = lines;
      lines = [];
      line = null;
      return result;
    }
  };
}

function pointEqual(a, b) {
  return abs(a[0] - b[0]) < epsilon && abs(a[1] - b[1]) < epsilon;
}

function Intersection(point, points, other, entry) {
  this.x = point;
  this.z = points;
  this.o = other; // another intersection
  this.e = entry; // is an entry?
  this.v = false; // visited
  this.n = this.p = null; // next & previous
}

// A generalized polygon clipping algorithm: given a polygon that has been cut
// into its visible line segments, and rejoins the segments by interpolating
// along the clip edge.
function clipRejoin(segments, compareIntersection, startInside, interpolate, stream) {
  var subject = [],
      clip = [],
      i,
      n;

  segments.forEach(function(segment) {
    if ((n = segment.length - 1) <= 0) return;
    var n, p0 = segment[0], p1 = segment[n], x;

    // If the first and last points of a segment are coincident, then treat as a
    // closed ring. TODO if all rings are closed, then the winding order of the
    // exterior ring should be checked.
    if (pointEqual(p0, p1)) {
      stream.lineStart();
      for (i = 0; i < n; ++i) stream.point((p0 = segment[i])[0], p0[1]);
      stream.lineEnd();
      return;
    }

    subject.push(x = new Intersection(p0, segment, null, true));
    clip.push(x.o = new Intersection(p0, null, x, false));
    subject.push(x = new Intersection(p1, segment, null, false));
    clip.push(x.o = new Intersection(p1, null, x, true));
  });

  if (!subject.length) return;

  clip.sort(compareIntersection);
  link(subject);
  link(clip);

  for (i = 0, n = clip.length; i < n; ++i) {
    clip[i].e = startInside = !startInside;
  }

  var start = subject[0],
      points,
      point;

  while (1) {
    // Find first unvisited intersection.
    var current = start,
        isSubject = true;
    while (current.v) if ((current = current.n) === start) return;
    points = current.z;
    stream.lineStart();
    do {
      current.v = current.o.v = true;
      if (current.e) {
        if (isSubject) {
          for (i = 0, n = points.length; i < n; ++i) stream.point((point = points[i])[0], point[1]);
        } else {
          interpolate(current.x, current.n.x, 1, stream);
        }
        current = current.n;
      } else {
        if (isSubject) {
          points = current.p.z;
          for (i = points.length - 1; i >= 0; --i) stream.point((point = points[i])[0], point[1]);
        } else {
          interpolate(current.x, current.p.x, -1, stream);
        }
        current = current.p;
      }
      current = current.o;
      points = current.z;
      isSubject = !isSubject;
    } while (!current.v);
    stream.lineEnd();
  }
}

function link(array) {
  if (!(n = array.length)) return;
  var n,
      i = 0,
      a = array[0],
      b;
  while (++i < n) {
    a.n = b = array[i];
    b.p = a;
    a = b;
  }
  a.n = b = array[0];
  b.p = a;
}

var sum = adder();

function polygonContains(polygon, point) {
  var lambda = point[0],
      phi = point[1],
      sinPhi = sin(phi),
      normal = [sin(lambda), -cos(lambda), 0],
      angle = 0,
      winding = 0;

  sum.reset();

  if (sinPhi === 1) phi = halfPi + epsilon;
  else if (sinPhi === -1) phi = -halfPi - epsilon;

  for (var i = 0, n = polygon.length; i < n; ++i) {
    if (!(m = (ring = polygon[i]).length)) continue;
    var ring,
        m,
        point0 = ring[m - 1],
        lambda0 = point0[0],
        phi0 = point0[1] / 2 + quarterPi,
        sinPhi0 = sin(phi0),
        cosPhi0 = cos(phi0);

    for (var j = 0; j < m; ++j, lambda0 = lambda1, sinPhi0 = sinPhi1, cosPhi0 = cosPhi1, point0 = point1) {
      var point1 = ring[j],
          lambda1 = point1[0],
          phi1 = point1[1] / 2 + quarterPi,
          sinPhi1 = sin(phi1),
          cosPhi1 = cos(phi1),
          delta = lambda1 - lambda0,
          sign$$1 = delta >= 0 ? 1 : -1,
          absDelta = sign$$1 * delta,
          antimeridian = absDelta > pi,
          k = sinPhi0 * sinPhi1;

      sum.add(atan2(k * sign$$1 * sin(absDelta), cosPhi0 * cosPhi1 + k * cos(absDelta)));
      angle += antimeridian ? delta + sign$$1 * tau : delta;

      // Are the longitudes either side of the point’s meridian (lambda),
      // and are the latitudes smaller than the parallel (phi)?
      if (antimeridian ^ lambda0 >= lambda ^ lambda1 >= lambda) {
        var arc = cartesianCross(cartesian(point0), cartesian(point1));
        cartesianNormalizeInPlace(arc);
        var intersection = cartesianCross(normal, arc);
        cartesianNormalizeInPlace(intersection);
        var phiArc = (antimeridian ^ delta >= 0 ? -1 : 1) * asin(intersection[2]);
        if (phi > phiArc || phi === phiArc && (arc[0] || arc[1])) {
          winding += antimeridian ^ delta >= 0 ? 1 : -1;
        }
      }
    }
  }

  // First, determine whether the South pole is inside or outside:
  //
  // It is inside if:
  // * the polygon winds around it in a clockwise direction.
  // * the polygon does not (cumulatively) wind around it, but has a negative
  //   (counter-clockwise) area.
  //
  // Second, count the (signed) number of times a segment crosses a lambda
  // from the point to the South pole.  If it is zero, then the point is the
  // same side as the South pole.

  return (angle < -epsilon || angle < epsilon && sum < -epsilon) ^ (winding & 1);
}

function clip(pointVisible, clipLine, interpolate, start) {
  return function(sink) {
    var line = clipLine(sink),
        ringBuffer = clipBuffer(),
        ringSink = clipLine(ringBuffer),
        polygonStarted = false,
        polygon,
        segments,
        ring;

    var clip = {
      point: point,
      lineStart: lineStart,
      lineEnd: lineEnd,
      polygonStart: function() {
        clip.point = pointRing;
        clip.lineStart = ringStart;
        clip.lineEnd = ringEnd;
        segments = [];
        polygon = [];
      },
      polygonEnd: function() {
        clip.point = point;
        clip.lineStart = lineStart;
        clip.lineEnd = lineEnd;
        segments = d3Array.merge(segments);
        var startInside = polygonContains(polygon, start);
        if (segments.length) {
          if (!polygonStarted) sink.polygonStart(), polygonStarted = true;
          clipRejoin(segments, compareIntersection, startInside, interpolate, sink);
        } else if (startInside) {
          if (!polygonStarted) sink.polygonStart(), polygonStarted = true;
          sink.lineStart();
          interpolate(null, null, 1, sink);
          sink.lineEnd();
        }
        if (polygonStarted) sink.polygonEnd(), polygonStarted = false;
        segments = polygon = null;
      },
      sphere: function() {
        sink.polygonStart();
        sink.lineStart();
        interpolate(null, null, 1, sink);
        sink.lineEnd();
        sink.polygonEnd();
      }
    };

    function point(lambda, phi) {
      if (pointVisible(lambda, phi)) sink.point(lambda, phi);
    }

    function pointLine(lambda, phi) {
      line.point(lambda, phi);
    }

    function lineStart() {
      clip.point = pointLine;
      line.lineStart();
    }

    function lineEnd() {
      clip.point = point;
      line.lineEnd();
    }

    function pointRing(lambda, phi) {
      ring.push([lambda, phi]);
      ringSink.point(lambda, phi);
    }

    function ringStart() {
      ringSink.lineStart();
      ring = [];
    }

    function ringEnd() {
      pointRing(ring[0][0], ring[0][1]);
      ringSink.lineEnd();

      var clean = ringSink.clean(),
          ringSegments = ringBuffer.result(),
          i, n = ringSegments.length, m,
          segment,
          point;

      ring.pop();
      polygon.push(ring);
      ring = null;

      if (!n) return;

      // No intersections.
      if (clean & 1) {
        segment = ringSegments[0];
        if ((m = segment.length - 1) > 0) {
          if (!polygonStarted) sink.polygonStart(), polygonStarted = true;
          sink.lineStart();
          for (i = 0; i < m; ++i) sink.point((point = segment[i])[0], point[1]);
          sink.lineEnd();
        }
        return;
      }

      // Rejoin connected segments.
      // TODO reuse ringBuffer.rejoin()?
      if (n > 1 && clean & 2) ringSegments.push(ringSegments.pop().concat(ringSegments.shift()));

      segments.push(ringSegments.filter(validSegment));
    }

    return clip;
  };
}

function validSegment(segment) {
  return segment.length > 1;
}

// Intersections are sorted along the clip edge. For both antimeridian cutting
// and circle clipping, the same comparison is used.
function compareIntersection(a, b) {
  return ((a = a.x)[0] < 0 ? a[1] - halfPi - epsilon : halfPi - a[1])
       - ((b = b.x)[0] < 0 ? b[1] - halfPi - epsilon : halfPi - b[1]);
}

var clipAntimeridian = clip(
  function() { return true; },
  clipAntimeridianLine,
  clipAntimeridianInterpolate,
  [-pi, -halfPi]
);

// Takes a line and cuts into visible segments. Return values: 0 - there were
// intersections or the line was empty; 1 - no intersections; 2 - there were
// intersections, and the first and last segments should be rejoined.
function clipAntimeridianLine(stream) {
  var lambda0 = NaN,
      phi0 = NaN,
      sign0 = NaN,
      clean; // no intersections

  return {
    lineStart: function() {
      stream.lineStart();
      clean = 1;
    },
    point: function(lambda1, phi1) {
      var sign1 = lambda1 > 0 ? pi : -pi,
          delta = abs(lambda1 - lambda0);
      if (abs(delta - pi) < epsilon) { // line crosses a pole
        stream.point(lambda0, phi0 = (phi0 + phi1) / 2 > 0 ? halfPi : -halfPi);
        stream.point(sign0, phi0);
        stream.lineEnd();
        stream.lineStart();
        stream.point(sign1, phi0);
        stream.point(lambda1, phi0);
        clean = 0;
      } else if (sign0 !== sign1 && delta >= pi) { // line crosses antimeridian
        if (abs(lambda0 - sign0) < epsilon) lambda0 -= sign0 * epsilon; // handle degeneracies
        if (abs(lambda1 - sign1) < epsilon) lambda1 -= sign1 * epsilon;
        phi0 = clipAntimeridianIntersect(lambda0, phi0, lambda1, phi1);
        stream.point(sign0, phi0);
        stream.lineEnd();
        stream.lineStart();
        stream.point(sign1, phi0);
        clean = 0;
      }
      stream.point(lambda0 = lambda1, phi0 = phi1);
      sign0 = sign1;
    },
    lineEnd: function() {
      stream.lineEnd();
      lambda0 = phi0 = NaN;
    },
    clean: function() {
      return 2 - clean; // if intersections, rejoin first and last segments
    }
  };
}

function clipAntimeridianIntersect(lambda0, phi0, lambda1, phi1) {
  var cosPhi0,
      cosPhi1,
      sinLambda0Lambda1 = sin(lambda0 - lambda1);
  return abs(sinLambda0Lambda1) > epsilon
      ? atan((sin(phi0) * (cosPhi1 = cos(phi1)) * sin(lambda1)
          - sin(phi1) * (cosPhi0 = cos(phi0)) * sin(lambda0))
          / (cosPhi0 * cosPhi1 * sinLambda0Lambda1))
      : (phi0 + phi1) / 2;
}

function clipAntimeridianInterpolate(from, to, direction, stream) {
  var phi;
  if (from == null) {
    phi = direction * halfPi;
    stream.point(-pi, phi);
    stream.point(0, phi);
    stream.point(pi, phi);
    stream.point(pi, 0);
    stream.point(pi, -phi);
    stream.point(0, -phi);
    stream.point(-pi, -phi);
    stream.point(-pi, 0);
    stream.point(-pi, phi);
  } else if (abs(from[0] - to[0]) > epsilon) {
    var lambda = from[0] < to[0] ? pi : -pi;
    phi = direction * lambda / 2;
    stream.point(-lambda, phi);
    stream.point(0, phi);
    stream.point(lambda, phi);
  } else {
    stream.point(to[0], to[1]);
  }
}

function clipCircle(radius) {
  var cr = cos(radius),
      delta = 6 * radians,
      smallRadius = cr > 0,
      notHemisphere = abs(cr) > epsilon; // TODO optimise for this common case

  function interpolate(from, to, direction, stream) {
    circleStream(stream, radius, delta, direction, from, to);
  }

  function visible(lambda, phi) {
    return cos(lambda) * cos(phi) > cr;
  }

  // Takes a line and cuts into visible segments. Return values used for polygon
  // clipping: 0 - there were intersections or the line was empty; 1 - no
  // intersections 2 - there were intersections, and the first and last segments
  // should be rejoined.
  function clipLine(stream) {
    var point0, // previous point
        c0, // code for previous point
        v0, // visibility of previous point
        v00, // visibility of first point
        clean; // no intersections
    return {
      lineStart: function() {
        v00 = v0 = false;
        clean = 1;
      },
      point: function(lambda, phi) {
        var point1 = [lambda, phi],
            point2,
            v = visible(lambda, phi),
            c = smallRadius
              ? v ? 0 : code(lambda, phi)
              : v ? code(lambda + (lambda < 0 ? pi : -pi), phi) : 0;
        if (!point0 && (v00 = v0 = v)) stream.lineStart();
        // Handle degeneracies.
        // TODO ignore if not clipping polygons.
        if (v !== v0) {
          point2 = intersect(point0, point1);
          if (!point2 || pointEqual(point0, point2) || pointEqual(point1, point2)) {
            point1[0] += epsilon;
            point1[1] += epsilon;
            v = visible(point1[0], point1[1]);
          }
        }
        if (v !== v0) {
          clean = 0;
          if (v) {
            // outside going in
            stream.lineStart();
            point2 = intersect(point1, point0);
            stream.point(point2[0], point2[1]);
          } else {
            // inside going out
            point2 = intersect(point0, point1);
            stream.point(point2[0], point2[1]);
            stream.lineEnd();
          }
          point0 = point2;
        } else if (notHemisphere && point0 && smallRadius ^ v) {
          var t;
          // If the codes for two points are different, or are both zero,
          // and there this segment intersects with the small circle.
          if (!(c & c0) && (t = intersect(point1, point0, true))) {
            clean = 0;
            if (smallRadius) {
              stream.lineStart();
              stream.point(t[0][0], t[0][1]);
              stream.point(t[1][0], t[1][1]);
              stream.lineEnd();
            } else {
              stream.point(t[1][0], t[1][1]);
              stream.lineEnd();
              stream.lineStart();
              stream.point(t[0][0], t[0][1]);
            }
          }
        }
        if (v && (!point0 || !pointEqual(point0, point1))) {
          stream.point(point1[0], point1[1]);
        }
        point0 = point1, v0 = v, c0 = c;
      },
      lineEnd: function() {
        if (v0) stream.lineEnd();
        point0 = null;
      },
      // Rejoin first and last segments if there were intersections and the first
      // and last points were visible.
      clean: function() {
        return clean | ((v00 && v0) << 1);
      }
    };
  }

  // Intersects the great circle between a and b with the clip circle.
  function intersect(a, b, two) {
    var pa = cartesian(a),
        pb = cartesian(b);

    // We have two planes, n1.p = d1 and n2.p = d2.
    // Find intersection line p(t) = c1 n1 + c2 n2 + t (n1 ⨯ n2).
    var n1 = [1, 0, 0], // normal
        n2 = cartesianCross(pa, pb),
        n2n2 = cartesianDot(n2, n2),
        n1n2 = n2[0], // cartesianDot(n1, n2),
        determinant = n2n2 - n1n2 * n1n2;

    // Two polar points.
    if (!determinant) return !two && a;

    var c1 =  cr * n2n2 / determinant,
        c2 = -cr * n1n2 / determinant,
        n1xn2 = cartesianCross(n1, n2),
        A = cartesianScale(n1, c1),
        B = cartesianScale(n2, c2);
    cartesianAddInPlace(A, B);

    // Solve |p(t)|^2 = 1.
    var u = n1xn2,
        w = cartesianDot(A, u),
        uu = cartesianDot(u, u),
        t2 = w * w - uu * (cartesianDot(A, A) - 1);

    if (t2 < 0) return;

    var t = sqrt(t2),
        q = cartesianScale(u, (-w - t) / uu);
    cartesianAddInPlace(q, A);
    q = spherical(q);

    if (!two) return q;

    // Two intersection points.
    var lambda0 = a[0],
        lambda1 = b[0],
        phi0 = a[1],
        phi1 = b[1],
        z;

    if (lambda1 < lambda0) z = lambda0, lambda0 = lambda1, lambda1 = z;

    var delta = lambda1 - lambda0,
        polar = abs(delta - pi) < epsilon,
        meridian = polar || delta < epsilon;

    if (!polar && phi1 < phi0) z = phi0, phi0 = phi1, phi1 = z;

    // Check that the first point is between a and b.
    if (meridian
        ? polar
          ? phi0 + phi1 > 0 ^ q[1] < (abs(q[0] - lambda0) < epsilon ? phi0 : phi1)
          : phi0 <= q[1] && q[1] <= phi1
        : delta > pi ^ (lambda0 <= q[0] && q[0] <= lambda1)) {
      var q1 = cartesianScale(u, (-w + t) / uu);
      cartesianAddInPlace(q1, A);
      return [q, spherical(q1)];
    }
  }

  // Generates a 4-bit vector representing the location of a point relative to
  // the small circle's bounding box.
  function code(lambda, phi) {
    var r = smallRadius ? radius : pi - radius,
        code = 0;
    if (lambda < -r) code |= 1; // left
    else if (lambda > r) code |= 2; // right
    if (phi < -r) code |= 4; // below
    else if (phi > r) code |= 8; // above
    return code;
  }

  return clip(visible, clipLine, interpolate, smallRadius ? [0, -radius] : [-pi, radius - pi]);
}

function clipLine(a, b, x0, y0, x1, y1) {
  var ax = a[0],
      ay = a[1],
      bx = b[0],
      by = b[1],
      t0 = 0,
      t1 = 1,
      dx = bx - ax,
      dy = by - ay,
      r;

  r = x0 - ax;
  if (!dx && r > 0) return;
  r /= dx;
  if (dx < 0) {
    if (r < t0) return;
    if (r < t1) t1 = r;
  } else if (dx > 0) {
    if (r > t1) return;
    if (r > t0) t0 = r;
  }

  r = x1 - ax;
  if (!dx && r < 0) return;
  r /= dx;
  if (dx < 0) {
    if (r > t1) return;
    if (r > t0) t0 = r;
  } else if (dx > 0) {
    if (r < t0) return;
    if (r < t1) t1 = r;
  }

  r = y0 - ay;
  if (!dy && r > 0) return;
  r /= dy;
  if (dy < 0) {
    if (r < t0) return;
    if (r < t1) t1 = r;
  } else if (dy > 0) {
    if (r > t1) return;
    if (r > t0) t0 = r;
  }

  r = y1 - ay;
  if (!dy && r < 0) return;
  r /= dy;
  if (dy < 0) {
    if (r > t1) return;
    if (r > t0) t0 = r;
  } else if (dy > 0) {
    if (r < t0) return;
    if (r < t1) t1 = r;
  }

  if (t0 > 0) a[0] = ax + t0 * dx, a[1] = ay + t0 * dy;
  if (t1 < 1) b[0] = ax + t1 * dx, b[1] = ay + t1 * dy;
  return true;
}

var clipMax = 1e9, clipMin = -clipMax;

// TODO Use d3-polygon’s polygonContains here for the ring check?
// TODO Eliminate duplicate buffering in clipBuffer and polygon.push?

function clipRectangle(x0, y0, x1, y1) {

  function visible(x, y) {
    return x0 <= x && x <= x1 && y0 <= y && y <= y1;
  }

  function interpolate(from, to, direction, stream) {
    var a = 0, a1 = 0;
    if (from == null
        || (a = corner(from, direction)) !== (a1 = corner(to, direction))
        || comparePoint(from, to) < 0 ^ direction > 0) {
      do stream.point(a === 0 || a === 3 ? x0 : x1, a > 1 ? y1 : y0);
      while ((a = (a + direction + 4) % 4) !== a1);
    } else {
      stream.point(to[0], to[1]);
    }
  }

  function corner(p, direction) {
    return abs(p[0] - x0) < epsilon ? direction > 0 ? 0 : 3
        : abs(p[0] - x1) < epsilon ? direction > 0 ? 2 : 1
        : abs(p[1] - y0) < epsilon ? direction > 0 ? 1 : 0
        : direction > 0 ? 3 : 2; // abs(p[1] - y1) < epsilon
  }

  function compareIntersection(a, b) {
    return comparePoint(a.x, b.x);
  }

  function comparePoint(a, b) {
    var ca = corner(a, 1),
        cb = corner(b, 1);
    return ca !== cb ? ca - cb
        : ca === 0 ? b[1] - a[1]
        : ca === 1 ? a[0] - b[0]
        : ca === 2 ? a[1] - b[1]
        : b[0] - a[0];
  }

  return function(stream) {
    var activeStream = stream,
        bufferStream = clipBuffer(),
        segments,
        polygon,
        ring,
        x__, y__, v__, // first point
        x_, y_, v_, // previous point
        first,
        clean;

    var clipStream = {
      point: point,
      lineStart: lineStart,
      lineEnd: lineEnd,
      polygonStart: polygonStart,
      polygonEnd: polygonEnd
    };

    function point(x, y) {
      if (visible(x, y)) activeStream.point(x, y);
    }

    function polygonInside() {
      var winding = 0;

      for (var i = 0, n = polygon.length; i < n; ++i) {
        for (var ring = polygon[i], j = 1, m = ring.length, point = ring[0], a0, a1, b0 = point[0], b1 = point[1]; j < m; ++j) {
          a0 = b0, a1 = b1, point = ring[j], b0 = point[0], b1 = point[1];
          if (a1 <= y1) { if (b1 > y1 && (b0 - a0) * (y1 - a1) > (b1 - a1) * (x0 - a0)) ++winding; }
          else { if (b1 <= y1 && (b0 - a0) * (y1 - a1) < (b1 - a1) * (x0 - a0)) --winding; }
        }
      }

      return winding;
    }

    // Buffer geometry within a polygon and then clip it en masse.
    function polygonStart() {
      activeStream = bufferStream, segments = [], polygon = [], clean = true;
    }

    function polygonEnd() {
      var startInside = polygonInside(),
          cleanInside = clean && startInside,
          visible = (segments = d3Array.merge(segments)).length;
      if (cleanInside || visible) {
        stream.polygonStart();
        if (cleanInside) {
          stream.lineStart();
          interpolate(null, null, 1, stream);
          stream.lineEnd();
        }
        if (visible) {
          clipRejoin(segments, compareIntersection, startInside, interpolate, stream);
        }
        stream.polygonEnd();
      }
      activeStream = stream, segments = polygon = ring = null;
    }

    function lineStart() {
      clipStream.point = linePoint;
      if (polygon) polygon.push(ring = []);
      first = true;
      v_ = false;
      x_ = y_ = NaN;
    }

    // TODO rather than special-case polygons, simply handle them separately.
    // Ideally, coincident intersection points should be jittered to avoid
    // clipping issues.
    function lineEnd() {
      if (segments) {
        linePoint(x__, y__);
        if (v__ && v_) bufferStream.rejoin();
        segments.push(bufferStream.result());
      }
      clipStream.point = point;
      if (v_) activeStream.lineEnd();
    }

    function linePoint(x, y) {
      var v = visible(x, y);
      if (polygon) ring.push([x, y]);
      if (first) {
        x__ = x, y__ = y, v__ = v;
        first = false;
        if (v) {
          activeStream.lineStart();
          activeStream.point(x, y);
        }
      } else {
        if (v && v_) activeStream.point(x, y);
        else {
          var a = [x_ = Math.max(clipMin, Math.min(clipMax, x_)), y_ = Math.max(clipMin, Math.min(clipMax, y_))],
              b = [x = Math.max(clipMin, Math.min(clipMax, x)), y = Math.max(clipMin, Math.min(clipMax, y))];
          if (clipLine(a, b, x0, y0, x1, y1)) {
            if (!v_) {
              activeStream.lineStart();
              activeStream.point(a[0], a[1]);
            }
            activeStream.point(b[0], b[1]);
            if (!v) activeStream.lineEnd();
            clean = false;
          } else if (v) {
            activeStream.lineStart();
            activeStream.point(x, y);
            clean = false;
          }
        }
      }
      x_ = x, y_ = y, v_ = v;
    }

    return clipStream;
  };
}

function extent() {
  var x0 = 0,
      y0 = 0,
      x1 = 960,
      y1 = 500,
      cache,
      cacheStream,
      clip;

  return clip = {
    stream: function(stream) {
      return cache && cacheStream === stream ? cache : cache = clipRectangle(x0, y0, x1, y1)(cacheStream = stream);
    },
    extent: function(_) {
      return arguments.length ? (x0 = +_[0][0], y0 = +_[0][1], x1 = +_[1][0], y1 = +_[1][1], cache = cacheStream = null, clip) : [[x0, y0], [x1, y1]];
    }
  };
}

var lengthSum = adder(),
    lambda0$2,
    sinPhi0$1,
    cosPhi0$1;

var lengthStream = {
  sphere: noop,
  point: noop,
  lineStart: lengthLineStart,
  lineEnd: noop,
  polygonStart: noop,
  polygonEnd: noop
};

function lengthLineStart() {
  lengthStream.point = lengthPointFirst;
  lengthStream.lineEnd = lengthLineEnd;
}

function lengthLineEnd() {
  lengthStream.point = lengthStream.lineEnd = noop;
}

function lengthPointFirst(lambda, phi) {
  lambda *= radians, phi *= radians;
  lambda0$2 = lambda, sinPhi0$1 = sin(phi), cosPhi0$1 = cos(phi);
  lengthStream.point = lengthPoint;
}

function lengthPoint(lambda, phi) {
  lambda *= radians, phi *= radians;
  var sinPhi = sin(phi),
      cosPhi = cos(phi),
      delta = abs(lambda - lambda0$2),
      cosDelta = cos(delta),
      sinDelta = sin(delta),
      x = cosPhi * sinDelta,
      y = cosPhi0$1 * sinPhi - sinPhi0$1 * cosPhi * cosDelta,
      z = sinPhi0$1 * sinPhi + cosPhi0$1 * cosPhi * cosDelta;
  lengthSum.add(atan2(sqrt(x * x + y * y), z));
  lambda0$2 = lambda, sinPhi0$1 = sinPhi, cosPhi0$1 = cosPhi;
}

function length(object) {
  lengthSum.reset();
  geoStream(object, lengthStream);
  return +lengthSum;
}

var coordinates = [null, null],
    object = {type: "LineString", coordinates: coordinates};

function distance(a, b) {
  coordinates[0] = a;
  coordinates[1] = b;
  return length(object);
}

var containsObjectType = {
  Feature: function(object, point) {
    return containsGeometry(object.geometry, point);
  },
  FeatureCollection: function(object, point) {
    var features = object.features, i = -1, n = features.length;
    while (++i < n) if (containsGeometry(features[i].geometry, point)) return true;
    return false;
  }
};

var containsGeometryType = {
  Sphere: function() {
    return true;
  },
  Point: function(object, point) {
    return containsPoint(object.coordinates, point);
  },
  MultiPoint: function(object, point) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n) if (containsPoint(coordinates[i], point)) return true;
    return false;
  },
  LineString: function(object, point) {
    return containsLine(object.coordinates, point);
  },
  MultiLineString: function(object, point) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n) if (containsLine(coordinates[i], point)) return true;
    return false;
  },
  Polygon: function(object, point) {
    return containsPolygon(object.coordinates, point);
  },
  MultiPolygon: function(object, point) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n) if (containsPolygon(coordinates[i], point)) return true;
    return false;
  },
  GeometryCollection: function(object, point) {
    var geometries = object.geometries, i = -1, n = geometries.length;
    while (++i < n) if (containsGeometry(geometries[i], point)) return true;
    return false;
  }
};

function containsGeometry(geometry, point) {
  return geometry && containsGeometryType.hasOwnProperty(geometry.type)
      ? containsGeometryType[geometry.type](geometry, point)
      : false;
}

function containsPoint(coordinates, point) {
  return distance(coordinates, point) === 0;
}

function containsLine(coordinates, point) {
  var ab = distance(coordinates[0], coordinates[1]),
      ao = distance(coordinates[0], point),
      ob = distance(point, coordinates[1]);
  return ao + ob <= ab + epsilon;
}

function containsPolygon(coordinates, point) {
  return !!polygonContains(coordinates.map(ringRadians), pointRadians(point));
}

function ringRadians(ring) {
  return ring = ring.map(pointRadians), ring.pop(), ring;
}

function pointRadians(point) {
  return [point[0] * radians, point[1] * radians];
}

function contains(object, point) {
  return (object && containsObjectType.hasOwnProperty(object.type)
      ? containsObjectType[object.type]
      : containsGeometry)(object, point);
}

function graticuleX(y0, y1, dy) {
  var y = d3Array.range(y0, y1 - epsilon, dy).concat(y1);
  return function(x) { return y.map(function(y) { return [x, y]; }); };
}

function graticuleY(x0, x1, dx) {
  var x = d3Array.range(x0, x1 - epsilon, dx).concat(x1);
  return function(y) { return x.map(function(x) { return [x, y]; }); };
}

function graticule() {
  var x1, x0, X1, X0,
      y1, y0, Y1, Y0,
      dx = 10, dy = dx, DX = 90, DY = 360,
      x, y, X, Y,
      precision = 2.5;

  function graticule() {
    return {type: "MultiLineString", coordinates: lines()};
  }

  function lines() {
    return d3Array.range(ceil(X0 / DX) * DX, X1, DX).map(X)
        .concat(d3Array.range(ceil(Y0 / DY) * DY, Y1, DY).map(Y))
        .concat(d3Array.range(ceil(x0 / dx) * dx, x1, dx).filter(function(x) { return abs(x % DX) > epsilon; }).map(x))
        .concat(d3Array.range(ceil(y0 / dy) * dy, y1, dy).filter(function(y) { return abs(y % DY) > epsilon; }).map(y));
  }

  graticule.lines = function() {
    return lines().map(function(coordinates) { return {type: "LineString", coordinates: coordinates}; });
  };

  graticule.outline = function() {
    return {
      type: "Polygon",
      coordinates: [
        X(X0).concat(
        Y(Y1).slice(1),
        X(X1).reverse().slice(1),
        Y(Y0).reverse().slice(1))
      ]
    };
  };

  graticule.extent = function(_) {
    if (!arguments.length) return graticule.extentMinor();
    return graticule.extentMajor(_).extentMinor(_);
  };

  graticule.extentMajor = function(_) {
    if (!arguments.length) return [[X0, Y0], [X1, Y1]];
    X0 = +_[0][0], X1 = +_[1][0];
    Y0 = +_[0][1], Y1 = +_[1][1];
    if (X0 > X1) _ = X0, X0 = X1, X1 = _;
    if (Y0 > Y1) _ = Y0, Y0 = Y1, Y1 = _;
    return graticule.precision(precision);
  };

  graticule.extentMinor = function(_) {
    if (!arguments.length) return [[x0, y0], [x1, y1]];
    x0 = +_[0][0], x1 = +_[1][0];
    y0 = +_[0][1], y1 = +_[1][1];
    if (x0 > x1) _ = x0, x0 = x1, x1 = _;
    if (y0 > y1) _ = y0, y0 = y1, y1 = _;
    return graticule.precision(precision);
  };

  graticule.step = function(_) {
    if (!arguments.length) return graticule.stepMinor();
    return graticule.stepMajor(_).stepMinor(_);
  };

  graticule.stepMajor = function(_) {
    if (!arguments.length) return [DX, DY];
    DX = +_[0], DY = +_[1];
    return graticule;
  };

  graticule.stepMinor = function(_) {
    if (!arguments.length) return [dx, dy];
    dx = +_[0], dy = +_[1];
    return graticule;
  };

  graticule.precision = function(_) {
    if (!arguments.length) return precision;
    precision = +_;
    x = graticuleX(y0, y1, 90);
    y = graticuleY(x0, x1, precision);
    X = graticuleX(Y0, Y1, 90);
    Y = graticuleY(X0, X1, precision);
    return graticule;
  };

  return graticule
      .extentMajor([[-180, -90 + epsilon], [180, 90 - epsilon]])
      .extentMinor([[-180, -80 - epsilon], [180, 80 + epsilon]]);
}

function graticule10() {
  return graticule()();
}

function interpolate(a, b) {
  var x0 = a[0] * radians,
      y0 = a[1] * radians,
      x1 = b[0] * radians,
      y1 = b[1] * radians,
      cy0 = cos(y0),
      sy0 = sin(y0),
      cy1 = cos(y1),
      sy1 = sin(y1),
      kx0 = cy0 * cos(x0),
      ky0 = cy0 * sin(x0),
      kx1 = cy1 * cos(x1),
      ky1 = cy1 * sin(x1),
      d = 2 * asin(sqrt(haversin(y1 - y0) + cy0 * cy1 * haversin(x1 - x0))),
      k = sin(d);

  var interpolate = d ? function(t) {
    var B = sin(t *= d) / k,
        A = sin(d - t) / k,
        x = A * kx0 + B * kx1,
        y = A * ky0 + B * ky1,
        z = A * sy0 + B * sy1;
    return [
      atan2(y, x) * degrees,
      atan2(z, sqrt(x * x + y * y)) * degrees
    ];
  } : function() {
    return [x0 * degrees, y0 * degrees];
  };

  interpolate.distance = d;

  return interpolate;
}

function identity(x) {
  return x;
}

var areaSum$1 = adder(),
    areaRingSum$1 = adder(),
    x00,
    y00,
    x0$1,
    y0$1;

var areaStream$1 = {
  point: noop,
  lineStart: noop,
  lineEnd: noop,
  polygonStart: function() {
    areaStream$1.lineStart = areaRingStart$1;
    areaStream$1.lineEnd = areaRingEnd$1;
  },
  polygonEnd: function() {
    areaStream$1.lineStart = areaStream$1.lineEnd = areaStream$1.point = noop;
    areaSum$1.add(abs(areaRingSum$1));
    areaRingSum$1.reset();
  },
  result: function() {
    var area = areaSum$1 / 2;
    areaSum$1.reset();
    return area;
  }
};

function areaRingStart$1() {
  areaStream$1.point = areaPointFirst$1;
}

function areaPointFirst$1(x, y) {
  areaStream$1.point = areaPoint$1;
  x00 = x0$1 = x, y00 = y0$1 = y;
}

function areaPoint$1(x, y) {
  areaRingSum$1.add(y0$1 * x - x0$1 * y);
  x0$1 = x, y0$1 = y;
}

function areaRingEnd$1() {
  areaPoint$1(x00, y00);
}

var x0$2 = Infinity,
    y0$2 = x0$2,
    x1 = -x0$2,
    y1 = x1;

var boundsStream$1 = {
  point: boundsPoint$1,
  lineStart: noop,
  lineEnd: noop,
  polygonStart: noop,
  polygonEnd: noop,
  result: function() {
    var bounds = [[x0$2, y0$2], [x1, y1]];
    x1 = y1 = -(y0$2 = x0$2 = Infinity);
    return bounds;
  }
};

function boundsPoint$1(x, y) {
  if (x < x0$2) x0$2 = x;
  if (x > x1) x1 = x;
  if (y < y0$2) y0$2 = y;
  if (y > y1) y1 = y;
}

// TODO Enforce positive area for exterior, negative area for interior?

var X0$1 = 0,
    Y0$1 = 0,
    Z0$1 = 0,
    X1$1 = 0,
    Y1$1 = 0,
    Z1$1 = 0,
    X2$1 = 0,
    Y2$1 = 0,
    Z2$1 = 0,
    x00$1,
    y00$1,
    x0$3,
    y0$3;

var centroidStream$1 = {
  point: centroidPoint$1,
  lineStart: centroidLineStart$1,
  lineEnd: centroidLineEnd$1,
  polygonStart: function() {
    centroidStream$1.lineStart = centroidRingStart$1;
    centroidStream$1.lineEnd = centroidRingEnd$1;
  },
  polygonEnd: function() {
    centroidStream$1.point = centroidPoint$1;
    centroidStream$1.lineStart = centroidLineStart$1;
    centroidStream$1.lineEnd = centroidLineEnd$1;
  },
  result: function() {
    var centroid = Z2$1 ? [X2$1 / Z2$1, Y2$1 / Z2$1]
        : Z1$1 ? [X1$1 / Z1$1, Y1$1 / Z1$1]
        : Z0$1 ? [X0$1 / Z0$1, Y0$1 / Z0$1]
        : [NaN, NaN];
    X0$1 = Y0$1 = Z0$1 =
    X1$1 = Y1$1 = Z1$1 =
    X2$1 = Y2$1 = Z2$1 = 0;
    return centroid;
  }
};

function centroidPoint$1(x, y) {
  X0$1 += x;
  Y0$1 += y;
  ++Z0$1;
}

function centroidLineStart$1() {
  centroidStream$1.point = centroidPointFirstLine;
}

function centroidPointFirstLine(x, y) {
  centroidStream$1.point = centroidPointLine;
  centroidPoint$1(x0$3 = x, y0$3 = y);
}

function centroidPointLine(x, y) {
  var dx = x - x0$3, dy = y - y0$3, z = sqrt(dx * dx + dy * dy);
  X1$1 += z * (x0$3 + x) / 2;
  Y1$1 += z * (y0$3 + y) / 2;
  Z1$1 += z;
  centroidPoint$1(x0$3 = x, y0$3 = y);
}

function centroidLineEnd$1() {
  centroidStream$1.point = centroidPoint$1;
}

function centroidRingStart$1() {
  centroidStream$1.point = centroidPointFirstRing;
}

function centroidRingEnd$1() {
  centroidPointRing(x00$1, y00$1);
}

function centroidPointFirstRing(x, y) {
  centroidStream$1.point = centroidPointRing;
  centroidPoint$1(x00$1 = x0$3 = x, y00$1 = y0$3 = y);
}

function centroidPointRing(x, y) {
  var dx = x - x0$3,
      dy = y - y0$3,
      z = sqrt(dx * dx + dy * dy);

  X1$1 += z * (x0$3 + x) / 2;
  Y1$1 += z * (y0$3 + y) / 2;
  Z1$1 += z;

  z = y0$3 * x - x0$3 * y;
  X2$1 += z * (x0$3 + x);
  Y2$1 += z * (y0$3 + y);
  Z2$1 += z * 3;
  centroidPoint$1(x0$3 = x, y0$3 = y);
}

function PathContext(context) {
  this._context = context;
}

PathContext.prototype = {
  _radius: 4.5,
  pointRadius: function(_) {
    return this._radius = _, this;
  },
  polygonStart: function() {
    this._line = 0;
  },
  polygonEnd: function() {
    this._line = NaN;
  },
  lineStart: function() {
    this._point = 0;
  },
  lineEnd: function() {
    if (this._line === 0) this._context.closePath();
    this._point = NaN;
  },
  point: function(x, y) {
    switch (this._point) {
      case 0: {
        this._context.moveTo(x, y);
        this._point = 1;
        break;
      }
      case 1: {
        this._context.lineTo(x, y);
        break;
      }
      default: {
        this._context.moveTo(x + this._radius, y);
        this._context.arc(x, y, this._radius, 0, tau);
        break;
      }
    }
  },
  result: noop
};

var lengthSum$1 = adder(),
    lengthRing,
    x00$2,
    y00$2,
    x0$4,
    y0$4;

var lengthStream$1 = {
  point: noop,
  lineStart: function() {
    lengthStream$1.point = lengthPointFirst$1;
  },
  lineEnd: function() {
    if (lengthRing) lengthPoint$1(x00$2, y00$2);
    lengthStream$1.point = noop;
  },
  polygonStart: function() {
    lengthRing = true;
  },
  polygonEnd: function() {
    lengthRing = null;
  },
  result: function() {
    var length = +lengthSum$1;
    lengthSum$1.reset();
    return length;
  }
};

function lengthPointFirst$1(x, y) {
  lengthStream$1.point = lengthPoint$1;
  x00$2 = x0$4 = x, y00$2 = y0$4 = y;
}

function lengthPoint$1(x, y) {
  x0$4 -= x, y0$4 -= y;
  lengthSum$1.add(sqrt(x0$4 * x0$4 + y0$4 * y0$4));
  x0$4 = x, y0$4 = y;
}

function PathString() {
  this._string = [];
}

PathString.prototype = {
  _radius: 4.5,
  _circle: circle$1(4.5),
  pointRadius: function(_) {
    if ((_ = +_) !== this._radius) this._radius = _, this._circle = null;
    return this;
  },
  polygonStart: function() {
    this._line = 0;
  },
  polygonEnd: function() {
    this._line = NaN;
  },
  lineStart: function() {
    this._point = 0;
  },
  lineEnd: function() {
    if (this._line === 0) this._string.push("Z");
    this._point = NaN;
  },
  point: function(x, y) {
    switch (this._point) {
      case 0: {
        this._string.push("M", x, ",", y);
        this._point = 1;
        break;
      }
      case 1: {
        this._string.push("L", x, ",", y);
        break;
      }
      default: {
        if (this._circle == null) this._circle = circle$1(this._radius);
        this._string.push("M", x, ",", y, this._circle);
        break;
      }
    }
  },
  result: function() {
    if (this._string.length) {
      var result = this._string.join("");
      this._string = [];
      return result;
    } else {
      return null;
    }
  }
};

function circle$1(radius) {
  return "m0," + radius
      + "a" + radius + "," + radius + " 0 1,1 0," + -2 * radius
      + "a" + radius + "," + radius + " 0 1,1 0," + 2 * radius
      + "z";
}

function index(projection, context) {
  var pointRadius = 4.5,
      projectionStream,
      contextStream;

  function path(object) {
    if (object) {
      if (typeof pointRadius === "function") contextStream.pointRadius(+pointRadius.apply(this, arguments));
      geoStream(object, projectionStream(contextStream));
    }
    return contextStream.result();
  }

  path.area = function(object) {
    geoStream(object, projectionStream(areaStream$1));
    return areaStream$1.result();
  };

  path.measure = function(object) {
    geoStream(object, projectionStream(lengthStream$1));
    return lengthStream$1.result();
  };

  path.bounds = function(object) {
    geoStream(object, projectionStream(boundsStream$1));
    return boundsStream$1.result();
  };

  path.centroid = function(object) {
    geoStream(object, projectionStream(centroidStream$1));
    return centroidStream$1.result();
  };

  path.projection = function(_) {
    return arguments.length ? (projectionStream = _ == null ? (projection = null, identity) : (projection = _).stream, path) : projection;
  };

  path.context = function(_) {
    if (!arguments.length) return context;
    contextStream = _ == null ? (context = null, new PathString) : new PathContext(context = _);
    if (typeof pointRadius !== "function") contextStream.pointRadius(pointRadius);
    return path;
  };

  path.pointRadius = function(_) {
    if (!arguments.length) return pointRadius;
    pointRadius = typeof _ === "function" ? _ : (contextStream.pointRadius(+_), +_);
    return path;
  };

  return path.projection(projection).context(context);
}

function transform(methods) {
  return {
    stream: transformer(methods)
  };
}

function transformer(methods) {
  return function(stream) {
    var s = new TransformStream;
    for (var key in methods) s[key] = methods[key];
    s.stream = stream;
    return s;
  };
}

function TransformStream() {}

TransformStream.prototype = {
  constructor: TransformStream,
  point: function(x, y) { this.stream.point(x, y); },
  sphere: function() { this.stream.sphere(); },
  lineStart: function() { this.stream.lineStart(); },
  lineEnd: function() { this.stream.lineEnd(); },
  polygonStart: function() { this.stream.polygonStart(); },
  polygonEnd: function() { this.stream.polygonEnd(); }
};

function fit(projection, fitBounds, object) {
  var clip = projection.clipExtent && projection.clipExtent();
  projection.scale(150).translate([0, 0]);
  if (clip != null) projection.clipExtent(null);
  geoStream(object, projection.stream(boundsStream$1));
  fitBounds(boundsStream$1.result());
  if (clip != null) projection.clipExtent(clip);
  return projection;
}

function fitExtent(projection, extent, object) {
  return fit(projection, function(b) {
    var w = extent[1][0] - extent[0][0],
        h = extent[1][1] - extent[0][1],
        k = Math.min(w / (b[1][0] - b[0][0]), h / (b[1][1] - b[0][1])),
        x = +extent[0][0] + (w - k * (b[1][0] + b[0][0])) / 2,
        y = +extent[0][1] + (h - k * (b[1][1] + b[0][1])) / 2;
    projection.scale(150 * k).translate([x, y]);
  }, object);
}

function fitSize(projection, size, object) {
  return fitExtent(projection, [[0, 0], size], object);
}

function fitWidth(projection, width, object) {
  return fit(projection, function(b) {
    var w = +width,
        k = w / (b[1][0] - b[0][0]),
        x = (w - k * (b[1][0] + b[0][0])) / 2,
        y = -k * b[0][1];
    projection.scale(150 * k).translate([x, y]);
  }, object);
}

function fitHeight(projection, height, object) {
  return fit(projection, function(b) {
    var h = +height,
        k = h / (b[1][1] - b[0][1]),
        x = -k * b[0][0],
        y = (h - k * (b[1][1] + b[0][1])) / 2;
    projection.scale(150 * k).translate([x, y]);
  }, object);
}

var maxDepth = 16, // maximum depth of subdivision
    cosMinDistance = cos(30 * radians); // cos(minimum angular distance)

function resample(project, delta2) {
  return +delta2 ? resample$1(project, delta2) : resampleNone(project);
}

function resampleNone(project) {
  return transformer({
    point: function(x, y) {
      x = project(x, y);
      this.stream.point(x[0], x[1]);
    }
  });
}

function resample$1(project, delta2) {

  function resampleLineTo(x0, y0, lambda0, a0, b0, c0, x1, y1, lambda1, a1, b1, c1, depth, stream) {
    var dx = x1 - x0,
        dy = y1 - y0,
        d2 = dx * dx + dy * dy;
    if (d2 > 4 * delta2 && depth--) {
      var a = a0 + a1,
          b = b0 + b1,
          c = c0 + c1,
          m = sqrt(a * a + b * b + c * c),
          phi2 = asin(c /= m),
          lambda2 = abs(abs(c) - 1) < epsilon || abs(lambda0 - lambda1) < epsilon ? (lambda0 + lambda1) / 2 : atan2(b, a),
          p = project(lambda2, phi2),
          x2 = p[0],
          y2 = p[1],
          dx2 = x2 - x0,
          dy2 = y2 - y0,
          dz = dy * dx2 - dx * dy2;
      if (dz * dz / d2 > delta2 // perpendicular projected distance
          || abs((dx * dx2 + dy * dy2) / d2 - 0.5) > 0.3 // midpoint close to an end
          || a0 * a1 + b0 * b1 + c0 * c1 < cosMinDistance) { // angular distance
        resampleLineTo(x0, y0, lambda0, a0, b0, c0, x2, y2, lambda2, a /= m, b /= m, c, depth, stream);
        stream.point(x2, y2);
        resampleLineTo(x2, y2, lambda2, a, b, c, x1, y1, lambda1, a1, b1, c1, depth, stream);
      }
    }
  }
  return function(stream) {
    var lambda00, x00, y00, a00, b00, c00, // first point
        lambda0, x0, y0, a0, b0, c0; // previous point

    var resampleStream = {
      point: point,
      lineStart: lineStart,
      lineEnd: lineEnd,
      polygonStart: function() { stream.polygonStart(); resampleStream.lineStart = ringStart; },
      polygonEnd: function() { stream.polygonEnd(); resampleStream.lineStart = lineStart; }
    };

    function point(x, y) {
      x = project(x, y);
      stream.point(x[0], x[1]);
    }

    function lineStart() {
      x0 = NaN;
      resampleStream.point = linePoint;
      stream.lineStart();
    }

    function linePoint(lambda, phi) {
      var c = cartesian([lambda, phi]), p = project(lambda, phi);
      resampleLineTo(x0, y0, lambda0, a0, b0, c0, x0 = p[0], y0 = p[1], lambda0 = lambda, a0 = c[0], b0 = c[1], c0 = c[2], maxDepth, stream);
      stream.point(x0, y0);
    }

    function lineEnd() {
      resampleStream.point = point;
      stream.lineEnd();
    }

    function ringStart() {
      lineStart();
      resampleStream.point = ringPoint;
      resampleStream.lineEnd = ringEnd;
    }

    function ringPoint(lambda, phi) {
      linePoint(lambda00 = lambda, phi), x00 = x0, y00 = y0, a00 = a0, b00 = b0, c00 = c0;
      resampleStream.point = linePoint;
    }

    function ringEnd() {
      resampleLineTo(x0, y0, lambda0, a0, b0, c0, x00, y00, lambda00, a00, b00, c00, maxDepth, stream);
      resampleStream.lineEnd = lineEnd;
      lineEnd();
    }

    return resampleStream;
  };
}

var transformRadians = transformer({
  point: function(x, y) {
    this.stream.point(x * radians, y * radians);
  }
});

function transformRotate(rotate) {
  return transformer({
    point: function(x, y) {
      var r = rotate(x, y);
      return this.stream.point(r[0], r[1]);
    }
  });
}

function scaleTranslate(k, dx, dy) {
  function transform$$1(x, y) {
    return [dx + k * x, dy - k * y];
  }
  transform$$1.invert = function(x, y) {
    return [(x - dx) / k, (dy - y) / k];
  };
  return transform$$1;
}

function scaleTranslateRotate(k, dx, dy, alpha) {
  var cosAlpha = cos(alpha),
      sinAlpha = sin(alpha),
      a = cosAlpha * k,
      b = sinAlpha * k,
      ai = cosAlpha / k,
      bi = sinAlpha / k,
      ci = (sinAlpha * dy - cosAlpha * dx) / k,
      fi = (sinAlpha * dx + cosAlpha * dy) / k;
  function transform$$1(x, y) {
    return [a * x - b * y + dx, dy - b * x - a * y];
  }
  transform$$1.invert = function(x, y) {
    return [ai * x - bi * y + ci, fi - bi * x - ai * y];
  };
  return transform$$1;
}

function projection(project) {
  return projectionMutator(function() { return project; })();
}

function projectionMutator(projectAt) {
  var project,
      k = 150, // scale
      x = 480, y = 250, // translate
      lambda = 0, phi = 0, // center
      deltaLambda = 0, deltaPhi = 0, deltaGamma = 0, rotate, // pre-rotate
      alpha = 0, // post-rotate
      theta = null, preclip = clipAntimeridian, // pre-clip angle
      x0 = null, y0, x1, y1, postclip = identity, // post-clip extent
      delta2 = 0.5, // precision
      projectResample,
      projectTransform,
      projectRotateTransform,
      cache,
      cacheStream;

  function projection(point) {
    return projectRotateTransform(point[0] * radians, point[1] * radians);
  }

  function invert(point) {
    point = projectRotateTransform.invert(point[0], point[1]);
    return point && [point[0] * degrees, point[1] * degrees];
  }

  projection.stream = function(stream) {
    return cache && cacheStream === stream ? cache : cache = transformRadians(transformRotate(rotate)(preclip(projectResample(postclip(cacheStream = stream)))));
  };

  projection.preclip = function(_) {
    return arguments.length ? (preclip = _, theta = undefined, reset()) : preclip;
  };

  projection.postclip = function(_) {
    return arguments.length ? (postclip = _, x0 = y0 = x1 = y1 = null, reset()) : postclip;
  };

  projection.clipAngle = function(_) {
    return arguments.length ? (preclip = +_ ? clipCircle(theta = _ * radians) : (theta = null, clipAntimeridian), reset()) : theta * degrees;
  };

  projection.clipExtent = function(_) {
    return arguments.length ? (postclip = _ == null ? (x0 = y0 = x1 = y1 = null, identity) : clipRectangle(x0 = +_[0][0], y0 = +_[0][1], x1 = +_[1][0], y1 = +_[1][1]), reset()) : x0 == null ? null : [[x0, y0], [x1, y1]];
  };

  projection.scale = function(_) {
    return arguments.length ? (k = +_, recenter()) : k;
  };

  projection.translate = function(_) {
    return arguments.length ? (x = +_[0], y = +_[1], recenter()) : [x, y];
  };

  projection.center = function(_) {
    return arguments.length ? (lambda = _[0] % 360 * radians, phi = _[1] % 360 * radians, recenter()) : [lambda * degrees, phi * degrees];
  };

  projection.rotate = function(_) {
    return arguments.length ? (deltaLambda = _[0] % 360 * radians, deltaPhi = _[1] % 360 * radians, deltaGamma = _.length > 2 ? _[2] % 360 * radians : 0, recenter()) : [deltaLambda * degrees, deltaPhi * degrees, deltaGamma * degrees];
  };

  projection.angle = function(_) {
    return arguments.length ? (alpha = _ % 360 * radians, recenter()) : alpha * degrees;
  };

  projection.precision = function(_) {
    return arguments.length ? (projectResample = resample(projectTransform, delta2 = _ * _), reset()) : sqrt(delta2);
  };

  projection.fitExtent = function(extent, object) {
    return fitExtent(projection, extent, object);
  };

  projection.fitSize = function(size, object) {
    return fitSize(projection, size, object);
  };

  projection.fitWidth = function(width, object) {
    return fitWidth(projection, width, object);
  };

  projection.fitHeight = function(height, object) {
    return fitHeight(projection, height, object);
  };

  function recenter() {
    var center = scaleTranslateRotate(k, 0, 0, alpha).apply(null, project(lambda, phi)),
        transform$$1 = (alpha ? scaleTranslateRotate : scaleTranslate)(k, x - center[0], y - center[1], alpha);
    rotate = rotateRadians(deltaLambda, deltaPhi, deltaGamma);
    projectTransform = compose(project, transform$$1);
    projectRotateTransform = compose(rotate, projectTransform);
    projectResample = resample(projectTransform, delta2);
    return reset();
  }

  function reset() {
    cache = cacheStream = null;
    return projection;
  }

  return function() {
    project = projectAt.apply(this, arguments);
    projection.invert = project.invert && invert;
    return recenter();
  };
}

function conicProjection(projectAt) {
  var phi0 = 0,
      phi1 = pi / 3,
      m = projectionMutator(projectAt),
      p = m(phi0, phi1);

  p.parallels = function(_) {
    return arguments.length ? m(phi0 = _[0] * radians, phi1 = _[1] * radians) : [phi0 * degrees, phi1 * degrees];
  };

  return p;
}

function cylindricalEqualAreaRaw(phi0) {
  var cosPhi0 = cos(phi0);

  function forward(lambda, phi) {
    return [lambda * cosPhi0, sin(phi) / cosPhi0];
  }

  forward.invert = function(x, y) {
    return [x / cosPhi0, asin(y * cosPhi0)];
  };

  return forward;
}

function conicEqualAreaRaw(y0, y1) {
  var sy0 = sin(y0), n = (sy0 + sin(y1)) / 2;

  // Are the parallels symmetrical around the Equator?
  if (abs(n) < epsilon) return cylindricalEqualAreaRaw(y0);

  var c = 1 + sy0 * (2 * n - sy0), r0 = sqrt(c) / n;

  function project(x, y) {
    var r = sqrt(c - 2 * n * sin(y)) / n;
    return [r * sin(x *= n), r0 - r * cos(x)];
  }

  project.invert = function(x, y) {
    var r0y = r0 - y;
    return [atan2(x, abs(r0y)) / n * sign(r0y), asin((c - (x * x + r0y * r0y) * n * n) / (2 * n))];
  };

  return project;
}

function conicEqualArea() {
  return conicProjection(conicEqualAreaRaw)
      .scale(155.424)
      .center([0, 33.6442]);
}

function albers() {
  return conicEqualArea()
      .parallels([29.5, 45.5])
      .scale(1070)
      .translate([480, 250])
      .rotate([96, 0])
      .center([-0.6, 38.7]);
}

// The projections must have mutually exclusive clip regions on the sphere,
// as this will avoid emitting interleaving lines and polygons.
function multiplex(streams) {
  var n = streams.length;
  return {
    point: function(x, y) { var i = -1; while (++i < n) streams[i].point(x, y); },
    sphere: function() { var i = -1; while (++i < n) streams[i].sphere(); },
    lineStart: function() { var i = -1; while (++i < n) streams[i].lineStart(); },
    lineEnd: function() { var i = -1; while (++i < n) streams[i].lineEnd(); },
    polygonStart: function() { var i = -1; while (++i < n) streams[i].polygonStart(); },
    polygonEnd: function() { var i = -1; while (++i < n) streams[i].polygonEnd(); }
  };
}

// A composite projection for the United States, configured by default for
// 960×500. The projection also works quite well at 960×600 if you change the
// scale to 1285 and adjust the translate accordingly. The set of standard
// parallels for each region comes from USGS, which is published here:
// http://egsc.usgs.gov/isb/pubs/MapProjections/projections.html#albers
function albersUsa() {
  var cache,
      cacheStream,
      lower48 = albers(), lower48Point,
      alaska = conicEqualArea().rotate([154, 0]).center([-2, 58.5]).parallels([55, 65]), alaskaPoint, // EPSG:3338
      hawaii = conicEqualArea().rotate([157, 0]).center([-3, 19.9]).parallels([8, 18]), hawaiiPoint, // ESRI:102007
      point, pointStream = {point: function(x, y) { point = [x, y]; }};

  function albersUsa(coordinates) {
    var x = coordinates[0], y = coordinates[1];
    return point = null,
        (lower48Point.point(x, y), point)
        || (alaskaPoint.point(x, y), point)
        || (hawaiiPoint.point(x, y), point);
  }

  albersUsa.invert = function(coordinates) {
    var k = lower48.scale(),
        t = lower48.translate(),
        x = (coordinates[0] - t[0]) / k,
        y = (coordinates[1] - t[1]) / k;
    return (y >= 0.120 && y < 0.234 && x >= -0.425 && x < -0.214 ? alaska
        : y >= 0.166 && y < 0.234 && x >= -0.214 && x < -0.115 ? hawaii
        : lower48).invert(coordinates);
  };

  albersUsa.stream = function(stream) {
    return cache && cacheStream === stream ? cache : cache = multiplex([lower48.stream(cacheStream = stream), alaska.stream(stream), hawaii.stream(stream)]);
  };

  albersUsa.precision = function(_) {
    if (!arguments.length) return lower48.precision();
    lower48.precision(_), alaska.precision(_), hawaii.precision(_);
    return reset();
  };

  albersUsa.scale = function(_) {
    if (!arguments.length) return lower48.scale();
    lower48.scale(_), alaska.scale(_ * 0.35), hawaii.scale(_);
    return albersUsa.translate(lower48.translate());
  };

  albersUsa.translate = function(_) {
    if (!arguments.length) return lower48.translate();
    var k = lower48.scale(), x = +_[0], y = +_[1];

    lower48Point = lower48
        .translate(_)
        .clipExtent([[x - 0.455 * k, y - 0.238 * k], [x + 0.455 * k, y + 0.238 * k]])
        .stream(pointStream);

    alaskaPoint = alaska
        .translate([x - 0.307 * k, y + 0.201 * k])
        .clipExtent([[x - 0.425 * k + epsilon, y + 0.120 * k + epsilon], [x - 0.214 * k - epsilon, y + 0.234 * k - epsilon]])
        .stream(pointStream);

    hawaiiPoint = hawaii
        .translate([x - 0.205 * k, y + 0.212 * k])
        .clipExtent([[x - 0.214 * k + epsilon, y + 0.166 * k + epsilon], [x - 0.115 * k - epsilon, y + 0.234 * k - epsilon]])
        .stream(pointStream);

    return reset();
  };

  albersUsa.fitExtent = function(extent, object) {
    return fitExtent(albersUsa, extent, object);
  };

  albersUsa.fitSize = function(size, object) {
    return fitSize(albersUsa, size, object);
  };

  albersUsa.fitWidth = function(width, object) {
    return fitWidth(albersUsa, width, object);
  };

  albersUsa.fitHeight = function(height, object) {
    return fitHeight(albersUsa, height, object);
  };

  function reset() {
    cache = cacheStream = null;
    return albersUsa;
  }

  return albersUsa.scale(1070);
}

function azimuthalRaw(scale) {
  return function(x, y) {
    var cx = cos(x),
        cy = cos(y),
        k = scale(cx * cy);
    return [
      k * cy * sin(x),
      k * sin(y)
    ];
  }
}

function azimuthalInvert(angle) {
  return function(x, y) {
    var z = sqrt(x * x + y * y),
        c = angle(z),
        sc = sin(c),
        cc = cos(c);
    return [
      atan2(x * sc, z * cc),
      asin(z && y * sc / z)
    ];
  }
}

var azimuthalEqualAreaRaw = azimuthalRaw(function(cxcy) {
  return sqrt(2 / (1 + cxcy));
});

azimuthalEqualAreaRaw.invert = azimuthalInvert(function(z) {
  return 2 * asin(z / 2);
});

function azimuthalEqualArea() {
  return projection(azimuthalEqualAreaRaw)
      .scale(124.75)
      .clipAngle(180 - 1e-3);
}

var azimuthalEquidistantRaw = azimuthalRaw(function(c) {
  return (c = acos(c)) && c / sin(c);
});

azimuthalEquidistantRaw.invert = azimuthalInvert(function(z) {
  return z;
});

function azimuthalEquidistant() {
  return projection(azimuthalEquidistantRaw)
      .scale(79.4188)
      .clipAngle(180 - 1e-3);
}

function mercatorRaw(lambda, phi) {
  return [lambda, log(tan((halfPi + phi) / 2))];
}

mercatorRaw.invert = function(x, y) {
  return [x, 2 * atan(exp(y)) - halfPi];
};

function mercator() {
  return mercatorProjection(mercatorRaw)
      .scale(961 / tau);
}

function mercatorProjection(project) {
  var m = projection(project),
      center = m.center,
      scale = m.scale,
      translate = m.translate,
      clipExtent = m.clipExtent,
      x0 = null, y0, x1, y1; // clip extent

  m.scale = function(_) {
    return arguments.length ? (scale(_), reclip()) : scale();
  };

  m.translate = function(_) {
    return arguments.length ? (translate(_), reclip()) : translate();
  };

  m.center = function(_) {
    return arguments.length ? (center(_), reclip()) : center();
  };

  m.clipExtent = function(_) {
    return arguments.length ? ((_ == null ? x0 = y0 = x1 = y1 = null : (x0 = +_[0][0], y0 = +_[0][1], x1 = +_[1][0], y1 = +_[1][1])), reclip()) : x0 == null ? null : [[x0, y0], [x1, y1]];
  };

  function reclip() {
    var k = pi * scale(),
        t = m(rotation(m.rotate()).invert([0, 0]));
    return clipExtent(x0 == null
        ? [[t[0] - k, t[1] - k], [t[0] + k, t[1] + k]] : project === mercatorRaw
        ? [[Math.max(t[0] - k, x0), y0], [Math.min(t[0] + k, x1), y1]]
        : [[x0, Math.max(t[1] - k, y0)], [x1, Math.min(t[1] + k, y1)]]);
  }

  return reclip();
}

function tany(y) {
  return tan((halfPi + y) / 2);
}

function conicConformalRaw(y0, y1) {
  var cy0 = cos(y0),
      n = y0 === y1 ? sin(y0) : log(cy0 / cos(y1)) / log(tany(y1) / tany(y0)),
      f = cy0 * pow(tany(y0), n) / n;

  if (!n) return mercatorRaw;

  function project(x, y) {
    if (f > 0) { if (y < -halfPi + epsilon) y = -halfPi + epsilon; }
    else { if (y > halfPi - epsilon) y = halfPi - epsilon; }
    var r = f / pow(tany(y), n);
    return [r * sin(n * x), f - r * cos(n * x)];
  }

  project.invert = function(x, y) {
    var fy = f - y, r = sign(n) * sqrt(x * x + fy * fy);
    return [atan2(x, abs(fy)) / n * sign(fy), 2 * atan(pow(f / r, 1 / n)) - halfPi];
  };

  return project;
}

function conicConformal() {
  return conicProjection(conicConformalRaw)
      .scale(109.5)
      .parallels([30, 30]);
}

function equirectangularRaw(lambda, phi) {
  return [lambda, phi];
}

equirectangularRaw.invert = equirectangularRaw;

function equirectangular() {
  return projection(equirectangularRaw)
      .scale(152.63);
}

function conicEquidistantRaw(y0, y1) {
  var cy0 = cos(y0),
      n = y0 === y1 ? sin(y0) : (cy0 - cos(y1)) / (y1 - y0),
      g = cy0 / n + y0;

  if (abs(n) < epsilon) return equirectangularRaw;

  function project(x, y) {
    var gy = g - y, nx = n * x;
    return [gy * sin(nx), g - gy * cos(nx)];
  }

  project.invert = function(x, y) {
    var gy = g - y;
    return [atan2(x, abs(gy)) / n * sign(gy), g - sign(n) * sqrt(x * x + gy * gy)];
  };

  return project;
}

function conicEquidistant() {
  return conicProjection(conicEquidistantRaw)
      .scale(131.154)
      .center([0, 13.9389]);
}

var A1 = 1.340264,
    A2 = -0.081106,
    A3 = 0.000893,
    A4 = 0.003796,
    M = sqrt(3) / 2,
    iterations = 12;

function equalEarthRaw(lambda, phi) {
  var l = asin(M * sin(phi)), l2 = l * l, l6 = l2 * l2 * l2;
  return [
    lambda * cos(l) / (M * (A1 + 3 * A2 * l2 + l6 * (7 * A3 + 9 * A4 * l2))),
    l * (A1 + A2 * l2 + l6 * (A3 + A4 * l2))
  ];
}

equalEarthRaw.invert = function(x, y) {
  var l = y, l2 = l * l, l6 = l2 * l2 * l2;
  for (var i = 0, delta, fy, fpy; i < iterations; ++i) {
    fy = l * (A1 + A2 * l2 + l6 * (A3 + A4 * l2)) - y;
    fpy = A1 + 3 * A2 * l2 + l6 * (7 * A3 + 9 * A4 * l2);
    l -= delta = fy / fpy, l2 = l * l, l6 = l2 * l2 * l2;
    if (abs(delta) < epsilon2) break;
  }
  return [
    M * x * (A1 + 3 * A2 * l2 + l6 * (7 * A3 + 9 * A4 * l2)) / cos(l),
    asin(sin(l) / M)
  ];
};

function equalEarth() {
  return projection(equalEarthRaw)
      .scale(177.158);
}

function gnomonicRaw(x, y) {
  var cy = cos(y), k = cos(x) * cy;
  return [cy * sin(x) / k, sin(y) / k];
}

gnomonicRaw.invert = azimuthalInvert(atan);

function gnomonic() {
  return projection(gnomonicRaw)
      .scale(144.049)
      .clipAngle(60);
}

function scaleTranslate$1(kx, ky, tx, ty) {
  return kx === 1 && ky === 1 && tx === 0 && ty === 0 ? identity : transformer({
    point: function(x, y) {
      this.stream.point(x * kx + tx, y * ky + ty);
    }
  });
}

function identity$1() {
  var k = 1, tx = 0, ty = 0, sx = 1, sy = 1, transform$$1 = identity, // scale, translate and reflect
      x0 = null, y0, x1, y1, // clip extent
      postclip = identity,
      cache,
      cacheStream,
      projection;

  function reset() {
    cache = cacheStream = null;
    return projection;
  }

  return projection = {
    stream: function(stream) {
      return cache && cacheStream === stream ? cache : cache = transform$$1(postclip(cacheStream = stream));
    },
    postclip: function(_) {
      return arguments.length ? (postclip = _, x0 = y0 = x1 = y1 = null, reset()) : postclip;
    },
    clipExtent: function(_) {
      return arguments.length ? (postclip = _ == null ? (x0 = y0 = x1 = y1 = null, identity) : clipRectangle(x0 = +_[0][0], y0 = +_[0][1], x1 = +_[1][0], y1 = +_[1][1]), reset()) : x0 == null ? null : [[x0, y0], [x1, y1]];
    },
    scale: function(_) {
      return arguments.length ? (transform$$1 = scaleTranslate$1((k = +_) * sx, k * sy, tx, ty), reset()) : k;
    },
    translate: function(_) {
      return arguments.length ? (transform$$1 = scaleTranslate$1(k * sx, k * sy, tx = +_[0], ty = +_[1]), reset()) : [tx, ty];
    },
    reflectX: function(_) {
      return arguments.length ? (transform$$1 = scaleTranslate$1(k * (sx = _ ? -1 : 1), k * sy, tx, ty), reset()) : sx < 0;
    },
    reflectY: function(_) {
      return arguments.length ? (transform$$1 = scaleTranslate$1(k * sx, k * (sy = _ ? -1 : 1), tx, ty), reset()) : sy < 0;
    },
    fitExtent: function(extent, object) {
      return fitExtent(projection, extent, object);
    },
    fitSize: function(size, object) {
      return fitSize(projection, size, object);
    },
    fitWidth: function(width, object) {
      return fitWidth(projection, width, object);
    },
    fitHeight: function(height, object) {
      return fitHeight(projection, height, object);
    }
  };
}

function naturalEarth1Raw(lambda, phi) {
  var phi2 = phi * phi, phi4 = phi2 * phi2;
  return [
    lambda * (0.8707 - 0.131979 * phi2 + phi4 * (-0.013791 + phi4 * (0.003971 * phi2 - 0.001529 * phi4))),
    phi * (1.007226 + phi2 * (0.015085 + phi4 * (-0.044475 + 0.028874 * phi2 - 0.005916 * phi4)))
  ];
}

naturalEarth1Raw.invert = function(x, y) {
  var phi = y, i = 25, delta;
  do {
    var phi2 = phi * phi, phi4 = phi2 * phi2;
    phi -= delta = (phi * (1.007226 + phi2 * (0.015085 + phi4 * (-0.044475 + 0.028874 * phi2 - 0.005916 * phi4))) - y) /
        (1.007226 + phi2 * (0.015085 * 3 + phi4 * (-0.044475 * 7 + 0.028874 * 9 * phi2 - 0.005916 * 11 * phi4)));
  } while (abs(delta) > epsilon && --i > 0);
  return [
    x / (0.8707 + (phi2 = phi * phi) * (-0.131979 + phi2 * (-0.013791 + phi2 * phi2 * phi2 * (0.003971 - 0.001529 * phi2)))),
    phi
  ];
};

function naturalEarth1() {
  return projection(naturalEarth1Raw)
      .scale(175.295);
}

function orthographicRaw(x, y) {
  return [cos(y) * sin(x), sin(y)];
}

orthographicRaw.invert = azimuthalInvert(asin);

function orthographic() {
  return projection(orthographicRaw)
      .scale(249.5)
      .clipAngle(90 + epsilon);
}

function stereographicRaw(x, y) {
  var cy = cos(y), k = 1 + cos(x) * cy;
  return [cy * sin(x) / k, sin(y) / k];
}

stereographicRaw.invert = azimuthalInvert(function(z) {
  return 2 * atan(z);
});

function stereographic() {
  return projection(stereographicRaw)
      .scale(250)
      .clipAngle(142);
}

function transverseMercatorRaw(lambda, phi) {
  return [log(tan((halfPi + phi) / 2)), -lambda];
}

transverseMercatorRaw.invert = function(x, y) {
  return [-y, 2 * atan(exp(x)) - halfPi];
};

function transverseMercator() {
  var m = mercatorProjection(transverseMercatorRaw),
      center = m.center,
      rotate = m.rotate;

  m.center = function(_) {
    return arguments.length ? center([-_[1], _[0]]) : (_ = center(), [_[1], -_[0]]);
  };

  m.rotate = function(_) {
    return arguments.length ? rotate([_[0], _[1], _.length > 2 ? _[2] + 90 : 90]) : (_ = rotate(), [_[0], _[1], _[2] - 90]);
  };

  return rotate([0, 0, 90])
      .scale(159.155);
}

exports.geoArea = area;
exports.geoBounds = bounds;
exports.geoCentroid = centroid;
exports.geoCircle = circle;
exports.geoClipAntimeridian = clipAntimeridian;
exports.geoClipCircle = clipCircle;
exports.geoClipExtent = extent;
exports.geoClipRectangle = clipRectangle;
exports.geoContains = contains;
exports.geoDistance = distance;
exports.geoGraticule = graticule;
exports.geoGraticule10 = graticule10;
exports.geoInterpolate = interpolate;
exports.geoLength = length;
exports.geoPath = index;
exports.geoAlbers = albers;
exports.geoAlbersUsa = albersUsa;
exports.geoAzimuthalEqualArea = azimuthalEqualArea;
exports.geoAzimuthalEqualAreaRaw = azimuthalEqualAreaRaw;
exports.geoAzimuthalEquidistant = azimuthalEquidistant;
exports.geoAzimuthalEquidistantRaw = azimuthalEquidistantRaw;
exports.geoConicConformal = conicConformal;
exports.geoConicConformalRaw = conicConformalRaw;
exports.geoConicEqualArea = conicEqualArea;
exports.geoConicEqualAreaRaw = conicEqualAreaRaw;
exports.geoConicEquidistant = conicEquidistant;
exports.geoConicEquidistantRaw = conicEquidistantRaw;
exports.geoEqualEarth = equalEarth;
exports.geoEqualEarthRaw = equalEarthRaw;
exports.geoEquirectangular = equirectangular;
exports.geoEquirectangularRaw = equirectangularRaw;
exports.geoGnomonic = gnomonic;
exports.geoGnomonicRaw = gnomonicRaw;
exports.geoIdentity = identity$1;
exports.geoProjection = projection;
exports.geoProjectionMutator = projectionMutator;
exports.geoMercator = mercator;
exports.geoMercatorRaw = mercatorRaw;
exports.geoNaturalEarth1 = naturalEarth1;
exports.geoNaturalEarth1Raw = naturalEarth1Raw;
exports.geoOrthographic = orthographic;
exports.geoOrthographicRaw = orthographicRaw;
exports.geoStereographic = stereographic;
exports.geoStereographicRaw = stereographicRaw;
exports.geoTransverseMercator = transverseMercator;
exports.geoTransverseMercatorRaw = transverseMercatorRaw;
exports.geoRotation = rotation;
exports.geoStream = geoStream;
exports.geoTransform = transform;

Object.defineProperty(exports, '__esModule', { value: true });

})));

},{"d3-array":2}],4:[function(_dereq_,module,exports){
(function (global){
!function(e){if("object"==typeof exports&&"undefined"!=typeof module)module.exports=e();else if("function"==typeof define&&define.amd)define([],e);else{("undefined"!=typeof window?window:"undefined"!=typeof global?global:"undefined"!=typeof self?self:this).fitAspect=e()}}(function(){return function a(o,s,d){function c(t,e){if(!s[t]){if(!o[t]){var i="function"==typeof _dereq_&&_dereq_;if(!e&&i)return i(t,!0);if(l)return l(t,!0);var n=new Error("Cannot find module '"+t+"'");throw n.code="MODULE_NOT_FOUND",n}var r=s[t]={exports:{}};o[t][0].call(r.exports,function(e){return c(o[t][1][e]||e)},r,r.exports,a,o,s,d)}return s[t].exports}for(var l="function"==typeof _dereq_&&_dereq_,e=0;e<d.length;e++)c(d[e]);return c}({1:[function(e,t,i){"use strict";var n=[{names:["square","1:1","instagram"],description:"Square",decimal:1,orientation:"landscape"},{names:["4:3","fullscreen","four three","1.33:1","ipad","pythagorean"],description:"Traditional TVs",decimal:1.333333,orientation:"landscape"},{names:["a4","√2:1","paper","lichtenberg","1:1.41"],description:"A4 paper",decimal:1.41},{names:["imax","1.43:1"],description:"IMAX film",decimal:1.43,orientation:"landscape"},{names:["3:2","35mm","photo","1.5:1","1.5"],description:"35mm photos",decimal:1.5,orientation:"landscape"},{names:["business card","bank card","1.58:1"],description:"Bank Cards",decimal:1.58577,orientation:"landscape"},{names:["golden","kepler","1.618","1.6:1"],description:"Golden ratio",decimal:1.61803,orientation:"landscape"},{names:["16:9","hd","hdtv","fhd","tv","computer","iphone","4k","8k","1.78:1"],description:"HD video",decimal:1.77777,orientation:"landscape"},{names:["widescreen","1.85:1"],description:"Movie-theatres",decimal:1.85,orientation:"landscape"},{names:["2:1","univisium","mobile","18:9"],description:"2:1",decimal:2,orientation:"landscape"},{names:["cinemascope","widescreen","wide","2.35:1","2.39:1"],description:"Widescreen",decimal:2.35,orientation:"landscape"},{names:["silver","1 + √2","2.41:1"],description:"Silver ratio",decimal:2.41,orientation:"landscape"}],r=n.map(function(e){return(e=Object.assign({},e)).decimal=1/e.decimal,e.orientation="portrait",e}),a={};n.forEach(function(t){t.names.forEach(function(e){a[e]=t})}),t.exports={lookup:a,portraits:r,list:n}},{}],2:[function(e,t,i){"use strict";var n=e("./aspects");t.exports=function(e,t){var i=e/t;return(i=parseInt(100*i,10)/100)<1?function(e,t){for(var i=0;i<t.length;i+=1)if(e>t[i].decimal){if(t[i-1]){var n=Math.abs(e-t[i].decimal);if(Math.abs(e-t[i-1].decimal)<n)return t[i-1]}return t[i]}return t[t.length-1]}(i,n.portraits):function(e,t){for(var i=0;i<t.length;i+=1)if(e<=t[i].decimal){if(t[i-1]){var n=Math.abs(e-t[i].decimal);if(Math.abs(e-t[i-1].decimal)<n)return t[i-1]}return t[i]}return t[t.length-1]}(i,n.list)}},{"./aspects":1}],3:[function(e,t,i){"use strict";var n=function(e,t){var i=1/t.decimal,n=e.orientation||"landscape";"portrait"===n&&(i=1/i);var r=e.width*i;return r=Math.round(r),{closest:t,width:e.width,height:r,orientation:n,original:e}},r=function(e,t){var i=t.decimal,n=e.orientation||"landscape";"portrait"===n&&(i=1/i);var r=e.height*i;return{closest:t,width:r=Math.round(r),height:e.height,orientation:n,original:e}};t.exports={both:function(e,t){var i=r(e,t);return i.width>e.width?n(e,t):i},width:r,height:n}},{}],4:[function(i,n,e){(function(e){"use strict";var o=i("./find-best-ratio"),s=i("./parse-ratio"),d=i("./fit"),t=function(){var e=0<arguments.length&&void 0!==arguments[0]?arguments[0]:{};if(!e.aspect&&!e.ratio){var t=o(e.width,e.height),i=1/t.decimal,n=e.width*i,r=(n-e.height)/e.height;return r=parseInt(1e3*r,10)/10,n=Math.round(n),{closest:t,percent_change:r,width:e.width,height:n}}var a=s(e.aspect||e.ratio||"");return null===a?(console.error("find-aspect-ratio error: Could not find a given aspect ratio."),e):"number"==typeof e.width&&"number"==typeof e.height?d.both(e,a):"number"==typeof e.width?d.height(e,a):"number"==typeof e.height?d.width(e,a):(console.error("find-aspect-ratio error: Please supply a height, width, or ratio value."),e)};"undefined"!=typeof self?self.nlp=t:"undefined"!=typeof window?window.nlp=t:void 0!==e&&(e.nlp=t),void 0!==n&&(n.exports=t)}).call(this,"undefined"!=typeof global?global:"undefined"!=typeof self?self:"undefined"!=typeof window?window:{})},{"./find-best-ratio":2,"./fit":3,"./parse-ratio":5}],5:[function(e,t,i){"use strict";var n=e("./aspects"),r=/^[0-9\.]+:[0-9\.]+$/;t.exports=function(e){if(e=(e=(e=(e=e.toLowerCase()).trim()).replace(" ratio","")).replace("-"," "),!0===n.lookup.hasOwnProperty(e))return n.lookup[e];if(!0!==r.test(e))return null;var t=e.split(":");return{description:"custom",decimal:parseFloat(t[0])/parseFloat(t[1])}}},{"./aspects":1}]},{},[4])(4)});

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{}],5:[function(_dereq_,module,exports){
!function(){var n=function(t,e,r,u){for(var o=1;o<e.length;o++){var f=e[o++],s="number"==typeof f?r[f]:f;1===e[o]?u[0]=s:2===e[o]?(u[1]=u[1]||{})[e[++o]]=s:3===e[o]?u[1]=Object.assign(u[1]||{},s):u.push(e[o]?t.apply(null,n(t,s,r,["",null])):s)}return u},t=function(n){for(var t,e,r=1,u="",o="",f=[0],s=function(n){1===r&&(n||(u=u.replace(/^\s*\n\s*|\s*\n\s*$/g,"")))?f.push(n||u,0):3===r&&(n||u)?(f.push(n||u,1),r=2):2===r&&"..."===u&&n?f.push(n,3):2===r&&u&&!n?f.push(!0,2,u):4===r&&e&&(f.push(n||u,2,e),e=""),u=""},p=0;p<n.length;p++){p&&(1===r&&s(),s(p));for(var h=0;h<n[p].length;h++)t=n[p][h],1===r?"<"===t?(s(),f=[f],r=3):u+=t:o?t===o?o="":u+=t:'"'===t||"'"===t?o=t:">"===t?(s(),r=1):r&&("="===t?(r=4,e=u,u=""):"/"===t?(s(),3===r&&(f=f[0]),r=f,(f=f[0]).push(r,4),r=0):" "===t||"\t"===t||"\n"===t||"\r"===t?(s(),r=2):u+=t)}return s(),f},e="function"==typeof Map,r=e?new Map:{},u=e?function(n){var e=r.get(n);return e||r.set(n,e=t(n)),e}:function(n){for(var e="",u=0;u<n.length;u++)e+=n[u].length+"-"+n[u];return r[e]||(r[e]=t(n))},o=function(t){var e=n(this,u(t),arguments,[]);return e.length>1?e:e[0]};"undefined"!=typeof module?module.exports=o:self.htm=o}();

},{}],6:[function(_dereq_,module,exports){
(function (global){
!function(e){if("object"==typeof exports&&"undefined"!=typeof module)module.exports=e();else if("function"==typeof define&&define.amd)define([],e);else{("undefined"!=typeof window?window:"undefined"!=typeof global?global:"undefined"!=typeof self?self:this).spencerColor=e()}}(function(){return function u(i,a,c){function f(r,e){if(!a[r]){if(!i[r]){var o="function"==typeof _dereq_&&_dereq_;if(!e&&o)return o(r,!0);if(d)return d(r,!0);var n=new Error("Cannot find module '"+r+"'");throw n.code="MODULE_NOT_FOUND",n}var t=a[r]={exports:{}};i[r][0].call(t.exports,function(e){return f(i[r][1][e]||e)},t,t.exports,u,i,a,c)}return a[r].exports}for(var d="function"==typeof _dereq_&&_dereq_,e=0;e<c.length;e++)f(c[e]);return f}({1:[function(e,r,o){"use strict";r.exports={blue:"#6699cc",green:"#6accb2",yellow:"#e1e6b3",red:"#cc7066",pink:"#F2C0BB",brown:"#705E5C",orange:"#cc8a66",purple:"#d8b3e6",navy:"#335799",olive:"#7f9c6c",fuscia:"#735873",beige:"#e6d7b3",slate:"#8C8C88",suede:"#9c896c",burnt:"#603a39",sea:"#50617A",sky:"#2D85A8",night:"#303b50",rouge:"#914045",grey:"#838B91",mud:"#C4ABAB",royal:"#275291",cherry:"#cc6966",tulip:"#e6b3bc",rose:"#D68881",fire:"#AB5850",greyblue:"#72697D",greygreen:"#8BA3A2",greypurple:"#978BA3",burn:"#6D5685",slategrey:"#bfb0b3",light:"#a3a5a5",lighter:"#d7d5d2",fudge:"#4d4d4d",lightgrey:"#949a9e",white:"#fbfbfb",dimgrey:"#606c74",softblack:"#463D4F",dark:"#443d3d",black:"#333333"}},{}],2:[function(e,r,o){"use strict";var n=e("./colors"),t={juno:["blue","mud","navy","slate","pink","burn"],barrow:["rouge","red","orange","burnt","brown","greygreen"],roma:["#8a849a","#b5b0bf","rose","lighter","greygreen","mud"],palmer:["red","navy","olive","pink","suede","sky"],mark:["#848f9a","#9aa4ac","slate","#b0b8bf","mud","grey"],salmon:["sky","sea","fuscia","slate","mud","fudge"],dupont:["green","brown","orange","red","olive","blue"],bloor:["night","navy","beige","rouge","mud","grey"],yukon:["mud","slate","brown","sky","beige","red"],david:["blue","green","yellow","red","pink","light"],neste:["mud","cherry","royal","rouge","greygreen","greypurple"],ken:["red","sky","#c67a53","greygreen","#dfb59f","mud"]};Object.keys(t).forEach(function(e){t[e]=t[e].map(function(e){return n[e]||e})}),r.exports=t},{"./colors":1}],3:[function(e,r,o){"use strict";var n=e("./colors"),t=e("./combos"),u={colors:n,list:Object.keys(n).map(function(e){return n[e]}),combos:t};r.exports=u},{"./colors":1,"./combos":2}]},{},[3])(3)});

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{}],7:[function(_dereq_,module,exports){
// https://github.com/topojson/topojson-client Version 3.0.0. Copyright 2017 Mike Bostock.
(function (global, factory) {
	typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
	typeof define === 'function' && define.amd ? define(['exports'], factory) :
	(factory((global.topojson = global.topojson || {})));
}(this, (function (exports) { 'use strict';

var identity = function(x) {
  return x;
};

var transform = function(transform) {
  if (transform == null) return identity;
  var x0,
      y0,
      kx = transform.scale[0],
      ky = transform.scale[1],
      dx = transform.translate[0],
      dy = transform.translate[1];
  return function(input, i) {
    if (!i) x0 = y0 = 0;
    var j = 2, n = input.length, output = new Array(n);
    output[0] = (x0 += input[0]) * kx + dx;
    output[1] = (y0 += input[1]) * ky + dy;
    while (j < n) output[j] = input[j], ++j;
    return output;
  };
};

var bbox = function(topology) {
  var t = transform(topology.transform), key,
      x0 = Infinity, y0 = x0, x1 = -x0, y1 = -x0;

  function bboxPoint(p) {
    p = t(p);
    if (p[0] < x0) x0 = p[0];
    if (p[0] > x1) x1 = p[0];
    if (p[1] < y0) y0 = p[1];
    if (p[1] > y1) y1 = p[1];
  }

  function bboxGeometry(o) {
    switch (o.type) {
      case "GeometryCollection": o.geometries.forEach(bboxGeometry); break;
      case "Point": bboxPoint(o.coordinates); break;
      case "MultiPoint": o.coordinates.forEach(bboxPoint); break;
    }
  }

  topology.arcs.forEach(function(arc) {
    var i = -1, n = arc.length, p;
    while (++i < n) {
      p = t(arc[i], i);
      if (p[0] < x0) x0 = p[0];
      if (p[0] > x1) x1 = p[0];
      if (p[1] < y0) y0 = p[1];
      if (p[1] > y1) y1 = p[1];
    }
  });

  for (key in topology.objects) {
    bboxGeometry(topology.objects[key]);
  }

  return [x0, y0, x1, y1];
};

var reverse = function(array, n) {
  var t, j = array.length, i = j - n;
  while (i < --j) t = array[i], array[i++] = array[j], array[j] = t;
};

var feature = function(topology, o) {
  return o.type === "GeometryCollection"
      ? {type: "FeatureCollection", features: o.geometries.map(function(o) { return feature$1(topology, o); })}
      : feature$1(topology, o);
};

function feature$1(topology, o) {
  var id = o.id,
      bbox = o.bbox,
      properties = o.properties == null ? {} : o.properties,
      geometry = object(topology, o);
  return id == null && bbox == null ? {type: "Feature", properties: properties, geometry: geometry}
      : bbox == null ? {type: "Feature", id: id, properties: properties, geometry: geometry}
      : {type: "Feature", id: id, bbox: bbox, properties: properties, geometry: geometry};
}

function object(topology, o) {
  var transformPoint = transform(topology.transform),
      arcs = topology.arcs;

  function arc(i, points) {
    if (points.length) points.pop();
    for (var a = arcs[i < 0 ? ~i : i], k = 0, n = a.length; k < n; ++k) {
      points.push(transformPoint(a[k], k));
    }
    if (i < 0) reverse(points, n);
  }

  function point(p) {
    return transformPoint(p);
  }

  function line(arcs) {
    var points = [];
    for (var i = 0, n = arcs.length; i < n; ++i) arc(arcs[i], points);
    if (points.length < 2) points.push(points[0]); // This should never happen per the specification.
    return points;
  }

  function ring(arcs) {
    var points = line(arcs);
    while (points.length < 4) points.push(points[0]); // This may happen if an arc has only two points.
    return points;
  }

  function polygon(arcs) {
    return arcs.map(ring);
  }

  function geometry(o) {
    var type = o.type, coordinates;
    switch (type) {
      case "GeometryCollection": return {type: type, geometries: o.geometries.map(geometry)};
      case "Point": coordinates = point(o.coordinates); break;
      case "MultiPoint": coordinates = o.coordinates.map(point); break;
      case "LineString": coordinates = line(o.arcs); break;
      case "MultiLineString": coordinates = o.arcs.map(line); break;
      case "Polygon": coordinates = polygon(o.arcs); break;
      case "MultiPolygon": coordinates = o.arcs.map(polygon); break;
      default: return null;
    }
    return {type: type, coordinates: coordinates};
  }

  return geometry(o);
}

var stitch = function(topology, arcs) {
  var stitchedArcs = {},
      fragmentByStart = {},
      fragmentByEnd = {},
      fragments = [],
      emptyIndex = -1;

  // Stitch empty arcs first, since they may be subsumed by other arcs.
  arcs.forEach(function(i, j) {
    var arc = topology.arcs[i < 0 ? ~i : i], t;
    if (arc.length < 3 && !arc[1][0] && !arc[1][1]) {
      t = arcs[++emptyIndex], arcs[emptyIndex] = i, arcs[j] = t;
    }
  });

  arcs.forEach(function(i) {
    var e = ends(i),
        start = e[0],
        end = e[1],
        f, g;

    if (f = fragmentByEnd[start]) {
      delete fragmentByEnd[f.end];
      f.push(i);
      f.end = end;
      if (g = fragmentByStart[end]) {
        delete fragmentByStart[g.start];
        var fg = g === f ? f : f.concat(g);
        fragmentByStart[fg.start = f.start] = fragmentByEnd[fg.end = g.end] = fg;
      } else {
        fragmentByStart[f.start] = fragmentByEnd[f.end] = f;
      }
    } else if (f = fragmentByStart[end]) {
      delete fragmentByStart[f.start];
      f.unshift(i);
      f.start = start;
      if (g = fragmentByEnd[start]) {
        delete fragmentByEnd[g.end];
        var gf = g === f ? f : g.concat(f);
        fragmentByStart[gf.start = g.start] = fragmentByEnd[gf.end = f.end] = gf;
      } else {
        fragmentByStart[f.start] = fragmentByEnd[f.end] = f;
      }
    } else {
      f = [i];
      fragmentByStart[f.start = start] = fragmentByEnd[f.end = end] = f;
    }
  });

  function ends(i) {
    var arc = topology.arcs[i < 0 ? ~i : i], p0 = arc[0], p1;
    if (topology.transform) p1 = [0, 0], arc.forEach(function(dp) { p1[0] += dp[0], p1[1] += dp[1]; });
    else p1 = arc[arc.length - 1];
    return i < 0 ? [p1, p0] : [p0, p1];
  }

  function flush(fragmentByEnd, fragmentByStart) {
    for (var k in fragmentByEnd) {
      var f = fragmentByEnd[k];
      delete fragmentByStart[f.start];
      delete f.start;
      delete f.end;
      f.forEach(function(i) { stitchedArcs[i < 0 ? ~i : i] = 1; });
      fragments.push(f);
    }
  }

  flush(fragmentByEnd, fragmentByStart);
  flush(fragmentByStart, fragmentByEnd);
  arcs.forEach(function(i) { if (!stitchedArcs[i < 0 ? ~i : i]) fragments.push([i]); });

  return fragments;
};

var mesh = function(topology) {
  return object(topology, meshArcs.apply(this, arguments));
};

function meshArcs(topology, object$$1, filter) {
  var arcs, i, n;
  if (arguments.length > 1) arcs = extractArcs(topology, object$$1, filter);
  else for (i = 0, arcs = new Array(n = topology.arcs.length); i < n; ++i) arcs[i] = i;
  return {type: "MultiLineString", arcs: stitch(topology, arcs)};
}

function extractArcs(topology, object$$1, filter) {
  var arcs = [],
      geomsByArc = [],
      geom;

  function extract0(i) {
    var j = i < 0 ? ~i : i;
    (geomsByArc[j] || (geomsByArc[j] = [])).push({i: i, g: geom});
  }

  function extract1(arcs) {
    arcs.forEach(extract0);
  }

  function extract2(arcs) {
    arcs.forEach(extract1);
  }

  function extract3(arcs) {
    arcs.forEach(extract2);
  }

  function geometry(o) {
    switch (geom = o, o.type) {
      case "GeometryCollection": o.geometries.forEach(geometry); break;
      case "LineString": extract1(o.arcs); break;
      case "MultiLineString": case "Polygon": extract2(o.arcs); break;
      case "MultiPolygon": extract3(o.arcs); break;
    }
  }

  geometry(object$$1);

  geomsByArc.forEach(filter == null
      ? function(geoms) { arcs.push(geoms[0].i); }
      : function(geoms) { if (filter(geoms[0].g, geoms[geoms.length - 1].g)) arcs.push(geoms[0].i); });

  return arcs;
}

function planarRingArea(ring) {
  var i = -1, n = ring.length, a, b = ring[n - 1], area = 0;
  while (++i < n) a = b, b = ring[i], area += a[0] * b[1] - a[1] * b[0];
  return Math.abs(area); // Note: doubled area!
}

var merge = function(topology) {
  return object(topology, mergeArcs.apply(this, arguments));
};

function mergeArcs(topology, objects) {
  var polygonsByArc = {},
      polygons = [],
      groups = [];

  objects.forEach(geometry);

  function geometry(o) {
    switch (o.type) {
      case "GeometryCollection": o.geometries.forEach(geometry); break;
      case "Polygon": extract(o.arcs); break;
      case "MultiPolygon": o.arcs.forEach(extract); break;
    }
  }

  function extract(polygon) {
    polygon.forEach(function(ring) {
      ring.forEach(function(arc) {
        (polygonsByArc[arc = arc < 0 ? ~arc : arc] || (polygonsByArc[arc] = [])).push(polygon);
      });
    });
    polygons.push(polygon);
  }

  function area(ring) {
    return planarRingArea(object(topology, {type: "Polygon", arcs: [ring]}).coordinates[0]);
  }

  polygons.forEach(function(polygon) {
    if (!polygon._) {
      var group = [],
          neighbors = [polygon];
      polygon._ = 1;
      groups.push(group);
      while (polygon = neighbors.pop()) {
        group.push(polygon);
        polygon.forEach(function(ring) {
          ring.forEach(function(arc) {
            polygonsByArc[arc < 0 ? ~arc : arc].forEach(function(polygon) {
              if (!polygon._) {
                polygon._ = 1;
                neighbors.push(polygon);
              }
            });
          });
        });
      }
    }
  });

  polygons.forEach(function(polygon) {
    delete polygon._;
  });

  return {
    type: "MultiPolygon",
    arcs: groups.map(function(polygons) {
      var arcs = [], n;

      // Extract the exterior (unique) arcs.
      polygons.forEach(function(polygon) {
        polygon.forEach(function(ring) {
          ring.forEach(function(arc) {
            if (polygonsByArc[arc < 0 ? ~arc : arc].length < 2) {
              arcs.push(arc);
            }
          });
        });
      });

      // Stitch the arcs into one or more rings.
      arcs = stitch(topology, arcs);

      // If more than one ring is returned,
      // at most one of these rings can be the exterior;
      // choose the one with the greatest absolute area.
      if ((n = arcs.length) > 1) {
        for (var i = 1, k = area(arcs[0]), ki, t; i < n; ++i) {
          if ((ki = area(arcs[i])) > k) {
            t = arcs[0], arcs[0] = arcs[i], arcs[i] = t, k = ki;
          }
        }
      }

      return arcs;
    })
  };
}

var bisect = function(a, x) {
  var lo = 0, hi = a.length;
  while (lo < hi) {
    var mid = lo + hi >>> 1;
    if (a[mid] < x) lo = mid + 1;
    else hi = mid;
  }
  return lo;
};

var neighbors = function(objects) {
  var indexesByArc = {}, // arc index -> array of object indexes
      neighbors = objects.map(function() { return []; });

  function line(arcs, i) {
    arcs.forEach(function(a) {
      if (a < 0) a = ~a;
      var o = indexesByArc[a];
      if (o) o.push(i);
      else indexesByArc[a] = [i];
    });
  }

  function polygon(arcs, i) {
    arcs.forEach(function(arc) { line(arc, i); });
  }

  function geometry(o, i) {
    if (o.type === "GeometryCollection") o.geometries.forEach(function(o) { geometry(o, i); });
    else if (o.type in geometryType) geometryType[o.type](o.arcs, i);
  }

  var geometryType = {
    LineString: line,
    MultiLineString: polygon,
    Polygon: polygon,
    MultiPolygon: function(arcs, i) { arcs.forEach(function(arc) { polygon(arc, i); }); }
  };

  objects.forEach(geometry);

  for (var i in indexesByArc) {
    for (var indexes = indexesByArc[i], m = indexes.length, j = 0; j < m; ++j) {
      for (var k = j + 1; k < m; ++k) {
        var ij = indexes[j], ik = indexes[k], n;
        if ((n = neighbors[ij])[i = bisect(n, ik)] !== ik) n.splice(i, 0, ik);
        if ((n = neighbors[ik])[i = bisect(n, ij)] !== ij) n.splice(i, 0, ij);
      }
    }
  }

  return neighbors;
};

var untransform = function(transform) {
  if (transform == null) return identity;
  var x0,
      y0,
      kx = transform.scale[0],
      ky = transform.scale[1],
      dx = transform.translate[0],
      dy = transform.translate[1];
  return function(input, i) {
    if (!i) x0 = y0 = 0;
    var j = 2,
        n = input.length,
        output = new Array(n),
        x1 = Math.round((input[0] - dx) / kx),
        y1 = Math.round((input[1] - dy) / ky);
    output[0] = x1 - x0, x0 = x1;
    output[1] = y1 - y0, y0 = y1;
    while (j < n) output[j] = input[j], ++j;
    return output;
  };
};

var quantize = function(topology, transform) {
  if (topology.transform) throw new Error("already quantized");

  if (!transform || !transform.scale) {
    if (!((n = Math.floor(transform)) >= 2)) throw new Error("n must be ≥2");
    box = topology.bbox || bbox(topology);
    var x0 = box[0], y0 = box[1], x1 = box[2], y1 = box[3], n;
    transform = {scale: [x1 - x0 ? (x1 - x0) / (n - 1) : 1, y1 - y0 ? (y1 - y0) / (n - 1) : 1], translate: [x0, y0]};
  } else {
    box = topology.bbox;
  }

  var t = untransform(transform), box, key, inputs = topology.objects, outputs = {};

  function quantizePoint(point) {
    return t(point);
  }

  function quantizeGeometry(input) {
    var output;
    switch (input.type) {
      case "GeometryCollection": output = {type: "GeometryCollection", geometries: input.geometries.map(quantizeGeometry)}; break;
      case "Point": output = {type: "Point", coordinates: quantizePoint(input.coordinates)}; break;
      case "MultiPoint": output = {type: "MultiPoint", coordinates: input.coordinates.map(quantizePoint)}; break;
      default: return input;
    }
    if (input.id != null) output.id = input.id;
    if (input.bbox != null) output.bbox = input.bbox;
    if (input.properties != null) output.properties = input.properties;
    return output;
  }

  function quantizeArc(input) {
    var i = 0, j = 1, n = input.length, p, output = new Array(n); // pessimistic
    output[0] = t(input[0], 0);
    while (++i < n) if ((p = t(input[i], i))[0] || p[1]) output[j++] = p; // non-coincident points
    if (j === 1) output[j++] = [0, 0]; // an arc must have at least two points
    output.length = j;
    return output;
  }

  for (key in inputs) outputs[key] = quantizeGeometry(inputs[key]);

  return {
    type: "Topology",
    bbox: box,
    transform: transform,
    objects: outputs,
    arcs: topology.arcs.map(quantizeArc)
  };
};

exports.bbox = bbox;
exports.feature = feature;
exports.mesh = mesh;
exports.meshArcs = meshArcs;
exports.merge = merge;
exports.mergeArcs = mergeArcs;
exports.neighbors = neighbors;
exports.quantize = quantize;
exports.transform = transform;
exports.untransform = untransform;

Object.defineProperty(exports, '__esModule', { value: true });

})));

},{}],8:[function(_dereq_,module,exports){
(function (global, factory) {
	typeof exports === 'object' && typeof module !== 'undefined' ? module.exports = factory() :
	typeof define === 'function' && define.amd ? define(factory) :
	(global.vhtml = factory());
}(this, (function () { 'use strict';

var emptyTags = ['area', 'base', 'br', 'col', 'command', 'embed', 'hr', 'img', 'input', 'keygen', 'link', 'meta', 'param', 'source', 'track', 'wbr'];

var esc = function esc(str) {
	return String(str).replace(/[&<>"']/g, function (s) {
		return '&' + map[s] + ';';
	});
};
var map = { '&': 'amp', '<': 'lt', '>': 'gt', '"': 'quot', "'": 'apos' };

var sanitized = {};

function h(name, attrs) {
	var stack = [];
	for (var i = arguments.length; i-- > 2;) {
		stack.push(arguments[i]);
	}

	if (typeof name === 'function') {
		(attrs || (attrs = {})).children = stack.reverse();
		return name(attrs);
	}

	var s = '<' + name;
	if (attrs) for (var _i in attrs) {
		if (attrs[_i] !== false && attrs[_i] != null) {
			s += ' ' + esc(_i) + '="' + esc(attrs[_i]) + '"';
		}
	}

	if (emptyTags.indexOf(name) === -1) {
		s += '>';

		while (stack.length) {
			var child = stack.pop();
			if (child) {
				if (child.pop) {
					for (var _i2 = child.length; _i2--;) {
						stack.push(child[_i2]);
					}
				} else {
					s += sanitized[child] === true ? child : esc(child);
				}
			}
		}

		s += '</' + name + '>';
	} else {
		s += '>';
	}

	sanitized[s] = true;
	return s;
}

return h;

})));


},{}],9:[function(_dereq_,module,exports){
"use strict";

function _templateObject() {
  var data = _taggedTemplateLiteral(["<svg ...", " class=\"outline\">\n      ", "\n    </svg>"]);

  _templateObject = function _templateObject() {
    return data;
  };

  return data;
}

function _taggedTemplateLiteral(strings, raw) { if (!raw) { raw = strings.slice(0); } return Object.freeze(Object.defineProperties(strings, { raw: { value: Object.freeze(raw) } })); }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

var fitAspect = _dereq_('fit-aspect-ratio');

var htm = _dereq_('htm');

var vhtml = _dereq_('vhtml');

var Dot = _dereq_('./shapes/Dot');

var Text = _dereq_('./shapes/Text');

var Shape = _dereq_('./shapes/Shape');

var Line = _dereq_('./shapes/Line');

var d3Geo = _dereq_('d3-geo');

var World =
/*#__PURE__*/
function () {
  function World() {
    var obj = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : {};

    _classCallCheck(this, World);

    this.width = obj.width || 600;
    this.height = obj.height || 400;

    if (obj.aspect) {
      this.aspect = obj.aspect;
      var res = fitAspect(obj);
      this.width = res.width || 600;
      this.height = res.height || 400;
    }

    this.shapes = [];
    this.html = htm.bind(vhtml);
    this._clip = true;
    this.projection = d3Geo.geoMercator().scale(1050).center([-79.3961, 43.6601]);
  }

  _createClass(World, [{
    key: "bind",
    value: function bind(fn) {
      this.html = htm.bind(fn);
    }
  }, {
    key: "dot",
    value: function dot(obj) {
      var dot = new Dot(obj, this);
      this.shapes.push(dot);
      return dot;
    }
  }, {
    key: "line",
    value: function line(obj) {
      var dot = new Line(obj, this);
      this.shapes.push(dot);
      return dot;
    }
  }, {
    key: "text",
    value: function text(obj) {
      var dot = new Text(obj, this);
      this.shapes.push(dot);
      return dot;
    }
  }, {
    key: "shape",
    value: function shape(obj) {
      var dot = new Shape(obj, this);
      this.shapes.push(dot);
      return dot;
    }
  }, {
    key: "clip",
    value: function clip(bool) {
      this._clip = bool;
      return this;
    }
  }, {
    key: "build",
    value: function build() {
      var h = this.html;
      var shapes = this.shapes.sort(function (a, b) {
        return a._order > b._order ? 1 : -1;
      });
      var elements = [];
      elements = elements.concat(shapes.map(function (shape) {
        return shape.build();
      }));
      var attrs = {
        width: this.width,
        height: this.height,
        viewBox: "0,0,".concat(this.width, ",").concat(this.height),
        preserveAspectRatio: 'xMidYMid meet',
        style: 'margin: 10px 20px 25px 25px;' // border:1px solid lightgrey;

      };

      if (this._clip) {
        attrs.style += 'overflow:hidden;';
      } else {
        attrs.style += 'overflow:visible;';
      }

      return h(_templateObject(), attrs, elements);
    }
  }]);

  return World;
}();

module.exports = World;

},{"./shapes/Dot":23,"./shapes/Line":24,"./shapes/Shape":25,"./shapes/Text":26,"d3-geo":3,"fit-aspect-ratio":4,"htm":5,"vhtml":8}],10:[function(_dereq_,module,exports){
"use strict";

module.exports = {
  shapes: {
    lakes: _dereq_('./shapes/lakes'),
    rivers: _dereq_('./shapes/rivers'),
    'great-lakes': _dereq_('./shapes/great-lakes'),
    world: _dereq_('./shapes/world'),
    'north-america': _dereq_('./shapes/north-america'),
    countries: _dereq_('./shapes/countries'),
    provinces: _dereq_('./shapes/provinces')
  },
  points: _dereq_('./points/index')
};

},{"./points/index":12,"./shapes/countries":15,"./shapes/great-lakes":16,"./shapes/lakes":17,"./shapes/north-america":18,"./shapes/provinces":19,"./shapes/rivers":20,"./shapes/world":21}],11:[function(_dereq_,module,exports){
"use strict";

module.exports = {
  tokyo: [35.68333, 139.68333],
  'sao paulo': [-23.55, -46.63333],
  jakarta: [-6.2, 106.81667],
  seoul: [37.56667, 126.96667],
  manila: [14.58, 121],
  'new york city': [40.7127, -74.0059],
  shanghai: [31.22861, 121.47472],
  cairo: [30.04444, 31.23583],
  delhi: [28.61, 77.23],
  'mexico city': [19.43333, -99.13333],
  beijing: [39.91667, 116.38333],
  guangzhou: [23.132, 113.266],
  mumbai: [18.975, 72.82583],
  dhaka: [23.71611, 90.39611],
  osaka: [34.69389, 135.50222],
  moscow: [55.75, 37.61667],
  bangkok: [13.7525, 100.49417],
  karachi: [24.86, 67.01],
  istanbul: [41.01361, 28.955],
  kolkata: [22.56667, 88.36667],
  tehran: [35.68917, 51.38889],
  'rio de janeiro': [-22.90833, -43.19639],
  'ho chi minh city': [10.8, 106.65],
  kinshasa: [-4.325, 15.32222],
  'los angeles': [34.05, -118.25],
  lagos: [6.45503, 3.38408],
  'buenos aires': [-34.60333, -58.38167],
  nanjing: [32.05, 118.76667],
  shenzhen: [22.55, 114.1],
  shantou: [23.354, 116.682],
  lahore: [31.54972, 74.34361],
  paris: [48.8567, 2.3508],
  chengdu: [30.657, 104.066],
  nagoya: [35.18333, 136.9],
  lima: [-12.04333, -77.02833],
  london: [51.50722, -0.1275],
  chicago: [41.83694, -87.68472],
  bogotá: [4.71111, -74.07222],
  kano: [12, 8.51667],
  surabaya: [-7.26528, 112.7425],
  chennai: [13.08333, 80.26667],
  bangalore: [12.98333, 77.58333],
  bandung: [-6.912, 107.6097],
  suzhou: [31.3, 120.6],
  dongguan: [23.021, 113.752],
  busan: [35.16667, 129.06667],
  chongqing: [29.55833, 106.56667],
  hyderabad: [17.37, 78.48],
  'hong kong': [22.3, 114.2],
  'kuala lumpur': [3.13333, 101.68333],
  philadelphia: [39.95278, -75.16361],
  baghdad: [33.33333, 44.38333],
  taipei: [25.06667, 121.51667],
  changsha: [28.228, 112.939],
  wuhan: [30.58333, 114.28333],
  casablanca: [33.53333, -7.58333],
  tianjin: [39.13333, 117.18333],
  hanoi: [21.02833, 105.85417],
  harbin: [45.75, 126.63333],
  santiago: [-33.45, -70.66667],
  wenzhou: [28, 120.7],
  luanda: [-8.83833, 13.23444],
  madrid: [40.38333, -3.71667],
  ahmedabad: [23.03, 72.58],
  shenyang: [41.8, 123.4],
  'saint petersburg': [59.9375, 30.30861],
  qingdao: [36.0669, 120.3827],
  foshan: [23.01667, 113.11667],
  toronto: [43.7, -79.4],
  zunyi: [27.68333, 106.9],
  vijayawada: [16.5193, 80.6305],
  berlin: [52.51667, 13.38889],
  jinan: [36.66667, 116.98333],
  riyadh: [24.63333, 46.71667],
  fukuoka: [33.58333, 130.4],
  singapore: [1.28333, 103.83333],
  khartoum: [15.50056, 32.56],
  visakhapatnam: [17.70417, 83.29778],
  yangon: [16.85, 96.18333],
  "xi'an": [34.265, 108.954],
  'belo horizonte': [-19.91667, -43.93333],
  dallas: [32.78333, -96.8],
  hamburg: [53.56528, 10.00139],
  sydney: [-33.865, 151.20944],
  guayaquil: [-2.18333, -79.88333],
  algiers: [36.75389, 3.05889],
  houston: [29.76278, -95.38306],
  ankara: [39.93333, 32.86667],
  abidjan: [5.31667, -4.03333],
  melbourne: [-37.81361, 144.96306],
  barcelona: [41.38333, 2.18333],
  quito: [-0.23333, -78.51667],
  alexandria: [31.2, 29.91667],
  surat: [21.17024, 72.83106],
  monterrey: [25.66667, -100.3],
  'guatemala city': [14.61333, -90.53528],
  fuzhou: [26.0753, 119.3082],
  johannesburg: [-26.20444, 28.04556],
  guadalajara: [20.67667, -103.3475],
  'porto alegre': [-30.03306, -51.23],
  'dar es salaam': [-6.8, 39.28333],
  rome: [41.9, 12.5],
  shijiazhuang: [38.06667, 114.48333],
  giza: [30.01667, 31.21667],
  changchun: [43.88694, 125.32472],
  zhengzhou: [34.764, 113.684],
  medan: [3.58333, 98.66667],
  fortaleza: [-3.71833, -38.54278],
  accra: [5.55, -0.2],
  chittagong: [22.36667, 91.8],
  isfahan: [32.64472, 51.6675],
  dalian: [38.92083, 121.63917],
  'new taipei city': [25.01111, 121.44583],
  'salvador, bahia': [-12.97472, -38.47667],
  brasília: [-15.79389, -47.88278],
  recife: [-8.05, -34.9],
  'cape town': [-33.92528, 18.42389],
  yokohama: [35.44417, 139.63806],
  'phoenix, arizona': [33.45, -112.06667],
  i̇zmir: [38.42, 27.14],
  cologne: [50.93639, 6.95278],
  hangzhou: [30.25, 120.16667],
  baku: [40.39528, 49.88222],
  ningbo: [29.869, 121.554],
  jeddah: [21.54333, 39.17278],
  oran: [35.69694, -0.63306],
  durban: [-29.88333, 31.05],
  kabul: [34.53333, 69.16667],
  montreal: [45.50889, -73.56167],
  cali: [3.42056, -76.52222],
  curitiba: [-25.41667, -49.25],
  'addis ababa': [9.03, 38.74],
  kiev: [50.45, 30.52333],
  ibadan: [7.39639, 3.91667],
  mashhad: [36.3, 59.6],
  hefei: [31.86667, 117.28333],
  budapest: [47.4925, 19.05139],
  pyongyang: [39.01944, 125.73806],
  kathmandu: [27.71139, 85.30861],
  milan: [45.46667, 9.18333],
  faisalabad: [31.41667, 73.09111],
  tangshan: [39.631, 118.18],
  semarang: [-6.96667, 110.41667],
  'ekurhuleni metropolitan municipality': [-26.21861, 28.16278],
  nairobi: [-1.28333, 36.81667],
  zhongshan: [22.517, 113.393],
  pune: [18.52028, 73.85667],
  warsaw: [52.23333, 21.01667],
  campinas: [-22.90083, -47.05722],
  jaipur: [26.9, 75.8],
  incheon: [37.48333, 126.63333],
  'san diego': [32.715, -117.1625],
  'quezon city': [14.63, 121.03],
  caracas: [10.48056, -66.90361],
  kanpur: [26.44992, 80.33187],
  dubai: [25.26306, 55.29722],
  lucknow: [26.85, 80.95],
  basra: [30.5, 47.81667],
  chaozhou: [23.66667, 116.63333],
  maracaibo: [10.65, -71.63333],
  munich: [48.13333, 11.56667],
  prague: [50.08333, 14.41667],
  vienna: [48.2, 16.36667],
  sapporo: [43.06667, 141.35],
  kaohsiung: [22.63333, 120.26667],
  taichung: [24.15, 120.66667],
  'davao city': [7.07, 125.6],
  ouagadougou: [12.35722, -1.53528],
  medellín: [6.23083, -75.59056],
  nagpur: [21.15, 79.08333],
  daegu: [35.86667, 128.6],
  lusaka: [-15.41667, 28.28333],
  almaty: [43.2775, 76.89583],
  'san antonio': [29.41667, -98.5],
  dakar: [14.69278, -17.44667],
  douala: [4.05, 9.68333],
  birmingham: [52.48333, -1.9],
  yaoundé: [3.86667, 11.51667],
  abuja: [9.06667, 7.48333],
  barranquilla: [10.96389, -74.79639],
  brisbane: [-27.46667, 153.03333],
  tripoli: [32.88722, 13.19139],
  tashkent: [41.3, 69.26667],
  tunis: [36.80639, 10.18167],
  'phnom penh': [11.55, 104.91667],
  managua: [12.13639, -86.25139],
  lanzhou: [36.062, 103.835],
  baoding: [38.86667, 115.46667],
  "sana'a": [15.34833, 44.20639],
  havana: [23.13333, -82.38333],
  rawalpindi: [33.6, 73.03333],
  patna: [25.6, 85.1],
  gujranwala: [32.15667, 74.19],
  makassar: [-5.13333, 119.41667],
  karaj: [35.83556, 51.01028],
  peshawar: [34.01667, 71.58333],
  minsk: [53.9, 27.56667],
  bucharest: [44.4325, 26.10389],
  tijuana: [32.525, -117.03333],
  shiraz: [29.61667, 52.53333],
  xiamen: [24.47984, 118.08942],
  islamabad: [33.71667, 73.06667],
  brazzaville: [-4.26778, 15.29194],
  maputo: [-25.96667, 32.58333],
  'santa cruz de la sierra': [-17.8, -63.18333],
  bhopal: [23.25, 77.41667],
  tabriz: [38.06667, 46.3],
  'hyderabad, pakistan': [25.37917, 68.36833],
  mandalay: [21.975, 96.08333],
  montevideo: [-34.88361, -56.18194],
  palembang: [-2.98611, 104.75556],
  belgrade: [44.81667, 20.46667],
  kharkiv: [50.00444, 36.23139],
  harare: [-17.82917, 31.05222],
  auckland: [-36.84056, 174.74],
  caloocan: [14.65, 120.97],
  novosibirsk: [55.05, 82.95],
  sofia: [42.7, 23.33],
  kobe: [34.69, 135.19556],
  daejeon: [36.351, 127.385],
  'córdoba, argentina': [-31.41667, -64.18333],
  kampala: [0.31361, 32.58111],
  'kawasaki, kanagawa': [35.51667, 139.7],
  tbilisi: [41.71667, 44.78333],
  gwangju: [35.16667, 126.91667],
  kyoto: [35.01167, 135.76833],
  quanzhou: [24.907, 118.587],
  yekaterinburg: [56.83333, 60.58333],
  'rosario, santa fe': [-32.95, -60.66667],
  ahvaz: [31.32028, 48.66917],
  'nizhny novgorod': [56.32694, 44.0075],
  saitama: [35.86139, 139.64556],
  suwon: [37.26667, 127.01667],
  qom: [34.64, 50.87639],
  astana: [51.16667, 71.43333],
  hiroshima: [34.38333, 132.45],
  tainan: [22.98333, 120.18333],
  ulsan: [35.55, 129.31667],
  omsk: [54.98333, 73.36667],
  'abu dhabi': [24.46667, 54.36667],
  'rostov-on-don': [47.23333, 39.7],
  allahabad: [25.45, 81.85],
  fez: [34.03333, -5],
  calgary: [51.05, -114.06667],
  yerevan: [40.18139, 44.51444],
  'cartagena, colombia': [10.4, -75.5]
};

},{}],12:[function(_dereq_,module,exports){
"use strict";

module.exports = Object.assign({}, _dereq_('./cities'), _dereq_('./ontario'), _dereq_('./manitoba'));

},{"./cities":11,"./manitoba":13,"./ontario":14}],13:[function(_dereq_,module,exports){
"use strict";

module.exports = {
  brandon: [49.84833, -99.95],
  'flin flon': [54.76667, -101.87778],
  'portage la prairie': [49.97278, -98.29194],
  dauphin: [51.14944, -100.04944],
  morden: [49.19194, -98.10056],
  winnipeg: [49.89944, -97.13917],
  thompson: [55.74333, -97.85528],
  steinbach: [49.52583, -96.68389],
  winkler: [49.18167, -97.93972],
  selkirk: [50.14361, -96.88389]
};

},{}],14:[function(_dereq_,module,exports){
"use strict";

module.exports = {
  brampton: [43.68333, -79.76667],
  barrie: [44.37111, -79.67694],
  belleville: [44.16667, -77.38333],
  brantford: [43.16667, -80.25],
  cornwall: [45.0275, -74.74],
  brockville: [44.58333, -75.68333],
  burlington: [43.31667, -79.8],
  cambridge: [43.36667, -80.31667],
  'clarence-rockland': [45.48333, -75.2],
  guelph: [43.55, -80.25],
  dryden: [49.78333, -92.83333],
  'elliot lake': [46.38333, -82.65],
  'greater sudbury': [46.49, -81.01],
  'haldimand county': [42.93333, -79.88333],
  hamilton: [43.25667, -79.86917],
  kitchener: [43.41861, -80.47278],
  kingston: [44.23333, -76.5],
  kenora: [49.76667, -94.48333],
  'kawartha lakes': [44.35, -78.75],
  london: [42.98361, -81.24972],
  mississauga: [43.6, -79.65],
  markham: [43.87667, -79.26333],
  'niagara falls': [43.06, -79.10667],
  'norfolk county': [42.85, -80.26667],
  ottawa: [45.42472, -75.695],
  'north bay': [46.3, -79.45],
  orillia: [44.6, -79.41667],
  oshawa: [43.9, -78.85],
  'owen sound': [44.56667, -80.93333],
  pickering: [43.83944, -79.08139],
  peterborough: [44.3, -78.31667],
  'port colborne': [42.88333, -79.25],
  pembroke: [45.81667, -77.1],
  sarnia: [42.99944, -82.30889],
  'st. catharines': [43.18333, -79.23333],
  'richmond hill': [43.86667, -79.43333],
  'quinte west': [44.18333, -77.56667],
  'sault ste. marie': [46.53333, -84.35],
  'thunder bay': [48.38222, -89.24611],
  stratford: [43.37083, -80.98194],
  'st. thomas': [42.775, -81.18333],
  thorold: [43.11667, -79.2],
  'temiskaming shores': [47.51667, -79.68333],
  toronto: [43.74167, -79.37333],
  waterloo: [43.46667, -80.51667],
  timmins: [48.46667, -81.33333],
  vaughan: [43.83333, -79.5],
  welland: [42.98333, -79.23333],
  windsor: [42.28333, -83],
  woodstock: [43.13056, -80.74667]
};

},{}],15:[function(_dereq_,module,exports){
module.exports={"type":"Topology","arcs":[[[25998,14837],[-26,20],[-75,-4],[-99,40],[-36,-11]],[[25762,14882],[-32,29],[-9,55],[-62,36],[-60,17],[-46,111],[4,118],[-52,11],[-7,73],[-67,56],[-63,77],[-50,195],[-33,66],[-33,96],[5,26]],[[25257,15848],[46,-6],[28,-28],[58,5],[55,-25],[29,18]],[[25473,15812],[33,62]],[[25506,15874],[8,16]],[[25514,15890],[28,33],[16,54],[54,81],[62,30],[14,-23],[27,43],[-28,17]],[[25687,16125],[19,120],[46,46],[55,29],[82,0]],[[25889,16320],[5,-77],[98,-6],[56,-42],[31,-47],[37,-5],[87,-57],[-9,-38],[11,-81],[-1,-182],[4,-61],[-17,-89],[-17,-52],[17,-22],[-9,-63],[28,-106],[-2,-25],[-39,-116],[-24,-43],[0,-53],[-15,-40],[2,-44],[-116,-212],[-18,-22]],[[25687,16125],[-57,-8],[-31,-54],[-31,-1],[-13,-54],[-31,-51],[-12,-55],[-20,-31],[-19,-59]],[[25257,15848],[-41,55],[-81,14],[-110,-35]],[[25025,15882],[-132,209],[-19,20],[-21,141],[0,262],[0,386],[229,0],[19,127],[-3,173],[9,50],[-9,118]],[[25098,17368],[43,-44],[24,-81],[19,19],[63,24],[21,-83],[83,-59],[70,-18],[44,16],[11,62],[47,-72],[13,-55],[102,-63],[24,-78],[15,-5],[42,-112],[23,-6],[43,33],[30,-32],[0,272],[-29,-11],[-17,-50],[-44,18],[-72,118],[-15,72],[23,150],[11,29],[-3,60],[2,124],[-28,133],[27,34],[37,81],[6,52],[219,59]],[[25932,17955],[10,-42],[30,-49],[34,1],[40,-67],[33,-31],[60,-23],[60,-55]],[[26199,17689],[10,-47],[24,-1],[19,-52],[-5,-38],[27,-44],[16,-78],[-49,-66],[15,-58],[-19,-83],[9,-78],[-6,-49],[18,-83],[-46,-30],[-7,-131],[-14,-51],[-22,-34],[39,-90],[26,2]],[[26234,16678],[-181,-106],[-24,-18],[-106,-53],[-55,-37],[21,-144]],[[25506,15874],[8,16]],[[28681,23408],[-62,-39],[-32,-39],[-19,-74],[6,-66],[-76,-70],[-72,-43],[-60,-22],[-45,-42],[-32,-5],[-68,-41],[-37,-40],[-47,-89],[-48,-12],[-35,11],[-44,-41],[-28,-44],[-76,-43],[-72,-9],[-67,-19],[-33,-59],[-28,-15],[-16,-40],[-51,0],[-32,-32],[-53,-12],[-56,49],[-30,94],[6,82],[-6,48],[-18,33],[-8,119],[-37,242],[15,83],[-5,74]],[[27415,23347],[45,70],[-6,55],[9,92],[28,34],[61,-42],[30,17],[122,6],[32,-24],[111,-19],[46,8],[36,-69],[20,-2],[37,37],[17,74],[73,155],[107,93],[162,43],[199,48]],[[28544,23923],[93,-372],[8,-7],[36,-136]],[[28764,22527],[52,6],[40,-24],[-47,-42],[-65,-4],[-35,42],[14,28],[41,-6]],[[35434,24475],[-69,-49],[-7,-51],[-83,-21],[9,-49],[-25,-46],[-4,-52],[-43,-65],[-22,-12],[-21,-77],[-3,-64],[-21,-72],[33,-102],[43,-62],[32,-60],[-3,-50],[32,-77],[47,-69],[8,-34],[67,-90],[13,-37],[24,2],[80,-186],[18,-23],[32,-145],[0,-36],[27,-153],[-4,-127],[19,-71],[-25,-68],[-3,-202],[-19,-29],[-7,-56],[-20,-5],[-68,-83],[-22,-8],[-30,-44],[-72,-66],[-22,22],[-33,-14],[-4,-70],[-23,-56],[-3,-48],[-50,-54],[-82,-66],[-22,-64],[-25,-38],[-28,-10],[-9,48],[3,176],[8,46],[23,19],[-8,36],[-28,30],[-17,-7],[-29,53]],[[34998,22039],[17,23],[35,4],[66,89],[45,11],[43,-42],[8,47],[-38,62],[-1,75],[19,27],[41,-17],[9,58],[64,29],[28,43],[43,20],[6,38],[-10,108],[16,89],[-34,151],[4,54],[19,73]],[[35378,22981],[-4,61],[21,60],[-57,108],[25,65],[-57,68],[-29,23],[-18,42],[-6,67],[-28,58],[-33,43],[-38,71],[-29,91],[-16,8],[-30,57],[4,51],[-49,27],[-101,118],[20,31],[-4,42],[64,-14],[24,32],[23,58],[-29,41],[-9,37],[-28,11],[23,49],[-48,58],[-17,6],[-51,-53],[-66,42],[-31,82],[12,67],[-31,20],[-14,54],[-56,104]],[[34715,24666],[37,77],[11,2],[60,-63],[43,66],[33,6],[36,-53],[17,47],[39,-11],[52,46],[9,49],[50,46],[34,-60],[36,-33],[38,10],[12,-22],[36,9],[30,-28],[-25,-46],[-6,-39],[24,-86],[34,-24],[43,-62],[55,8],[21,-30]],[[14764,21630],[22,-53],[-24,-56],[-43,-49],[-17,-3],[-25,-69],[40,-122],[-40,-39],[-52,-16],[-7,-23],[9,-82],[-31,-68],[78,-154]],[[14674,20896],[17,-46],[-17,-48],[-32,-53],[-34,-4],[-34,-50],[-73,-42],[-32,13],[-37,-30],[-7,-76],[-25,-17],[-40,72],[-89,1],[-29,47],[-39,-1],[1,-52],[42,-69],[1,-84],[22,-88],[-2,-61],[16,-16],[66,-9],[-1,-41],[-79,-70],[-10,-63],[-23,-36],[-22,0],[-31,-44],[-45,-32],[-32,-53],[-39,26],[-47,-44],[-45,8],[-55,88]],[[13920,20022],[-1,30],[-41,226],[-47,84],[-29,27],[59,94],[5,21],[-43,98],[-24,141],[4,168],[22,63],[2,33],[23,69],[-16,48],[-35,10],[-76,-29],[-76,6],[-43,-11],[-80,177],[-36,7],[-33,24],[-44,-19],[-75,8],[-34,-9],[-23,23],[-15,60],[-29,19],[-1,112],[11,27],[-7,65],[-31,53],[-16,106],[-32,33],[-37,-7],[43,123],[17,155],[22,75],[30,61],[25,18],[35,103],[79,43]],[[13373,22357],[-21,-31],[-47,-20],[-9,-47],[43,-166],[-9,-47],[-35,-74],[-20,-64],[17,-57],[27,-48],[9,-58],[47,12],[27,49],[4,78],[-26,96],[-15,27],[-20,113],[9,40],[26,8],[54,46],[73,36],[22,35],[20,-19],[9,50],[-47,-11],[-11,58],[10,46],[25,18],[21,-40],[15,-95],[10,-18],[49,8],[49,-19],[53,-59],[20,-130],[45,-22],[109,31],[91,4],[49,-82],[89,-41],[34,6],[81,79],[44,-2],[-5,47],[46,-4],[98,23],[57,-11],[-17,-33],[-39,3],[-18,-32],[22,-48],[18,-5],[24,-86],[25,44],[43,-37],[32,7],[25,-47],[36,-17],[27,-43],[-32,-58],[-31,-150],[62,40],[39,-10],[18,18],[40,-18]],[[42666,16499],[8,-72],[40,5],[8,-77],[-46,-33],[-24,50],[2,43],[-15,79],[27,5]],[[42748,16221],[52,-78],[-47,-23],[-17,97],[12,4]],[[42857,16168],[-34,5],[18,51],[16,-56]],[[42875,15903],[10,-56],[-34,4],[3,50],[21,2]],[[42985,15597],[-43,15],[19,52],[24,-67]],[[30878,29027],[-13,-50],[-81,-99],[74,-68],[77,-27],[2,45],[45,16],[11,-48],[35,-5],[56,-64],[61,-21],[-93,-84],[-28,8],[-22,-40],[-32,-23],[-48,29],[-42,-10]],[[30880,28586],[-38,-8],[-33,55],[45,63],[-38,62],[-18,-29],[-72,-51],[-36,30],[-26,-50],[12,-53],[-9,-30],[-71,-2],[27,-60],[-29,-15],[-21,-66],[-92,18],[-36,-34],[-9,-49],[33,-20],[7,-31],[52,-10],[-6,-65],[11,-63],[24,-31],[-5,-45],[-59,-120],[-7,-69]],[[30486,27913],[-30,20],[-24,-12],[-31,27],[-67,3]],[[30334,27951],[0,96],[13,32],[-95,74],[-30,-6],[-117,109],[-43,53],[-18,-5],[-55,51],[-26,42],[-106,105],[-25,39],[-8,65],[-35,77],[-17,77],[-56,53],[-21,-17],[-108,5],[-45,40],[2,111],[-12,62],[-18,23],[-50,1],[-41,49],[-30,9],[-35,47],[-41,-20],[-55,-109],[-24,-6]],[[29238,29008],[-28,38],[-25,-33],[-12,-53]],[[29173,28960],[-15,-19],[6,-130],[-128,13]],[[29036,28824],[0,453],[0,353],[305,120],[17,-1]],[[29358,29749],[-44,-151],[2,-81],[26,-4],[38,66],[-15,57],[18,39],[-6,59]],[[29377,29734],[85,-72]],[[29462,29662],[-29,-87],[13,-35],[62,-89],[28,50],[2,96]],[[29538,29597],[117,-99],[19,-50],[102,-148],[150,30],[152,-17],[57,36],[73,-89],[37,-95],[37,25],[-11,-217],[60,-2],[21,-159],[18,-28],[138,12],[22,-37],[-8,-48],[30,-34],[38,1],[-2,48],[49,72],[10,35],[73,67],[54,33],[39,57],[24,-6],[41,43]],[[15581,12348],[-51,-140],[-47,-64],[-90,-58],[-95,35],[-55,-29],[-90,51],[-39,49],[-81,-6],[-70,124],[-5,43],[10,101],[34,80],[-15,93]],[[14987,12627],[10,32],[-8,88],[18,34],[-2,60],[30,140],[-8,58],[33,89]],[[15060,13128],[48,-21],[12,30],[35,9],[81,-115],[37,-97],[39,43],[38,-77],[82,-57],[45,-81],[56,-44],[20,-76],[59,-74],[-50,-95],[-1,-99],[20,-26]],[[14011,23802],[62,-14],[1,-30],[-43,-59],[-54,-5],[-97,9],[5,111],[119,-5],[7,-7]],[[1050,33722],[53,-12],[75,25],[50,-45],[73,-32],[87,-8],[-16,-31],[-114,-17],[-92,63],[-66,14],[-42,-23],[-48,19],[40,47]],[[1706,33007],[50,-15],[17,-88],[-73,-31],[-117,52],[-37,43],[74,3],[44,36],[42,0]],[[3335,32446],[47,6],[37,-61],[-57,-51],[-73,-26],[-64,-52],[-34,32],[-33,-47],[-46,91],[21,49],[32,18],[57,-10],[76,50],[37,1]],[[6149,32021],[-15,-52],[14,-78],[-22,-60]],[[6126,31831],[-44,-56]],[[6082,31775],[-34,8],[-25,77],[24,44],[-8,83],[-45,61],[-61,-30],[-11,-50],[-43,16],[46,96],[-41,58],[-78,136],[-80,28],[-11,115],[-59,83],[-4,30],[-87,54],[-90,164],[56,-197],[-31,-11],[-69,71],[-23,46],[-64,0],[57,-62],[-3,-34],[-56,-23],[-118,74],[-59,79],[-40,32],[-178,101],[26,66],[-78,-27],[-77,5],[-102,49],[-157,27],[-105,-18],[-135,71],[-123,31],[-31,40],[-89,76],[-88,-6],[-93,-28],[18,-150],[-28,-39],[-51,-9],[-73,11],[-94,-95],[-47,0],[-47,-71],[-95,-12],[-14,48],[81,51],[-71,13],[4,64],[35,47],[16,118],[108,62],[55,-20],[6,79],[-82,2],[-138,-86],[-2,-34],[-64,-50],[-47,-64],[-13,-69],[-43,-10],[-121,-113],[-5,-53],[59,-18],[39,-40],[-84,-80],[-28,-73],[-94,-31],[-181,-151],[1,-49],[-136,-104],[-103,-43],[17,-53],[-63,-46],[-90,-39],[-54,-3],[-62,-52],[-78,-33],[-79,-3],[-41,-56],[-98,-40],[-12,50],[100,118],[115,66],[45,-55],[39,52],[29,67],[63,54],[62,28],[77,77],[42,57],[87,67],[0,98],[21,27],[-10,55],[67,69],[-20,31],[-109,-52],[-34,0],[-25,46],[-44,-29],[11,-44],[-38,-11],[-72,97],[-31,-20],[-55,51],[-123,-84],[-48,-13],[-5,110],[-28,38],[24,67],[-51,129],[-40,-41],[-80,-32],[-85,-8],[-93,109],[-85,52],[68,77],[66,-37],[104,7],[-40,47],[-31,-23],[-145,22],[-48,30],[-28,85],[-38,34],[3,35],[57,14],[-11,51],[63,82],[49,33],[-6,39],[53,91],[38,10],[92,-47],[90,48],[42,58],[38,-17],[108,23],[39,57],[-26,95],[-49,42],[56,30],[6,48],[-37,28],[-71,-24],[-128,-97],[-112,47],[-73,-1],[-73,-27],[-153,27],[-42,32],[9,43],[-51,38],[20,53],[-91,18],[-84,52],[124,48],[41,38],[119,59],[136,52],[76,12],[34,-17],[-31,-70],[41,-29],[231,-6],[32,43],[163,37],[-3,34],[-74,21],[-66,-26],[-57,33],[29,61],[-42,15],[-127,-4],[-89,36],[-50,91],[-155,96],[-70,24],[-52,59],[21,101],[86,-4],[149,16],[53,22],[87,77],[26,81],[131,126],[109,-6],[207,126],[39,-17],[123,10],[99,54],[51,50],[199,-49],[47,-45],[72,-20],[119,29],[123,-27],[-30,-35],[65,-44],[89,-7],[69,20],[63,-15],[108,13],[173,-43],[20,-19],[231,-13],[77,-33],[244,24],[223,-102],[50,-1]],[[4799,35041],[0,-320],[0,-385],[0,-449],[0,-321],[0,-256],[0,-321],[126,-26],[37,32],[73,4],[-13,-58],[196,-188],[19,-71],[77,55],[29,0],[31,99],[105,42],[52,-47],[13,-64],[62,-45],[27,-51],[49,-32],[67,-106],[174,-314],[36,-42],[102,-57],[41,-48],[48,-24],[-1,-27]],[[2033,31821],[15,-53],[-28,-24],[-60,-1],[-48,-43],[-45,-2],[-7,42],[44,60],[81,34],[48,-13]],[[5892,32069],[-2,-36],[-65,29],[35,93],[38,-56],[-6,-30]],[[6032,31933],[-13,-49],[-89,-9],[-5,42],[25,37],[2,54],[44,27],[38,-62],[-2,-40]],[[5738,32265],[51,-27],[-1,-44],[-31,-47],[-48,0],[-10,73],[-22,60],[61,-15]],[[3970,32882],[-10,20],[86,101],[16,-24],[-63,-80],[-29,-17]],[[5541,32342],[43,-139],[-1,-91],[-87,114],[17,45],[-73,63],[55,48],[46,-40]],[[5576,32520],[55,-4],[-4,-57],[41,-67],[1,-53],[-62,-62],[-20,8],[13,85],[-21,35],[-12,78],[9,37]],[[5447,32538],[18,-22],[72,-21],[-10,-118],[-51,33],[-15,-40],[-36,-7],[-20,50],[-61,65],[7,27],[96,33]],[[5714,32120],[44,-5],[13,-50],[42,-21],[40,-68],[2,-44],[54,-61],[0,-88],[-99,53],[-82,165],[26,27],[-36,49],[-4,43]],[[5661,32231],[17,-58],[-42,-37],[-22,93],[47,2]],[[3394,32563],[38,-43],[-108,-36],[-49,19],[87,80],[32,-20]],[[1647,31584],[21,-54],[-40,-39],[-45,11],[19,96],[45,-14]],[[43413,31390],[76,-36],[-61,-22],[-15,58]],[[3005,23927],[-37,13],[-1,68],[-20,81],[28,58],[-1,57],[78,-62],[16,-54],[32,-49],[-30,-45],[-32,-13],[-33,-54]],[[2893,24348],[26,4],[35,-42],[-14,-25],[-41,-6],[-6,69]],[[2732,24463],[-39,-30],[-20,59],[39,25],[20,-54]],[[2538,24568],[-29,-5],[-22,29],[26,40],[28,-1],[-3,-63]],[[12956,29632],[188,0],[204,0],[24,62],[53,-6],[24,64],[50,78],[6,75],[22,42],[7,59],[94,165],[51,-57],[64,34],[62,-60],[3,-303],[43,-22],[-5,-72],[43,-23]],[[13889,29668],[2,-62],[-31,-42],[-67,-53],[-46,12],[-41,-44],[-43,28],[-13,-74],[-20,-25],[-69,-39],[-48,-9],[-3,-31],[-39,-61],[-38,-115],[12,-34],[-39,-74],[38,-22],[41,-133],[-31,-20],[-94,22],[-12,-56],[-92,-19],[-71,-6],[-91,-53],[-45,-41],[-34,-50],[-2,-34],[35,-28],[-27,-132],[-31,-74],[-57,-53],[-57,42],[-4,-54],[44,-129],[-40,-96],[-29,-96],[-36,-2],[28,73],[-68,90],[-7,52],[12,56],[-13,103],[-16,-4],[-14,-89],[2,-81],[29,-140],[0,-118],[-16,-31],[14,-36],[35,-30],[13,-56],[-31,-83],[-41,-55],[87,-29],[0,-55],[-49,-64],[-42,-19],[-29,-72],[2,-51],[-50,0],[-76,-93],[-32,-82],[-69,-8],[-42,-47],[-44,-115],[-92,-113],[-26,-10],[-76,-102],[-9,-37],[-34,-52],[3,-38],[-42,-145],[33,-237],[43,-163],[46,-124],[-15,-67],[52,-214],[12,-27],[11,-113],[-11,-161],[-21,-47],[-8,-63],[-46,-39],[-46,-4],[2,40],[-33,112],[-43,34],[-19,99],[-21,26],[-18,63],[-32,48],[-24,102],[-3,43],[-22,28],[22,140],[1,88],[-14,36],[-64,88],[-50,104],[-43,39],[-33,-8],[-9,-35],[-51,-29],[-64,-21],[-4,43],[-50,64],[-51,36],[-10,36],[-52,-20],[-38,6],[-23,24],[-13,-44],[-63,-10],[-30,85],[-12,-69],[-69,-3],[-26,14],[-51,-16],[-33,-39],[-78,47],[-24,-53],[30,-24],[22,6],[50,-31],[7,-32],[-23,-33],[25,-44],[61,-47],[-17,-41],[-24,11],[-65,86],[-21,-6],[-6,-55],[-34,25],[-46,-36],[-66,34],[-6,50],[-22,9],[-20,42],[-27,19],[-45,-61],[-50,9],[-63,40],[-63,-2],[-50,-22],[-97,-66],[-48,-89],[-71,-71],[-52,3],[-21,-13],[-14,-50],[-30,-32],[-43,-17],[-2,-58],[-26,-101],[-18,-23],[-3,-76],[18,-125],[32,-96]],[[10195,25452],[-26,-20],[-55,38],[-34,5],[-126,84],[-15,69],[-27,60],[-7,112],[-31,40],[-67,131],[-4,39],[-40,124],[-96,155],[-68,5],[-34,16],[-58,-50],[-19,-94],[-45,-47],[-90,71],[-51,55],[-34,91],[0,32],[-37,112],[-36,36],[-66,95],[-55,67],[-24,49],[-217,2],[0,-99],[-130,0],[-218,-1],[-233,130],[-234,130],[14,45],[-296,-39]],[[7736,26895],[-14,28],[-9,96],[-19,43],[-75,93],[-41,5],[-11,60],[-40,2],[-39,19],[-15,32],[-42,35],[-107,12],[-20,24],[1,109],[-28,29],[0,34],[-60,94],[-65,118],[-5,53],[16,35],[-12,45],[-34,12],[-29,47],[-13,74],[4,65],[-61,57],[-3,39],[-88,139],[-14,101],[4,55],[-12,53],[-54,86],[-6,52],[29,105],[8,91],[-21,89],[4,43],[-18,30],[-6,101],[-17,51],[24,116],[24,77],[11,238],[16,175],[-2,124],[-15,31],[22,83],[-22,19],[-38,200],[-30,54],[-10,55],[9,49],[81,-46],[137,-17],[25,-46],[13,-104],[-32,-17],[13,-48],[43,57],[-5,99],[19,42],[-67,215]],[[7040,30507],[170,0],[210,0],[262,0],[157,0],[210,0],[262,0],[158,0],[210,0],[210,0],[262,0],[262,0],[263,0],[262,0],[315,0],[186,0],[1,83],[37,-14],[6,-67],[23,-57],[69,-18],[43,-29],[55,22],[33,-3],[101,-60],[80,-62],[84,31],[12,-21],[80,3],[12,-22],[51,-3]],[[11126,30290],[-125,-69],[-45,-36],[-52,-68],[-84,-87],[28,-26],[110,56],[26,-59],[33,-22],[74,45],[63,20],[64,58],[76,-43],[71,-27],[34,-67],[40,5],[49,-20],[50,44],[76,8],[74,20],[-3,-61],[63,-15],[24,12],[43,-116],[-67,7],[-34,-26],[-21,25],[-62,26],[-60,-37],[-74,-12],[-64,-34],[-92,-173],[61,-36],[-24,-86],[0,-46],[-24,-59],[3,-40],[-25,-105],[15,-96],[-2,-131],[35,-108],[32,-18],[56,33],[21,27],[41,111],[9,106],[-39,162],[12,34],[-9,61],[27,57],[7,85],[47,53],[46,-26],[13,96],[51,25],[-15,47],[39,41],[71,-40],[8,-17],[80,-43],[19,-44],[-16,-25],[19,-75],[-4,-67],[-35,-70],[-33,-21],[-7,-49],[25,-25],[33,23],[16,44],[48,30],[35,-47],[30,-183]],[[12007,29196],[-8,-61],[-20,-40]],[[11979,29095],[-56,-60],[-20,-74],[-24,-45],[115,-74],[99,23],[61,57],[110,59],[77,57],[61,56],[41,54],[-4,59],[-19,8]],[[12420,29215],[0,39]],[[12420,29254],[109,24],[80,-27],[68,2],[95,65],[-4,49],[17,37],[-17,48],[54,67]],[[12822,29519],[76,90],[58,23]],[[13226,28750],[21,-20],[-105,-53],[-70,-14],[-14,44],[37,28],[131,15]],[[11346,30160],[-79,-63],[5,58],[35,18],[39,-13]],[[17583,7887],[49,-9],[47,-32],[65,-112],[-36,-23],[-52,78],[-28,10],[-117,86],[38,20],[34,-18]],[[14907,8501],[42,-9],[4,-44],[62,27],[17,-71],[-104,-50],[3,-35],[-51,-17],[-71,9],[9,46],[63,52],[-4,43],[30,49]],[[14730,8458],[126,8],[-81,-119],[-40,-4],[-13,-33],[-41,-11],[-34,29],[89,62],[-6,68]],[[21820,31084],[-77,-51],[-33,6],[-101,71],[-64,-25],[-27,18],[4,63],[87,45],[50,75],[-15,82],[-66,-5],[65,72],[128,45],[18,68],[-30,86],[-50,75],[17,109],[-65,-40],[-68,-5],[-76,21],[51,124],[-21,61],[5,77],[-49,-26],[-38,-119],[-26,-6],[19,152],[27,101],[-67,20],[46,138],[-31,44],[25,98],[34,79],[40,5],[-7,56],[43,-2],[205,27],[-20,-68],[-95,-80],[15,-79],[57,24],[39,-8],[124,7],[37,-46],[-60,-138],[-41,-67],[-7,-53],[-65,-52],[16,-28],[56,16],[55,-27],[61,-73],[45,-175],[76,-59],[72,-85],[-15,-21],[69,-189],[-38,-56],[29,-21],[34,34],[61,-1],[74,-45],[11,-63],[-23,-84],[-35,-53],[-39,-8],[-8,-52],[-25,-40],[90,-6],[-3,-39],[-53,-57],[-93,-35],[-51,11],[-71,-11],[-62,20],[-28,-24],[-63,-5],[-1,-26],[-63,3],[-56,22],[-49,-19],[-23,-68],[-25,-20],[-47,35],[-68,-22],[-61,-59],[-15,49],[94,117],[7,44],[41,46],[43,9],[87,-5],[57,91]],[[21438,32121],[-49,-12],[9,78],[42,-34],[-2,-32]],[[21386,32564],[-28,-75],[-54,-20],[-27,55],[104,71],[5,-31]],[[21392,32376],[2,-42],[-77,28],[63,46],[12,-32]],[[21383,31626],[-53,-7],[-35,69],[-52,-55],[-61,20],[-32,44],[28,55],[41,23],[41,71]],[[21260,31846],[92,33],[42,-6],[81,-157],[-16,-50],[-51,-49],[-25,9]],[[29076,25384],[8,-18],[3,-130]],[[29087,25236],[-40,-52],[-33,28],[-4,-137],[27,-40],[-53,-15],[-5,-59],[-38,-151],[-2,-73]],[[28939,24737],[-10,-18],[-157,34],[-157,34],[-118,252],[-3,45]],[[28494,25084],[24,-7],[17,-59],[43,2],[49,35],[153,-17],[62,45],[43,116],[43,51],[52,100],[52,65],[17,59]],[[29049,25474],[8,-81],[19,-9]],[[26851,30090],[-83,-4],[-25,-34],[-130,-57],[-56,-8],[-53,-53],[-21,13],[-46,-55],[18,-100]],[[26455,29792],[-38,53],[-127,53],[-8,-27]],[[26282,29871],[-48,18],[-89,-20],[-55,39],[-3,37],[43,10],[-47,51],[-145,-35],[-35,-98],[-36,-52],[-76,-40],[13,-93]],[[25804,29688],[-17,24],[-43,7],[-57,-37],[-38,13],[-29,35]],[[25620,29730],[90,131],[-2,82],[27,23],[85,-39],[11,42],[-6,63],[-38,30],[-3,62],[-47,46],[4,93],[-52,53],[-55,9],[-62,56],[-34,13],[-39,-23],[-47,3],[-28,-28]],[[25424,30346],[-39,-12],[-17,-46],[-86,-18],[-70,-43],[-39,47],[-49,-5],[-120,39],[-41,-30]],[[24963,30278],[-91,100]],[[24872,30378],[50,147]],[[24922,30525],[20,21],[-7,81],[48,79],[83,105],[32,7],[17,96],[-55,102],[-7,45]],[[25053,31061],[43,15],[50,62],[112,16],[81,-6],[65,-22],[84,-13],[16,-31],[112,-1],[52,-14],[66,16],[26,-51],[100,21],[21,-34],[37,7],[-12,53],[16,47],[39,51],[96,12]],[[26057,31189],[64,3],[28,43],[48,-12],[30,19],[73,1],[81,-124],[-34,-22],[20,-81],[53,-30],[43,6],[31,-35],[12,-110],[32,-43],[27,20],[89,-50],[18,18],[81,27],[35,-67],[42,-41],[26,29],[81,-50],[32,7],[74,-63],[31,6],[5,-87],[-47,-43],[34,-169],[-28,-86],[-103,2],[-22,-30],[-53,-34],[-9,-103]],[[26059,19534],[-115,0],[-41,-14]],[[25903,19520],[-19,-2],[-32,-65],[-34,8],[-30,-12]],[[25788,19449],[-2,90],[10,104]],[[25796,19643],[28,44],[-23,40]],[[25801,19727],[4,48],[27,88],[1,71],[29,33],[17,47],[20,12]],[[25899,20026],[14,-44],[46,107],[37,35],[16,46],[-17,58]],[[25995,20228],[-37,53],[-28,12],[14,86],[-11,42],[18,81],[-8,18]],[[25943,20520],[39,65],[40,-23],[39,26],[37,-60],[29,39],[62,21],[20,17],[19,-23],[41,-4],[60,102]],[[26329,20680],[23,-89],[28,-27],[7,-116],[36,-76],[28,-166],[0,-92],[-22,-57],[0,-33],[-38,-41],[-8,-39],[-31,-57],[-23,-86]],[[26329,19801],[-49,2],[-37,52],[-17,-46],[-114,-53],[-32,-31],[5,-42],[-33,-103],[7,-46]],[[29173,28960],[32,-81],[58,27],[-25,102]],[[30334,27951],[-51,15],[-42,34],[-26,-70],[-57,-3],[-34,-23],[-38,-174],[-79,-72],[-93,-42],[7,-31],[-14,-51],[-54,-47],[-44,4],[-34,42],[-54,3],[-34,36]],[[29687,27572],[-18,224],[-96,-1],[-34,72],[-62,47],[-39,75],[-52,36],[-47,-10],[-89,59],[-36,9],[-28,59],[-93,7],[-26,-39],[-80,6],[-62,-43],[-28,-39],[-12,-49],[-68,-37],[-35,2]],[[28782,27950],[-11,129],[6,224],[-31,58],[-61,48],[14,38],[29,14],[-8,64],[-52,6],[-29,51],[-2,39],[26,140],[25,-47],[58,-1],[31,-37],[55,41],[48,12],[-14,66],[-61,72],[-31,125],[-27,10],[-57,-9],[-24,-26],[-17,-139],[-27,41],[-14,55]],[[28608,28924],[46,62],[47,32],[107,28],[68,-56],[28,-35],[43,-102],[36,-42],[53,13]],[[27256,28866],[118,-17],[39,28],[81,-101]],[[27494,28776],[34,-85],[-18,-52],[27,-90],[61,-7],[59,-74]],[[27657,28468],[6,-11]],[[27663,28457],[-28,26],[-25,-76],[-45,-10],[31,-119],[2,-91],[18,-14],[-27,-99],[43,-36],[1,-68],[27,-31],[-3,-33]],[[27657,27906],[-20,8],[-40,-44],[-9,49],[-24,25],[-22,-20],[-39,5],[-52,27],[-39,1],[-51,-58]],[[27361,27899],[-11,37],[-93,-42],[-100,2],[-85,-59],[-81,-32],[-72,3],[-40,37],[-31,8],[-93,-56],[-46,2],[-50,33],[-15,-76],[13,-43],[-33,-20],[-3,-37],[-27,-38],[-29,19]],[[26565,27637],[8,18],[-18,68],[46,77],[-17,55],[-63,-69],[-44,9],[-45,36],[-26,-3],[-37,-40],[-34,-58],[-40,-35],[-111,-32],[-52,33],[-43,77],[-83,58],[-87,14],[-20,-121],[-31,-1],[-66,-33],[-68,53],[-12,65],[-40,1],[-29,23],[-83,-15],[40,64],[-96,-2],[21,51],[-37,30],[-18,66],[19,64],[-67,48],[-48,17],[17,36],[53,-15],[-12,74],[30,39],[-40,89],[17,59],[-87,-21],[8,115],[16,8],[53,82],[67,12],[23,-29],[100,18],[62,29],[71,60],[-40,44],[6,37],[29,11],[73,-17],[52,10],[58,-25],[54,5],[25,47],[78,59],[27,31],[132,63],[168,-13],[32,23],[36,-73],[32,-21],[61,11],[15,-57],[41,-36],[33,23],[35,-40],[86,-23],[76,-34],[129,40],[47,-27],[56,-5],[85,55],[68,67]],[[25596,28966],[-4,-26],[44,-85],[88,-52],[-13,-48],[-21,-8],[-75,24],[-30,-20],[-53,-4],[-29,-63],[-60,-41],[-55,-83],[-9,42],[66,69],[-84,-3],[-8,25]],[[25353,28693],[35,50],[0,62],[37,36],[-6,44],[-32,25]],[[25387,28910],[36,55],[84,25],[30,-35],[59,11]],[[23564,27037],[-6,-119],[6,-50],[-83,-73],[-43,-75],[-21,-7],[-24,-53],[16,-95],[-3,-55],[-40,-87],[-46,-34]],[[23320,26389],[-59,404],[-87,104],[-15,84],[-17,35],[-38,32],[-33,132],[2,54],[40,72],[35,34],[15,38],[18,103],[-18,146],[13,109],[-15,39],[43,86]],[[23204,27861],[30,13],[39,44],[68,32],[62,-30],[-1,-38],[28,-66],[19,32],[60,43],[9,-44],[-41,-83],[-33,-38],[-6,-32],[14,-63],[51,-56],[14,-86],[-31,-78],[-41,-75],[-41,-44],[-19,-63],[14,-45],[69,-74],[44,10],[15,-68],[37,-15]],[[14641,21978],[-72,-15],[16,104],[-15,46],[83,20],[-15,-37],[9,-76],[-6,-42]],[[22259,22167],[-14,-62],[-1,-77],[68,-80],[2,-94],[10,-62],[21,-52],[3,-171],[0,-279],[-5,-68],[24,-87],[-19,-17]],[[22348,21118],[-53,-28]],[[22295,21090],[-25,51],[-31,29],[-26,87],[14,111],[-17,42],[13,40],[-3,92],[13,34],[-9,39],[-24,38],[9,44],[-3,58],[8,53],[-33,30],[7,42],[5,121],[-52,86],[9,86],[-10,21]],[[22140,22194],[68,-28],[51,1]],[[37538,17666],[-13,48],[-5,77]],[[37520,17791],[31,65],[77,34],[101,7],[42,32],[41,-24],[-47,-64],[-31,-15],[-34,-38],[-27,-9],[-93,-61],[-42,-52]],[[37411,17703],[50,34]],[[37461,17737],[-20,-53],[-30,19]],[[34468,24213],[24,15],[28,-65],[-18,-73],[28,-56],[42,24],[23,-7],[17,-132],[-30,-125],[11,-26],[-3,-56],[-25,-99],[24,-23],[72,90],[11,32],[40,39],[76,-70],[41,30],[39,87],[71,-23],[59,-136],[38,-53],[10,-35],[-7,-49],[-1,-94],[36,-107],[44,-38],[-4,-22],[31,-41],[-13,-75],[4,-94],[-9,-88],[-36,-40]],[[35091,22903],[-50,18],[-45,-13],[-73,-2],[-27,13],[-49,-19],[-36,-43],[-22,-65],[-50,-66],[1,-56],[16,-60],[5,-75],[31,-54],[-2,-74],[24,-84]],[[34814,22323],[-21,67],[-27,22],[-62,92],[-39,35],[-34,-15],[-72,20],[13,158],[-38,20],[-52,-9],[-27,-28],[9,-68],[-15,-78],[3,-114],[-45,-155],[-17,-126],[-40,-125],[0,-129],[29,-114],[40,22],[22,-44],[7,-97],[40,-88],[21,-182],[-4,-56],[20,-2],[53,-69],[60,1],[37,-86],[37,-51]],[[34712,21124],[-28,-92],[-28,-6],[-11,24],[-37,-26],[-18,-33],[-16,29],[12,41],[-3,63],[-22,1],[-10,39],[-55,27],[-28,-23]],[[34468,21168],[-52,95],[3,51],[-24,49],[-21,9],[-7,41],[-31,72],[-32,38],[-26,62],[-19,-36],[-23,54],[1,75],[16,115],[41,197]],[[34294,21990],[6,104],[54,97],[52,149],[-19,89],[-7,79],[-23,42],[-14,80],[9,28],[-5,106],[-25,74],[-45,68],[-40,99],[-6,86],[46,44],[3,142],[12,56],[-26,93],[-5,54],[-79,155],[-9,78],[-29,100],[33,28],[-6,68],[12,44],[2,72],[36,66],[16,-15],[83,18],[13,59],[39,8],[37,53],[31,23],[28,-24]],[[27008,18399],[-31,-23],[-6,76],[22,-4],[15,-49]],[[25932,17955],[-54,218],[-25,63],[-49,65],[-21,67],[-7,63],[14,83],[-31,146],[7,117]],[[25766,18777],[39,-1],[53,81],[34,109],[45,69],[-1,64],[-33,14],[-11,61],[16,53]],[[25908,19227],[34,14],[-2,170],[-42,95],[5,14]],[[26059,19534],[10,-20],[-24,-184],[1,-76],[30,-87],[22,61],[32,0],[56,-26],[31,23],[43,-27],[49,75],[-66,15],[49,106],[-3,15],[50,116]],[[26339,19525],[219,-218],[222,-222],[6,-44],[-10,-47],[23,-47],[176,-223]],[[26975,18724],[-30,-181],[-20,-79],[7,-100],[73,-120],[6,-54],[-28,-86],[19,-109],[-17,-95],[18,-109],[23,-55],[-2,-47],[20,-111],[83,-121]],[[27127,17457],[-58,-78],[-52,-35],[-30,-31],[-65,-24],[-37,-40],[-39,29],[-32,-3],[-11,-41],[-35,-43],[-40,-2],[-29,26],[-57,-33],[-40,10],[-14,29],[-36,19],[-19,-28],[-83,1]],[[26450,17213],[-1,6]],[[26449,17219],[-43,112],[8,72],[-10,41],[-7,108],[-25,66],[-40,52],[-6,-36]],[[26326,17634],[-65,11],[-11,19],[-51,25]],[[30880,28586],[-63,-42],[-59,34],[-54,-23],[-37,-74],[8,-52],[37,11],[125,-2],[15,-33],[83,40],[40,-38],[38,-16],[147,5],[49,16]],[[31209,28412],[-3,-48],[23,-50],[-12,-32],[13,-54],[27,-15],[31,26],[58,-33],[11,-23],[-6,-66],[18,-73],[-5,-37],[28,-47],[-28,-34]],[[31364,27926],[-29,35],[-49,5],[-63,-40],[-3,41],[-30,11],[-72,-44],[-29,-53],[-37,-10],[-68,-63],[-26,8],[-20,87],[15,177],[-34,-4],[-3,86],[-46,33],[-32,-27],[-50,-90],[5,-57],[-94,-25],[-9,-76],[-19,-21],[-37,47],[-36,-14],[-74,-68],[-38,49]],[[37038,24719],[-13,-53],[-7,-100],[-19,24],[-13,71],[-32,41],[-31,133],[11,123],[58,169],[50,121],[68,53],[42,-66],[-14,-32],[1,-64],[-23,-89],[-20,-155],[-22,-100],[-36,-76]],[[27361,27899],[-71,-113],[-45,-18],[-15,-28],[-6,-68],[13,-96],[-17,-77],[-3,-114],[-25,-74],[-37,-22],[-236,-210]],[[26919,27079],[-240,-232],[-55,15],[-72,77]],[[26552,26939],[16,47],[-6,106]],[[26562,27092],[9,51],[42,50],[37,72],[-19,86],[-56,4]],[[26575,27355],[-9,48],[5,82],[-22,76],[16,76]],[[23324,30182],[-3,-53]],[[23321,30129],[6,-47]],[[23327,30082],[52,-37],[56,-5]],[[23435,30040],[-3,-69],[-42,11],[-19,-52],[-84,-17],[-28,-51],[-9,-49],[-61,91],[-46,2],[-11,-52],[-120,-20]],[[23012,29834],[-30,53],[-2,55],[-86,2],[116,197],[76,59]],[[23086,30200],[71,3],[52,35],[32,-24],[37,3]],[[23278,30217],[46,-35]],[[24496,32448],[-33,-28],[-14,-102],[-51,-35],[-22,41],[4,63],[49,60],[67,1]],[[23550,32712],[31,-25],[19,80],[-14,67],[31,59],[68,54],[13,75],[-37,121],[48,10],[24,67],[-89,81],[18,124],[-23,67],[13,90],[-27,64],[22,66],[60,76],[78,30],[81,-14],[23,35],[-8,64],[-53,26],[102,158],[19,119],[-11,62],[116,39],[-8,41],[121,124],[-34,81],[56,44],[25,59],[66,46],[73,-31],[57,131],[196,-45],[0,41],[80,109]],[[24686,34907],[169,-114],[97,-28],[105,-96],[-21,-111],[64,-140],[-35,-72],[56,-149]],[[25121,34197],[-57,5],[-73,-20],[-90,24],[-14,-54],[-85,-42],[1,-62],[-54,-69],[41,-95],[-55,-44],[-32,-68],[-128,-89],[-45,3],[-39,-52],[-53,-13],[-18,-62],[-69,-15],[-8,-72],[-57,-16],[24,-47],[-18,-104],[-36,-36],[15,-192],[87,-24],[74,-74],[46,-69],[5,-40],[-57,-77],[-28,-65],[-101,-55],[-100,-94],[15,-48],[-26,-88],[12,-68],[-22,-52],[6,-43],[-22,-79],[-43,-107],[-83,-16],[-67,3],[-69,-73],[15,-66],[-20,-29],[-46,7],[-59,-18],[-54,14],[11,62],[-47,97],[26,28],[10,78],[-57,63],[-85,178],[-3,56],[-32,32],[-22,79],[-15,112],[30,11]],[[26080,14056],[13,-13],[-2,-58],[9,-122]],[[26100,13863],[-15,5],[-4,-108],[-60,3],[-50,40],[-33,76],[1,77],[49,125],[26,22],[56,-52],[10,5]],[[15485,20930],[-37,-76],[1,-116],[12,-94],[44,-102],[-26,-99],[1,-71],[-26,-78],[-26,-30]],[[15428,20264],[-48,59],[-47,-34],[-42,-7],[-25,24],[-24,-47],[26,-58],[-13,-43],[-57,22]],[[15198,20180],[-27,21],[-62,186],[-12,108],[-30,-5],[-35,71],[-24,71],[-3,37],[26,109],[-9,34],[25,37],[47,7],[15,42],[-13,27],[15,47]],[[15111,20972],[17,85],[11,12],[62,-12],[64,-31],[14,36],[122,6],[59,-17],[37,-22],[1,-44],[-13,-55]],[[25943,20520],[-34,36],[-7,40],[-38,32],[-64,133],[-25,5],[-39,-49],[-53,26],[-37,-40],[-29,6],[-49,64],[-48,102]],[[25520,20875],[-18,40],[-5,72],[-16,35],[-44,48],[-26,16],[-25,85],[6,39],[-58,94],[-62,59],[-24,41],[7,45],[-48,91],[-70,34],[-17,82]],[[25120,21656],[8,22],[39,26],[16,99],[15,42],[0,54],[27,62],[8,52],[90,28],[11,-48],[35,-52],[47,-96],[26,-6],[38,25],[100,-3],[20,-60],[98,0],[19,62],[58,35],[16,65],[49,47],[58,-66],[33,-52],[60,13],[53,92],[34,97],[58,87],[-9,144],[-32,64],[82,1],[-2,47],[59,-2],[-16,-137],[12,-163],[89,-143],[8,-77],[-8,-84],[23,-1]],[[26342,21830],[2,-193],[-18,-30],[-67,1],[-21,-11],[-27,-98],[26,-42],[54,-20],[29,-35],[20,-62],[26,-49],[26,-23],[28,-52],[45,-217],[23,-40]],[[26488,20959],[-159,-279]],[[25120,21656],[-58,14],[-17,19],[-9,65],[19,50],[3,106],[-41,124],[-56,117]],[[24961,22151],[8,93],[-35,37],[-10,34],[3,70],[-11,12],[-17,136],[-43,2],[-22,26],[20,71],[30,47],[-15,88],[8,40],[41,48],[-16,92],[36,38],[11,61],[21,36],[0,81],[22,37],[108,17],[-1,327],[0,195],[0,294]],[[25099,24033],[0,110],[123,1],[0,438]],[[25222,24582],[141,0],[235,0],[188,-1],[230,1],[145,0],[161,0],[162,0],[201,0]],[[26685,24582],[7,-90],[35,-88],[-7,-18],[9,-120],[-4,-96],[7,-118],[14,-108],[13,-59],[56,-58],[34,-68],[50,-53]],[[26899,23706],[-44,-93],[-57,-27],[-34,-38],[-12,-49],[-50,-1],[-15,-95],[4,-72],[-49,-230],[-11,-26],[12,-192]],[[26643,22883],[-17,-138],[-22,-78],[-13,-120],[-35,-9],[-27,-32],[-43,-128],[-17,-31],[-19,-118],[-3,-91],[-20,-26],[-25,30],[-36,-77],[2,-89],[-21,-87],[-5,-59]],[[31990,21908],[33,-4],[57,-94],[81,-206],[8,-68],[28,-74],[26,-109],[-2,-84],[-27,-105],[-32,-40],[-81,-58],[-56,7],[-21,31],[-29,149],[-8,166],[-11,105],[13,-3],[16,140],[-1,46],[23,94],[-17,107]],[[22536,28487],[38,-20],[-26,-68],[-21,-19],[-38,24],[-16,32],[-33,12],[51,53],[45,-14]],[[21928,29281],[47,-36],[15,-40],[65,-26],[21,-31],[36,5],[31,-29],[83,0],[8,34],[81,-34],[9,-21]],[[22324,29103],[34,-20]],[[22358,29083],[41,-33],[20,15],[45,-17],[50,27],[30,-8]],[[22544,29067],[4,-107],[-30,-39],[-85,-66],[-28,-39],[-129,-50],[-40,-52],[-43,-111],[-54,-97],[-31,-78],[15,-100],[50,-67],[-71,-71],[-33,-62],[-22,-126],[-62,-4],[-58,-72],[-38,-100],[-60,6],[-24,-20],[-58,9],[-136,-8],[-38,-46],[-61,-19],[-23,-63],[-33,-24],[-51,36],[-28,89],[-14,9],[-13,70],[-49,53],[-64,-4]],[[21237,27914],[-12,89],[39,93],[20,8],[-39,91],[8,56],[29,42],[-50,139],[-1,31],[52,8],[17,67],[-15,41],[26,30],[-1,68],[-14,78],[88,115],[-40,30],[-10,60],[-51,2],[-45,-26],[-64,11],[-29,-16],[-7,71],[-70,-43]],[[21068,28959],[-13,36],[24,38],[-36,68],[-7,50],[-24,36],[7,43],[38,35],[41,1],[35,22],[0,31],[68,41],[54,-38],[196,-2],[141,-37],[113,23],[69,-33],[21,19],[66,-28],[67,17]],[[20139,25983],[-11,-50],[-29,-32],[-31,73],[48,16],[23,-7]],[[20402,25937],[23,118],[22,-27],[-12,-73],[-33,-18]],[[20254,25932],[1,-60],[-40,-20],[-12,46],[16,35],[35,-1]],[[19954,26008],[-20,58],[31,6],[-11,-64]],[[37730,28046],[4,30],[43,71],[19,14],[107,1],[42,69]],[[37945,28231],[30,-98],[53,-109],[35,-89],[17,-117],[-9,-151],[20,-23],[-18,-99],[-25,-69],[-29,-18],[-58,0],[-8,-51],[-90,19],[-40,-46],[-8,-62],[-36,32],[-64,-57],[-14,128],[-13,49],[25,76],[0,32],[28,27],[-3,54],[-19,33],[1,68],[-41,87],[45,32],[20,38],[-22,115],[8,14]],[[37693,27046],[1,52],[44,20],[22,-46],[-67,-26]],[[25998,14837],[30,-193],[2,-44],[31,-90],[23,-125],[0,-257],[-4,-72]],[[26100,13863],[95,-2]],[[26195,13861],[-5,-51],[-23,-116],[-15,-130],[-31,-92],[-32,-48],[-30,-22],[-55,-96],[-38,-115],[-91,-235],[-39,-77],[-29,-33],[-74,-117],[-34,-63],[-122,-167],[-97,-103],[-79,-52],[-55,11],[-41,-31],[-2,-35],[-78,8],[-22,-43],[-28,-1],[-51,25],[-73,16],[-39,-21],[-88,16],[-38,-13],[-56,-67],[-67,-7],[-23,9],[-65,-21],[-63,-71],[-47,7],[-39,40],[-28,48],[-32,-3],[-3,56],[-34,5],[-22,33],[10,48],[-69,166],[-3,30],[49,38],[9,32],[-1,84],[-13,84],[-65,158],[-60,202],[-30,153],[-26,87],[-36,86]],[[24172,13473],[38,36],[22,84],[15,8],[41,-70],[7,-76],[107,-41],[35,6],[53,-16],[79,93],[38,12],[0,421],[0,386]],[[24607,14316],[45,-56],[16,-42],[39,-152],[3,-55],[-23,-61],[1,-66],[12,-23],[112,-1],[23,38],[32,22],[15,42],[52,69],[30,135],[47,42],[32,-17],[46,-56],[36,-7],[17,-24],[52,-17],[85,23],[18,23],[40,189],[60,29],[53,82],[19,118],[20,39],[91,99],[40,84],[21,22],[56,24],[36,59],[29,4]],[[25684,13147],[45,40],[5,48],[19,29],[12,65],[-11,40],[-43,46],[-40,65],[-49,-26],[-61,-52],[-54,-127],[-30,-32],[38,-109],[4,-34],[44,-63],[33,-9],[46,109],[42,10]],[[27259,19382],[-1,27],[-67,154],[-1,308],[0,185],[-1,316],[46,85],[67,170]],[[27302,20627],[17,35],[95,34],[28,59],[64,61],[60,24],[112,-8],[123,238],[127,226],[124,213]],[[28052,21509],[118,319],[0,397]],[[28170,22225],[56,18],[31,24],[57,17],[44,44],[21,47],[19,9],[57,-33],[-21,-109],[7,-149],[-9,-49],[-21,-40],[-9,-181],[-48,-128],[-40,-141],[-31,-53],[-12,-66],[-23,-80],[-27,-67],[-32,-128],[-5,-51],[-49,-149],[-51,-119],[-32,-100],[-135,-266],[-102,-178],[-28,-36],[-111,-110],[-72,-92],[-107,-169],[-112,-207],[-71,-143],[-30,-100],[-25,-58]],[[28052,21509],[-123,0],[-84,52],[-172,102],[-108,63],[-49,77],[-17,9],[-37,110],[-42,71],[-22,92],[32,83]],[[27430,22168],[40,110]],[[27470,22278],[24,-34],[50,-123],[66,-78],[69,2],[48,46],[59,41],[92,-20],[50,40],[54,54],[37,-13],[38,6],[81,40],[32,-14]],[[41907,17164],[-71,48],[10,19],[61,-67]],[[41561,17945],[16,-58],[-38,13],[-3,39],[-26,-2],[10,64],[41,-56]],[[41805,17719],[27,-36],[48,3],[33,-37],[24,-60],[-21,-14],[-40,23],[-58,7],[-28,57],[15,57]],[[41821,17881],[-29,35],[-86,73],[-43,62],[17,34],[46,-66],[40,-27],[51,-66],[4,-45]],[[41527,18145],[-6,-21],[-42,22],[-24,32],[-26,59],[-11,59],[53,-55],[14,-49],[42,-47]],[[41928,17929],[31,-65],[-7,-42],[27,-35],[12,-78],[-27,2],[-20,33],[-35,172],[19,13]],[[42047,17474],[48,-15],[22,-56],[-20,-25],[-41,24],[-31,33],[22,39]],[[24872,30378],[-45,-15],[-39,47],[-47,-10],[-71,5],[-19,-51],[-54,-36],[-34,20],[-19,-24],[-68,-14],[-24,-57],[-118,-4],[-76,52]],[[24258,30291],[-34,83],[11,47]],[[24235,30421],[28,57],[71,6],[40,39],[9,42],[54,52],[29,4]],[[24466,30621],[39,-24],[36,43],[44,-89],[100,44],[126,4],[47,-44],[64,-30]],[[24181,29960],[-109,-66],[-1,-80],[-43,-22],[8,-59],[-75,9],[-105,-17],[-37,19]],[[23819,29744],[18,16]],[[23837,29760],[-29,99],[18,30],[-29,31],[37,44]],[[23834,29964],[99,-22],[61,46],[118,11],[17,41]],[[24129,30040],[23,-2],[29,-78]],[[20883,21616],[-13,-72],[-22,-19],[-10,-68],[-76,-116],[-29,-71]],[[20733,21270],[-28,40],[-93,65],[-3,61],[-45,42],[-10,64],[-23,15],[7,46],[-17,59],[9,45],[-17,33]],[[20513,21740],[32,12],[9,35],[25,24],[31,107],[44,3],[29,26],[87,-3],[63,-146],[23,-137],[-26,-78],[19,-4],[34,37]],[[24942,29464],[-35,-51],[-3,-66],[20,-55],[52,-59],[-34,-67],[-29,-9],[7,-88],[-22,-28]],[[24898,29041],[-13,10],[-83,-24]],[[24802,29027],[23,92],[-44,18],[-41,75],[-54,23],[0,-36],[-34,-45]],[[24652,29154],[-49,56],[-49,28],[-44,71]],[[24510,29309],[36,37],[-28,53],[40,22],[-53,57],[0,49],[25,73],[-43,2]],[[24487,29602],[7,59],[33,29],[-40,28],[6,29],[-26,67],[8,21]],[[24475,29835],[37,12],[40,37],[87,-10]],[[24639,29874],[66,-79],[0,-58],[88,-79],[-13,-56],[31,-35],[56,-37],[50,36],[-1,-59],[26,-43]],[[20637,23004],[12,-35],[-5,-42],[26,-55],[7,-73],[-12,-53],[20,-41],[32,-17],[30,-94],[-7,-90],[7,-28]],[[20747,22476],[-80,-1],[-31,-15],[-62,27],[-35,40],[-80,9]],[[20459,22536],[-180,1],[-47,-42],[-70,-7],[-33,-20],[-37,-2]],[[20092,22466],[-6,155]],[[20086,22621],[14,20],[100,0],[3,38],[65,15],[16,35],[59,-45],[52,-25],[50,22],[-16,45],[-53,-8],[-65,61],[-22,6],[-39,-18],[-10,-31],[-129,0]],[[20111,22736],[-23,56],[-6,55],[-46,120],[-30,34],[33,28],[37,82],[38,119]],[[20114,23230],[11,81],[34,75],[67,-9],[63,30],[71,3],[29,-17],[40,-59],[13,-36],[39,-1],[17,-18],[19,-79],[33,-43],[16,-58],[71,-95]],[[27415,23347],[-13,80],[-38,85],[-11,69],[-67,99],[-39,124],[-25,50],[-14,88],[-44,148],[-84,112],[-24,6],[-33,49],[-42,100],[-22,74],[6,46],[-19,80],[13,112],[-19,107],[-59,182],[-21,44],[-46,60],[-46,24],[-44,116],[10,30],[-14,73],[-29,77],[-19,24],[-11,63],[-19,15],[-33,107],[-82,184],[-20,66],[-30,66],[-69,25],[20,79],[21,186]],[[26449,26197],[131,-36],[57,67],[34,81],[88,28],[20,70],[43,41],[-126,218],[247,110],[22,29]],[[26965,26805],[151,-41],[210,-188],[96,-129],[226,-284],[205,-30],[21,7]],[[27874,26140],[111,-23],[29,-100],[95,2]],[[28109,26019],[23,-90],[18,-38],[16,-73],[33,-42],[29,-56],[71,-78],[5,-66],[23,-48],[-22,-43],[30,-120],[22,-31],[13,-74],[30,-65]],[[28400,25195],[20,-47],[37,7]],[[28457,25155],[5,-59],[32,-12]],[[28939,24737],[56,-154],[-82,-440],[-161,-96],[-208,-124]],[[943,16798],[20,-48],[-6,-26],[-39,3],[-25,46],[50,25]],[[25908,19227],[-18,19],[-33,-22],[-26,16],[-7,-82],[-59,-21],[-12,30],[-34,-10]],[[25719,19157],[-20,60]],[[25699,19217],[36,37],[26,57],[-17,73]],[[25744,19384],[44,65]],[[25,34885],[114,-40],[22,-46],[54,-1],[104,-65],[254,-125],[70,-196],[67,-33],[28,32],[-55,71],[95,9],[128,-48],[104,2],[53,-56],[144,-105],[51,-13],[-94,-60],[-16,-54],[-190,-43],[0,-83],[-87,-75],[53,-35],[-64,-55],[-88,15],[-73,60],[-88,39],[-44,-3],[-57,46],[-29,101],[-119,31],[-53,-24],[-101,-4],[-24,56],[-46,43],[-5,54],[-86,-12],[-14,-71],[47,-63],[-44,-72],[-36,-27],[0,54],[0,54],[0,54],[0,53],[0,54],[0,54],[0,53],[0,54],[0,54],[0,54],[0,53],[0,54],[0,54],[0,53],[0,54],[0,54],[25,-10]],[[38229,29039],[-20,51]],[[38209,29090],[7,60],[60,20],[24,117],[-11,60],[10,80],[-24,129],[3,55],[45,16],[50,76],[21,-21]],[[38394,29682],[5,-108],[37,-46],[30,39],[24,78]],[[38490,29645],[38,14],[9,80],[30,25],[10,60],[43,81],[0,55],[42,138],[-5,38],[39,32],[33,59],[-23,67],[0,66],[-33,11],[-56,-22],[-46,-38],[-39,1],[-54,-34],[-18,-40],[-96,-19],[-96,3],[-36,95],[12,47],[-31,58],[0,56],[-44,7],[-86,109],[-47,-6],[-51,52],[-87,-6],[-35,22],[-26,44],[11,74],[-35,54],[0,56],[-47,86],[-14,89],[-37,92],[-24,104],[-37,46],[4,29],[-49,66],[-83,34],[-20,-14],[-148,91],[-106,-17],[-50,3],[-115,-37],[-52,-7],[-34,-25],[-75,-84],[-4,-34],[73,-14],[11,-104],[-29,-54],[-55,-54],[-107,-268],[19,-23],[-11,-46],[-99,-49],[-72,-73],[-77,25],[-69,44]],[[36506,30690],[-67,41],[-52,-29],[-54,15],[-66,63],[-55,8],[-89,-58],[-50,-46],[-13,-40],[-120,-42],[-92,-15],[-62,-42],[-78,1],[-118,36],[-77,-2],[-71,59],[-11,73],[-88,14],[-64,71],[-61,-2],[-102,37],[-113,-37],[-35,-30],[-68,-6],[-117,54],[-48,44],[-22,168],[-90,22],[-104,57],[-75,10],[-25,25],[-107,44],[-26,-59],[-35,-19],[-17,-46],[-33,-35],[-15,-79],[25,-75],[31,-25],[-22,-100],[-110,-76],[-28,34],[-82,4],[-57,24],[-41,-22],[-111,24],[-32,44],[-13,73],[-141,11],[-20,38],[-39,-20],[-34,39],[-19,-36],[-48,-2],[-57,-49],[-40,-12],[-169,-112],[1,-41],[-101,-59],[-79,1],[-9,-43],[-37,-21]],[[32954,30544],[-61,-16]],[[32893,30528],[-87,115],[-53,-25],[-120,37],[-27,87],[-90,50],[-38,106],[-73,49],[-73,-23],[-33,-36],[-126,3],[-10,47],[-32,-2],[0,54],[-31,20],[-53,-22],[-3,-52],[-53,-38],[-64,158],[-122,251],[-95,163],[-139,124],[26,83],[-64,-2],[-75,-46],[-75,-61],[-73,-39],[-144,-11],[14,47],[-60,65],[-64,3],[-48,18],[-82,-10],[-30,12],[12,85],[-34,116],[-21,39],[-69,-32],[-85,43],[-63,7],[-56,-40],[-39,-10],[-16,-44],[-47,-19],[-181,-45],[-92,-10],[-39,-57],[-86,5],[-56,-26],[-100,-21],[-14,-15],[-141,-35],[-74,23],[-23,-19],[-20,-75],[33,-15],[-6,-58],[100,-40],[-5,-36],[-114,7],[-30,-50],[17,-102],[-108,-67],[50,-71],[110,-46],[27,-47],[-24,-81],[-55,-36],[-37,-7],[-84,30],[-18,-48],[-114,25],[-65,81],[-64,6],[-48,-45],[-53,39],[-64,-10],[-86,-92],[-31,0],[-100,87],[-38,-47],[-47,73],[-43,30],[-43,51],[-95,-1],[-29,44],[-40,-2],[-49,-41],[-33,-2],[-22,38],[-45,18],[-68,-97],[-92,-45],[-21,-51],[-64,-55],[-13,-54],[18,-75],[-51,-40],[-26,22],[-63,98],[-47,-23],[-4,-47],[-37,-45],[-24,-107],[26,-59],[-38,-64],[-5,-86],[42,-28],[11,-74],[25,-46],[23,14],[77,-13],[55,-93],[50,-120],[-43,-45],[76,-51]],[[28206,29924],[-18,-32],[-49,-23],[-65,-76],[-78,-30],[-4,-32],[-46,-121],[-44,-55],[-1,-51],[59,-50],[34,-150],[-5,-104],[94,-188],[42,-74]],[[28125,28938],[-22,-53],[-41,-31],[-24,-54],[-74,22],[-7,31],[-56,78],[-39,17]],[[27862,28948],[-59,32],[-39,37],[11,60],[-47,12],[-58,49],[-45,-1],[-68,-40],[-21,39],[-85,53],[-13,23],[-70,29],[-103,-1],[-115,69],[-61,8],[-21,-33]],[[27068,29284],[-80,105],[-75,86],[-66,28],[-41,62],[-54,8],[-26,52],[-72,47],[40,22],[87,20],[-4,41],[39,96],[31,20],[-27,115],[163,107],[-99,15],[-33,-18]],[[26057,31189],[-23,47],[-5,106],[-34,48],[19,40],[89,-23],[70,54],[-62,97],[-55,7],[4,53],[-47,31],[-35,100],[-40,29],[22,59],[-20,50],[12,64],[-83,60],[-37,0],[-56,-35],[-13,56],[-35,18],[-36,-17],[-29,32],[-51,12]],[[25612,32077],[-5,88],[-37,71],[-1,98],[-55,47]],[[25514,32381],[36,107]],[[25550,32488],[56,6],[-64,57],[30,40],[-19,109]],[[25553,32700],[42,110]],[[25595,32810],[1,53],[114,23],[25,37],[120,-28],[-49,71],[-81,-1],[-52,41],[1,51],[-105,-16]],[[25569,33041],[75,79],[104,86],[54,56],[196,225],[31,69],[-44,71],[-138,106],[55,73],[-2,47],[-47,28],[-4,88],[-58,45],[14,144],[47,12],[-24,90],[-100,166],[0,27],[105,127],[4,31],[-78,82],[-81,28],[-26,65],[61,117]],[[25713,34903],[52,61],[74,21],[29,53],[78,-21],[1,54]],[[25947,35071],[72,-21],[65,58],[118,-44],[-70,-60],[133,-11],[28,-26],[268,-26],[230,-110],[86,-74],[50,-7],[90,-55],[58,-13],[42,-40],[72,-26],[21,-98],[28,-12],[-21,-85],[-82,-83],[-152,-69],[-110,-15],[-135,35],[-220,38],[-84,47],[-166,34],[-78,56],[-47,-23],[136,-130],[102,-43],[48,-58],[-21,-77],[-26,-26],[49,-89],[8,-94],[21,-26],[75,-14],[61,-41],[27,-41],[133,-42],[77,47],[-14,64],[-95,20],[-81,96],[44,71],[79,-14],[59,-50],[216,-66],[64,92],[-59,77],[2,55],[63,34],[44,46],[97,35],[90,87],[49,-21],[77,-2],[107,-39],[48,96],[-25,93],[-62,34],[54,163],[-22,122],[127,6],[100,-22],[98,-111],[1,-27],[-143,-20],[-73,-62],[0,-27],[77,-37],[52,-73],[87,-6],[99,23],[26,25],[20,119],[124,32],[241,136],[211,53],[16,31],[138,58],[20,-34],[-26,-46],[26,-89],[73,10],[106,65],[77,17],[133,-20],[129,73],[91,25],[56,-58],[-34,-64],[99,-5],[0,62],[73,5],[54,57],[-33,27],[-61,105],[104,57],[93,-19],[196,-19],[102,-31],[103,-58],[171,-71],[72,-16],[90,-79],[95,-36],[40,48],[38,84],[-74,3],[-52,59],[-14,54],[-115,47],[-17,70],[27,45],[17,114],[-76,36],[31,83],[76,29],[90,59],[41,73],[54,144],[74,50],[228,4],[162,-49],[-23,-135],[-43,-84],[-45,-32],[21,-53],[61,-34],[19,-72],[-25,-143],[11,-151],[-11,-89],[35,-58],[90,-53],[-50,-86],[2,-70],[-69,-61],[-92,-127],[-59,-10],[23,-44],[-76,-47],[-68,12],[-62,54],[-107,-13],[23,-46],[141,-52],[99,6],[95,-27],[49,19],[12,50],[114,54],[55,42],[35,91],[86,78],[-4,67],[-43,76],[23,73],[67,24],[165,26],[95,-112],[-7,-152],[101,106],[-43,141],[-203,73],[-71,1],[-75,-33],[-128,24],[15,48],[-41,64],[33,103],[59,104],[-99,134],[-51,40],[72,88],[162,66],[13,60],[-38,86],[71,-3],[46,-110],[-61,-106],[32,-35],[-21,-68],[96,-28],[182,-11],[-14,32],[-144,62],[-34,69],[103,27],[84,-42],[84,25],[-94,52],[125,45],[115,-3],[165,-64],[92,-76],[182,1],[-9,-72],[-73,-35],[2,-100],[65,33],[20,-83],[-43,-79],[101,28],[29,44],[-72,122],[52,108],[-46,62],[-68,11],[-57,68],[-167,55],[-21,49],[23,52],[-41,30],[9,106],[152,20],[212,2],[375,50],[76,-4],[-78,80],[-70,16],[86,46],[-6,35],[125,98],[202,75],[108,26],[159,13],[255,45],[-46,23],[172,42],[165,-8],[118,-40],[220,68],[17,52],[187,0],[95,49],[-13,63],[264,140],[108,23],[232,-53],[20,-22],[-141,-47],[98,-38],[151,7],[60,-24],[-97,-77],[121,-14],[50,46],[380,1],[185,-75],[22,-56],[77,-3],[55,-57],[-32,-123],[-84,-61],[-130,-61],[-250,-91],[-4,-34],[-90,-30],[-108,-74],[-114,-16],[-98,-106],[157,8],[259,65],[-13,55],[63,65],[98,-17],[134,-51],[80,8],[155,-35],[157,26],[243,-23],[140,-2],[-2,-75],[162,-59],[154,-9],[109,5],[99,-15],[58,14],[57,52],[-39,74],[133,49],[151,-51],[174,6],[88,-15],[141,-64],[44,-106],[-16,-48],[36,-38],[-16,-49],[-113,-1],[172,-212],[172,-83],[50,35],[82,163],[41,55],[106,-89],[83,-24],[220,55],[176,-60],[41,-33],[58,79],[109,-2],[25,-25],[95,10],[-35,46],[18,99],[-87,39],[57,37],[105,-1],[77,21],[-33,67],[154,-38],[179,-5],[242,-34],[137,-53],[-104,-79],[-77,-3],[8,-51],[75,11],[70,39],[93,81],[157,4],[135,-32],[57,-38],[-10,-32],[-118,-29],[206,-42],[121,-51],[114,-100],[158,10],[259,47],[264,-12],[161,-54],[47,-31],[34,-75],[-34,-96],[12,-19],[133,-39],[15,-112],[66,-42],[-4,104],[103,59],[220,15],[43,-24],[154,-5],[130,-19],[100,53],[64,-36],[19,-67],[124,-42],[37,-64],[114,8],[56,48],[-51,119],[-12,115],[255,-31],[89,-32],[82,11],[266,-3],[101,-47],[236,-56],[142,-89],[0,-860],[-68,-53],[-114,-49],[-95,25],[-65,32],[-120,-4],[60,-48],[90,25],[-4,-64],[64,-49],[53,8],[33,-65],[7,-97],[77,-73],[-9,-40],[38,-70],[-55,-81],[-136,50],[-82,9],[-80,-18],[-46,-35],[-130,-53],[-93,-62],[-110,-23],[-60,-68],[-34,14],[-57,-90],[-178,-118],[-42,-19],[-31,-103],[-46,22],[-44,81],[-49,35],[-134,-5],[-112,-37],[-32,-21],[-78,-96],[-24,24],[27,111],[-156,-85],[-8,-52],[-79,43],[-74,-4],[-45,-45],[-17,-117],[-33,-62],[-68,-68],[-38,-58],[-23,-81],[10,-35],[53,-33],[31,39],[52,-24],[9,-33],[-54,-72],[3,-120],[56,-27],[9,-100],[-35,-41],[-64,47],[-54,-37],[-38,-96],[-6,-65],[46,-133],[-47,-48],[-117,2],[-86,-77],[-28,-89],[6,-89],[-116,-73],[-42,-38],[-16,-54],[-2,-71],[-44,-109],[-71,-72],[-40,-59],[-56,-53],[-31,111],[-1,96],[-15,131],[-17,26],[-16,83],[-24,202],[-35,206],[-8,106],[20,159],[33,137],[91,100],[31,69],[-23,62],[81,11],[27,49],[67,1],[102,87],[68,84],[31,71],[86,92],[38,17],[73,88],[38,28],[31,58],[119,80],[91,29],[-12,43],[48,51],[-25,26],[27,55],[23,127],[48,39],[-42,50],[-114,-32],[-40,-167],[10,-48],[-86,20],[-166,-154],[-34,-46],[-66,18],[19,42],[-53,23],[-8,36],[53,110],[-16,24],[-78,-40],[-59,44],[-123,-37],[-74,10],[-97,-70],[-7,-46],[-70,-60],[-25,-48],[-109,-88],[-83,-119],[-18,-67],[53,2],[72,-39],[0,-38],[-57,-10],[-30,16],[-55,-31],[-39,33],[-41,-2],[-59,-61],[-59,19],[-39,-25],[-93,-8],[-25,46],[106,17],[-62,80],[-126,14],[-100,40],[-71,-31],[17,-31],[-164,-22],[-44,-26],[-132,35],[-21,-46],[-42,-5],[-46,47],[-132,-8],[-76,8],[-83,-9],[-75,-29],[-69,-52],[-33,-56],[-94,-72],[-38,-45],[-29,-88],[-55,-27],[-68,-79],[-33,-15],[-64,-65],[-119,-181],[-110,-98],[-41,-26],[-36,-49],[-78,-52],[-34,-38],[0,-46],[73,-33],[116,9],[-14,-152],[58,-24],[27,101],[33,-91],[-39,-72],[119,12],[45,33],[5,126],[76,-27],[48,18],[58,-49],[21,-53],[81,-70],[49,-68],[-39,-63],[21,-12],[-11,-104],[39,-42],[-15,-56],[-53,-66],[-30,-85],[-26,-151],[7,-91],[-9,-48],[7,-69],[-26,-119],[-17,-117],[-21,-44],[-77,-95],[-46,-111],[-51,-72],[-31,-112],[-80,-160],[-108,-142],[-82,-149],[-33,-26],[-48,-118],[-43,-68],[-54,-51],[-121,-102],[-68,-29],[-55,40],[-50,1],[4,78],[-38,-26],[-37,18],[-70,-126],[-40,-23],[-18,-48]],[[27986,37500],[-17,-54],[-143,-35],[-125,45],[285,44]],[[28335,37516],[175,-52],[-91,-32],[-254,-38],[14,-36],[-169,-21],[-122,42],[124,44],[189,13],[7,57],[127,23]],[[30487,36487],[-183,-56],[-322,-69],[-284,-77],[-139,-114],[-121,-23],[-87,-43],[-11,-81],[-83,-25],[13,-36],[-112,-106],[-66,-15],[-175,34],[-87,-22],[-66,91],[108,42],[116,147],[113,72],[72,94],[86,27],[135,76],[267,58],[18,36],[218,-10],[184,31],[172,52],[24,29],[181,57],[117,-16],[51,-79],[-139,-74]],[[28955,35844],[137,-23],[-43,-91],[-44,-5],[-45,-66],[20,-54],[-27,-62],[91,-130],[105,-92],[90,-43],[-59,-31],[-179,23],[-134,-3],[-108,30],[-17,52],[37,32],[-59,31],[-5,44],[-191,-11],[-46,63],[18,65],[82,12],[45,38],[37,81],[-48,15],[88,52],[-6,38],[69,32],[129,21],[63,-18]],[[34026,37549],[160,-61],[-99,-58],[19,-64],[-322,-25],[-85,-26],[-205,37],[-84,40],[130,38],[84,57],[-24,43],[337,64],[89,-45]],[[34167,37348],[199,-31],[95,-53],[-47,-100],[31,-79],[-61,-28],[-232,-2],[-148,42],[-209,24],[-53,61],[-141,20],[235,132],[331,14]],[[34808,37149],[97,-23],[96,-59],[85,-13],[21,-70],[-59,-33],[-137,-20],[-310,-15],[-138,-48],[-97,14],[90,66],[32,75],[161,147],[159,-21]],[[39381,36397],[110,-28],[11,64],[56,32],[120,-51],[151,-9],[193,-61],[-46,-69],[-106,-50],[-213,-42],[-37,31],[-212,-32],[-81,22],[-63,-64],[-124,31],[-133,96],[35,25],[-3,88],[29,35],[161,64],[152,-82]],[[40211,36297],[202,9],[5,-31],[200,-11],[67,-61],[-129,-37],[-185,11],[-239,82],[79,38]],[[44157,35299],[-29,49],[33,40],[136,68],[0,-120],[-140,-37]],[[39715,31693],[26,-56],[-7,-76],[38,-109],[12,-73],[0,-77],[-19,-58],[-2,-58],[21,-110],[16,-25],[45,-261],[28,-85],[28,-128],[-18,-22],[-49,22],[-77,-25],[-65,-247],[-2,-73],[14,-44],[41,-69],[34,-146],[-56,10],[-31,21],[-21,-32],[-24,-102],[-31,-17],[-16,96],[26,152],[-9,98],[27,93],[-6,61],[-33,101],[19,72],[15,108],[1,142],[-10,65],[17,156],[-25,65],[-35,47],[-8,118],[18,62],[9,106],[-7,66],[40,34],[47,-10],[18,114],[-42,69],[53,25]],[[138,35464],[65,-10],[101,-59],[-36,-43],[-207,-31],[-61,15],[0,120],[138,8]],[[28333,34939],[-78,-71],[-118,-28],[-48,67],[5,51],[38,36],[74,17],[127,-72]],[[29946,37467],[-31,58],[279,20],[-103,-66],[-145,-12]],[[29280,37340],[-70,3],[-39,78],[173,-6],[-64,-75]],[[29798,37496],[-11,-48],[-126,-43],[-215,22],[35,65],[247,17],[70,-13]],[[29172,37390],[0,-56],[-108,-4],[-22,58],[130,2]],[[29262,37652],[25,-64],[-283,-14],[169,77],[89,1]],[[28882,37558],[178,-19],[188,-52],[-244,-34],[-203,39],[81,66]],[[33553,37244],[-192,48],[129,31],[201,-31],[-138,-48]],[[39499,35996],[-61,-18],[-50,58],[92,20],[19,-60]],[[39644,35973],[143,-72],[-19,-76],[-196,20],[-115,31],[52,91],[135,6]],[[39122,31846],[6,-60],[-41,-36],[-27,79],[17,40],[45,-23]],[[42968,35026],[-105,18],[-69,38],[70,39],[126,-29],[-22,-66]],[[42283,32617],[16,90],[107,21],[-2,-49],[-41,-11],[-80,-51]],[[41334,30795],[-16,-22],[-67,-24],[-3,45],[57,22],[25,62],[29,-2],[-25,-81]],[[40567,29772],[-3,43],[84,82],[-18,-61],[-63,-64]],[[40433,29701],[-41,-23],[-43,-49],[-18,44],[102,28]],[[40139,29521],[-12,-56],[-42,-67],[-13,13],[55,110],[12,0]],[[36100,36084],[-74,-67],[-140,39],[54,60],[160,-32]],[[30845,35797],[-93,-2],[9,60],[116,34],[85,-74],[-117,-18]],[[31701,35621],[-90,16],[104,59],[76,-33],[-90,-42]],[[31931,35715],[-104,38],[62,44],[47,-38],[-5,-44]],[[29587,35104],[-1,-46],[-107,3],[-77,33],[-59,83],[66,43],[178,-116]],[[24727,31887],[-44,-65],[73,-11],[5,73]],[[24761,31884],[103,-44],[61,-1],[33,-49],[-19,-76],[11,-29]],[[24950,31685],[-74,0],[-241,14],[-74,8]],[[24561,31707],[42,64],[3,37],[67,17],[47,64]],[[24720,31889],[7,-2]],[[26455,29792],[46,-84],[23,-9],[75,31],[50,-13],[-22,-72],[-65,-13],[-49,20],[-47,-65],[-46,1],[-53,-59],[-46,-33],[-56,37],[19,78],[-7,41],[-45,22],[-33,33],[-50,12],[78,76],[64,44],[-9,32]],[[25804,29688],[-19,-92],[-62,-18],[-50,-102],[-7,-121]],[[25666,29355],[-66,17],[-20,37],[-57,7],[-42,32],[-107,-35],[-49,-53],[-39,-21],[-131,27],[-25,-7],[-124,25],[-24,44],[-40,36]],[[24639,29874],[46,6],[66,37],[21,67],[57,116],[27,81],[35,49],[48,15],[24,33]],[[25424,30346],[44,-23],[44,-113],[64,-115],[27,-30],[21,-116],[-18,-105],[14,-114]],[[28400,25195],[5,22],[-10,122],[29,117],[32,38],[35,-55],[-7,-83],[14,-83],[-2,-40],[-19,-65],[-20,-13]],[[20033,26968],[61,-24],[-17,-24],[-41,16],[-3,32]],[[21237,27914],[-52,-38],[-94,26],[-50,-20],[23,87],[3,67],[-11,49],[9,50],[-13,72],[-37,-14],[-5,46],[-27,16],[12,133],[28,45],[38,126],[-6,14],[28,186],[-19,137],[4,63]],[[24466,30621],[-33,81],[-66,34],[-19,-14],[-89,89],[-58,-61],[-83,112],[-80,44],[-45,45],[-22,-34]],[[23971,30917],[25,86],[-51,127],[19,55],[-25,61],[8,37],[-60,77],[35,89],[-19,98]],[[23903,31547],[-6,31]],[[23897,31578],[0,17]],[[23897,31595],[62,15],[181,60],[46,58],[87,38],[130,24],[43,-89],[37,-18],[78,24]],[[24950,31685],[30,-2],[54,-45],[4,-45]],[[25038,31593],[50,-200],[-5,-80],[-45,-25],[-35,-64],[56,-48],[-13,-72],[7,-43]],[[37050,23840],[91,-71],[56,24],[-15,-74],[-3,-68],[10,-70],[35,-68],[-11,-67],[-36,-140],[-43,-23],[-24,-32],[2,-58],[-27,-76],[38,-129],[-6,-56],[19,-80],[50,-41],[-1,48],[36,38],[45,-15],[46,-111],[25,48],[36,-18],[-15,-78],[26,-57],[-4,-35],[44,-16],[-10,-103],[-23,27],[5,55],[-36,-6],[-38,29],[-19,87],[-32,33],[-37,69],[-11,-55],[19,-57],[-18,-30],[-14,51],[-41,65],[-36,33],[-34,-21],[-7,-29],[-30,-15],[-44,53],[-25,-17],[-3,84],[38,67],[-5,49],[-42,11],[6,-62],[-19,-7],[-23,75],[-21,12],[-18,128],[-5,89],[-15,37],[7,70],[41,-61],[28,38],[-10,68],[12,93],[2,106],[-8,44],[30,191],[26,21],[36,3]],[[36583,21606],[5,60],[66,116],[28,23],[112,221],[32,49],[-3,59],[25,99],[8,-75],[19,-88],[-11,-32],[-38,-34],[-11,-42],[-51,-32],[-1,-33],[-41,-112],[-37,-34],[-18,-49],[-58,-74],[-26,-22]],[[37221,22303],[54,-19],[28,2],[-18,-92],[-26,-28],[-4,-36],[-30,-29],[-40,-16],[-26,-36],[-3,91],[11,50],[6,119],[23,29],[25,-35]],[[37299,21743],[-16,-1],[-16,57],[-38,36],[-18,46],[7,59],[48,28],[-5,91],[20,84],[34,24],[31,-16],[7,-31],[-50,-201],[-1,-57],[19,-54],[-22,-65]],[[37477,22243],[44,7],[12,-35],[-2,-94],[21,-44],[9,-70],[-28,-51],[-31,30],[-5,59],[6,75],[-15,40],[-27,-9],[-14,135],[30,-43]],[[37559,22504],[10,-46],[26,-28],[-9,-53],[6,-89],[15,-91],[-49,3],[-24,43],[-19,95],[-25,54],[-36,49],[-11,72],[33,-10],[72,10],[11,-9]],[[37001,22712],[61,-10],[40,-66],[-6,-65],[8,-43],[-18,-74],[-20,-18],[-38,64],[-20,101],[-31,63],[-13,57],[37,-9]],[[37653,21800],[23,-10],[16,-95],[-22,-48],[27,-31],[12,-62],[-3,-81],[17,-34],[1,-94],[-18,-52],[-44,-61],[-31,131],[-52,-118],[-3,-23],[26,-49],[10,-107],[-27,-69],[-20,-7],[-14,91],[-17,-39],[-49,28],[-52,52],[-17,37],[-11,130],[27,88],[-30,59],[-37,34],[-21,-2],[-13,-88],[-26,26],[-41,0],[-6,43],[-36,-11],[-40,-147],[-26,-8],[-5,51],[15,36],[10,98],[26,48],[70,28],[13,53],[40,38],[12,29],[30,-18],[21,-41],[3,-54],[43,18],[28,74],[34,-10],[17,90],[34,-23],[8,35],[36,-3],[4,28],[-15,116],[13,20],[44,-54],[16,-42]],[[37329,21828],[-7,9],[9,104],[40,112],[27,107],[15,-8],[0,-75],[-13,-59],[-38,-65],[-18,-94],[-15,-31]],[[37479,21902],[-28,-34],[-53,-2],[-14,43],[44,69],[20,6],[29,-29],[2,-53]],[[37172,21165],[-17,-3],[-18,43],[30,28],[29,-22],[-24,-46]],[[37150,22726],[24,-17],[-13,-57],[-22,27],[11,47]],[[37371,22451],[24,-26],[16,-45],[-7,-32],[-69,82],[-19,44],[24,24],[31,-47]],[[37450,22746],[-22,-22],[-17,29],[10,69],[13,22],[24,-46],[-8,-52]],[[13539,18824],[-33,-14],[-21,33],[-48,-5],[-22,-39],[-107,-34],[-50,-54],[-55,-69],[-23,-13],[-11,-112],[-32,-102],[12,-80],[-44,-47],[-32,-50],[1,-94],[-27,-10],[-4,-39],[34,-49],[-6,-34],[27,-90],[24,-29],[7,-39],[40,-74],[1,-28],[-30,-64],[49,0],[53,-22],[25,-108],[70,0],[46,8],[74,107],[-1,-336],[52,-8],[32,26],[47,-5]],[[13587,17350],[110,-340],[-9,-41],[-27,-42],[-1,-136],[-11,-40],[24,-114],[-44,-87],[-17,-80],[25,-60],[-24,-60]],[[13613,16350],[-39,53],[-31,-29],[-9,-94],[24,0],[45,-45],[-4,-33],[43,-9]],[[13642,16193],[19,-4]],[[13661,16189],[-7,-64]],[[13654,16125],[-50,-122],[-23,-24],[14,-68]],[[13595,15911],[-42,-43],[7,-43],[-16,-67],[-60,-31]],[[13484,15727],[-65,91],[-48,55],[-8,57],[-16,28],[-72,64],[-43,65],[-41,20],[-74,68],[-52,33],[-40,53],[-118,110],[-53,113],[-49,58],[-44,110],[-1,40],[14,95],[-79,256],[-40,64],[-8,87],[-51,82],[-13,99],[-55,164],[-32,158],[-17,47],[-22,118],[-31,90],[-45,82],[-47,169],[-43,91],[-86,81],[-41,49],[-5,26],[37,40],[-7,62],[-28,68],[7,30],[-28,79],[7,76],[59,130],[59,75]],[[12265,19010],[18,-107],[-41,-42],[3,-72],[13,-15],[30,37],[42,-40],[19,5],[38,-104],[31,-14],[27,72],[21,22],[1,52],[32,121],[9,76],[20,17],[40,81],[145,92],[73,95],[64,132],[18,125],[21,2],[-1,87],[-20,32],[-7,55],[24,11]],[[12885,19730],[18,13],[41,-33],[48,-84],[21,-85],[47,-50],[24,-11],[21,-63],[0,-35],[37,-30],[5,-98],[33,-29],[32,13],[28,-17],[51,22],[29,38],[43,-40],[60,26],[100,-97],[3,-20],[-81,-226],[25,-19],[23,12],[46,-93]],[[14992,15328],[19,-66],[22,-117],[-19,-230],[20,-21],[71,-12],[38,-17],[37,8],[13,30],[41,-40],[42,-5],[29,-80],[-5,-31],[36,-259],[47,4],[44,35],[47,-51],[-28,-258],[-18,-78]],[[15428,14140],[-7,-161],[-19,-75],[-38,-61],[-36,-17],[-44,-90],[-46,21],[-34,-51],[-20,19],[-63,-1],[-86,34],[-44,10],[-54,-9],[-4,26],[39,74],[16,50],[-1,54],[15,54],[23,29],[38,132],[-29,59],[-39,35],[-28,3],[-124,111],[-64,80],[-117,51],[-45,66],[-73,83],[-70,163],[-35,45]],[[14439,14874],[47,256],[0,111],[44,111],[19,90],[216,76],[113,3],[112,-117],[2,-76]],[[40971,18710],[-36,56],[2,48],[-12,63],[-50,107],[102,-131],[13,-33],[-19,-110]],[[40841,18811],[25,18],[36,-28],[-1,-86],[-18,-49],[-30,-10],[-2,-32],[18,-46],[-34,-46],[-43,3],[-35,-81],[-52,-33],[-47,-45],[-95,-3],[-33,47],[-32,-11],[-39,46],[-50,33],[-7,49],[80,13],[31,-19],[72,10],[34,21],[23,-23],[41,3],[40,21],[23,69],[23,33],[42,17],[-1,54],[-13,64],[5,32],[39,-21]],[[39495,17752],[0,271],[0,215],[-14,36],[14,87],[0,465],[0,355]],[[39495,19181],[26,-4],[126,-100],[86,-52],[74,-24],[62,-78],[57,-9],[75,-115],[31,-8],[53,-96],[10,-142],[68,-32],[78,-66],[42,-7],[23,-24],[29,-56],[5,-77],[-16,-13],[-73,0],[-19,-45],[28,-99],[66,-109],[49,-50],[15,-99],[25,-31],[16,-78],[26,-8],[38,16],[18,-12],[-5,-74],[33,-40],[61,-16],[-26,-31],[14,-46],[98,-54],[-22,-37],[-1,-46],[-57,13],[-33,49],[-100,22],[-28,19],[-40,-4],[-76,26],[-26,34],[-18,52],[-48,63],[-13,65],[-27,15],[-63,171],[-19,37],[-37,27],[-32,5],[-61,28],[-25,35],[-32,17],[-35,-44],[-40,20],[5,-59],[-40,-56],[-64,-24],[2,-36],[35,-72],[-4,-35],[-88,-81],[-52,35],[-74,-9],[-27,13],[-33,-15],[-20,22]],[[40245,19324],[17,-49],[-57,7],[-11,39],[51,3]],[[40659,19170],[-33,0],[-15,32],[23,28],[24,-18],[1,-42]],[[40671,17702],[43,-36],[2,-34],[-57,9],[12,61]],[[41339,18286],[-30,-39],[-46,31],[-17,43],[0,48],[-24,22],[-28,52],[-5,94],[31,1],[58,-133],[44,-52],[17,-67]],[[12628,21654],[-4,-51],[24,-86],[-17,-72],[-29,-36],[-20,-1],[-19,-67]],[[12563,21341],[-33,69],[-31,113],[38,70],[-44,15],[-5,40],[-71,81],[-44,2],[-30,-35],[-8,-56],[-46,-54],[-30,-13],[-11,-46],[47,-90],[8,-37],[-34,-16],[-18,-33],[-50,-12],[-24,107],[-28,-18],[-29,21],[-28,91],[-117,40],[-23,-16],[-2,-39]],[[11950,21525],[-14,65],[17,19],[-7,63],[21,47],[-24,23],[0,86],[17,31],[29,-3]],[[11989,21856],[24,-33],[16,-87],[20,-21],[37,5],[52,-39],[64,23],[35,43],[52,28],[68,85],[69,-19],[10,-18],[53,-5],[52,-37],[31,-37],[16,-39],[40,-51]],[[31594,27581],[35,-121]],[[31629,27460],[-6,-26],[-50,-56],[-68,-15],[-41,-36],[-64,30],[-109,27],[-42,-24],[-18,-72],[20,-30],[3,-88],[18,-62],[-19,-58],[38,-55],[7,-49],[38,-2],[-1,-53],[40,-12],[43,-40],[-10,-31],[-63,-42],[-28,-52],[10,-54],[-9,-61],[14,-33],[-51,-59],[-39,-73],[4,-47],[-15,-28],[-53,-35],[-19,-84],[-40,-115],[-69,-60],[-20,-73],[-28,-54],[-10,-47],[-40,-20],[-99,-31],[-14,46],[-28,19],[-31,-39],[-31,-82],[-40,-66],[-12,-81],[32,-39],[40,-10],[12,-24],[-10,-88],[23,-80],[48,-9],[0,-53],[49,-162],[-1,-63],[-40,-35],[-21,39],[-55,-28],[-10,-21],[-37,-5],[-19,22],[-103,-1],[0,-66],[-55,-8],[-14,-16]],[[30536,24990],[-61,-10],[-44,80],[-17,127],[-58,23],[0,81],[-34,76],[-92,-48],[-33,4],[-78,-15],[-14,-27],[-63,33],[-53,12],[-28,-39],[-113,10],[-32,-26],[-70,0],[-19,14]],[[29727,25285],[10,125],[10,16],[14,88],[46,25],[24,45],[90,23],[10,46],[-10,82],[-50,0],[6,54],[-9,111],[3,44],[-24,7],[-26,40],[-57,28],[-33,54],[-37,128],[-59,106]],[[29635,26307],[201,-98],[134,19],[66,-23],[52,38],[70,-1],[133,60],[16,107],[-2,63],[13,67],[65,86],[45,-19],[55,27],[-19,36],[78,66],[47,-1],[34,-37],[50,67],[-5,109],[21,54],[11,74],[52,21],[45,56],[-51,116],[15,34],[52,-22],[78,21],[6,71],[-16,35],[80,143],[-2,49],[-22,94],[-29,46],[53,86],[66,58],[58,28],[141,13],[28,-14],[67,44]],[[31321,27880],[99,-24],[9,-39],[51,-19],[17,-59],[-4,-80],[38,-45],[63,-33]],[[29087,25236],[31,-111],[33,-71],[38,-49],[48,-27],[63,-19],[25,-20],[31,4],[24,-27],[32,-85],[62,-121],[36,-16],[-3,-63],[-53,-158],[-59,-85],[-51,-155],[-33,4],[-13,32],[-30,-72],[-18,-140],[12,-129],[-78,-25],[-43,-33],[-21,-36],[-13,-93],[-35,-47],[-97,-24],[-27,-56],[4,-45],[-28,-75],[-36,-17],[-25,15],[-62,-6],[-56,-54],[-64,-24]],[[29049,25474],[10,32],[33,26],[-7,-114],[-9,-34]],[[23550,32712],[-68,32],[-24,45],[-57,-83],[-69,-14],[-135,-142],[-44,-34],[-86,-27],[-57,0],[-42,53],[-33,1],[-84,56],[-24,45],[5,54],[37,-3],[14,50],[-15,36],[-60,-41],[-28,13],[21,92],[58,4],[-5,44],[-76,-6],[16,161],[-24,30],[-15,221],[74,131],[72,30],[29,36],[115,-2],[-34,78],[69,31],[65,-1],[37,97],[55,24],[59,-16],[46,12],[47,-42],[84,6],[-5,36],[-98,-23],[-27,75],[106,122],[81,57],[15,50],[100,60],[14,85],[86,99],[16,114],[64,57],[139,154],[5,68],[64,57],[92,29],[64,58],[8,35],[91,39],[39,66],[49,12],[23,70],[86,19],[26,41],[107,29],[135,2],[29,62],[158,44],[44,-60],[169,132],[30,66],[136,-33],[-61,-66],[25,-46],[127,125],[90,8],[44,32],[185,-51],[137,-70],[47,-4],[43,-55],[-165,-65],[23,-55],[133,13]],[[25713,34903],[-15,34],[60,65],[-24,44],[-90,33],[-64,53],[-96,-35],[-72,2],[-63,-57],[-32,-93],[0,-53],[-100,-87],[-133,47],[-66,-35],[-112,16],[-97,121],[-61,-42],[-62,-9]],[[24020,34886],[2,-72],[-105,21],[103,51]],[[24518,35133],[41,-22],[-99,-85],[-54,41],[112,66]],[[24302,35029],[62,-20],[-7,-67],[-120,-13],[4,49],[61,51]],[[24088,34802],[18,-38],[-121,-18],[60,68],[43,-12]],[[24807,37005],[54,-4],[31,-77],[142,-17],[-41,-35],[138,-20],[-49,-88],[-128,-49],[-91,49],[-124,-13],[68,105],[-70,31],[-100,92],[170,26]],[[24720,37368],[182,-35],[9,68],[543,-53],[41,-56],[-147,-63],[-44,-47],[-209,-46],[-154,30],[-226,15],[-262,79],[-47,66],[248,75],[66,-33]],[[24214,37292],[97,-4],[124,-69],[62,-91],[169,-6],[5,-33],[109,-53],[-199,-25],[-164,-132],[-26,-110],[-74,-27],[-114,-180],[-287,144],[-44,46],[112,54],[-139,36],[-13,38],[162,27],[178,77],[-222,-19],[-65,-32],[-118,-7],[-170,105],[34,59],[-103,32],[-58,86],[16,60],[232,-4],[77,-65],[206,43],[94,1],[119,49]],[[37730,28046],[-32,22],[-53,-1],[-21,24],[-39,-56],[-32,68],[-62,20],[47,94],[30,23],[-19,46],[31,91],[-6,52],[-72,51],[-27,7],[-24,47]],[[37451,28534],[3,22],[62,78],[65,44],[38,46],[32,8],[62,99],[13,63],[38,30],[34,-55],[114,-31],[22,38],[-30,94],[108,11],[34,39],[14,43],[50,14],[22,114],[42,-23],[1,-32],[34,-46]],[[38229,29039],[-28,0],[-48,-57],[-38,-73],[-9,-48],[10,-41],[-7,-98],[-45,-29],[-29,-52],[-50,-38],[-49,-62],[-42,-9],[-49,-46],[-2,-71],[-19,-55],[48,-28],[73,-101]],[[23823,22762],[56,-138]],[[23879,22624],[16,-152],[40,-19],[13,-41],[-8,-136],[-44,-49],[-38,-28],[-44,-118],[-15,-95],[-18,-29],[-11,-109],[-31,-25],[-15,-119],[-27,-57],[-22,-7],[-21,-68],[-27,-152],[-31,-70],[12,-34],[-35,-50],[-3,-42],[-28,-47],[-21,-10],[-15,57],[-52,80],[-16,-38],[-34,0],[-7,23],[-45,-51],[-15,-51],[-35,-46],[-39,-68],[-15,-50],[-17,-128],[-30,-97]],[[23201,20798],[-32,-44],[-140,-9],[-18,-26],[-31,-12],[-84,-11],[-60,78],[-25,106],[0,59],[-19,18],[-45,120],[-53,70],[-38,14],[-78,-1],[-96,-8]],[[22482,21152],[8,75],[-7,59],[4,91],[-8,105],[6,199],[5,59],[33,8],[11,80],[24,72],[28,28],[11,56],[-8,24],[31,74],[-14,104],[-28,69],[13,66]],[[22591,22321],[6,183],[37,54],[17,61],[8,89],[34,41],[49,25],[52,-1],[30,26],[110,-59],[52,-109],[31,-24],[37,24],[53,50],[38,-10],[80,-84],[56,-19],[51,-2],[38,71],[32,29],[95,23],[56,-4],[71,-36],[16,-22],[42,1],[24,51],[82,75],[35,8]],[[22591,22321],[-29,40],[-59,107],[-63,-32],[2,-71]],[[22442,22365],[-36,84],[15,57],[-17,38],[-35,-22],[-28,5],[-71,89],[-1,112],[-44,34],[-24,59],[-5,60],[-27,55],[6,91]],[[22175,23027],[9,15],[83,1],[41,63],[210,15],[7,19],[55,-16],[2,28],[44,60],[17,96],[13,36],[14,140],[-1,471]],[[22669,23955],[187,64],[11,10],[142,218],[60,88],[106,112],[128,134],[191,201],[127,134]],[[23621,24916],[186,-74],[93,-124],[92,83]],[[23992,24801],[23,-236],[2,-87],[52,-125],[-3,-48],[43,-74],[2,-21],[-26,-87],[-12,-240],[-20,-418],[-137,-254],[-69,-172],[-20,-74],[-24,-55],[20,-148]],[[11916,23045],[-17,-54],[13,-90],[-27,-75],[-19,-149],[6,-82],[1,-117],[-18,-27],[-14,-78],[14,-64],[-26,-75],[4,-37],[24,-47]],[[11857,22150],[-21,-38],[-44,8],[-22,44],[-44,16],[-25,-24],[-88,53],[-15,-27]],[[11598,22182],[-27,59],[-62,89],[-36,92],[-45,61],[-67,103],[16,30],[25,-13]],[[11402,22603],[34,2],[17,51],[26,20],[-3,99],[50,1],[22,49],[45,-33],[13,31],[61,79],[3,40],[21,50],[24,8],[18,-28],[96,29],[56,47],[31,-3]],[[43450,10693],[14,-1],[62,65],[70,-7],[-8,-64],[-18,-42],[24,-52],[-91,-170],[-40,-101],[-62,-62],[11,-80],[33,-12],[-19,-47],[-54,36],[-12,-22],[-89,-63],[-27,-5],[-26,-84],[-14,-114],[-22,-39],[-29,-102],[9,-41],[-44,-16],[-90,-134],[-72,-17],[-41,14],[-48,-9],[-23,54],[-43,-1],[-20,38],[-38,-11],[-79,10],[13,93],[-12,50],[56,128],[73,81],[72,112],[54,20],[44,47],[61,41],[79,106],[54,40],[63,101],[28,147],[43,31],[22,46],[16,109],[16,42],[45,52],[16,-56],[27,-16],[16,-95]],[[42838,9468],[12,-46],[-41,-24],[-34,9],[40,98],[23,-37]],[[43469,12086],[71,-20],[58,-48],[28,-74],[-19,-42],[50,-118],[3,-66],[-12,-50],[71,-34],[30,-46],[-11,160],[52,-106],[14,-105],[14,-46],[82,-54],[69,-22],[83,96],[32,-3],[33,-27],[-18,-59],[-15,-130],[-36,-37],[-1,-93],[-69,14],[-34,-23],[-22,-40],[19,-67],[-33,-106],[-56,-112],[-49,-120],[-83,-87],[-18,42],[-32,-3],[-33,31],[65,147],[11,73],[-30,74],[-54,30],[-27,38],[-51,29],[-21,42],[10,39],[68,39],[24,39],[15,123],[26,92],[-25,78],[7,111],[-38,0],[-11,48],[8,60],[-24,45],[-38,16],[-102,213],[9,41],[-60,123],[32,-17],[38,-88]],[[22886,30893],[-37,5],[17,77],[-40,33],[-123,44],[-34,-19]],[[22669,31033],[-79,14],[110,188],[24,110],[84,63],[21,38],[65,30],[93,8],[47,-35]],[[23034,31449],[-1,-62],[-22,-81],[3,-55],[-35,-62],[-1,-41],[-98,-33],[31,-60],[-37,-111],[12,-51]],[[22870,31291],[-31,26],[-67,-13],[-3,-61],[56,-20],[45,68]],[[22669,31033],[-40,-40],[-68,38]],[[22561,31031],[108,2]],[[22862,31269],[-37,-36],[-41,17],[61,46],[17,-27]],[[32990,25871],[5,-27],[-20,-135],[22,-89],[-14,-65],[-94,-15],[-33,42],[-39,-26],[-41,31],[-78,14],[-21,44],[-33,-22],[-75,66],[-9,57],[-64,42],[-33,-25],[-47,19],[-19,-21],[-68,33],[-7,34],[-28,-1],[-57,54],[-17,-11],[-66,68],[-18,35],[-71,69],[-21,-8],[-43,48],[7,59],[34,138],[62,97],[13,-1]],[[32117,26375],[30,-16],[20,54],[27,11],[50,-14],[14,-46],[85,-95],[38,-16],[52,-94],[44,21],[20,-13],[16,-67],[59,-70],[55,-1],[-4,-60],[68,-9],[34,-76],[15,34],[41,-34],[18,34],[72,-54],[67,-4],[52,11]],[[25257,15848],[-44,-6],[-46,-51],[-21,17],[-57,-55],[-39,-49],[-35,93],[-24,4],[-134,-42],[-128,-26],[1,-313],[0,-268],[-1,-227],[-122,0],[0,-242],[0,-367]],[[24172,13473],[-89,143],[-47,127],[-25,131],[3,41],[-24,62],[-15,130],[-1,152],[-41,183],[-4,202],[-8,69],[15,58],[-25,112],[-43,93],[-17,65],[-47,122],[-35,160],[-16,35],[-72,242],[-44,84],[-40,119],[-5,55],[1,110]],[[23593,15968],[34,17],[66,-9],[68,54],[37,-9],[66,-84],[224,0],[95,0],[229,-2],[40,-67],[44,-25],[69,-11],[93,-4],[28,-24],[60,9],[38,-10],[241,79]],[[26234,16678],[53,-122],[86,35],[21,-41],[5,-154],[-36,-129],[18,-69],[110,-197],[-1,64],[-14,54],[23,88],[49,22],[9,140],[2,165],[-73,170],[-48,82]],[[26438,16786],[-11,123],[2,70],[-12,48],[4,59],[26,62],[3,65]],[[27127,17457],[19,-43],[-16,-23],[7,-66],[-17,-58],[22,-287],[0,-162],[-5,-19],[12,-204],[22,-17],[2,-71],[-27,-71],[-7,-78],[-55,-111],[-33,-101],[-74,-78],[-19,-40],[-56,-23],[-60,-36],[-111,-109],[-37,-56],[-5,-29],[-44,-86],[-17,-55],[-34,-16],[-59,-70],[-35,-73],[-51,-70],[-24,-2],[-7,-128],[35,-88],[17,-85],[1,-44],[17,-56],[8,-85],[-2,-79],[28,4],[-7,-61],[11,-67],[-25,-184],[21,-5],[-13,-77],[-35,-81],[-67,-61],[-95,-54],[-60,-43],[-69,-84],[-24,-79],[42,-54],[-6,-131]],[[21875,27459],[53,-78],[-7,-31],[17,-95],[4,-188],[21,-97],[54,-90],[-21,-28],[-4,-55],[-145,9],[-59,-19],[-11,-46],[-52,-28],[-47,-10],[-1,-102],[16,-55],[-35,-7],[-41,-47],[-80,-51],[-26,-66],[-33,-46],[-40,-19],[-89,-13],[-17,-53],[-57,12],[-48,-51],[-24,-9],[-120,-139],[-3,-233]],[[21080,25824],[-17,0],[8,-102],[-30,-22],[-51,-1],[-39,-50],[-37,11],[-27,-11],[-37,29],[-46,4],[-57,-28],[6,-54],[-36,-75],[-10,-42],[-39,-16],[-34,-202],[-15,-61],[-78,-94],[-30,-93],[-57,-41],[-31,-84],[-18,-143],[-6,-98],[-51,-99],[-5,-57],[-21,-33],[-76,0],[-36,9],[-120,-5],[-34,-11]],[[20056,24455],[9,106],[17,57],[34,38],[20,57],[18,111],[50,143],[-12,20],[39,50],[49,90],[18,15],[23,71],[7,117],[34,114],[13,72],[57,52],[46,54],[49,202],[28,57],[119,47],[68,55],[43,73],[73,77],[78,164],[23,65],[2,75],[-28,59],[9,155],[16,63],[40,82],[13,107],[50,76],[30,59],[36,41],[91,58],[81,73],[18,36],[50,141],[52,221],[65,32],[36,-101],[33,-41],[63,-27],[78,26],[60,-9],[29,37],[16,-61],[76,-5]],[[21080,25824],[0,-81]],[[21080,25743],[0,-283],[-208,0],[-202,0],[0,-296],[-1,-259],[-73,-43],[-51,-60],[-17,-54],[9,-57],[10,-255],[-273,0],[-213,-1],[-10,-114]],[[20051,24321],[-6,11],[11,123]],[[24652,29154],[-39,-21],[4,-41]],[[24617,29092],[-34,-15],[-16,33],[-46,-100],[8,-66]],[[24529,28944],[-56,83],[-46,40]],[[24427,29067],[-10,28]],[[24417,29095],[13,25],[-12,65],[61,116],[31,8]],[[36506,30690],[-90,-249],[-16,-25],[-4,-72],[-32,-26],[4,-41],[41,-56],[41,37],[66,3],[72,-48],[52,74],[89,-1],[33,-50],[117,-133],[23,-64],[-24,-56],[-172,25],[-78,-43],[-33,11],[-9,-46],[-58,5],[-37,-21],[-43,-89],[4,-19],[-69,-75],[-143,-20],[-13,-36],[-47,-57],[-55,-43],[-122,38],[-23,31],[-63,1],[-34,-52],[-27,-101],[54,-95],[12,-49],[-48,-47],[-66,-34],[-75,-125],[-118,-70],[-93,-8],[-63,7],[-173,-35],[-193,-121],[-18,-27],[-68,9],[0,48],[-97,-27],[-79,56],[-135,46],[-32,55],[-213,46],[-63,-24],[-279,49],[-101,-16],[-10,46],[-48,61],[-47,164],[-21,12],[-1,59],[-37,-4],[-41,20],[-64,65],[-29,6],[-37,50],[-107,30],[-45,-6],[-103,15],[-66,31],[-21,-5],[-27,72],[6,45],[36,67],[-11,52],[15,56],[-20,94],[-46,72],[-23,86],[-35,44],[-31,-11],[-26,39],[-63,0],[-58,48],[-8,36],[-66,37],[6,40],[-35,32],[9,61]],[[23064,29357],[-8,-5]],[[23056,29352],[8,5]],[[10195,25452],[-10,-82],[-34,-126],[-20,-137],[-12,-238],[2,-80],[-14,-69],[13,-132],[14,-92],[17,-46],[47,-171],[51,-94],[31,-70],[20,-115],[38,-64],[20,-66],[79,-12],[47,-40],[31,-75],[44,4],[78,52],[82,8],[22,32],[90,23],[21,-54],[33,-3],[32,37],[-8,60],[48,55],[26,44],[6,83],[23,40],[3,142],[16,99],[65,58],[116,31],[51,34],[41,10],[109,-37],[26,32],[26,-37],[6,-60],[-11,-58],[-44,-83],[-30,-89],[4,-45],[-31,-57],[32,-12],[-41,-249],[-12,-39],[-42,99],[-12,-55]],[[11284,23808],[-28,-6],[-35,-105],[-30,8],[-14,-41]],[[11177,23664],[-224,0],[-1,-123],[-51,0],[53,-85],[33,-35],[9,-43],[27,-27],[-4,-69],[-158,-1],[-56,-164],[14,-55],[-20,-116]],[[10799,22946],[-70,131],[-137,200],[-56,51],[-35,-18],[-30,47],[-29,-53],[-40,-44],[-92,-62],[-37,-9],[-36,17],[-47,40],[-70,12],[-47,53],[-47,22],[-30,50],[-114,41],[-41,44],[-102,61],[-92,99],[-30,60],[-46,7],[-59,23],[-92,58],[-58,110],[-60,58],[-66,48],[-21,55],[-46,91],[-23,90],[51,43],[-29,43],[31,75],[4,82],[-28,28],[-26,81],[0,74],[-18,66],[-75,125],[-66,150],[-72,106],[1,28],[-54,27],[-52,128],[-46,50],[-34,12],[-44,55],[-5,67],[28,60],[-10,50],[-25,38],[-34,-1],[-23,82],[-41,19],[-24,35],[-17,73],[10,46],[-48,5],[-25,17],[-43,92],[-25,19],[-33,77],[-27,43],[-7,55],[-64,157],[-10,69],[-35,109],[7,84],[-71,37],[-1,27],[-38,35],[-25,-27],[-49,50],[-36,14],[-5,-141],[14,-43],[16,-99],[-2,-59],[34,-90],[24,-21],[52,-80],[32,-97],[36,-28],[14,-63],[27,-19],[17,-132],[50,-66],[17,-74],[22,-48],[54,-57],[29,-128],[5,-74],[32,-57],[17,-84],[27,-78],[-7,-44],[23,-82],[64,-9],[39,-82],[4,-31],[32,-40],[-5,-58],[-57,-72],[-20,26],[-12,74],[-22,58],[-33,29],[-50,81],[-47,49],[-81,112],[-7,44],[9,98],[-14,93],[-25,66],[-35,23],[-44,59],[-23,60],[-48,-30],[-30,54],[-75,55],[-11,47],[-56,67],[52,10],[33,20],[14,30],[17,91],[-12,39],[-104,171],[-21,10],[-63,72],[-17,120],[-22,24],[-9,86],[-28,36],[-5,51],[-40,80],[5,63],[-28,32],[-35,118]],[[8342,26120],[-9,-52],[-29,18],[5,70],[23,16],[10,-52]],[[29242,15258],[-33,-5],[1,61],[23,51],[27,-48],[-18,-59]],[[20114,23230],[0,99],[23,140],[33,136],[6,75],[-7,139],[-15,106],[-38,79],[29,93],[9,97],[-27,93],[-24,-4],[-31,99],[-21,-61]],[[21080,25743],[162,-170],[139,-146],[174,-187]],[[21555,25240],[-218,0],[28,-432],[21,-328],[14,-218],[24,-382],[18,-273],[14,-217],[33,-62],[-19,-173],[-252,0],[-199,0],[-94,-27],[-58,14],[-35,-3],[-24,-60],[-75,107],[-32,-46],[-13,-95],[-29,-55],[-22,14]],[[21555,25240],[193,-223],[121,-141],[120,-141],[121,-141],[179,-208],[3,-63],[55,-57],[9,-39],[30,-24],[36,-5],[23,-40],[56,-24],[42,-37],[6,-83],[-18,-57],[36,-35],[102,33]],[[22175,23027],[-56,32],[-64,-2],[-36,-50],[-88,-74],[-25,-6],[-18,-63],[-42,23],[-51,-71],[-5,-60],[-35,-1],[-13,-86],[-33,-19],[-35,39],[-24,2],[-34,-58],[12,-71],[-31,-27],[7,-73],[-46,-67],[-38,-14],[-22,-31],[7,-67],[-9,-70],[-23,-36],[4,-59],[-8,-76]],[[21469,22042],[-40,-8],[-23,-43],[-19,46],[-9,71],[-37,-31],[-16,-52],[-32,-1],[-19,-31],[-32,13],[-16,39],[-20,-2],[-41,-59]],[[21165,21984],[-2,35],[-27,26],[-13,121],[-41,10],[32,77],[-51,63],[0,55],[-27,105],[-40,17],[1,-49],[-57,-47],[-56,37],[-23,-16],[-26,-52],[-32,67],[-46,-42],[-24,40],[14,45]],[[34712,21124],[29,-16],[24,-67],[55,-75],[27,-57],[27,-91],[6,-100],[-13,-137],[11,-54],[-2,-129],[12,-35],[34,-43],[50,-188],[9,-53],[-14,-26],[-23,20],[-36,-1],[-27,-26],[-15,47],[-77,68],[-22,41],[-50,45],[-45,72],[-32,24],[-27,44],[0,81],[-72,156],[10,13],[-22,77],[0,61],[-18,86],[-14,119],[-2,87],[-27,101]],[[36616,20669],[-59,37],[-72,7],[-12,-14],[-64,9],[-27,-34],[-13,-56],[0,-108],[-14,-90],[-26,-2],[-19,-41],[11,-69],[-48,-60],[5,-60],[-16,-28],[-19,-84],[-19,7],[-60,-14],[-35,-44],[-75,44],[-8,29],[-58,-2],[-35,-26],[-14,-65],[-38,-32],[-59,10],[-43,-6],[-53,-34],[-24,29],[-81,136],[-14,62],[11,29]],[[35638,20199],[11,-37],[33,-31],[45,0],[40,-37],[32,-7],[27,48],[14,88],[-3,65],[24,-11],[0,55],[35,48],[48,14],[77,34],[37,29],[50,118],[59,110],[17,77]],[[36184,20762],[19,-25],[48,-100],[20,32],[8,49],[-11,71],[34,40]],[[36302,20829],[10,-112],[20,54],[-16,58]],[[36316,20829],[47,33],[-13,80],[17,34],[30,-7],[32,76],[10,54],[44,86],[31,100],[7,-63],[13,-6],[27,67],[15,-11],[5,-53],[43,-44],[-3,-117],[25,0],[23,24],[13,-41],[42,-22],[17,-38],[45,-36],[33,-3],[3,-47],[-15,-22],[-56,-30],[-36,12],[-24,-41],[46,-72],[-8,-31],[-60,-24],[-33,19],[-20,-37]],[[36654,20673],[-29,-4]],[[36625,20669],[29,4]],[[26326,17634],[-5,-32],[18,-78],[25,-63],[-8,-26],[4,-146],[12,-93],[-33,-87],[2,-67],[32,-129],[-1,-96],[37,-81],[-11,-52],[24,-62],[28,26],[40,-47],[-23,118],[-28,28],[-1,39]],[[26450,17213],[-1,6]],[[28244,17025],[49,-140],[29,-213],[8,-152],[25,-91],[5,-52],[-33,-126],[-24,34],[-15,76],[-30,-24],[8,-112],[14,-39],[-8,-123],[-28,-48],[-12,-69],[5,-121],[-15,-96],[-37,-172],[-32,-183],[-24,-110],[-31,-196],[-55,-245],[-20,-169],[-23,-140],[-22,-76],[-24,-125],[-18,-42],[-38,-37],[-70,-18],[-80,-73],[-48,4],[-37,46],[-58,24],[-38,51],[-23,100],[-20,39],[-5,136],[9,46],[-18,99],[-22,42],[-16,111],[0,73],[21,89],[8,63],[37,39],[14,69],[40,107],[20,100],[6,108],[-26,78],[-1,73],[-23,100],[-8,197],[54,151],[7,106],[53,10],[32,42],[53,-2],[13,39],[77,22],[18,44],[50,62],[14,-3],[35,66],[66,87],[-5,39],[27,90],[-12,50],[18,30],[26,-27],[19,39],[48,60],[14,73],[-2,104],[38,83],[41,-77]],[[24898,29041],[61,-70],[20,-61],[-11,-83]],[[24968,28827],[-52,-48],[-69,-3],[-44,-52],[-75,-4]],[[24728,28720],[-27,12],[-31,81],[-5,54],[14,78]],[[24679,28945],[26,43],[35,22],[24,-16],[38,33]],[[22901,30755],[17,-55],[29,-16],[-18,-76]],[[22929,30608],[-68,19]],[[22861,30627],[11,23],[-17,61],[46,44]],[[24720,31889],[7,-2]],[[24761,31884],[-23,177]],[[24738,32061],[75,53],[53,20],[213,-17],[37,-14],[97,30],[20,-44],[64,-15],[85,-84],[39,-18]],[[25421,31972],[-17,-71],[-35,-75],[-39,-18],[-38,-129],[-96,-78],[-158,-8]],[[23321,30129],[6,-47]],[[23564,27037],[38,-20],[58,-51],[58,-13],[65,25],[45,-25],[72,-26],[34,-37],[82,-27],[22,-51],[17,-110],[26,-50],[51,-36],[81,-11],[70,-29],[104,-67],[59,-79],[33,-28],[43,0],[52,44],[37,68],[17,61],[-23,105],[-5,57],[24,88],[62,80],[54,43],[45,5],[26,30],[68,-4],[111,-66],[2,-63],[22,-26],[63,-12],[41,-33],[68,2],[42,-29],[15,-50]],[[25243,26702],[-36,-71],[14,-122],[-11,-70],[-19,-46],[31,-234],[0,-247],[0,-345],[0,-197],[0,-345],[0,-443]],[[25099,24033],[-118,104],[-177,156],[-235,207],[-177,156],[-177,156],[-100,88],[-123,-99]],[[23621,24916],[-53,169],[-105,58],[-30,-16],[-23,24],[-28,101],[-3,62],[-68,161],[5,59],[46,48],[4,65],[-18,106],[21,99],[-13,170],[4,90],[-21,131],[-45,121],[26,25]],[[20883,21616],[27,-12],[37,23],[30,-41],[19,-155],[-11,-44],[42,-49],[28,10],[22,87],[27,-22]],[[21104,21413],[24,-107],[-5,-59],[-33,-65],[49,-47],[39,-13],[11,-57],[42,-29],[3,-114],[-17,-53],[-3,-57],[6,-103]],[[21220,20709],[-88,52],[-107,102],[-141,225],[-63,51],[-7,34],[-55,49],[-26,48]],[[26562,27092],[-41,-40],[-15,-38],[-37,1]],[[26469,27015],[61,208],[5,48],[40,84]],[[24738,32061],[3,165],[35,43],[13,66],[33,55],[102,34],[71,-88],[19,-52],[44,-26],[78,44],[15,34],[-10,120]],[[25141,32456],[97,42],[109,-49],[66,-68],[45,17],[56,-17]],[[25612,32077],[-31,-15],[-39,-61],[-65,7],[-56,-36]],[[35378,22981],[-50,-63],[-36,-22],[-39,58],[-18,-28],[-46,-24],[18,-65],[-27,-27],[-46,51],[-22,-10],[-21,52]],[[34468,24213],[16,91],[39,21],[-4,36],[27,71],[47,56]],[[34593,24488],[11,-76],[39,2],[20,-18],[4,148],[-26,94],[18,46],[56,-18]],[[30878,29027],[10,63],[28,43],[69,20],[181,-70],[17,20],[16,89],[81,53],[103,-74],[72,-20],[19,26],[89,-3],[67,11],[22,-12],[97,-12],[62,1],[79,-23],[40,-66],[53,-10],[35,-49]],[[32018,29014],[1,-35],[-114,-54],[-114,-91],[-30,-65],[-66,-18],[-74,10],[-20,-12],[-42,-117],[-75,-32],[-26,6],[-9,60],[-71,-34],[-73,-69],[-52,-20],[-19,-54],[9,-42],[-34,-35]],[[27874,26140],[46,97],[26,93],[26,26],[42,3],[38,-24]],[[28052,26335],[20,-80],[-21,-1],[-30,-44],[40,-13],[16,-83],[32,-95]],[[24679,28945],[-10,76],[-30,26],[-22,45]],[[26488,20959],[59,-32],[7,-41],[-5,-77],[39,-79],[102,-8],[145,-167],[64,-11],[109,-32],[20,27],[23,60],[114,93],[39,-62],[17,-11],[81,8]],[[27259,19382],[-18,-38],[-48,-18],[-12,-70],[-31,-60],[-30,-3],[-22,-29],[-13,-124],[-32,-71],[-5,-46],[-40,-152],[-33,-47]],[[26339,19525],[13,36],[-10,66],[17,30],[73,36],[-12,41],[-36,-23],[-17,-34],[-31,63],[-7,61]],[[29538,29597],[-2,148],[-37,31],[24,57],[60,-16],[44,56],[55,37],[5,38],[-111,42],[-44,-17],[21,-37],[46,-15],[12,-55],[-67,13],[-29,-12],[-30,-62],[-45,-14],[0,60],[-66,-34],[-16,-68]],[[28608,28924],[-4,70],[19,72],[-3,73],[-92,33],[-28,55],[-40,3],[1,68],[-20,36],[-38,120],[-61,29],[9,65],[56,1],[39,-28],[-17,105],[46,83],[39,9],[98,0],[68,-20],[-38,59],[45,135],[-7,78],[11,27],[-31,63],[-54,8],[-36,-33],[-66,39],[-58,20],[-73,-50],[-21,-1],[-51,-53],[-80,-26],[-15,-40]],[[32893,30528],[-70,-17],[-4,-69],[-21,-37],[-98,-31],[-28,-103],[16,-145],[-21,-42],[-31,-6],[-55,-45],[-15,31],[-80,-1],[-101,48],[-23,-34],[-55,-192],[-28,-140],[31,-17],[-13,-70],[-26,21],[-34,-13],[-48,33],[-92,-40],[-97,-27],[-8,-46],[56,-11],[-14,-67],[2,-75],[17,-45],[36,-161],[-49,-26],[19,-29],[-46,-54],[5,-104]],[[29377,29734],[51,69],[27,-47],[7,-94]],[[26919,27079],[35,-193],[11,-81]],[[26449,26197],[3,44]],[[26452,26241],[21,129],[4,83],[32,134],[2,76]],[[26511,26663],[13,63],[-1,138]],[[26523,26864],[3,54],[26,21]],[[39846,29437],[34,-3],[47,-34],[39,-2],[68,64],[-31,-102],[30,-129],[-39,-32],[-49,-19],[-41,8],[-40,-22],[-67,-102],[-24,-92],[-89,57],[-81,70],[-55,-7],[-52,-45],[-34,47],[-28,1],[-20,-49],[26,-45],[25,-3],[52,-68],[-19,-15],[-42,17],[-34,-65],[-29,-21],[-18,33],[13,74],[-35,104],[9,58],[41,32],[32,56],[0,63],[36,-27],[64,-3],[14,40],[-2,57],[30,83],[17,153],[-24,96],[10,54],[34,24],[58,-85],[36,-67],[72,-93],[66,-61]],[[38289,27129],[50,11],[16,-33],[-21,-50],[44,-4],[10,-90],[-30,-55],[-34,-156],[0,-47],[-15,-58],[-29,-33],[-38,3],[-25,-20],[-48,25],[15,68],[-16,37],[1,70],[32,47],[22,69],[-31,104],[-32,3],[6,-56],[-43,-27],[9,65],[-39,47],[10,28],[87,59],[14,44],[29,21],[29,-13],[27,-59]],[[38681,27273],[34,-7],[13,-89],[-45,-47],[-24,-79],[-27,44],[-40,14],[-43,-33],[-38,-114],[-59,17],[2,82],[-18,48],[34,39],[17,67],[26,21],[17,-36],[55,20],[16,49],[45,26],[35,-22]],[[39526,28835],[28,7],[-7,-68],[8,-107],[41,-70],[22,-98],[0,-91],[-9,-70],[-30,-30],[-23,-125],[-45,-14],[-22,-86],[14,-106],[-9,-102],[-21,-34],[-21,-75],[0,-98],[31,-73],[-52,-47],[-5,-53],[-44,-51],[-32,-18],[4,75],[32,54],[-31,25],[-23,-54],[-49,-29],[-17,-40],[-16,-87],[-27,0],[5,60],[-28,25],[-48,-108],[-79,15],[-33,24],[-48,7],[-24,39],[-19,-67],[42,-53],[-3,-24],[-64,-33],[-51,-135],[-28,-16],[-29,14],[-35,76],[-9,86],[35,46],[-3,34],[-39,-5],[-37,29],[-66,-14],[-29,-38],[-64,-21],[-38,-28],[-59,-13],[-29,24],[-23,-28],[-18,-81],[-43,43],[-63,-23],[-38,6],[-4,63],[14,29],[43,4],[62,69],[131,172],[43,10],[13,-21],[103,17],[30,20],[88,26],[19,-49],[44,-5],[51,58],[-11,49],[85,165],[2,100],[16,41],[21,-120],[40,-15],[21,40],[99,59],[48,75],[22,62],[63,66],[22,100],[28,62],[30,128],[1,62],[-18,60],[13,64],[-11,63],[44,54],[13,84],[29,-7],[7,-67],[54,-2],[17,48],[-7,31],[-49,-23],[16,81],[36,-29]],[[12642,23805],[43,-15],[26,-29],[43,-23],[17,-52],[-39,-11],[-40,24],[-27,-16],[-16,-41],[-20,26],[-50,10],[-37,69],[-27,5],[9,51],[42,16],[76,-14]],[[23435,30040],[70,-15],[27,39],[112,-8],[29,-51],[161,-41]],[[23837,29760],[-63,40],[-22,-29],[-66,-39],[-27,-3],[-6,-45],[36,-60],[-33,-54],[18,-110],[36,-50],[108,-93],[29,-86],[25,-107],[22,-41],[84,-99],[37,-26],[98,1],[23,-40],[-29,-30],[12,-41],[134,-82],[46,-48],[59,-41],[46,-62],[19,-59],[-17,-62],[-33,26],[-26,75],[-48,7],[-32,38],[-36,-6],[-49,-131],[9,-49],[63,-57],[8,-84],[-52,-23],[-24,-39],[-2,-67],[-32,-35],[-28,-67],[-41,-1],[-10,52],[22,28],[19,90],[27,10],[-21,130],[-41,140],[-49,18],[-42,37],[5,30],[-28,64],[-37,-8],[-19,40],[-25,3],[-46,90],[-60,11],[-19,-10],[-56,50],[-69,103],[-33,31],[-21,45],[-57,55],[-57,88],[-48,126],[-9,75],[-24,37],[-39,18],[-54,47],[-65,23],[-26,-17],[-58,-93],[-72,-34]],[[23071,29360],[19,59],[-84,54],[-15,50],[23,46],[-36,45],[-14,43],[55,26],[10,36],[-43,91],[26,24]],[[24065,28143],[-42,-96],[-17,-71],[24,-98],[-22,-72],[-75,25],[-44,67],[-29,-1],[-121,103],[-28,0],[-32,55],[13,51],[23,29],[21,-33],[32,34],[23,-2],[41,-39],[101,9],[42,28],[34,-4],[56,15]],[[23334,28727],[21,-84],[-20,-51],[7,-38],[-17,-204],[-62,16],[-22,-71],[-28,3],[-29,61],[4,113],[-5,43],[8,83],[-35,78],[2,48],[33,-8],[43,26],[57,55],[41,-40],[2,-30]],[[26523,26864],[-44,31],[-16,-16],[-13,-66],[4,-55],[26,-30],[-31,-38],[-10,-45],[29,-7],[43,25]],[[26452,26241],[-9,-17]],[[26443,26224],[-30,156],[-51,224]],[[26362,26604],[29,82]],[[26391,26686],[25,69],[53,260]],[[26362,26604],[-6,25]],[[26356,26629],[35,57]],[[21383,31626],[-15,-23],[25,-90],[14,-142],[-24,-85],[-42,-77],[-27,4],[-92,-24],[-17,-26],[-141,-87],[-59,-19],[-67,-3],[32,45],[-99,39],[52,59],[7,59],[39,34],[35,118],[-61,65],[-10,43],[3,91],[-17,60],[83,5],[90,-15],[-5,25],[49,36],[-66,38],[48,45],[12,57],[58,12],[69,33],[3,-57]],[[27657,27906],[32,-135],[26,-37],[15,-75],[45,-39],[54,-4],[-24,-65],[20,-87],[-56,-72],[-5,-50],[-25,-34],[12,-58],[-17,-40],[42,-81],[18,2],[33,-84],[-4,-59],[32,-7],[92,-101],[30,-10],[18,-59],[39,-79],[-19,-86],[0,-88],[41,-2],[1,-115],[39,-40],[26,-71]],[[28122,26330],[-50,17],[-20,-12]],[[29062,25663],[-28,-49],[-24,25],[38,37],[14,-13]],[[29687,27572],[-23,-168],[-16,-44],[-27,-30],[-30,-93],[3,-100],[50,-29],[-44,-81],[33,-195],[-4,-129],[7,-39],[100,-22],[19,-68],[-4,-53],[-116,-214]],[[29727,25285],[-22,-22],[-92,40],[-79,22],[-50,5],[-20,17],[-50,-14],[-31,30],[-123,21],[-57,31],[-28,127],[-15,117],[-21,41],[-56,24],[-87,-49],[-28,-45],[-16,3],[-49,-50],[-31,-11],[-49,41],[-66,7],[-31,48],[-58,43],[-36,40],[-26,64],[-55,46],[-45,4],[-48,63],[-23,84],[-3,47],[-24,31],[1,43],[-25,18],[-3,60],[-59,110],[-23,63],[-53,-39],[-78,81],[0,-61],[-46,-35]],[[27663,28457],[40,-87],[42,-55],[78,-28]],[[27823,28287],[46,6]],[[27869,28293],[45,53],[113,110],[27,8],[40,-63],[-22,-19],[13,-73],[-27,-35],[48,-51],[22,-38],[34,5]],[[28162,28190],[11,-120],[15,-49],[48,-37],[81,-20],[49,-86],[72,-60],[80,-28],[52,2],[146,54],[79,18],[-13,86]],[[34143,20075],[56,-107],[-7,-76],[-24,-7],[-10,52],[-24,25],[-20,104],[29,9]],[[34350,19364],[-35,25],[-34,102],[9,50],[24,12],[49,-156],[-13,-33]],[[36501,17863],[-33,-69],[-37,40],[-1,68],[41,51],[41,-40],[-11,-50]],[[36354,17964],[32,-55],[-46,-46],[-12,-40],[-48,81],[-29,11],[-17,57],[65,-12],[19,23],[36,-19]],[[35197,19387],[14,-43],[6,-71],[20,-60],[55,-24],[-17,-29],[-16,-77],[-68,51],[-7,73],[-19,68],[-29,22],[-25,-10],[-26,18],[30,51],[-1,34],[26,29],[40,4],[17,-36]],[[37305,18755],[-15,-43],[-10,-86],[28,-30],[-26,-26],[-22,-61],[-28,28],[23,73],[10,130],[23,49],[17,-34]],[[37240,18597],[-10,-26],[-35,15],[14,55],[-3,66],[41,33],[7,-69],[-18,-45],[4,-29]],[[35476,20564],[-28,35],[-10,42],[24,35],[24,-47],[-10,-65]],[[35463,19096],[-19,-50],[-14,21],[-40,-18],[-2,90],[9,52],[50,-7],[26,-51],[-10,-37]],[[38807,19611],[63,-17],[48,-74],[-21,-33],[-34,21],[-56,103]],[[38818,19404],[90,-13],[-1,-45],[-87,42],[-2,16]],[[38245,19753],[57,-32],[-22,-40],[-43,-25],[-8,33],[-55,19],[7,24],[64,21]],[[37954,20204],[-19,-4],[-10,58],[14,38],[34,28],[10,-27],[-8,-55],[-21,-38]],[[37917,19389],[-72,-15],[-21,19],[8,42],[35,20],[50,-66]],[[37655,19361],[-8,-28],[-59,-5],[5,31],[62,2]],[[37526,19379],[43,-22],[-39,-30],[-21,11],[-51,-24],[-5,69],[35,12],[38,-16]],[[38308,17998],[-30,29],[22,87],[37,74],[16,-67],[-8,-41],[-37,-82]],[[37751,18071],[-40,-62],[-78,29],[16,34],[30,-10],[48,30],[24,-21]],[[37477,17967],[68,-14],[-4,-33],[-80,-20],[-2,58],[18,9]],[[38268,19465],[-27,13],[-14,65],[28,15],[22,-17],[-9,-76]],[[34022,20902],[43,9],[43,-14],[43,-1],[45,-71],[14,-54],[28,-48],[7,-71],[43,-36],[15,-37],[46,-39],[71,-88],[49,-117],[40,-86],[29,-30],[17,27],[27,3],[31,-54],[22,-70],[38,-15],[46,-80],[10,-59],[28,-46],[47,-14],[27,-49],[50,-3],[24,-42],[14,-53],[-44,-52],[0,-75],[13,-49],[23,-29],[59,-37],[20,4],[25,-194],[34,-38],[-21,-64],[8,-40],[36,46],[45,-5],[23,-24],[57,-135],[-25,-111],[11,-59],[-12,-62],[7,-86],[0,-100],[-9,-146],[-24,-27],[-34,55],[-33,-43],[-54,49],[-5,-84],[-65,114],[-29,67],[-113,134],[-47,69],[-50,122],[-68,95],[-33,96],[-23,31],[-33,97],[0,46],[-45,140],[-22,103],[-55,113],[-32,91],[-54,55],[-45,251],[-28,89],[-69,74],[-38,26],[-13,108],[-25,28],[-52,131],[-20,30],[-44,23],[-117,208],[-36,115],[3,61],[21,14],[55,-25],[35,-48],[45,-14]],[[37257,17864],[-39,-25],[-46,-5],[-55,-33],[-29,18],[-10,-22],[-36,-4],[-60,29],[-53,5],[-26,-17],[-12,35],[13,55],[39,34],[46,11],[82,-52],[21,-22],[64,27],[37,-37],[21,5],[21,44],[30,21],[-8,-67]],[[36916,17696],[34,-60],[26,-6],[41,-79],[-17,-37],[-32,-20],[-36,21],[-26,51],[-41,43],[-63,14],[-16,41],[28,30],[53,7],[49,-5]],[[36698,17928],[30,10],[62,-14],[2,-91],[-65,-23],[-12,40],[-24,-37],[-138,-56],[-34,20],[6,104],[40,36],[50,-13],[29,-62],[20,-4],[33,30],[-52,55],[-7,42],[45,6],[15,-43]],[[37516,19972],[-24,-37],[-33,-78],[-26,-20],[-57,-16],[-60,4],[-23,35],[-98,-1],[-54,-9],[-51,12],[-51,-11],[-39,16],[-43,-15],[-27,-62],[-14,-79],[10,-100],[21,-55],[31,-30],[18,-73],[45,-8],[60,121],[56,-17],[38,39],[75,0],[34,41],[26,-18],[0,-78],[-41,29],[-29,-20],[-37,-84],[-43,-54],[-37,-22],[-17,-37],[-46,-18],[62,-85],[28,-92],[26,-34],[13,-65],[-10,-16],[-8,-76],[44,-66],[10,-36],[23,-5],[3,-53],[-85,-32],[-36,-77],[-40,19],[-13,39],[16,107],[-37,39],[-49,79],[-1,34],[18,52],[0,91],[-8,18],[-41,0],[-39,-44],[-11,-40],[17,-65],[6,-79],[-9,-83],[7,-117],[-17,-116],[13,-54],[-9,-33],[-45,-7],[-29,-26],[-42,59],[-2,24],[20,96],[11,100],[2,85],[-25,123],[-53,-14],[-15,31],[-7,53],[6,50],[-9,36],[37,62],[10,75],[19,46],[-2,115],[25,110],[25,49],[7,44],[-6,86],[18,29],[-7,43],[30,100],[26,61],[31,-34],[29,48],[19,56],[21,8],[44,-25],[19,-32],[34,5],[70,-16],[56,-37],[48,17],[69,-19],[53,40],[32,46],[8,35],[23,17],[25,51],[30,-43],[-42,-112]],[[35361,18435],[36,-46],[58,-17],[19,7],[27,-46],[20,-70],[76,-16],[25,12],[40,-20],[74,-10],[20,31],[18,73],[29,8],[23,-51],[47,5],[68,-54],[55,-7],[14,-65],[18,-18],[0,-55],[55,-36],[62,3],[40,15],[38,-30],[8,-27],[-8,-112],[25,-80],[-110,64],[-54,42],[-71,-27],[-69,19],[-75,4],[-111,34],[-70,57],[-93,41],[-66,8],[-36,-29],[-66,16],[-45,40],[-32,16],[-81,12],[-26,39],[12,41],[-70,42],[-57,17],[21,69],[16,3],[10,74],[26,45],[92,-41],[27,43],[41,-23]],[[36616,20669],[-14,-21],[33,-61],[-12,-35],[-57,-10],[26,-50],[-4,-37],[32,-29],[4,-52],[18,-17],[34,-94],[-25,-77],[27,-59],[69,-84],[42,-74],[-11,-21],[-44,-16],[-42,14],[-22,36],[-11,-46],[-22,-22],[-28,-108],[-7,-123],[12,-98],[-39,-34],[-30,-58],[-28,-17],[-27,-46],[-16,-68],[0,-59],[18,-54],[-5,-46],[-19,-15],[-5,-71],[-46,-152],[-156,-126],[-8,12],[-3,90],[-8,46],[-56,47],[-35,-39],[-20,8],[-13,52],[-23,-14],[-38,69],[-8,-56],[-45,-47],[-39,18],[-41,-46],[-16,-1],[0,105],[-56,27],[-40,-27],[-42,8],[-16,28],[-42,-7],[-1,52],[-18,163],[-14,18],[7,102],[-28,85],[-42,31],[-9,51],[-27,31],[-4,51],[16,65],[-38,71],[-5,96],[21,154],[37,94],[31,23]],[[37866,19940],[18,-4],[11,47],[23,25],[0,35],[33,44],[32,12],[2,-102],[-50,-51],[-4,-31],[43,-41],[10,-41],[-88,24],[-11,-38],[0,-54],[11,-66],[34,-106],[-26,6],[-19,62],[-24,40],[2,116],[-19,44],[-4,96],[-13,72],[15,47],[11,82],[33,65],[-1,-67],[15,-29],[0,-81],[-44,-70],[10,-36]],[[38114,19125],[29,-25],[48,-2],[24,-31],[12,-58],[23,-39],[-6,-63],[-97,84],[-22,33],[-40,-1],[-6,-27],[-61,28],[-13,20],[-28,-44],[-27,4],[-17,39],[-24,12],[14,70],[106,6],[7,-22],[38,32],[40,-16]],[[37758,19076],[45,-67],[1,-53],[-67,-41],[-58,47],[-20,41],[-3,55],[65,23],[37,-5]],[[37538,17666],[-13,-34],[-66,-106],[-31,-7],[-39,-36],[-26,8],[8,51],[-15,24],[15,77],[40,60]],[[37461,17737],[59,54]],[[38729,18501],[1,-102],[-15,-42],[-24,6],[-35,60],[18,12],[5,65],[34,60],[16,-59]],[[38703,18340],[-22,-82],[-20,-21],[-17,31],[14,130],[45,-58]],[[39195,17938],[-30,-29],[-52,6],[-23,25],[40,137],[35,44],[59,11],[27,-67],[-25,-80],[-31,-47]],[[39495,17752],[-38,60],[-22,50],[-47,70],[-21,44],[-67,-20],[-24,23],[-32,-39],[-10,26],[29,122],[-42,74],[-18,69],[32,17],[-52,113],[-23,144],[-20,-4],[-3,53],[-38,46],[-82,73],[-58,23],[-79,63],[-65,20],[-31,-2],[-54,56],[-10,25],[-65,62],[-21,-4],[-31,52],[-9,47],[-50,-153],[-34,-7],[-27,86],[15,33],[-15,57],[-49,70],[-37,13],[-10,28],[32,23],[61,-23],[38,66],[19,11],[52,-24],[36,34],[2,63],[-53,-27],[-65,-10],[-56,12],[-24,-5],[-35,55],[-12,94],[-78,37],[-22,-14],[17,132],[68,33],[40,55],[32,22],[29,-1],[55,-34],[49,-47],[62,-4],[35,-135],[-19,-79],[6,-104],[37,-140],[3,55],[18,10],[9,-87],[41,-88],[56,-2],[46,76],[19,58],[28,33],[18,68],[56,16],[40,38],[-6,41],[84,78],[104,-67],[140,-123],[45,0],[58,-21],[15,-35],[28,-1]],[[34663,20210],[8,-30],[-21,-50],[-16,10],[-7,58],[36,12]],[[34826,19918],[-60,18],[1,66],[56,-60],[3,-24]],[[35017,20021],[9,-56],[-10,-29],[-40,41],[14,36],[27,8]],[[35004,19680],[1,-61],[-28,33],[27,28]],[[35648,19495],[-35,-14],[6,57],[33,-12],[-4,-31]],[[36157,18194],[-46,-25],[-53,2],[-34,16],[13,52],[136,6],[-16,-51]],[[37845,19684],[14,-33],[-26,-38],[-21,31],[33,40]],[[38188,19383],[8,-62],[-21,-17],[-44,14],[13,50],[44,15]],[[34018,20272],[-54,52],[-27,13],[0,57],[75,-88],[6,-34]],[[37309,19497],[42,-37],[-23,-37],[-26,44],[-32,-62],[-10,34],[12,55],[37,3]],[[36459,18905],[-26,-41],[-9,78],[12,79],[18,-12],[5,-104]],[[37277,17359],[14,47],[38,49],[6,-39],[-58,-57]],[[37442,17926],[-28,-49],[-11,43],[39,6]],[[37397,17938],[-63,-69],[-11,33],[26,38],[30,16],[18,-18]],[[36654,20673],[-29,-4]],[[31629,27460],[93,85]],[[31722,27545],[30,-10],[-4,-44],[33,-131],[48,-30],[32,-36],[-22,-58],[6,-129],[38,-60],[11,-62],[2,-98],[-37,-31],[-27,52],[-38,-16],[13,-67],[30,-57],[-5,-48],[12,-98],[38,24],[29,-65],[28,-34],[28,4],[49,-46],[0,-42],[60,-33],[41,-55]],[[32990,25871],[63,49],[28,-49],[-12,-76],[17,-46]],[[33086,25749],[-18,-30],[14,-47],[36,-32],[23,7],[53,-33],[44,12],[27,30],[49,-26],[67,4],[17,17],[30,-14],[41,11],[10,96],[-11,35],[-31,-2],[-18,25],[5,45]],[[33424,25847],[42,-7],[54,21],[31,27],[4,44],[52,56],[10,41],[57,22],[29,25],[48,75],[41,37],[18,-30],[80,-27],[11,33],[64,53],[25,-44],[-9,-74],[34,31],[18,-63],[-38,-74],[40,7],[21,-20],[37,0],[31,-33]],[[34124,25947],[1,-61],[-53,-65],[-3,-102],[-26,35],[-74,-26],[-20,-38],[-54,-64],[-41,-34],[-9,-27],[9,-95],[-17,-60],[-45,-79],[-9,-37],[19,-42],[-15,-62],[-36,-98],[-20,-98],[-55,29],[-46,3],[13,-74],[-6,-121],[-35,-91],[10,-78],[-24,-79],[-34,28],[-14,-33]],[[33540,24578],[-11,155],[-16,53],[-3,87],[-11,79],[-39,1],[1,-40],[-22,-48],[2,-40],[-19,-27],[-37,27],[-19,122],[25,95],[63,22],[23,40],[23,114],[-25,60],[-70,-6],[-23,8],[-105,-4],[-77,32],[1,140],[-66,21],[-22,36],[-82,37],[-36,-62],[-8,-44],[45,-69],[39,-18],[15,-65],[-53,-2],[-18,-67],[-20,7],[-15,-63],[15,-31],[61,-35],[11,-21],[-21,-122],[22,-52],[-2,-40],[15,-47],[25,-208]],[[33106,24603],[0,-96],[-98,-7],[-38,44],[-33,-37],[-71,-34],[-30,-58],[-2,-29],[17,-89],[-28,-85],[-31,-31],[-27,-55],[-87,-50],[-14,45],[-26,-31],[-3,-54],[-56,-84],[-38,-96],[-44,-87],[-55,-49],[-56,-101],[-75,-74],[-28,-39],[0,-69],[-13,-48],[-61,-51],[-44,8],[-20,-22],[-13,-67],[-19,-44],[-41,30],[-44,-41],[-23,-85],[-6,-55],[14,-109],[-7,-80],[16,-96],[-22,-37],[34,-54],[-14,-147],[-10,-52],[-35,-102],[-13,-91],[12,-83],[-2,-192],[-55,-3],[-16,-60],[-32,-77],[-10,-50],[8,-41],[-69,-36],[-28,-47],[-16,-111],[-67,-67],[-27,15],[-41,57],[-51,109],[-28,120],[17,19],[-14,85],[-52,189],[-25,127],[-24,74],[-41,78],[-31,112],[-21,112],[-12,130],[-22,83],[-14,98],[-53,128],[-2,70],[-20,39],[-36,106],[-18,89],[-22,255],[-20,104],[-15,129],[-7,135],[-18,116],[26,161],[-8,122],[-15,13],[-18,114],[-3,61],[9,67],[-32,-2],[-11,-53],[-26,-45],[27,-36],[0,-28],[-29,-83],[-76,-63],[-46,-28],[-38,0],[-29,22],[-44,56],[-46,90],[-69,107],[-22,45],[66,46],[78,36],[18,54],[-9,34],[-72,-47],[-53,20],[-73,75],[-28,83],[-22,6],[-9,57]],[[33558,22286],[-23,74],[19,75],[6,124],[8,22],[6,105],[27,-24],[-23,-68],[13,-88],[-16,-23],[2,-45],[-16,-52],[6,-26],[-9,-74]],[[20236,34290],[99,-39],[-2,-65],[140,-58],[-23,-30],[29,-62],[-120,-137],[-59,-33],[-128,-40],[-71,-54],[-145,-34],[-16,-39],[-87,-28],[-190,33],[-129,83],[-25,-17],[-142,-5],[12,34],[83,47],[-14,68],[-42,16],[-22,46],[-192,15],[165,35],[84,103],[-111,19],[-122,-38],[-40,18],[85,148],[112,-36],[-43,71],[67,23],[99,-74],[39,-47],[-20,-53],[29,-47],[61,52],[47,13],[0,69],[40,4],[54,-59],[25,65],[75,20],[107,-5],[40,-38],[53,51],[50,-16],[39,32],[-14,37],[69,15],[54,-63]],[[24475,29835],[-58,-36],[-58,1],[-80,49],[-70,93],[-28,18]],[[24129,30040],[44,31],[-2,79],[23,18],[15,52],[39,5],[10,66]],[[11402,22603],[-19,82],[-40,10]],[[11343,22695],[12,91],[-53,38],[-37,-30],[-8,28],[-41,32],[-26,45],[-37,19]],[[11153,22918],[23,42],[-6,57],[10,45],[20,15],[92,129]],[[11292,23206],[44,33],[31,7],[16,-26],[25,9],[47,-16],[68,5],[70,48],[37,-22],[63,19],[87,-36],[62,-133],[26,17],[22,-12],[26,-54]],[[13316,24082],[15,-129],[-18,-32],[8,-56],[-32,-29],[29,-56],[0,-67]],[[13318,23713],[-36,42],[-55,-2],[-46,-15],[-62,22],[-45,-14],[-17,-32],[-25,38],[-41,29],[4,61],[20,8],[78,-31],[99,-19],[54,53],[-57,87],[14,81],[-43,37],[-43,11],[0,32],[35,21],[59,0],[51,-34],[54,-6]],[[14764,21630],[22,-38],[21,-3],[57,-63],[71,-116],[18,-48],[0,-63],[22,-35],[39,-21],[97,-151],[0,-120]],[[15198,20180],[-43,-13],[-35,29],[-36,-23],[-22,-45],[-48,-12],[-7,-29],[-37,15],[-21,-32],[0,-34],[-43,-18],[-46,38],[-64,115],[0,82],[-17,19],[-13,72],[3,66],[17,79],[-2,52],[21,25],[16,51],[-35,119],[-34,8],[15,125],[-19,35],[-74,-8]],[[20459,22536],[3,-80],[-30,-29],[26,-37],[1,-60],[-66,-17],[-51,-33],[-45,-125]],[[20297,22155],[-21,47],[-22,14],[-10,42],[-3,69],[8,48],[-47,-39],[-17,43],[-22,-9],[-37,62],[-34,34]],[[20513,21740],[3,27],[-37,81],[-15,-1],[-1,69],[-90,87],[-43,150],[-33,2]],[[21165,21984],[-18,-30],[0,-116],[30,-18],[-7,-138],[32,-49],[-16,-21],[-51,-2],[-2,-45],[29,-23],[-26,-130],[-32,1]],[[11177,23664],[-9,-422],[42,0]],[[11210,23242],[36,-28],[8,19],[38,-27]],[[11153,22918],[-63,-81],[-27,-68]],[[11063,22769],[-63,42],[-67,-1],[-28,14],[-55,53],[-51,69]],[[24619,28409],[-24,6],[-29,58],[36,5],[17,-69]],[[25030,28305],[7,-24],[50,-36],[27,-3],[22,-99],[35,-16],[-8,-39],[-35,32],[-21,60],[-35,2],[-62,88],[20,35]],[[25398,28386],[21,-48],[-14,-30],[-69,36],[1,33],[31,19],[30,-10]],[[25279,28530],[-10,-39],[-37,10],[0,32],[47,-3]],[[25574,27640],[-15,53],[24,38],[29,-30],[-38,-61]],[[25083,27553],[39,14],[18,-51],[50,13],[93,-26],[32,9],[1,-45],[68,28],[-6,-49],[-178,-24],[-11,33],[-30,16],[-111,30],[1,52],[34,0]],[[25353,28693],[-23,26],[-92,33],[-39,-30],[-50,12],[-62,-81],[-27,-71],[-43,-1],[-84,61],[-5,-101],[41,-89],[9,-65],[-19,-16],[15,-50],[-23,-51],[59,-30],[53,-68],[35,-17],[10,-124],[-58,66],[-48,-12],[1,-74],[36,-39],[-41,-24],[-48,13],[38,-136],[-58,2],[-14,-47],[-51,101],[-23,-64],[-38,76],[14,50],[-15,51],[-30,28],[-23,55],[32,61],[32,-5],[20,33],[90,-47],[30,-29],[40,20],[-87,81],[-24,-19],[-32,13],[-60,-20],[-45,14],[-14,59],[-27,34],[0,44],[-59,70],[-36,84]],[[24610,28470],[25,-13],[21,42],[-8,30],[42,30],[19,72],[27,47],[-8,42]],[[24968,28827],[89,11],[52,31],[56,5],[10,-24],[81,-44],[111,42],[-11,53],[31,9]],[[22295,21090],[-30,-61],[-24,-11],[-61,-1],[-74,-56],[-56,-60],[-32,-10],[-71,-44],[-45,-48],[-49,37],[-88,35]],[[21765,20871],[4,8]],[[21769,20879],[8,1]],[[21777,20880],[28,12],[0,91],[-21,9],[-29,155],[-5,101],[31,87],[3,55],[21,105],[35,61],[-11,129],[-18,68],[6,82]],[[21817,21835],[-9,39],[-2,163],[-16,42],[10,89],[262,-2],[50,39],[28,-11]],[[23086,30200],[-11,18],[7,85],[31,126],[36,74],[-84,39],[-88,2],[-48,64]],[[22901,30755],[31,43],[-16,61],[-30,34]],[[23034,31449],[-18,20],[29,67],[89,2],[65,-29],[10,70],[34,-9],[1,93],[-31,30],[38,31],[-36,80]],[[23215,31804],[83,-21],[49,4]],[[23347,31787],[76,-85],[80,-56],[48,-52],[49,44],[39,5],[57,66],[56,-13],[51,-59],[34,3],[18,-66],[48,-27]],[[23971,30917],[-54,9],[-104,-46],[-7,-16],[-125,-59],[-22,-37],[22,-50],[-8,-47],[30,-61],[113,-110],[32,-42]],[[23848,30458],[-3,-40],[-37,-1],[-10,-41],[-73,-51],[17,-60],[8,-90],[-41,42],[-119,-19],[-51,-35],[-39,-6],[-14,27],[-54,5],[-15,-50],[-34,35],[-42,15]],[[23341,30189],[-63,28]],[[23835,31690],[-42,-30],[-23,66],[30,34],[35,-70]],[[23897,31578],[0,17]],[[27256,28866],[24,42],[7,58],[-33,151],[-45,37],[-36,52],[-38,12],[-67,66]],[[27862,28948],[-31,-51],[54,-69],[-20,-60],[-35,28],[-47,6],[-9,25],[-54,25],[-34,-35]],[[27686,28817],[-20,-18],[-75,1],[-97,-24]],[[20086,22621],[-8,61],[33,54]],[[23784,20228],[-16,-112],[14,-90],[71,39],[40,-10],[26,-81],[6,-42],[-14,-41],[-29,-19],[-25,-76],[-3,-87],[26,-14],[50,-77],[-8,-78],[5,-96],[-9,-105],[-16,-25],[-20,-102],[-25,1],[-12,66],[-40,-51],[-58,18],[-24,84],[-25,23],[-20,-23],[2,-88],[-47,-18],[-14,15],[-46,-11],[-5,-97],[28,-32],[-6,-43],[27,-31],[-31,-79],[-22,35],[-27,-26],[-19,-61]],[[23518,18894],[-22,56],[-38,58],[-36,84],[-73,109],[-60,152],[8,64],[-33,52],[-27,84],[-3,45],[39,30],[19,49],[15,121],[40,-26],[4,24],[-37,37],[-17,55],[20,12],[12,80]],[[23329,19980],[21,10],[32,-16],[161,-1],[-1,257]],[[23542,20230],[3,29],[125,-1],[62,-11],[43,2],[9,-21]],[[23315,29149],[9,-148],[-19,-88],[-26,-76],[-47,45],[-13,111],[-16,27],[18,73],[94,56]],[[23071,29360],[-7,-3]],[[23056,29352],[-24,-16],[-84,-107],[-47,-21],[-38,5],[-49,29],[-43,47],[-53,-12],[-71,45],[-98,-88],[-25,-61],[5,-71],[15,-35]],[[22358,29083],[-34,20]],[[21928,29281],[38,35],[17,100],[12,118],[7,132],[15,142],[-10,111],[-78,44],[-34,65],[5,50],[-92,105],[-48,39],[-91,32],[-51,-4],[-39,58],[29,62],[-40,-1],[0,39],[81,37],[42,0],[60,29],[28,-11],[38,-56],[31,25],[124,-2],[-16,36],[-2,87],[-34,106],[74,-1],[14,-64],[120,-20],[55,36],[-12,53],[53,35],[78,30],[37,65],[4,98],[41,55],[75,23]],[[22459,30969],[9,-48],[30,-37],[33,15],[20,-54],[59,-41],[52,-22],[0,-62],[59,9],[26,-47],[51,-24],[28,-37],[35,6]],[[29014,15070],[-53,14],[-16,48],[9,34],[44,-1],[21,-51],[-5,-44]],[[14664,22935],[-29,-6],[-19,74],[37,-15],[11,-53]],[[14570,23267],[-21,12],[-4,53],[30,-7],[-5,-58]],[[15485,20930],[20,70],[18,23],[48,-48],[68,-30],[55,-89],[20,-17],[46,-92],[16,10],[17,-112]],[[15793,20645],[-42,-79],[-45,-142],[-42,-151],[-33,-40],[-32,5],[-13,29],[-22,-19],[-31,23],[-45,-52],[-60,45]],[[14395,23720],[-14,0]],[[14381,23720],[14,0]],[[3775,15871],[-32,-10],[-4,45],[29,2],[7,-37]],[[42353,15310],[29,-8],[66,-87],[27,-20],[58,-109],[79,-81],[47,-70],[31,-30],[4,-51],[-24,-12],[-74,63],[-10,30],[-105,94],[-39,52],[-58,101],[-35,77],[4,51]],[[42747,15109],[-33,22],[19,37],[14,-59]],[[30662,8975],[43,29],[0,-61],[39,7],[57,36],[29,-31],[-32,-49],[-47,16],[-7,-43],[-30,-23],[-91,45],[-13,83],[52,-9]],[[25569,33041],[-68,0],[-87,-27],[-109,-18],[-61,-30],[-177,-50],[-114,24],[-27,62],[-140,52],[-9,81],[30,137],[-43,88],[12,51],[-31,88],[12,37],[56,54],[-13,36],[70,9],[15,43],[124,100],[127,136],[34,63],[101,46],[-4,103],[-94,61],[-52,10]],[[44297,16205],[-31,-61],[22,-65],[-62,-14],[-27,21],[-24,-42],[-37,-16],[-15,77],[27,-2],[35,43],[33,13],[24,28],[55,21],[0,-3]],[[44085,15941],[39,-62],[9,-94],[-41,-12],[-21,-25],[-39,-1],[-65,39],[-7,47],[17,51],[27,37],[70,33],[11,-13]],[[26643,22883],[35,13],[26,-10],[29,40],[39,-67],[38,155],[76,-94],[38,13],[32,32],[16,-33],[41,12],[44,-21],[41,-2],[37,-45],[30,-18],[73,-141],[50,-69],[23,-67],[22,-24],[30,-67]],[[27363,22490],[-75,-193],[4,-133],[19,-9],[74,31],[45,-18]],[[25141,32456],[26,106],[-51,-19],[-43,20],[-32,66],[-2,118],[73,21],[36,39],[131,11],[8,26],[138,-19],[43,-22],[113,-8],[14,15]],[[25553,32700],[-69,-6],[-17,-33],[83,-173]],[[24932,32621],[83,-31],[-66,-49],[-75,-32],[-33,33],[25,47],[66,32]],[[24969,32666],[-47,-30],[-16,38],[22,47],[39,-19],[2,-36]],[[26899,23706],[37,-127],[29,-153],[19,-178],[26,-85],[44,-14],[21,-29],[21,-71],[52,-29],[21,-30],[46,-27],[59,-140],[73,-87],[19,-82],[15,2],[33,-79],[36,-8],[4,-26]],[[27454,22543],[-29,-19],[-26,-54],[-36,20]],[[23223,20579],[26,-29],[-29,-89],[-29,9],[32,109]],[[23329,19980],[-26,24],[14,65],[17,28],[21,80],[-1,83]],[[23354,20260],[22,-30],[166,0]],[[11343,22695],[-14,-48],[-72,1],[-43,21],[-51,43],[-65,18],[-35,39]],[[25243,26702],[29,-31],[63,23],[166,-53],[36,-36],[53,-26],[67,-10],[69,-48],[64,39],[41,48],[36,7],[22,43],[37,-11],[37,13],[24,27],[40,-28],[45,18],[-14,-55],[41,-44],[27,36],[34,-41],[65,22],[66,-8],[65,42]],[[26443,26224],[-20,-46],[-15,-112],[-21,-88],[-6,-75],[-22,-55],[-56,62],[-63,114],[-22,111],[-43,99],[-19,99],[-26,-76],[26,-53],[8,-87],[33,-93],[37,-79],[42,-68],[1,-64],[42,-122],[8,-88],[46,-137],[29,-73],[77,-267],[53,-90],[-15,-63],[2,-74],[22,-109],[27,-45],[38,-24],[23,-52],[51,-65],[5,-22]],[[12265,19010],[36,35],[13,31],[24,112],[-49,11],[-19,-39],[-21,17],[-63,96],[24,25],[-8,97],[4,55],[-12,67],[43,51],[12,57],[-4,48],[43,79],[14,92],[-8,82],[43,43],[103,50],[5,54]],[[12445,20073],[67,-90],[104,-89],[6,-47],[62,-34],[52,-8],[23,47],[30,-22],[23,-43],[73,-57]],[[11033,19584],[-25,21],[8,35],[25,7],[10,-38],[-18,-25]],[[10918,19759],[12,-54],[24,-43],[2,-39],[20,-34],[-13,-42],[-58,-16],[-10,47],[38,39],[-27,74],[-8,58],[20,10]],[[13316,24082],[28,39],[39,-11],[35,15],[39,-31],[24,-1],[38,-29],[22,8],[26,-82],[51,6],[1,-28],[78,-65],[37,-51],[3,-29],[-37,-70],[-34,41],[-103,8],[-50,-42],[-37,-8],[-20,26],[-47,-13],[-10,-45],[-36,-93],[-29,28],[-16,58]],[[14608,23101],[-25,60],[20,14],[5,-74]],[[27454,22543],[22,-54],[14,-60],[-17,-48],[-58,-51],[-20,-39],[64,1],[11,-14]],[[18463,38095],[493,-60],[181,-91],[240,-19],[123,-42],[-117,-47],[-171,-21],[-703,-28],[4,-48],[210,24],[331,-11],[186,-57],[66,59],[217,13],[23,-82],[-71,-76],[120,-8],[91,55],[343,-29],[205,88],[228,-9],[186,-31],[82,-50],[-237,-84],[-127,-14],[-2,-47],[-278,-41],[102,-39],[-115,-44],[-148,-5],[-167,18],[-89,-54],[2,-45],[105,-27],[-16,-73],[41,-43],[-231,-131],[-97,-168],[108,27],[169,-43],[-146,-22],[55,-55],[115,-32],[118,-1],[-21,-96],[-191,33],[-108,-9],[-115,-69],[37,-62],[181,-16],[73,-101],[10,-116],[-130,18],[-62,-52],[123,-22],[94,-109],[-18,-43],[-216,-38],[-147,42],[-2,-52],[242,-47],[-17,-78],[-128,-14],[-78,-35],[-211,28],[129,-69],[100,-36],[-5,-114],[-27,-62],[-218,83],[-74,-13],[222,-102],[105,-62],[32,-46],[10,-137],[18,-72],[-52,-19],[-170,2],[-57,23],[-72,130],[-157,84],[-11,-76],[-120,-51],[-120,9],[-104,-113],[64,-17],[156,16],[73,-50],[77,24],[290,-53],[11,-51],[-77,-56],[-226,-127],[-19,-35],[-79,-41],[-358,-88],[-76,2],[-136,-54],[-55,14],[-88,57],[-28,-38],[11,-64],[-108,-65],[-117,-198],[-65,-63],[-123,-65],[-92,-68],[-68,-11],[-101,-41],[-30,25],[-79,-18],[-163,-15],[48,-48],[-133,-93],[58,-96],[-74,-56],[21,-20],[8,-89],[-61,-47],[-8,-52],[-65,-66],[-63,-89],[-21,-91],[26,-65],[-38,-71],[-37,-168],[-54,-96],[4,-59],[-24,-30],[-92,0],[-92,27],[-69,34],[-26,69],[-56,21],[-12,36],[16,70],[-106,-62],[-114,3],[-26,-23],[-46,62],[-14,48],[-52,12],[-45,68],[-12,66],[-30,24],[9,51],[-94,53],[-2,78],[-50,49],[-90,130],[3,57],[-38,83],[-41,29],[-21,163],[-62,95],[-53,1],[11,84],[-63,39],[-2,57],[73,91],[-80,35],[-24,32],[11,62],[46,35],[-19,57],[52,82],[192,-34],[48,60],[-122,-23],[-72,19],[-74,-3],[42,69],[91,20],[64,-34],[65,42],[22,104],[71,159],[-216,27],[-95,53],[-122,27],[-60,51],[42,36],[216,-29],[109,-55],[48,117],[-30,33],[-85,0],[-129,42],[45,42],[-90,30],[-59,-40],[-90,-20],[-111,41],[34,122],[-42,25],[7,51],[100,41],[6,51],[-114,29],[45,56],[-99,65],[22,74],[-46,59],[-50,12],[-12,99],[-221,157],[10,67],[-329,103],[-276,39],[-117,-7],[-132,-39],[-61,25],[-107,-60],[-163,22],[-151,58],[15,78],[-189,40],[-9,67],[232,4],[104,41],[207,3],[-30,63],[-158,-36],[-79,19],[-90,-30],[-283,95],[-144,65],[52,68],[419,77],[186,57],[186,3],[68,46],[77,145],[-252,17],[-19,74],[166,59],[139,69],[132,3],[55,43],[192,-16],[41,57],[-12,77],[94,38],[142,-8],[69,46],[492,57],[240,-6],[231,-94],[148,17],[-126,58],[-34,49],[262,-12],[440,-127],[93,5],[36,114],[-42,38],[-234,90],[405,66],[233,-33],[116,50],[198,0],[113,34],[596,25],[373,-8]],[[15660,35106],[84,-30],[-8,-70],[-180,-51],[-25,55],[-126,27],[17,73],[-20,43],[56,35],[123,-21],[79,-61]],[[15379,35731],[-63,-49],[-57,47],[120,2]],[[16628,37770],[-232,58],[-64,68],[107,8],[203,-64],[-14,-70]],[[19019,35320],[4,-59],[-101,-43],[-47,21],[-159,-21],[23,97],[90,-8],[143,40],[47,-27]],[[19934,36305],[51,-91],[-134,2],[-23,70],[106,19]],[[23695,31998],[-43,-70],[24,-39],[-40,-22],[-5,-82],[-38,22],[-11,60],[-45,4],[-38,113],[43,7],[18,42],[49,5],[43,33],[44,-12],[-1,-61]],[[23215,31804],[-6,113],[-60,40],[9,84],[-10,35],[5,102],[15,54],[41,57],[84,8],[31,19],[50,76],[37,8],[34,-37],[-2,-45],[-28,-54],[-1,-83],[69,-21],[1,-50],[-66,-18],[-17,-76],[-72,-82],[18,-147]],[[23458,31960],[22,-64],[-4,-41],[-43,-19],[-55,26],[-16,77],[96,21]],[[23546,31802],[47,-19],[-35,-39],[-52,32],[40,26]],[[24235,30421],[-59,44],[-52,-10],[-92,48],[-41,2],[-6,-46],[-29,-38],[-79,1],[-29,36]],[[26333,27450],[-59,-6],[-18,28],[-57,-17],[-25,18]],[[26174,27473],[28,49],[64,-12],[74,30],[-15,-40],[8,-50]],[[26333,27450],[-38,-21],[-23,-36],[-49,-24],[-21,-26],[-61,34],[-16,49],[49,47]],[[12079,24838],[69,-8],[12,-15],[65,10],[35,-35],[35,-1],[50,-43],[15,-37],[34,-37],[72,-9],[67,-56],[34,-46],[46,-7],[104,-128],[48,-13],[23,-21],[43,-5],[16,-25],[-6,-57],[93,-19],[46,-58],[42,-20],[-10,-47],[-47,-5],[-59,-34],[-54,-2],[-75,21],[-103,-21],[-88,-8],[20,50],[42,48],[13,41],[-15,34],[-95,16],[-61,74],[-17,101],[-23,23],[-56,-12],[-78,39],[-40,28],[-33,42],[-57,-2],[-59,29],[-46,3],[-32,42],[40,17],[-18,47],[-103,2],[-79,-103],[-64,-11],[-16,-50],[-52,-33],[-6,30],[21,32],[-4,67],[39,63],[97,66],[73,17],[69,32],[33,-6]],[[11990,24489],[-36,-29],[-29,77],[12,33],[34,-12],[19,-69]],[[24487,29602],[-42,46],[-82,14],[-23,-14],[-84,21],[-26,23],[-47,-14],[-30,-45],[-32,39],[-30,-2],[-6,-91],[48,-64],[11,-56],[61,-96],[46,-58],[66,-104],[-5,-23]],[[24312,29178],[-84,100],[-62,33],[-51,-5],[-5,30],[-54,55],[-47,77],[18,27],[-46,71],[-5,83],[-57,57],[-52,-111],[-41,60],[-7,89]],[[24427,29067],[-105,102]],[[24322,29169],[95,-74]],[[21765,20871],[4,8]],[[21469,22042],[52,-40],[16,-69],[42,-47],[36,-15],[8,21],[48,34],[81,0],[29,-46],[21,-57],[15,12]],[[21777,20880],[-40,0],[-85,22],[-188,-31],[-43,-17],[-115,-75],[-86,-70]],[[11950,21525],[-30,62],[-4,52],[-32,7],[-15,-38],[-24,37],[15,41],[-3,51],[-32,53],[-94,81],[-6,56],[-22,-18],[-28,-52],[-24,50],[-38,20],[-22,51],[-4,58],[20,76],[-29,34],[20,36]],[[11857,22150],[23,-99],[40,-93],[43,-82],[26,-20]],[[25098,17368],[-17,-32],[-53,9],[-30,-22],[-73,4],[-35,-30],[-13,66],[13,24],[3,85],[-13,78],[-47,126],[11,170],[-1,77],[-15,104],[8,98],[-11,25],[-141,6],[-15,80],[-75,-16],[-43,-34],[-6,-72],[-13,-40],[-4,-69],[-54,-7],[-30,14],[-91,-29],[-47,-2],[-47,135],[-21,45],[-8,71],[-21,69],[-3,83],[-23,60],[-27,22],[-72,1],[-164,-7],[-51,8],[-112,-1]],[[23757,18467],[-35,-3],[-46,-24],[-25,50]],[[23651,18490],[36,14],[2,120],[7,33],[31,57],[30,23]],[[23757,18737],[42,-45],[33,33],[4,51],[28,-7],[51,41],[7,-117],[36,-11],[49,93],[61,94],[34,21],[33,103],[9,96],[0,187],[40,74],[42,135],[49,49],[58,99],[-4,60],[21,112],[4,67],[-3,127],[19,92],[2,105],[17,88],[35,112],[16,83],[-2,38]],[[24438,20517],[3,105],[-8,66],[33,59],[29,80],[53,52],[38,-8],[51,-57],[41,-81],[49,-11],[72,-37],[45,2],[63,-26],[42,113],[44,19],[37,-16],[111,73],[55,-14],[59,21],[34,63],[80,-31],[80,-24],[36,30],[35,-20]],[[25995,20228],[-48,-66],[-48,-86],[0,-50]],[[25801,19727],[-23,-13],[-21,-76],[39,5]],[[25744,19384],[-18,17],[-20,-88],[-7,-96]],[[25719,19157],[25,-73],[0,-189],[14,-40],[8,-78]],[[23784,20228],[147,-1],[11,9],[40,-41],[37,3],[17,-17],[49,-7],[40,-52],[9,10],[-7,84],[13,36]],[[24140,20252],[35,123],[3,83],[22,72],[41,4],[53,29],[89,-41],[39,27],[16,-32]],[[23757,18737],[-24,41],[-61,-38],[-45,-85]],[[23627,18655],[-29,72],[0,25],[-80,142]],[[27497,17142],[-30,32],[19,31],[11,-63]],[[13920,20022],[-27,-3],[0,89],[-32,104],[-18,4],[-39,-69],[-44,-16],[-141,1],[-32,10],[-33,-13],[-1,-143],[61,1],[24,-44],[-1,-50],[-37,21],[-71,-33],[-2,-158],[49,-75],[15,-96],[18,-60],[-33,-324],[-37,-344]],[[12445,20073],[-20,37],[28,50],[27,-17],[-5,117],[39,45],[22,-5],[55,86],[53,155],[9,61],[-32,43],[18,145],[-7,26],[-8,127],[-16,26],[35,53],[-27,88],[14,71],[-67,160]],[[12628,21654],[64,-124],[10,48],[-14,68],[24,16],[51,65],[17,60],[29,37],[32,4],[-5,61],[11,58],[-10,45],[28,91],[74,109],[48,-26],[13,-39],[18,99],[18,17],[32,-14],[59,4],[73,91],[55,39],[17,66],[52,51],[27,3],[29,-22],[15,-63],[-22,-41]],[[34593,24488],[-2,39],[-35,-20],[-37,-43],[-39,1],[-36,78],[-1,46],[-89,21],[18,82],[6,72],[15,29],[-21,30],[-59,21],[-3,95],[-19,48],[10,60],[-68,-1],[-65,-47],[-2,106],[-17,24],[7,63],[17,14],[12,90],[23,9],[20,67],[20,-2],[40,72],[9,147],[-10,217],[-18,16],[-25,-21],[-25,130],[-28,48],[-41,32],[-26,-64]],[[33424,25847],[-4,42],[-46,26],[-19,-21],[-44,21],[-45,2],[0,36],[-46,15],[-54,-45],[-73,-141],[-7,-33]],[[31722,27545],[-44,-5],[-84,41]],[[31321,27880],[-21,30],[43,29],[21,-13]],[[38490,29645],[-28,51],[-43,7],[-25,-21]],[[37451,28534],[-32,-36],[-93,-17],[-39,-19],[-86,-68],[-26,-47],[-58,-63],[-40,-12],[-27,24],[71,42],[13,64],[-63,-4],[-1,35],[31,21],[0,45],[35,24],[48,89],[10,40],[-51,66],[-84,13],[-50,-68],[-36,-79],[-109,-72],[-25,-33],[-20,-76],[-39,-53],[-35,3],[-40,-24],[-32,35],[-31,-21],[-28,-111],[25,-69],[31,-28],[114,-31],[18,-75],[-17,-81],[20,-29],[41,-17],[39,7],[15,43],[60,77],[46,29],[58,-51],[52,-31],[86,-12],[32,4],[-10,-105],[-22,-25],[-50,28],[-116,-79],[-38,-57],[0,-38],[-35,-25],[-37,14],[-5,-56],[-77,-123],[-26,-64],[16,-57],[29,-37],[76,-56],[15,-34],[22,-117],[45,-136],[-2,-78],[68,-64],[0,-36],[53,-68],[-19,-41],[-40,32],[-37,-30],[75,-88],[26,-88],[-43,-17],[-65,-62],[-22,-45],[-45,8],[22,-56],[33,6],[32,31],[33,-17],[31,-54],[32,-19],[3,-88],[-6,-78],[-25,26],[-10,-107],[-19,-28],[15,-60],[-30,-30],[-32,13],[-17,-55],[-46,-101],[-5,-57],[-33,-49],[-23,-93],[-37,-28],[12,-46],[-52,-61],[23,-30],[-14,-72],[-34,-27],[-3,-58],[-26,4],[-14,-65],[-43,-76],[-58,10],[-10,-32],[6,-51],[-53,-90],[-20,1],[-74,-86],[-40,-60],[-8,-51],[-76,-32],[-26,12],[-18,-30],[-37,22],[-42,-44],[-25,30],[-19,-50],[-28,3]],[[36209,24701],[-31,-6]],[[36178,24695],[-23,21],[-26,56],[-35,11],[27,-70],[-3,-76]],[[36118,24637],[-6,-12]],[[36112,24625],[-26,-23],[-32,7],[-52,-75],[-35,-8],[-44,24],[-42,-64],[-62,-17],[-40,-21],[-71,-77],[-2,-39],[41,-75],[-8,-35],[-40,-21],[-57,144],[2,47],[31,77],[-47,12],[-49,-25],[-15,45],[-86,14],[-10,-29],[-34,-11]],[[35793,24142],[15,-74],[-23,-21],[-31,-85],[-15,-93],[-114,-116],[-61,32],[-40,37],[-8,82],[7,94],[60,74],[5,37],[53,31],[31,-5],[33,20],[51,-17],[11,36],[26,-32]],[[36118,24637],[-6,-12]],[[36209,24701],[-31,-6]],[[13595,15911],[25,-96],[27,-23],[0,-51],[15,-150],[62,-103],[-29,-63],[15,-30],[-18,-46],[-2,-86],[67,-184],[1,-70],[25,-97],[13,-95],[0,-73],[37,-15],[47,16]],[[13880,14745],[23,-40],[-42,-226],[-111,-79],[-31,-52],[-7,-46],[22,-56],[-19,-31],[-8,-55],[23,-147],[-20,-80],[32,-79],[-3,-33],[-62,-27],[-40,-169],[-60,-108],[-21,-151],[-24,-49],[12,-97],[-3,-130],[-26,-18],[-23,-101],[-21,-54],[-8,-92],[51,-198],[19,-176],[27,-28],[-12,-73],[5,-108],[-24,-17],[-29,-95],[-23,-130],[11,-98],[-3,-63],[-80,-102],[-17,-70],[10,-59],[-6,-143],[17,-61],[9,-153],[-41,-23],[-17,-31],[0,-66],[-14,-70],[-22,-7],[7,-84],[-20,-35],[16,-34],[-30,-79],[8,-44],[-5,-166],[19,-99],[-33,-10],[-12,-33],[-3,-152],[45,-39],[-15,-60],[19,-33],[8,-74],[-16,-39],[-1,-61],[75,-13],[-6,-70],[-86,-6],[45,-42],[30,-55],[-19,-62],[-33,-46],[17,-51],[-30,-45],[23,-92],[-31,-55],[5,-81],[-54,-64],[-22,-85],[24,-51],[-4,-56],[-31,-34],[-1,-60],[-45,-40],[-19,-71],[-40,-3],[-14,-59],[13,-47],[-4,-72],[14,-23],[29,-112],[79,29],[29,-66],[-16,-139],[56,-74],[4,-24],[241,-4],[93,-29],[94,-48]],[[13727,8262],[-98,33],[-40,-47],[-152,-76],[-19,-101],[-5,-121],[-37,-23],[-108,55],[-25,35],[8,48],[56,-7],[56,39],[30,58],[-29,18],[-91,-66],[-40,-42],[-11,-45],[-62,47],[20,82],[53,13],[79,60],[-102,2],[-46,-79],[-40,-36],[-83,91],[-34,114],[56,-25],[39,42],[-29,54],[-27,3],[10,99],[-113,59],[-35,84],[51,4],[-11,42],[18,61],[36,48],[5,81],[-9,45],[3,177],[-29,70],[-1,61],[90,-10],[40,-34],[-23,119],[-28,-60],[-35,-9],[-53,59],[22,59],[40,55],[0,45],[-105,52],[-17,31],[-63,-3],[90,99],[-17,63],[74,7],[38,16],[14,77],[60,-14],[23,111],[35,12],[42,41],[9,65],[-74,59],[5,60],[28,58],[-12,39],[27,89],[11,164],[36,70],[-40,16],[19,55],[-43,25],[-28,-52],[-34,-3],[-36,61],[-19,88],[37,227],[32,65],[22,123],[-31,132],[-6,57],[8,71],[-24,75],[7,112],[47,5],[13,105],[29,65],[12,91],[24,48],[-5,38],[50,108],[20,105],[7,99],[21,76],[20,36],[-9,123],[18,27],[17,68],[4,60],[-11,40],[-2,88],[-16,139],[-6,119],[5,66],[33,41],[10,108],[-4,66],[-17,33],[-4,60],[26,55],[15,65],[17,143],[15,30],[4,88],[31,189],[1,74],[-10,45],[28,90],[5,45],[-15,116],[10,189],[12,47],[-24,44],[3,68],[28,46],[30,297],[0,53],[-13,116],[6,202],[-16,118],[-10,191],[-7,11]],[[13704,8197],[-3,-483]],[[13701,7714],[-53,-12],[-50,11],[-29,32],[-38,-23],[-57,2],[-116,41],[67,62],[24,-32],[50,0],[-70,82],[-1,50],[42,56],[18,-79],[-30,-4],[56,-82],[53,16],[-30,43],[-20,77],[98,62],[-35,32],[-56,-19],[-45,47],[8,39],[77,65],[33,52],[41,-38],[66,3]],[[13895,7648],[-85,-23],[-37,8],[-5,64],[106,-10],[21,-39]],[[13071,10240],[-42,-3],[-30,21],[18,84],[6,146],[15,92],[63,-22],[7,-125],[-39,-26],[43,-77],[-38,-49],[-3,-41]],[[12984,8966],[-14,-188],[-87,35],[-1,77],[27,50],[-12,40],[29,28],[19,56],[28,-11],[11,-87]],[[12857,9051],[-17,39],[16,56],[-5,57],[21,11],[29,-89],[-1,-43],[-43,-31]],[[12973,9088],[-44,-7],[-35,106],[29,48],[30,-44],[20,-103]],[[13175,8015],[55,-23],[34,-48],[-20,-38],[-62,-32],[13,57],[-66,-17],[-47,113],[93,-12]],[[12942,8422],[5,-43],[31,-37],[-40,-61],[-31,106],[35,35]],[[13572,7700],[98,-22],[77,-52],[-2,-59],[-70,16],[-3,41],[-50,21],[-14,-61],[-70,66],[34,50]],[[13364,7894],[45,-17],[-14,-58],[-41,32],[-61,5],[-30,35],[27,36],[74,-33]],[[13076,10010],[-31,-22],[32,-137],[-15,-49],[-31,3],[-31,91],[-33,58],[14,38],[50,19],[24,56],[21,-57]],[[13168,9925],[-45,-12],[-12,43],[29,67],[55,-47],[-27,-51]],[[12907,9035],[-19,-51],[-28,-3],[-16,43],[63,11]],[[12913,8714],[-48,-10],[15,73],[55,-22],[-22,-41]],[[13005,9725],[-47,6],[47,108],[0,-114]],[[24961,22151],[-45,16],[-55,-38],[-42,-54],[3,-46],[-41,-80],[-16,-6],[-60,-125],[-53,-61],[-34,1],[-49,-24],[-65,-1],[-31,-28],[27,-51],[-25,-55],[-56,-82],[-99,-11],[-49,-37],[-57,-58],[-32,68],[-18,-39],[-66,-45],[-45,10]],[[24053,21405],[9,58],[-13,14],[-40,155],[-97,142],[-40,84],[-4,23],[33,63],[73,-8],[87,4],[-1,28],[-32,59],[-25,108],[-5,58],[11,94],[-5,67],[-24,93],[-15,85],[-27,36],[-3,35],[-56,21]],[[24140,20252],[-10,45],[-5,95],[-35,48],[-80,154],[1,53],[-14,63],[-39,68],[-8,113],[-3,150],[-23,38],[38,53],[36,111],[26,105],[29,57]],[[19264,23048],[-32,-10],[0,68],[32,-58]],[[5826,31634],[117,-31],[-38,-158],[-63,-16],[-28,26],[-37,103],[-2,80],[51,-4]],[[5937,31430],[14,-60],[-11,-37],[-54,6],[-22,78],[73,13]],[[6497,30869],[61,-28],[150,-43],[104,-195],[79,-45],[22,-31],[40,-110],[-12,-55],[-135,58],[-53,36],[24,56],[-66,-16],[-45,33],[-1,40],[-50,37],[-38,-5],[-17,100],[-53,1],[-43,62],[-49,-8],[-5,70],[-44,36],[3,48],[40,7],[88,-48]],[[4799,35041],[126,-6],[98,-23],[129,-80],[174,-58],[103,6],[23,87],[68,34],[155,9],[88,-16],[40,47],[111,23],[177,85],[105,-31],[-162,-81],[-111,-18],[-101,-63],[-85,-82],[82,-16],[-1,58],[98,46],[16,31],[96,0],[81,49],[134,44],[36,-31],[148,112],[-49,43],[51,22],[75,-61],[65,-112],[68,-58],[83,-26],[4,61],[72,80],[78,-74],[-36,-60],[101,-1],[48,36],[24,59],[117,2],[136,-35],[83,-52],[175,-36],[95,-47],[119,-30],[131,-11],[53,26],[150,-69],[55,-57],[-16,-29],[-82,0],[-81,-76],[36,-23],[152,-24],[180,-5],[177,23],[117,41],[51,-54],[127,-31],[75,-73],[-58,-43],[160,-39],[-100,186],[26,66],[114,21],[129,88],[-83,-3],[-18,-35],[-125,-5],[-14,-35],[-65,-5],[-52,28],[46,73],[108,17],[156,51],[77,-44],[14,-57],[152,-94],[88,18],[97,-66],[139,-26],[136,32],[161,-26],[90,23],[118,-46],[31,52],[-152,46],[27,56],[117,29],[36,-43],[116,-21],[-42,-133],[73,-79],[64,15],[-48,104],[23,63],[88,11],[135,103],[-4,76],[-98,-32],[-10,37],[59,47],[-23,72],[-228,92],[-52,100],[39,47],[-40,59],[18,102],[148,138],[80,16],[30,-47],[130,-66],[37,-44],[-4,-90],[73,-44],[102,-104],[-146,-102],[103,-40],[167,-6],[-33,-47],[44,-93],[-13,-86],[40,-45],[45,56],[26,108],[75,55],[125,-100],[29,-89],[-66,-26],[18,-114],[115,-128],[89,73],[21,68],[51,54],[27,122],[57,25],[50,73],[-64,34],[-15,135],[156,-2],[71,-30],[129,-2],[-3,-51],[161,-74],[-76,-46],[70,-13],[12,-43],[-89,-42],[-67,-4],[71,-128],[86,-88],[-24,-86],[-50,-18],[-94,-88],[-146,-65],[-80,-25],[-76,32],[-43,47],[-146,-1],[5,-46],[72,-30],[-5,-36],[-153,-147],[-84,-1],[-219,129],[-45,-12],[140,-117],[240,-33],[-31,-82],[-102,-141],[-67,-38],[-92,7],[-91,-13],[16,-39],[-162,-65],[71,-34],[7,-53],[-22,-36],[-71,-31],[-112,3],[19,-51],[-43,-9],[-57,-100],[-44,-35],[-6,-49],[-78,-86],[-1,-39],[-72,-158],[-17,-102],[0,-151],[9,-80],[53,-41],[124,9],[103,-325],[-21,-52],[45,-6],[140,51],[63,-4],[99,-53],[104,-29],[107,-84],[63,-90],[232,-100],[76,-70],[44,-5],[98,12],[164,-37],[44,-73],[-25,-101],[34,-118],[-2,-121],[-12,-67],[81,-116],[-8,-30],[86,-71],[38,-47],[35,-94],[66,-34],[42,87],[94,-16],[-30,52],[59,114],[-28,82],[0,50],[-22,41],[-26,148],[14,63],[-30,21],[-52,131],[93,41],[125,78],[70,70],[83,121],[16,131],[-7,104],[-39,127],[-32,57],[-106,84],[-62,62],[11,48],[82,105],[5,65],[48,27],[2,53],[-45,85],[-23,79],[-36,14],[41,51],[11,78],[28,26],[-54,45],[-23,76],[8,54],[86,48],[61,-11],[130,-46],[59,0],[80,-29],[131,60],[190,-160],[55,-24],[-25,-33],[51,-67],[86,-23],[55,3],[78,-83],[-14,-67],[16,-124],[-7,-105],[35,-1],[-2,-55],[30,-42],[58,2],[46,-68],[93,-84],[119,74],[30,53],[39,6],[54,64],[21,146],[32,29],[30,75],[53,4],[70,-137],[48,-68],[44,-105],[73,-84],[-10,-34],[45,-76],[50,-24],[-8,-55],[14,-51],[64,-81],[-5,-72],[-60,-22],[42,-42],[53,-115],[31,20],[79,-105],[-32,-61],[64,13],[115,-24],[27,-69],[47,-14],[25,20],[97,-65],[-36,-40],[-56,-8],[-59,-66],[-26,-1],[-80,-46],[-41,0],[-39,-52],[39,-21],[62,31],[50,51],[73,39],[16,38],[96,-14],[32,-69],[-41,-47],[19,-36],[60,59],[47,6],[43,-39],[39,-83],[-5,-198],[17,-38],[-39,-45],[-118,-104],[-98,-7],[-31,-23],[-60,-5],[-81,-114],[-112,-115],[-149,-11],[-54,-22],[-24,28],[-97,15],[-65,-13],[-61,14],[-139,-6],[-50,9],[-105,-26],[-46,3],[-54,-48],[-53,-142],[-112,-33],[-80,-81],[-55,-97],[-49,-63],[-16,-59],[-67,-89],[-32,-63],[-69,-75],[-75,-24],[-99,-88],[-37,-17],[-61,-103],[-29,-6],[-40,-46],[8,-34],[-90,-75]],[[12822,29519],[-32,-27],[-90,-41],[-33,-65],[-65,34],[-182,-45],[-62,-56],[-25,-57],[53,-27],[34,19]],[[12420,29215],[1,-53],[-85,-1],[-58,-19],[-29,-35],[-98,9],[-30,-14],[-48,-70],[-71,-47],[-59,-11],[-17,59],[60,8],[-7,54]],[[12007,29196],[29,15],[49,53],[9,32],[-2,129],[16,55],[40,73],[-1,32],[-38,82],[36,18],[5,-48],[29,-50],[46,-17],[23,-32],[46,-19],[-3,69],[35,29],[-32,38],[-90,185],[-96,9],[-2,27],[-71,4],[-19,13],[-84,-3],[-131,39],[-9,45],[-33,35],[-41,69],[21,62],[-51,70],[19,68],[-93,-2],[-37,20],[-23,39],[-36,115],[-31,14],[-71,-4],[-91,40],[-33,-80],[-41,10],[-78,-42],[-18,-72],[-29,-26]],[[7040,30507],[-32,0],[-36,103],[-64,5],[-35,68],[-42,3],[-37,47],[-33,87],[-166,25],[1,60],[-29,24],[-53,-10],[-80,62],[8,72],[-56,68],[-32,81],[40,35],[-29,26],[20,93],[-42,4],[-25,80],[-78,33],[-34,-31],[-95,104],[36,90],[-48,63],[67,163],[-25,102],[8,57]],[[8636,37017],[52,-69],[-126,-9],[-93,20],[-246,-23],[3,28],[337,76],[73,-23]],[[8557,36897],[84,-32],[-52,-95],[-216,-40],[-149,41],[-6,82],[339,44]],[[7930,36734],[-90,-57],[53,-29],[-45,-79],[-222,-36],[-165,-49],[-67,-80],[-61,-7],[-76,42],[-109,4],[-107,42],[51,38],[92,9],[311,189],[133,17],[72,-14],[130,50],[100,-40]],[[8823,36448],[143,-40],[91,40],[97,-28],[15,-41],[-47,-121],[-158,-59],[-148,4],[-59,27],[-83,-40],[-82,-10],[-94,-44],[-195,-49],[-122,3],[-93,39],[17,35],[159,46],[-198,25],[-158,-27],[-92,45],[-138,22],[55,44],[22,60],[72,14],[-27,58],[127,80],[170,3],[46,-54],[138,-1],[146,-86],[57,-57],[242,-9],[17,42],[-105,36],[38,62],[-96,59],[115,76],[87,-38],[62,-78],[-21,-38]],[[8057,35687],[126,22],[36,63],[55,1],[183,-59],[60,-40],[71,28],[-56,76],[80,-5],[109,-57],[52,-49],[61,-165],[35,-24],[76,56],[-48,56],[-67,210],[64,49],[74,-30],[78,1],[131,-90],[23,-82],[109,-215],[-27,-73],[75,-75],[105,-55],[43,3],[94,-52],[103,-31],[27,-94],[-88,-4],[-66,32],[-54,-64],[80,-32],[14,-84],[-127,-44],[-71,-3],[-109,26],[-92,-2],[10,35],[-138,18],[-64,61],[-96,-96],[-114,-15],[-71,-38],[-125,-29],[-166,-19],[-224,-10],[-60,75],[-9,78],[-78,17],[-102,-1],[-167,35],[-73,83],[-4,65],[313,47],[241,-20],[105,11],[-41,38],[-202,54],[-275,-23],[-196,9],[-84,57],[56,59],[195,44],[30,25],[-271,-8],[55,49],[-139,6],[-9,66],[90,61],[-20,59],[101,65],[120,50],[239,69],[72,-66],[-58,-105]],[[7415,36021],[218,33],[95,-20],[207,-122],[8,-39],[-432,-166],[-55,-71],[-95,-32],[-33,-128],[-21,-27],[-319,-84],[-98,121],[-159,65],[-57,36],[95,99],[19,107],[128,152],[-57,41],[-54,86],[416,40],[94,-39],[122,-26],[-22,-26]],[[13598,37975],[580,-31],[160,-23],[48,-45],[208,-27],[15,-36],[-111,-52],[-295,-68],[-112,-98],[-387,-133],[-95,-56],[-91,-6],[-154,-139],[-226,-26],[-29,-33],[-221,-16],[20,-48],[-143,-43],[225,-62],[-177,-157],[-136,-18],[-128,4],[-8,-94],[-97,-41],[-210,-3],[37,-36],[115,1],[15,-49],[114,10],[38,-49],[-37,-41],[-151,-57],[-145,-29],[-59,74],[-334,-14],[-155,-32],[-142,41],[-47,-25],[-211,6],[-145,19],[9,73],[135,61],[188,21],[-112,65],[-38,49],[132,38],[213,-73],[-32,81],[-82,34],[-165,21],[8,53],[84,78],[194,28],[241,-29],[26,37],[-154,8],[14,55],[-110,85],[-163,51],[14,104],[318,-20],[229,-111],[104,-6],[-240,128],[53,21],[305,34],[205,56],[-30,57],[-292,-85],[-295,13],[-83,-35],[-151,-7],[-445,62],[-126,80],[119,81],[-328,32],[142,42],[481,68],[207,58],[291,-45],[51,78],[275,74],[296,-12],[213,36],[297,-11],[64,18],[350,8],[59,-22]],[[10399,36828],[269,-11],[-27,-58],[-303,2],[-20,50],[81,17]],[[10638,36222],[0,-74],[-142,-11],[-193,61],[-37,39],[110,95],[153,14],[84,-54],[25,-70]],[[9844,35984],[103,-47],[194,34],[72,-48],[-34,-62],[-64,-22],[9,-54],[48,-22],[35,-66],[75,-33],[-44,-38],[22,-105],[-119,-45],[-73,7],[-3,-49],[-57,-30],[-62,14],[-70,86],[-105,86],[-76,37],[-63,-1],[-122,103],[26,49],[71,11],[71,-68],[106,6],[25,76],[-144,68],[59,37],[7,45],[113,31]],[[11699,34078],[38,41],[65,-51],[83,-25],[91,-71],[76,-29],[33,-48],[9,-89],[109,15],[65,-72],[-92,-66],[-164,54],[-11,48],[-126,38],[-29,-63],[-71,-50],[-40,-61],[-106,-37],[-27,114],[-182,3],[37,55],[78,47],[-15,94],[51,250],[50,47],[39,-27],[0,-62],[39,-55]],[[10127,36537],[12,-135],[24,-65],[-58,-30],[10,-66],[-95,-23],[-205,-1],[-59,88],[98,60],[-323,-37],[36,58],[101,47],[-87,64],[215,73],[168,0],[163,-33]],[[9422,37163],[95,-49],[96,-1],[279,-111],[-18,-62],[77,-40],[-5,-57],[-174,7],[-60,66],[-204,39],[-196,-22],[99,110],[-44,29],[-204,29],[18,63],[241,-1]],[[10842,37562],[153,-119],[246,-45],[119,-3],[86,-105],[-20,-51],[202,-29],[26,-71],[-213,-68],[-218,-155],[-224,-9],[-158,19],[-96,33],[-183,133],[199,63],[-315,2],[-77,28],[53,52],[-143,42],[-38,64],[142,57],[-63,65],[106,68],[250,28],[166,1]],[[10546,36635],[59,1],[135,-69],[174,17],[55,-53],[193,-30],[-119,-54],[253,-119],[116,23],[146,-25],[180,33],[83,36],[85,-15],[134,18],[225,-44],[82,-40],[18,-125],[-86,-29],[-8,-36],[-189,-23],[-212,24],[-110,-17],[-207,6],[-183,-15],[-102,5],[-43,49],[-132,-37],[-115,33],[-89,-10],[-56,31],[-14,73],[-43,49],[44,68],[-15,41],[-111,112],[-71,-18],[-198,-2],[-63,60],[-105,36],[-14,59],[137,20],[166,-33]],[[10311,36991],[162,-41],[-54,-27],[29,-51],[-234,-41],[-155,137],[5,83],[200,-26],[47,-34]],[[10684,36031],[117,-41],[139,8],[91,-31],[-156,-185],[-61,-60],[-151,11],[-53,-29],[27,-55],[-60,-86],[-142,0],[-7,104],[-36,62],[-11,200],[73,72],[69,20],[161,10]],[[10159,35040],[140,-66],[76,-133],[-88,-59],[-132,16],[-51,29],[-168,42],[-39,37],[137,66],[-12,46],[35,46],[59,24],[43,-48]],[[14630,29838],[36,-45],[42,46],[-4,56],[78,-10],[-18,-61],[-46,-50],[-35,-14],[-75,-3],[-26,80],[77,188],[31,39],[35,-45],[-30,-131],[-65,-50]],[[12956,29632],[141,92],[49,129],[53,37],[29,38],[73,61],[79,27],[91,61],[116,174],[13,31],[60,68],[92,77],[83,50],[170,79],[81,11],[85,-16],[70,-59],[1,-82],[-63,-71],[-60,-45],[-83,37],[-34,-44],[44,-17],[22,-47],[54,25],[65,-19],[-26,-78],[-50,-59],[60,-9],[-9,-38],[22,-47],[23,-95],[77,-16],[-17,-32],[94,-59],[72,-3],[27,-26],[65,54],[10,-33],[47,-6],[26,-54],[0,-46],[-152,-86],[-97,-44],[-29,3],[-27,-37],[-29,35],[-36,-21],[-11,-54],[-60,-96],[-57,-44],[-19,-32],[-31,9],[-19,51],[-29,4],[-8,73],[11,49],[55,86],[92,79],[56,30],[27,-31],[18,65],[-63,0],[-37,47],[-121,-88],[-28,7],[-41,-35],[-53,-6],[-31,22]],[[15324,31065],[-43,-75],[-7,-63],[-102,-186],[-19,-64],[12,-25],[62,96],[60,-28],[-16,-84],[26,-32],[38,12],[43,-40],[55,40],[76,-10],[42,-27],[5,-39],[-34,-73],[7,-62],[95,10],[-14,-38],[-48,-37],[-21,-68],[12,-57],[39,74],[41,7],[-28,-83],[65,-23],[-32,-98],[-20,-92],[-67,0],[4,60],[-72,-16],[37,111],[-15,81],[-25,23],[-30,-88],[-51,-17],[-57,-105],[-58,-8],[-6,45],[42,19],[37,63],[-21,47],[-26,-43],[-35,14],[2,57],[-46,-26],[-90,-21],[-120,22],[-50,0],[-97,-23],[-30,65],[11,28],[69,73],[27,41],[-27,19],[38,106],[38,-1],[-14,68],[24,34],[88,254],[21,19],[47,107],[89,64],[41,-21],[28,14]],[[11494,35340],[116,31],[105,-13],[16,108],[-72,20],[-82,68],[131,92],[-60,4],[-38,77],[47,39],[-23,34],[53,50],[290,87],[155,-21],[31,-71],[67,-42],[41,-82],[-49,-41],[-33,-77],[101,1],[98,39],[40,-26],[87,19],[-23,33],[59,39],[136,7],[210,-63],[35,-73],[78,-14],[-6,-43],[74,-26],[29,-53],[71,43],[178,-47],[82,-79],[28,-80],[75,24],[96,-18],[223,-160],[12,-69],[-105,7],[-49,-39],[56,-22],[86,-1],[73,-30],[-4,-40],[-84,3],[-33,-28],[-35,-107],[157,-28],[-34,-21],[100,-81],[113,-4],[35,-34],[61,6],[45,-42],[18,-60],[52,6],[68,-31],[-27,-46],[108,-26],[50,25],[83,-82],[-70,-73],[-82,-19],[67,-42],[-83,-87],[-86,-21],[-5,-93],[-92,-12],[-63,23],[-91,130],[12,54],[-68,17],[-76,42],[-39,61],[-64,-55],[16,-59],[-45,-25],[-87,4],[117,-100],[-3,-28],[57,-58],[79,-52],[56,-1],[61,-48],[32,-137],[41,4],[20,-167],[-68,-1],[49,-75],[-61,3],[-12,-48],[-88,63],[-81,17],[-67,50],[-65,72],[-16,-41],[-80,61],[-52,-5],[95,-117],[57,-17],[86,-86],[41,-22],[45,-65],[27,-86],[-24,-10],[-146,62],[-115,19],[-87,37],[-57,75],[-85,4],[-125,61],[-89,85],[64,18],[-33,43],[-61,0],[-151,157],[-131,54],[-68,-47],[-52,19],[-34,-33],[-118,-33],[-131,28],[-48,55],[10,70],[80,49],[15,64],[104,-18],[80,-32],[104,34],[175,36],[-23,49],[-84,82],[141,117],[31,12],[30,66],[70,49],[-104,188],[-32,35],[-163,96],[-33,61],[-87,-21],[-96,-43],[-19,70],[75,4],[20,63],[-165,74],[21,38],[-83,16],[-26,76],[-54,-4],[-96,83],[-41,-37],[72,-57],[-10,-44],[-89,-18],[-182,43],[-13,-24],[-115,-31],[-118,36],[-94,-9],[-194,32],[-109,8],[-37,56],[-62,2],[-95,-35],[-117,60],[-44,52],[-20,67],[187,-27],[-3,58],[-148,18],[-95,45],[-22,100],[45,46],[-25,55],[62,87],[12,58],[69,73],[121,70],[117,25],[204,-6],[17,-25],[-137,-95],[-67,-86],[38,-90],[-3,-74],[38,-77],[125,-92],[-47,-28],[-146,-47]],[[14544,30529],[-51,-3],[-102,32],[-64,38],[-39,57],[-69,37],[38,25],[104,-27],[53,-24],[128,-93],[2,-42]],[[14297,29953],[16,21],[68,-31],[118,14],[18,-14],[-65,-56],[2,-42],[-43,6],[-29,55],[-20,-22],[-44,17],[-21,52]],[[12059,33572],[-14,-67],[-111,-97],[-84,-11],[-25,70],[64,94],[170,11]],[[12361,33452],[32,-36],[-5,-48],[-48,-91],[-69,54],[5,65],[45,55],[40,1]],[[12837,34750],[64,-20],[3,-153],[-80,-55],[-151,-4],[-37,97],[70,115],[131,20]],[[12362,35920],[154,3],[132,-37],[113,-92],[-13,-57],[-177,17],[-205,-31],[-54,24],[-30,79],[-79,34],[-3,75],[54,10],[108,-25]],[[12215,31331],[-35,3],[-102,43],[-1,51],[88,4],[46,-62],[4,-39]],[[12436,32104],[-30,-84],[-28,11],[-44,-28],[-29,27],[49,69],[25,7],[15,71],[38,-39],[4,-34]],[[11095,36713],[-162,26],[0,49],[120,-3],[42,-72]],[[9283,36686],[-81,8],[-105,120],[59,2],[129,-87],[-2,-43]],[[9193,35759],[-110,84],[-91,39],[38,47],[136,16],[106,-34],[8,-57],[-87,-95]],[[9570,36438],[-44,-51],[-90,-4],[-128,92],[155,28],[63,-6],[44,-59]],[[9349,36563],[119,-29],[-150,-27],[31,56]],[[7589,36343],[-132,8],[95,63],[123,13],[-86,-84]],[[6126,31831],[-44,-56]],[[6330,31373],[6,-92],[-71,67],[-11,43],[17,57],[50,-37],[9,-38]],[[14171,33223],[-69,43],[32,27],[57,-20],[-20,-50]],[[13090,34632],[-118,9],[24,59],[116,-28],[-22,-40]],[[12566,33685],[42,-40],[-50,-38],[-65,53],[73,25]],[[12375,35071],[-15,-34],[-75,-21],[-61,47],[151,8]],[[12402,34865],[75,115],[46,-15],[-77,-85],[-44,-15]],[[11493,34621],[-33,16],[-12,64],[67,18],[4,-56],[-26,-42]],[[9337,36223],[-95,25],[67,62],[87,-59],[-59,-28]],[[11940,29814],[78,32],[48,0],[43,-41],[-36,-57],[-133,66]],[[23354,20260],[11,134],[8,36],[-9,50],[-25,50],[-31,85],[-52,37],[-12,75],[-43,71]],[[34998,22039],[-20,28],[-40,10],[-34,-17],[-16,21],[23,63],[-23,56],[-22,-49],[-25,-2],[-7,66],[4,54],[-24,54]],[[34294,21990],[-23,36],[-7,71],[26,68],[9,56],[-1,77],[-13,32],[8,107],[-15,206],[-41,125],[-16,-4],[-2,98],[-35,153],[-13,223],[-15,32],[4,114],[-30,-4],[-24,119],[-26,56],[-19,-119],[-23,-47],[-75,-68],[-31,-18],[-10,-42],[-36,-56],[-26,22],[-29,-1],[-6,35],[-38,68],[-35,-62],[-11,26],[7,86],[10,27],[29,204],[-19,139],[-22,89],[-22,63],[-45,29],[36,91],[-48,73],[-31,62],[-44,4],[-26,30],[12,28],[-19,49],[-15,-20],[-35,70]],[[33509,24317],[-18,110],[19,33],[32,6],[-2,112]],[[22442,22365],[-12,-45],[-38,-60],[-68,6],[-35,-43],[-10,-42],[-20,-14]],[[25666,29355],[-15,-78],[-41,2],[-25,-46],[-5,-96],[-50,-62],[28,-26],[38,-83]],[[36302,20829],[14,0]],[[36184,20762],[29,3],[85,91],[4,-27]],[[15793,20645],[23,56],[30,-49],[18,-92],[3,-86],[48,-252],[25,-67],[33,-9],[28,-28],[10,-53],[-2,-56],[-49,-72],[-35,-91],[-29,-54],[-58,-57],[-15,-67],[-36,-82],[-3,-56],[-24,-35],[-12,-48],[-26,8],[4,-61],[30,12],[81,79],[46,24],[28,-142],[36,-56],[49,41],[35,-21],[50,43],[-30,-173],[12,3],[30,129],[27,19],[35,75],[32,-7],[0,82],[41,90],[38,16],[19,-12],[31,20],[34,-26],[39,-7],[23,-42],[49,-14],[71,-69],[23,-2],[22,-75],[25,51],[35,-57],[17,-4],[33,-136],[-17,-9],[-24,-162],[35,44],[16,89],[24,10],[22,-20],[58,18],[10,27],[54,-19],[42,-43],[43,-29],[46,10],[29,-29],[39,-13],[104,31],[63,-14],[117,-116],[40,-58],[26,-14],[53,-110],[52,-82],[40,-27],[15,-43],[39,-12],[33,-29],[75,9],[53,-16],[20,-27],[19,-69],[17,-136],[13,-46],[23,-196],[-9,-101],[5,-49],[-40,-210],[-22,-66],[-52,-104],[-78,-171],[-46,-42],[-21,-32],[-58,-149],[-34,-132],[-41,-108],[-27,-55],[-32,-25],[-3,43],[-32,-5],[-6,-80],[-24,-48],[-7,-49],[6,-88],[12,-9],[-15,-137],[8,-131],[14,-134],[-30,-198],[-11,-121],[8,-85],[-41,-63],[-20,-57],[-11,-85],[5,-140],[-18,-82],[-20,-20],[-36,-121],[-12,-61],[-49,-74],[-31,-131],[7,-91],[-17,-36],[-71,-50],[-34,-59],[5,-46],[-13,-35],[-112,-4],[-71,-18],[-41,30],[-35,-22],[-60,-10],[-4,-61],[-35,-10],[-58,-67],[-67,-24],[-110,-97],[-34,-56],[-91,-111],[-6,-37],[-32,-31],[-26,8],[0,-72],[-18,-48],[-16,-194],[15,-108],[-11,-80],[3,-113],[-22,-110],[-58,-65],[-59,-108],[-35,-96],[-33,-137],[-39,-104],[-66,-128],[-1,80],[58,91],[-3,61],[-54,13],[-25,-68],[-13,-83],[-63,-73],[-27,-110],[8,-62],[-26,-59],[-38,-154],[-33,-58],[-56,-74]],[[15060,13128],[47,89],[35,41],[76,163],[74,142],[46,55],[30,19],[34,70],[62,28],[15,29],[45,37],[20,197],[-18,58],[-9,64],[-32,32],[-36,-23],[-21,11]],[[14992,15328],[37,40],[-33,52],[40,146],[2,36],[26,140],[-1,34],[-31,134],[-19,0],[-50,61],[-9,128],[15,35],[-23,36],[-202,13],[-8,173],[-36,74],[32,12],[-3,103],[-22,95],[9,37],[-13,50],[-70,66],[-88,-8],[-47,86],[-22,0],[-51,30],[-30,42],[-7,30],[-34,-2],[-30,35],[-43,-2],[-59,19],[-12,42],[-59,60],[-12,52],[-28,66],[-9,42],[7,78],[-12,97],[15,43],[0,84],[-10,35],[-107,-25],[-57,-33],[-47,-64],[-58,-52],[-17,-39],[-44,-4],[-40,-70],[-40,-22],[-15,24],[-90,13]],[[16042,19703],[63,15],[65,-15],[26,-27],[-20,-98],[-33,-115],[-15,-35],[-38,-25],[-33,19],[-10,-44],[-27,-17],[-25,18],[-61,-18],[-12,20],[-19,101],[-2,130],[16,82],[49,34],[76,-25]],[[16028,19813],[-12,-58],[-20,-8],[-42,16],[9,42],[48,16],[17,-8]],[[15916,19725],[-34,-43],[-8,49],[42,-6]],[[16065,19729],[-33,-7],[0,35],[25,15],[8,-43]],[[15771,19439],[3,51],[15,25],[1,51],[16,45],[36,24],[11,-28],[-14,-77],[-44,-71],[-24,-20]],[[24322,29169],[-10,9]],[[13654,16125],[22,-3],[32,68],[-29,10],[-26,87],[-40,63]],[[14439,14874],[-23,52],[-120,-2],[-20,-15],[-49,-151],[-28,117],[-56,27],[-87,1],[-41,58],[-18,1],[-15,-62],[-42,-22],[-11,-43],[-24,-22],[-25,-68]],[[13642,16193],[19,-4]],[[22482,21152],[-134,-34]],[[11284,23808],[20,-26],[4,-51],[-21,-112],[8,-21],[-11,-71],[-2,-122],[-31,-76],[-16,-9],[-25,-78]],[[22459,30969],[102,62]],[[33509,24317],[-33,84],[-6,112],[-39,180],[-26,65],[-41,-41],[-25,-3],[-38,109],[-21,-9],[-4,-66],[20,-70],[-1,-39],[-53,-97],[-8,44],[-34,-7],[-30,-47],[-41,-10],[-18,33],[-5,48]],[[12593,25076],[17,-83],[-31,-26],[-10,63],[-18,40],[42,6]],[[13163,24405],[-17,-47],[-61,-9],[17,56],[61,0]],[[12582,25177],[0,-54],[-37,-38],[-31,102],[33,85],[35,-95]],[[27657,28468],[33,14],[52,-60],[38,-25],[43,-110]],[[27869,28293],[7,65],[-42,86],[-51,63],[21,78],[-48,40],[-24,48],[26,46],[-72,98]],[[28125,28938],[66,-119],[15,-60],[67,-97],[27,-2],[39,-36],[-93,-48],[-28,-129],[4,-56],[-20,-61],[-23,13],[-20,-70],[3,-83]],[[23324,30182],[17,7]],[[39766,17130],[-8,-77],[36,-74],[22,-176],[-4,-65],[25,-134],[26,-25],[30,36],[33,15],[21,-57],[79,-99],[-2,-117],[23,-127],[-4,-77],[58,-146],[28,-124],[-12,-140],[37,-62],[-4,-63],[23,-52],[81,-69],[34,4],[21,-54],[40,-54],[56,-48],[28,-13],[20,-32],[-9,-63],[13,-34],[51,-86],[31,-140],[17,-109],[13,-40],[46,61],[19,-46],[38,-41],[28,-3],[2,-132],[7,-62],[81,-122],[24,-5],[26,-36],[28,-87],[40,-45],[24,-88],[32,-50],[1,-57],[30,-60],[-10,-75],[10,-149],[-6,-46],[57,-230],[3,-134],[-31,-96],[-10,-132],[-30,-148],[3,-75],[-13,-116],[-20,-77],[-28,-57],[-10,-86],[-38,-70],[-43,-32],[-41,-98],[-31,-149],[-33,-59],[-33,-193],[-16,-8],[-45,-136],[-17,-158],[-12,-64],[3,-91],[-7,-59],[-78,-60],[-127,-7],[-48,-22],[-59,-63],[-66,-97],[-79,-27],[-35,-26],[-17,52],[-49,29],[18,31],[-78,11],[26,55],[-28,42],[-43,-39],[15,-29],[-69,-55],[-36,-52],[-33,-27],[-39,17],[-75,64],[-53,12],[-57,28],[-37,-21],[-50,63],[-48,11],[-29,29],[-22,56],[-52,87],[-6,41],[15,87],[-38,124],[-63,89],[35,58],[-16,17],[-83,-64],[-41,6],[26,63],[14,66],[-3,58],[-49,130],[-21,-63],[-5,-60],[-23,-91],[-28,3],[-71,-24],[16,71],[46,1],[13,69],[0,96],[19,66],[35,62],[-10,83],[17,23],[-10,71],[-58,-93],[-25,-96],[-74,-58],[-25,-29],[-96,-200],[-25,79],[-32,148],[-37,61],[-12,66],[-22,30],[-38,5],[-25,92],[17,44],[-84,80],[-42,0],[-56,50],[-67,-11],[-61,68],[-71,44],[-44,-24],[-80,6],[-146,-28],[-108,-79],[-92,-44],[-67,-10],[-79,12],[-27,-9],[-56,-57],[-86,-72],[-48,-16],[-30,-38],[-32,-102],[-27,-53],[-54,-33],[-53,21],[-77,-22],[-11,26],[-81,10],[-73,-10],[-48,-20],[-70,-2],[-27,-28],[-22,-59],[-69,-24],[-93,-111],[-68,-25],[-131,24],[-66,43],[-32,59],[-55,48],[-33,11],[-2,163],[23,-29],[41,25],[21,74],[-8,116],[11,21],[4,150],[-5,42],[-64,195],[-23,130],[-4,171],[-13,64],[-28,60],[-11,72],[-46,101],[-8,118],[-9,43],[-39,110],[-59,133],[20,30],[23,-99],[19,-9],[15,58],[-33,51],[-24,85],[18,19],[42,-88],[8,-51],[33,-6],[0,96],[-27,67],[-40,125],[-31,118],[1,67],[16,88],[25,69],[5,86],[-14,85],[34,153],[25,-84],[26,-4],[28,88],[32,45],[73,54],[68,101],[86,83],[53,2],[33,-17],[34,17],[64,59],[68,25],[43,58],[59,-9],[23,17],[52,11],[84,54],[37,41],[39,81],[42,138],[21,19],[44,78],[-16,15],[-11,92],[2,52],[33,75],[36,40],[31,77],[21,-93],[47,-137],[17,101],[14,34],[-35,83],[29,66],[16,-42],[88,59],[-24,78],[57,131],[49,48],[-8,50],[62,72],[42,-24],[11,84],[36,21],[21,-31],[41,91],[48,-42],[46,-57],[27,-64],[38,-58],[-3,-64],[37,56],[121,-32],[36,32],[-18,50],[-28,37],[10,39],[31,51],[16,90],[36,27],[16,32],[-16,34],[4,42],[29,59],[27,9],[6,53],[44,14],[24,36],[27,-22],[55,10],[60,-2],[26,29],[11,101],[-25,35],[-42,-2],[-39,50],[18,20],[37,-5],[52,-67],[21,26],[21,-14],[20,-55],[42,-25],[46,-3],[41,-39],[37,-11],[24,16],[37,-46],[23,-6],[70,69],[27,-62],[9,-51],[22,-2],[1,66],[33,38],[22,-59],[29,-27],[-51,-95],[7,-48],[-16,-49],[-21,19],[-45,-36],[8,-111],[-13,-76],[-59,-132],[15,-54],[83,-88],[11,-37],[36,-30],[24,-41],[27,4],[37,-42],[50,-36],[25,-54],[41,-53],[47,-13],[47,-26],[30,-95],[23,-11],[74,-71],[58,17],[39,47],[17,87],[31,81],[24,127],[5,102],[21,120],[-12,129],[8,69],[-15,78],[22,119],[-4,69],[29,80],[-21,19],[15,90],[18,41],[23,136],[4,72],[35,53],[33,-67],[14,-65],[5,-113],[37,-29]],[[40984,14102],[-12,44],[7,78],[39,29],[-34,-151]],[[38971,16724],[4,-59],[18,-48],[-65,14],[7,80],[36,13]],[[38201,17190],[-32,5],[-6,40],[23,30],[15,-75]],[[38221,17257],[50,9],[30,32],[33,-54],[-30,-61],[-42,-47],[-38,40],[-27,66],[24,15]],[[39079,11910],[56,-4],[-5,-34],[-41,3],[-28,-39],[-85,9],[-27,31],[13,31],[85,35],[32,-32]],[[39996,10802],[29,3],[66,-56],[61,-30],[32,4],[51,33],[57,-2],[21,35],[43,14],[40,-37],[1,-240],[-35,-48],[-12,-68],[7,-129],[-24,-14],[-26,83],[-20,-11],[-20,-71],[-21,-13],[-25,-74],[-40,23],[-62,-8],[-9,59],[-21,7],[-39,70],[-27,84],[-3,115],[-57,132],[-16,96],[21,57],[28,-14]],[[39858,10949],[-8,83],[29,22],[6,-70],[-27,-35]],[[40360,11027],[36,-50],[3,-34],[-27,-26],[-26,54],[14,56]],[[14987,12627],[4,-107],[-31,-34],[-16,-121],[19,-116],[-16,-23],[30,-85],[64,-46],[56,-65],[18,-69],[-24,-47],[11,-93],[70,-62],[3,-94],[-51,-130],[-38,-66],[-18,-75],[-78,-76],[-102,-53],[-101,-36],[-158,-34],[-61,-1],[-57,17],[-33,-50],[35,-49],[-10,-99],[-19,-16],[-14,-66],[1,-57],[18,-48],[-18,-47],[-69,-48],[-101,-9],[-43,31],[-89,43],[-35,-15],[-3,-49],[17,-102],[-5,-88],[20,-42],[46,-30],[13,-30],[39,9],[38,61],[24,-64],[-12,-88],[-54,-12],[-23,64],[-41,9],[-46,-50],[49,-33],[24,-33],[-68,-52],[-36,-74],[5,-91],[-15,-95],[-35,-40],[1,-76],[-68,10],[-49,-48],[-43,-17],[-76,-156],[-1,-83],[22,-45],[75,-99],[96,-20],[32,-55],[-25,-106],[16,-25],[-72,-88],[-79,-62],[-53,-72],[-27,-64],[-12,-135],[-45,-51],[-97,-64],[-37,-124],[26,-109],[-1,-42],[79,-147],[-6,-11]],[[13704,8197],[43,-73],[-26,-60],[43,-13],[16,-54],[88,-107],[130,-106],[61,-26],[69,-5],[-20,-44],[-61,-7],[-83,-27],[-51,23],[-159,24],[-53,-8]],[[23651,18490],[-4,95],[-20,70]],[[23593,15968],[10,119],[-3,158],[-6,34],[27,43],[39,219],[12,131],[15,63],[6,69],[43,90],[10,56],[53,56],[33,87],[13,68],[7,166],[-15,93],[-23,47],[-40,157],[-7,69],[-19,75],[44,79],[3,70],[-36,129],[-28,120],[-5,61],[-37,80],[-27,110],[60,19],[35,31]],[[21875,27459],[38,-2],[71,59],[52,67],[60,42],[29,-9],[29,18],[13,35],[45,44],[91,56],[165,18],[46,41],[68,2],[31,22],[121,0],[54,-48],[28,0],[79,41],[52,49],[54,-37],[84,18],[37,-32],[82,18]],[[24610,28470],[-19,73],[-56,53],[-6,74],[14,68],[-2,108],[16,80],[-28,18]],[[15132,5844],[-153,-86],[-70,-93],[-7,-44],[-84,-31],[-48,22],[-86,-66],[-49,-55],[-83,-11],[-50,-66],[19,-56],[-39,-55],[57,-62],[71,31],[80,10],[-35,-56],[-120,-23],[-99,12],[-3,-97],[-68,92],[-61,-3],[-40,-83],[4,-60],[-79,27],[-24,-48],[8,-62],[-77,-5],[3,-73],[-27,-103],[79,-61],[12,-58],[191,-30],[-14,-51],[124,-119],[5,-54],[55,-64],[-51,-54],[46,-29],[5,-88],[80,7],[36,-66],[-24,-87],[22,-56],[-111,-18],[77,-46],[73,9],[-1,-78],[42,-128],[53,-5],[-178,-95],[65,-32],[-21,-112],[-65,-16],[69,-62],[-149,-6],[41,-46],[-101,-5],[-149,-59],[85,-33],[-16,-65],[-370,-123],[-141,-21],[-220,-53],[-98,-67],[-209,-22],[-256,13],[-171,25],[-237,-11],[3,-44],[113,-97],[106,-44],[312,-25],[-49,-67],[-192,-62],[-360,52],[-352,41],[-104,-29],[441,-102],[53,-19],[-54,-68],[-298,-17],[-364,102],[-63,-23],[62,-68],[196,-72],[95,-88],[91,51],[449,-12],[58,-66],[-65,-61],[105,-105],[58,-93],[271,-25],[248,66],[96,-52],[675,-146],[278,-8],[-6,-27],[-265,6],[69,-48],[-219,-7],[0,-62],[122,-41],[280,24],[-11,-58],[101,-99],[40,-80],[255,-25],[266,135],[191,80],[224,61],[410,49],[282,18],[285,-58],[-58,-65],[181,4],[169,36],[275,214],[328,91],[192,-40],[153,38],[79,53],[378,42],[360,63],[626,24],[102,45],[-228,21],[-577,36],[-45,95],[-135,4],[-194,-18],[-311,65],[-89,35],[8,67],[134,138],[66,24],[59,67],[179,62],[68,-4],[227,110],[158,63],[222,47],[132,44],[220,42],[222,23],[127,-5],[117,42],[122,-10],[120,33],[-29,32],[179,199],[107,2],[93,-18],[135,103],[-144,-1],[-23,41],[-72,25],[10,48],[103,72],[73,11],[84,-17],[40,27],[-37,59],[134,-12],[57,35],[154,47],[66,112],[-146,63],[-25,50],[84,12],[103,-48],[38,56],[78,59],[81,-32],[72,-110],[114,28],[12,93],[-31,40],[59,34],[68,-13],[111,29],[43,-31],[-61,-66],[16,-41],[185,3],[91,-10],[110,12],[99,-25],[137,18],[44,-79],[85,69],[84,43],[219,66],[110,12],[198,42],[281,35],[34,30],[104,-25],[76,63],[133,-71],[121,-48],[61,-9],[94,87],[57,35],[110,-29],[83,10],[161,-6],[96,27],[16,-43],[174,-31],[38,55],[94,-1],[-21,-83],[29,-51],[106,-3],[116,16],[33,77],[45,52],[82,-48],[-17,-37],[92,-36],[53,19],[48,70],[54,-10],[64,-97],[150,-32],[151,28],[61,30],[149,42],[134,62],[170,16],[152,49],[46,83],[-52,121],[9,45],[59,39],[132,-3],[-47,-87],[96,1],[78,-119],[184,-3],[48,-35],[55,20],[72,-23],[44,-51],[44,11],[104,128],[20,100],[156,86],[157,49],[50,50],[232,49],[35,29],[174,33],[-10,48],[69,28],[64,-30],[-35,-36],[46,-34],[58,15],[85,-29],[8,126],[-31,38],[96,23],[94,-51],[73,5],[-10,73],[-26,14],[3,68],[191,93],[220,36],[157,-13],[94,-36],[80,-64],[61,-11],[40,-42],[-90,-37],[-37,-106],[91,46],[90,9],[137,-47],[63,-56],[152,21],[101,-34],[170,-23],[124,30],[108,-24],[236,-32],[84,-1],[142,-26],[134,32],[70,-158],[-56,-60],[12,-109],[-89,-31],[35,-62],[-13,-46],[-111,6],[-112,-95],[82,-33],[162,-2],[-47,-133],[-103,-78],[-74,-134],[-67,-208],[-39,-57],[92,-20],[80,47],[0,74],[55,52],[140,30],[131,117],[81,48],[23,103],[54,98],[105,69],[8,65],[62,55],[76,23],[77,-30],[105,-3],[58,70],[160,81],[77,30],[82,97],[30,77],[69,34],[335,95],[159,19],[23,35],[122,72],[151,-10],[50,23],[78,4],[151,54],[212,-7],[104,19],[58,29],[130,21],[81,-26],[77,9],[75,-23],[107,41],[87,-28],[108,7],[83,22],[77,-24],[66,31],[74,-45],[63,5],[187,64],[72,89],[86,0],[61,19],[157,-27],[300,-93],[97,-13],[121,-37],[27,-28],[134,-31],[142,84],[-4,47],[39,52],[67,25],[203,35],[75,-28],[77,-95],[120,-44],[40,-46],[-44,-56],[-158,-40],[4,-51],[331,87],[51,-15],[124,7],[141,-26],[-19,-37],[290,61],[86,51],[220,55],[77,-22],[61,15],[37,47],[58,14],[98,-30],[26,-57],[82,-64],[110,-15],[99,17],[62,126],[46,39],[103,21],[108,-9],[128,12],[70,22],[91,-43],[13,-44],[84,32],[47,45],[104,-37],[44,-32],[41,22],[107,-15],[47,-30],[136,-6],[81,-31],[123,-9],[47,-17],[85,5],[43,-31],[70,-11],[68,27],[111,-31],[34,-27],[-84,-143],[57,0],[89,37],[104,0],[105,-75],[-3,-68],[204,-42],[321,18],[9,-83],[52,15],[121,-8],[66,-28],[77,29],[9,83],[51,-15],[97,-93],[66,-40],[119,-35],[69,1],[53,-29],[83,22],[200,-69],[43,-46],[10,-53],[54,-23],[48,-52],[47,-112],[106,-31],[-30,75],[22,64],[56,7],[83,-70],[80,-3],[50,28],[178,-30],[95,-4],[116,-32],[75,-82],[77,-20],[167,-78],[86,-53],[-78,-17],[-33,-100],[-140,-33],[124,-38],[-34,-71],[-294,-21],[56,-42],[-76,-45],[-103,2],[-7,-44],[-163,-65],[65,-124],[-139,16],[-124,-18],[-70,-48],[-8,-76],[-105,-11],[140,-138],[-47,-67],[36,-16],[4,-114],[-38,-46],[78,-17],[64,-78],[1,-43],[100,-101],[219,-90],[-11,-26],[-183,-7],[-330,-93],[-93,11],[-88,-49],[-33,-78],[43,-100],[-221,-62],[-23,-28],[242,-1],[45,-211],[235,-105],[106,-70],[-309,-49],[258,-12],[182,24],[306,-131],[169,-25],[-61,-65],[383,-25],[85,-37],[922,-124],[73,-31],[0,-1240],[-44297,0],[0,91],[0,64],[0,155],[0,155],[0,155],[0,155],[0,155],[0,155],[0,155],[198,3],[823,-45],[590,-61],[489,-18],[714,-63],[18,83],[-810,62],[-138,82],[75,46],[-107,30],[-320,-1],[-199,79],[-603,121],[276,12],[419,-84],[271,6],[382,-28],[113,-42],[338,-5],[248,43],[-31,48],[274,32],[301,117],[-210,112],[94,52],[-117,42],[-262,42],[62,35],[553,29],[481,28],[-301,120],[43,45],[234,17],[18,65],[-285,50],[-197,67],[-19,31],[-156,-4],[-209,34],[-155,72],[216,52],[-77,45],[-237,-1],[-125,54],[-27,38],[36,137],[190,-13],[56,24],[172,-5],[111,-22],[56,-41],[214,-2],[246,-82],[170,54],[-50,46],[223,18],[67,-47],[97,4],[-45,85],[26,70],[-195,68],[-155,-12],[-130,28],[164,58],[300,-71],[30,22],[-92,50],[49,46],[109,2],[188,69],[153,16],[102,-43],[223,104],[255,31],[131,-15],[40,67],[107,32],[232,-35],[175,19],[278,-28],[237,38],[309,1],[180,-13],[289,5],[230,21],[164,60],[157,-20],[241,4],[39,102],[55,14],[103,-36],[-49,-122],[-3,-74],[275,42],[-14,114],[77,18],[108,-39],[0,-75],[-143,-93],[44,-13],[226,1],[180,-29],[146,-5],[212,50],[393,-3],[82,-64],[305,50],[-210,81],[18,68],[-139,6],[-57,101],[-130,31],[44,59],[176,-31],[148,8],[71,33],[-101,37],[-144,18],[-182,0],[-50,30],[-37,80],[110,18],[9,-55],[196,1],[345,-11],[229,-61],[113,20],[122,-21],[66,22],[180,8],[163,-31],[60,17],[17,58],[117,19],[27,43],[74,-2],[40,-34],[-27,-81],[125,8],[45,-36],[130,32],[123,-64],[344,-78],[108,26],[5,93],[101,81],[114,-32],[137,-102],[191,21],[-35,-78],[104,2],[123,37],[159,-14],[132,56],[320,38],[185,37],[142,59],[45,49],[59,115],[-78,123],[-6,103],[-47,140],[-63,88],[-9,83],[-29,46],[70,19],[95,-14],[48,55],[-51,65],[74,137],[-44,47],[55,109],[-109,21],[12,70],[66,28],[42,-49],[30,132],[74,-4],[18,108],[119,59],[72,68],[4,83],[34,32],[79,20],[67,41],[116,33],[68,63],[26,53],[73,33],[94,23],[78,59],[81,22],[102,49],[45,-32]],[[16584,2456],[185,-2],[93,-145],[-3,-71],[-66,-85],[-697,-99],[-72,-31],[-540,-19],[-23,66],[107,83],[145,11],[241,145],[-48,43],[63,149],[140,122],[171,46],[107,12],[159,-21],[147,-47],[65,-42],[-17,-69],[-138,-14],[-19,-32]],[[14798,2119],[-11,-85],[-93,-47],[-290,40],[-164,3],[-210,80],[426,3],[148,14],[-47,38],[43,54],[94,33],[108,-30],[-4,-103]],[[13529,4568],[78,-87],[40,-111],[93,-168],[12,-212],[-49,-85],[-71,-72],[-388,-29],[-89,57],[147,63],[97,16],[35,31],[-122,19],[-74,35],[-63,-57],[-120,-51],[-109,20],[-69,40],[7,58],[136,51],[193,-1],[-54,57],[86,12],[79,-16],[64,35],[114,5],[-4,28],[-92,16],[-1,50],[90,38],[0,45],[-98,-9],[-90,50],[14,120],[-49,68],[33,38],[179,33],[38,-37],[7,-50]],[[10079,3970],[76,-19],[74,32],[92,-9],[63,-38],[4,-75],[-58,-37],[-148,12],[-142,-6],[-220,61],[-255,30],[17,33],[251,32],[31,-46],[72,17],[78,45],[65,-32]],[[7314,3565],[35,-51],[-88,-37],[-9,-39],[-151,-15],[-70,13],[29,83],[-60,62],[314,-16]],[[6604,3668],[73,-26],[91,-93],[107,-1],[36,-86],[-128,4],[-196,87],[-3,24],[-113,41],[18,46],[115,4]],[[2577,2237],[-153,-8],[-193,31],[-290,84],[12,71],[107,60],[199,-39],[241,-110],[77,-89]],[[42708,2783],[168,-15],[111,-30],[-127,-35],[-114,8],[-82,-45],[-61,67],[105,50]],[[17973,2344],[-313,6],[116,44],[197,-50]],[[13494,2265],[-147,12],[-31,40],[65,84],[87,39],[377,114],[93,-5],[-311,-181],[-59,-90],[-74,-13]],[[15031,5695],[73,-10],[-10,-31],[-120,-17],[20,48],[37,10]],[[15316,5888],[-105,-58],[-1,45],[106,13]],[[14374,5603],[19,-31],[-76,-42],[-65,49],[103,69],[19,-45]],[[13783,4944],[-73,-57],[-64,29],[2,44],[92,131],[80,-75],[-60,-44],[23,-28]],[[13079,4250],[-61,-63],[-243,-46],[-40,42],[159,37],[25,36],[119,12],[41,-18]],[[12922,4449],[67,-48],[-50,-51],[-108,18],[5,48],[86,33]],[[12999,3709],[-38,-38],[11,-74],[-163,61],[-14,54],[71,20],[8,38],[111,-16],[14,-45]],[[10931,3691],[-62,123],[106,3],[4,-84],[-48,-42]],[[34574,5338],[-53,-5],[-17,52],[58,19],[41,-21],[-29,-45]],[[25453,4306],[-105,2],[-2,38],[86,18],[21,-58]],[[14006,2178],[-90,37],[-9,74],[119,-12],[62,-72],[-82,-27]],[[14520,4451],[-61,95],[83,-20],[-22,-75]],[[7438,3487],[-31,65],[104,1],[-73,-66]],[[2403,1846],[-405,32],[161,28],[244,-60]],[[24145,4388],[-75,15],[10,36],[68,15],[-3,-66]],[[21837,4221],[54,-12],[-34,-49],[-53,14],[33,47]],[[19613,3536],[22,-94],[-51,-26],[-55,98],[84,22]],[[6021,3383],[-95,10],[9,48],[97,-20],[-11,-38]],[[3786,2827],[-151,10],[72,44],[79,-54]],[[14395,23720],[-14,0]]],"transform":{"scale":[0.0081269611937603,0.004556033256075637],"translate":[-180,-89.99892578125002]},"objects":{"ne_50m_admin_0_countries_lakes":{"type":"GeometryCollection","geometries":[{"arcs":[[0,1,2,3,4,5,6,7]],"type":"Polygon"},{"arcs":[[-7,8,-3,9,10,11,12,13,14]],"type":"Polygon"},{"arcs":[[[16,17,18]],[[19]]],"type":"MultiPolygon"},{"arcs":[[20,21,22,23]],"type":"Polygon"},{"arcs":[[24,25,26,27]],"type":"Polygon"},{"type":null},{"arcs":[[[28]],[[29]],[[30]],[[31]],[[32]]],"type":"MultiPolygon"},{"arcs":[[33,34,35,36,37,38,39,40,41,42,43]],"type":"Polygon"},{"arcs":[[44,45,46]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[47]],"type":"Polygon"},{"arcs":[[[48]],[[49]],[[50]],[[51,52,53,54]],[[55]],[[56]],[[57]],[[58]],[[59]],[[60]],[[61]],[[62]],[[63]],[[64]],[[65]],[[66]],[[67]],[[68]],[[69]],[[70]],[[71]],[[72,73,74,75,76,77,78,79,80,81,82]],[[83]],[[84]]],"type":"MultiPolygon"},{"arcs":[[85]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[[86]],[[87]]],"type":"MultiPolygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[[88]],[[89]],[[90]],[[91]],[[92,93]]],"type":"MultiPolygon"},{"arcs":[[94,95,96,97,98]],"type":"Polygon"},{"arcs":[[99,100,101,102,103,104,105,106,107,108,109]],"type":"Polygon"},{"arcs":[[110,111,112,113,114,115,116,117,118,119]],"type":"Polygon"},{"arcs":[[-39,120,-37,121,122,123,124]],"type":"Polygon"},{"arcs":[[[125,126,127,128,129,130,131]],[[132,133,134]]],"type":"MultiPolygon"},{"arcs":[[135,136,137]],"type":"Polygon"},{"arcs":[[138]],"type":"Polygon"},{"type":null},{"arcs":[[139,140,141,142]],"type":"Polygon"},{"arcs":[[[143,144]],[[145,146]]],"type":"MultiPolygon"},{"arcs":[[147,148,149,150,151,152]],"type":"Polygon"},{"arcs":[[[153]],[[-13,154,155,156,-111,157,158,159,160,161,162,163]]],"type":"MultiPolygon"},{"arcs":[[-35,164,165,166]],"type":"Polygon"},{"arcs":[[167]],"type":"Polygon"},{"arcs":[[-131,168,169,170,171,172]],"type":"Polygon"},{"arcs":[[173,174,175,176,177,178,179]],"type":"Polygon"},{"arcs":[[[180]],[[181,182,183]]],"type":"MultiPolygon"},{"arcs":[[184,185]],"type":"Polygon"},{"arcs":[[186,187,188,189]],"type":"Polygon"},{"arcs":[[-118,190,191,192,193,194]],"type":"Polygon"},{"arcs":[[-193,195,196,197,198,199,200,201]],"type":"Polygon"},{"arcs":[[202]],"type":"Polygon"},{"arcs":[[[203]],[[204,205,206,207,208,209]],[[210]],[[211]],[[212]],[[213]]],"type":"MultiPolygon"},{"arcs":[[[214,215]],[[216]]],"type":"MultiPolygon"},{"arcs":[[-1,217,-186,218,219,220,221],[222]],"type":"Polygon"},{"arcs":[[223,224,225,226]],"type":"Polygon"},{"arcs":[[-226,227,228,229]],"type":"Polygon"},{"arcs":[[[230]],[[231]],[[232]],[[233]],[[234]],[[235]],[[236]]],"type":"MultiPolygon"},{"arcs":[[-107,237,238,239,240]],"type":"Polygon"},{"arcs":[[241,242,243,244,245]],"type":"Polygon"},{"type":null},{"arcs":[[246,247,248]],"type":"Polygon"},{"type":null},{"arcs":[[249,250,251,252,253,254,255,256]],"type":"Polygon"},{"arcs":[[257,258,259,260,261,262,263]],"type":"Polygon"},{"arcs":[[-18,264,265,266,267,268,269,270,-97,271]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[272]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"arcs":[[-112,-157,273,274,275,276]],"type":"Polygon"},{"arcs":[[[277]],[[278,279,280,281,282,283,284,285,286,287,288,-110,289,290,291,292,293,294,295,296,297]],[[298]],[[299]],[[300]],[[301]],[[302]],[[303]],[[304]],[[305]],[[306]],[[307]],[[308]],[[309]],[[310]],[[311]],[[312]],[[313]],[[314]],[[315]],[[316]],[[317]],[[318]],[[319]],[[320]],[[321]],[[322]],[[323]],[[324]],[[325]],[[326]],[[327]],[[328]],[[329]],[[330]],[[331]],[[332,333,334,335,336]],[[-101,337]]],"type":"MultiPolygon"},{"arcs":[[-103,338,339,-257,340,-105,341]],"type":"Polygon"},{"arcs":[[-270,342]],"type":"Polygon"},{"arcs":[[[343]],[[-209,344]]],"type":"MultiPolygon"},{"arcs":[[-108,-241,345,346,347,348,349,-335,350,351]],"type":"Polygon"},{"arcs":[[[352]],[[353]],[[354]],[[355]],[[356]],[[357]],[[358]],[[359]],[[360]],[[361]],[[362]],[[363]],[[364]],[[365]]],"type":"MultiPolygon"},{"arcs":[[366,367,368,369,370,371,372,373,374,375]],"type":"Polygon"},{"arcs":[[376,377,378]],"type":"Polygon"},{"arcs":[[[379]],[[380]],[[381,382]],[[383]],[[384]],[[385]],[[386]]],"type":"MultiPolygon"},{"arcs":[[387,388,389,390]],"type":"Polygon"},{"type":null},{"arcs":[[391,392,393,394,395,396]],"type":"Polygon"},{"arcs":[[[-19,-272,-96,397]],[[-99,398]]],"type":"MultiPolygon"},{"arcs":[[[-182,399,-297,400]],[[401]],[[402]],[[403]],[[404]],[[405]],[[406]],[[407]]],"type":"MultiPolygon"},{"arcs":[[-215,408,409,-279,410]],"type":"Polygon"},{"arcs":[[411,412,413,414,415]],"type":"Polygon"},{"arcs":[[-416,416,417,418,419,420,421]],"type":"Polygon"},{"arcs":[[422,423,424,425]],"type":"Polygon"},{"arcs":[[[426]],[[427]],[[428]]],"type":"MultiPolygon"},{"type":null},{"type":null},{"arcs":[[[429,430,431],[432]],[[433,434]],[[435]]],"type":"MultiPolygon"},{"type":null},{"type":null},{"arcs":[[436,437]],"type":"Polygon"},{"type":null},{"arcs":[[-10,438,-221,439,440]],"type":"Polygon"},{"arcs":[[-8,-15,441,442,-161,443,-219,-185,-218]],"type":"Polygon"},{"arcs":[[444,445,446]],"type":"Polygon"},{"arcs":[[-446,447,448,449]],"type":"Polygon"},{"arcs":[[-253,450,451,452,453,454]],"type":"Polygon"},{"arcs":[[-283,455]],"type":"Polygon"},{"arcs":[[-104,-342]],"type":"Polygon"},{"arcs":[[456,457]],"type":"Polygon"},{"arcs":[[[-75,458,459,460,461]],[[462]]],"type":"MultiPolygon"},{"arcs":[[463]],"type":"Polygon"},{"arcs":[[-264,464,-449,465,466]],"type":"Polygon"},{"type":null},{"arcs":[[-258,-467,467,-419,468,469,470]],"type":"Polygon"},{"type":null},{"arcs":[[[-151,471]],[[472,473,474,475,476]]],"type":"MultiPolygon"},{"arcs":[[-14,-164,479,-442]],"type":"Polygon"},{"arcs":[[481]],"type":"Polygon"},{"arcs":[[-251,482,483,484,485]],"type":"Polygon"},{"arcs":[[486,487,488]],"type":"Polygon"},{"arcs":[[-334,490,491,492,-351]],"type":"Polygon"},{"arcs":[[-175,493]],"type":"Polygon"},{"arcs":[[-136,494,495,-198,496,-421,497]],"type":"Polygon"},{"arcs":[[-247,498,499,500]],"type":"Polygon"},{"arcs":[[-223]],"type":"Polygon"},{"arcs":[[-172,501,502]],"type":"Polygon"},{"arcs":[[-492,503,504,-291,505]],"type":"Polygon"},{"arcs":[[-23,506,-148,507,508]],"type":"Polygon"},{"arcs":[[-34,509,510,-165]],"type":"Polygon"},{"arcs":[[-268,511,512]],"type":"Polygon"},{"arcs":[[-252,-486,513,-451]],"type":"Polygon"},{"type":null},{"arcs":[[-119,-195,514,-224,515,-159,516]],"type":"Polygon"},{"arcs":[[[-44,517,-40,-125,518,-285,519,-510]],[[-42,520]]],"type":"MultiPolygon"},{"arcs":[[-170,521,-266,522,523,524,525]],"type":"Polygon"},{"arcs":[[[526]],[[527]],[[528]],[[529]]],"type":"MultiPolygon"},{"arcs":[[530]],"type":"Polygon"},{"arcs":[[[-177,531,-244,532,533]],[[534]],[[535]]],"type":"MultiPolygon"},{"arcs":[[-171,-526,536,-524,537,538,539,540,-502]],"type":"Polygon"},{"arcs":[[[-540,541,542]],[[-525,-537]]],"type":"MultiPolygon"},{"arcs":[[-93,543]],"type":"Polygon"},{"arcs":[[-130,544,545,-512,-267,-522,-169]],"type":"Polygon"},{"arcs":[[[546]],[[-123,547,-395,548,-545,-129,549,550,551,552]]],"type":"MultiPolygon"},{"arcs":[[[553]],[[554]],[[555]],[[556]],[[557]],[[558]],[[559]],[[560]],[[561]],[[562]],[[563]],[[564]],[[565]],[[566]],[[567]],[[568]],[[569]],[[570]],[[571]],[[572]],[[573]],[[574]],[[575]],[[576]],[[577]],[[578]],[[-473,579]],[[580]],[[581]],[[582]],[[-144,583,-147,584]],[[585]],[[586]],[[587]],[[-382,588]],[[589]],[[590]],[[591]],[[592]],[[593]],[[594]],[[595]],[[596]],[[597]],[[598]],[[599]],[[600]],[[601]],[[602]]],"type":"MultiPolygon"},{"arcs":[[[-393,604,605,-437,606,607,608,609,610,611]],[[612]]],"type":"MultiPolygon"},{"arcs":[[613]],"type":"Polygon"},{"arcs":[[-106,-341,-256,614,-246,615,-238]],"type":"Polygon"},{"arcs":[[-426,616,617,618,619]],"type":"Polygon"},{"arcs":[[620,621]],"type":"Polygon"},{"arcs":[[-25,622,-189,623]],"type":"Polygon"},{"arcs":[[-260,624,625]],"type":"Polygon"},{"arcs":[[-249,626,-625,-259,-471,627,-499]],"type":"Polygon"},{"arcs":[[-461,628,629,-619,630,631]],"type":"Polygon"},{"type":null},{"arcs":[[[632]],[[633]],[[634]],[[635]],[[636]],[[637]],[[-134,638,639,-484,640]]],"type":"MultiPolygon"},{"arcs":[[-142,641,642,643,644,645]],"type":"Polygon"},{"arcs":[[[-179,646,-487,647,-432,648,649,650,-347,651,652,653]],[[654]]],"type":"MultiPolygon"},{"arcs":[[-126,656,-288,657,658]],"type":"Polygon"},{"arcs":[[-262,659]],"type":"Polygon"},{"arcs":[[660,661,662,663]],"type":"Polygon"},{"arcs":[[[664]],[[-178,-534,665,-458,666,-207,667,-205,668,669,-488,-647]],[[670]],[[671]],[[672]],[[-187,673,674]]],"type":"MultiPolygon"},{"type":null},{"type":null},{"arcs":[[675,676]],"type":"Polygon"},{"type":null},{"arcs":[[677]],"type":"Polygon"},{"arcs":[[[678]],[[679]]],"type":"MultiPolygon"},{"arcs":[[680]],"type":"Polygon"},{"type":null},{"arcs":[[-183,-401,-296,681]],"type":"Polygon"},{"arcs":[[[682]],[[683]]],"type":"MultiPolygon"},{"arcs":[[-194,-202,684,685,-228,-225,-515]],"type":"Polygon"},{"arcs":[[[-505,686,-294,687,-292]],[[688]],[[689]]],"type":"MultiPolygon"},{"arcs":[[-201,690,691,-685]],"type":"Polygon"},{"arcs":[[[692]],[[-663,693,694]]],"type":"MultiPolygon"},{"arcs":[[-618,695,-631]],"type":"Polygon"},{"arcs":[[-199,-496,696,-542,-539,697]],"type":"Polygon"},{"arcs":[[[-375,698,699]],[[700]],[[701]]],"type":"MultiPolygon"},{"arcs":[[-621,702]],"type":"Polygon"},{"arcs":[[703]],"type":"Polygon"},{"arcs":[[-229,-686,-692,704]],"type":"Polygon"},{"arcs":[[[705]],[[706]],[[707]],[[708]],[[709]],[[710]]],"type":"MultiPolygon"},{"type":null},{"arcs":[[[711]],[[-650,712]],[[713]],[[714]]],"type":"MultiPolygon"},{"arcs":[[-240,715,-652,-346]],"type":"Polygon"},{"arcs":[[716,717]],"type":"Polygon"},{"arcs":[[-717,718]],"type":"Polygon"},{"arcs":[[[719]],[[720]]],"type":"MultiPolygon"},{"arcs":[[[-242,-615,-255,721,722]],[[-454,723,724]]],"type":"MultiPolygon"},{"arcs":[[-470,726,-645,727,-500,-628]],"type":"Polygon"},{"arcs":[[-390,728,-424,729]],"type":"Polygon"},{"arcs":[[-12,730,731,732,733,734,-191,-117,735,-115,736,-113,-277,737,-275,738,-155]],"type":"Polygon"},{"arcs":[[-661,739,740,-734,741,742]],"type":"Polygon"},{"arcs":[[743]],"type":"Polygon"},{"arcs":[[-27,744,-376,-700,745,-388,746]],"type":"Polygon"},{"arcs":[[[-24,-509,747,-609,748,-607,-438,-606,749,-397,750,-166,-511,-520,-284,-456,-282,751,-280,-410,752,753,754,755,756]],[[757]]],"type":"MultiPolygon"},{"arcs":[[-756,758]],"type":"Polygon"},{"arcs":[[-754,759]],"type":"Polygon"},{"arcs":[[[-373,760,761,762]],[[763,764]],[[765]],[[766]],[[767]],[[768]],[[769]],[[770]],[[771]],[[772]],[[773]],[[774]],[[775]],[[776]],[[777]],[[778]]],"type":"MultiPolygon"},{"arcs":[[-197,779,780,-412,-422,-497]],"type":"Polygon"},{"arcs":[[-192,-735,-741,781,-780,-196]],"type":"Polygon"},{"arcs":[[782]],"type":"Polygon"},{"arcs":[[[783]],[[784]],[[785]],[[-55,786,-83,787,-81,788,-79,789,-77,790]],[[791]],[[792]],[[793]],[[794]],[[795]],[[796]],[[797]],[[798]],[[799]],[[800]],[[801]],[[802]],[[803]],[[804]],[[805]],[[806]],[[807]],[[808]],[[809]],[[-73,810]],[[811]],[[812]],[[813]],[[814]],[[815]],[[816]],[[817]],[[818]],[[819]],[[820]],[[821]],[[822]],[[823]],[[824]],[[825]],[[826]],[[828]],[[829]],[[830]],[[831]],[[832]],[[833]],[[834]],[[835]],[[836]]],"type":"MultiPolygon"},{"arcs":[[-413,-781,-782,-740,-664,-695,837]],"type":"Polygon"},{"arcs":[[-22,838,-149,-507]],"type":"Polygon"},{"arcs":[[-153,839,840,-610,-748,-508]],"type":"Polygon"},{"arcs":[[-156,-739,-274]],"type":"Polygon"},{"arcs":[[-143,-646,-727,-469,-418,841]],"type":"Polygon"},{"arcs":[[-135,-641,-483,-250,-340,842]],"type":"Polygon"},{"arcs":[[[-476,843]],[[-475,844]]],"type":"MultiPolygon"},{"arcs":[[[-26,-624,-188,-675,845,-47,846,-377,847,-367,-745]],[[848]],[[849]],[[850]],[[851]],[[852]]],"type":"MultiPolygon"},{"arcs":[[-2,-222,-439]],"type":"Polygon"},{"arcs":[[-254,-455,-725,853,-722]],"type":"Polygon"},{"arcs":[[-372,854,-368,-848,-379,855,-761]],"type":"Polygon"},{"arcs":[[-608,-749]],"type":"Polygon"},{"arcs":[[-140,-842,-417,-415,857]],"type":"Polygon"},{"arcs":[[-460,858,-629]],"type":"Polygon"},{"arcs":[[-430,-648,-489,-670,859,-434]],"type":"Polygon"},{"arcs":[[-109,-352,-493,-506,-290]],"type":"Polygon"},{"type":null},{"arcs":[[-611,-841,860]],"type":"Polygon"},{"type":null},{"arcs":[[[861]],[[862]],[[863]]],"type":"MultiPolygon"},{"arcs":[[[-128,864,-550]],[[-552,865,-658,-287,866]]],"type":"MultiPolygon"},{"arcs":[[-174,867,-653,-716,-239,-616,-245,-532,-176,-494]],"type":"Polygon"},{"arcs":[[[868]],[[869]],[[870]],[[871]],[[872]],[[873]],[[874]],[[875]],[[876]]],"type":"MultiPolygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[-127,-659,-866,-551,-865]],"type":"Polygon"},{"arcs":[[[-46,877,-762,-856,-378,-847]],[[-764,878]]],"type":"MultiPolygon"},{"type":null},{"arcs":[[[-733,879,-742]],[[-11,-441,880,-731]]],"type":"MultiPolygon"},{"arcs":[[-206,-668]],"type":"Polygon"},{"arcs":[[-137,-498,-420,-468,-466,-448,-445,881]],"type":"Polygon"},{"arcs":[[-452,-514,-485,-640,882]],"type":"Polygon"},{"arcs":[[-36,-167,-751,-396,-548,-122]],"type":"Polygon"},{"arcs":[[-392,-750,-605]],"type":"Polygon"},{"arcs":[[[883]],[[884]],[[885]],[[886]],[[887]],[[888]],[[889]],[[890]],[[891]],[[892]],[[893]],[[894]],[[895]],[[896]],[[897]],[[898]],[[899]],[[900]],[[901]],[[902]],[[903]],[[904]],[[905]],[[906]],[[907]],[[908]],[[909]],[[910]],[[911]],[[912]]],"type":"MultiPolygon"},{"arcs":[[-676,913]],"type":"Polygon"}]}}}
},{}],16:[function(_dereq_,module,exports){
module.exports={
  "type": "Topology",
  "arcs": [
    [[9011, 1404], [28, -17]],
    [
      [9039, 1387],
      [-28, -96],
      [-13, -150],
      [10, -78],
      [46, -32],
      [125, -70],
      [103, -14],
      [144, 61],
      [45, 2],
      [41, -19],
      [37, -40],
      [56, -107],
      [15, 19],
      [6, 96],
      [29, 83],
      [49, 72],
      [101, 87],
      [151, 102],
      [116, 55],
      [107, 18],
      [8, 17],
      [55, 56],
      [9, -11],
      [-21, -54],
      [0, -27],
      [20, 3],
      [26, 58],
      [32, 115],
      [79, 124],
      [235, 228],
      [91, 62],
      [152, 37],
      [213, 14],
      [213, -30],
      [210, -75],
      [144, -35],
      [77, 6],
      [303, -49],
      [31, 12],
      [-42, 27],
      [-201, 29],
      [-62, 29],
      [-79, -6],
      [-5, 15],
      [56, 98],
      [31, 12],
      [38, 38],
      [41, 49],
      [65, 57],
      [121, 44],
      [177, 32],
      [176, 57],
      [48, 5],
      [37, -18],
      [68, 6],
      [208, 42],
      [122, -1],
      [88, -12],
      [55, -24],
      [65, 8],
      [120, 71],
      [14, 26],
      [-24, 44],
      [-64, 62],
      [-29, 48],
      [5, 36],
      [-8, 31],
      [-31, 47],
      [7, 20],
      [159, -56],
      [13, -20],
      [-5, -38],
      [-25, -57],
      [5, -78],
      [33, -99],
      [14, -79],
      [-5, -58],
      [-41, -55],
      [-75, -53],
      [-69, -80],
      [-62, -109],
      [-85, -88],
      [-106, -69],
      [-101, -87],
      [-96, -104],
      [-192, -144],
      [-373, -243],
      [-48, -28],
      [-44, -31],
      [-84, -53],
      [-82, -54],
      [-535, -233],
      [-360, -195],
      [-219, -156],
      [-190, -194],
      [-86, -60],
      [-95, -16],
      [-149, 11],
      [-155, -23],
      [-161, -57],
      [-105, -47],
      [-47, -35],
      [-46, -11],
      [-44, 14],
      [-119, 95],
      [-8, -7],
      [-22, -14],
      [-70, 16],
      [-79, 2],
      [-69, -35],
      [-69, -5],
      [-68, 23],
      [57, 34],
      [271, 72],
      [-3, 13],
      [-98, 41],
      [-30, 38],
      [-25, -3],
      [-19, -44],
      [-36, -19],
      [-52, 6],
      [-53, 33],
      [-53, 60],
      [-174, 97],
      [-48, 49],
      [-62, 33],
      [-91, 17],
      [-13, 28],
      [15, 38],
      [79, 179],
      [46, 69],
      [19, 19],
      [25, -1],
      [11, 18],
      [-4, 39],
      [10, 23],
      [51, 66],
      [28, 35],
      [-12, 54],
      [13, 99],
      [33, 129],
      [38, 93]
    ],
    [
      [13129, 2518],
      [30, -5],
      [31, 54],
      [16, 31],
      [-3, 18],
      [-71, 16],
      [-27, -19],
      [-5, -43],
      [8, -34],
      [21, -18]
    ],
    [[9478, 551], [11, 39], [-16, 70], [-11, 19], [-20, -34], [-11, -78], [1, -33], [46, 17]],
    [
      [16336, 4877],
      [11, -40],
      [-15, -46],
      [-42, -51],
      [-406, -339],
      [-114, -121],
      [43, -69],
      [3, -31],
      [26, 5],
      [83, 59],
      [24, -6],
      [12, -16],
      [7, -38],
      [-38, -67],
      [17, -15],
      [73, 15],
      [32, -7],
      [-28, -50],
      [-38, -32],
      [-6, -32],
      [-56, -82],
      [-24, -17],
      [-14, 11],
      [4, 37],
      [-52, -40],
      [-4, -32],
      [30, -31],
      [23, -48],
      [17, -65],
      [9, -138],
      [10, -27],
      [4, -13],
      [-5, -16],
      [-5, -33],
      [12, -29],
      [-9, -33],
      [-30, -39],
      [-57, -20],
      [-87, -4],
      [-97, -44],
      [-107, -84],
      [-73, -77],
      [-40, -70],
      [-30, -29],
      [-19, 12],
      [-144, -57],
      [-51, -39],
      [-13, -34],
      [-14, -7],
      [-13, 22],
      [-103, 16],
      [-192, 11],
      [-150, -11],
      [-108, -32],
      [-96, 19],
      [-83, 69],
      [-131, 54],
      [-337, 59],
      [-282, -11],
      [-233, -41],
      [-364, -121],
      [-160, -68],
      [-78, -55],
      [-46, -15],
      [-64, 17],
      [-127, 2],
      [-110, 40],
      [-90, 28],
      [-45, 76],
      [-1, 30],
      [118, 210],
      [46, 96],
      [14, 61],
      [27, 41],
      [70, 54],
      [19, 38],
      [34, 17],
      [50, -2],
      [94, 70],
      [137, 142],
      [105, 87],
      [72, 33],
      [196, 32],
      [124, 42],
      [88, 15],
      [53, -11],
      [82, 17],
      [111, 45],
      [149, 33],
      [187, 22],
      [119, 23],
      [52, 27],
      [102, -2],
      [18, 38],
      [27, 18],
      [99, -10],
      [34, -22],
      [5, -16],
      [-24, -11],
      [14, -27],
      [51, -44],
      [82, -14],
      [97, 28],
      [65, -12],
      [-16, -19],
      [-35, 0],
      [-5, -9],
      [26, -37],
      [22, -2],
      [77, 41],
      [5, -16],
      [-34, -60],
      [11, -37],
      [57, -15],
      [97, 26],
      [120, 74],
      [56, 21],
      [8, 35],
      [-39, 0],
      [-104, -14],
      [-29, -21],
      [-1, 14],
      [19, 56],
      [94, 81],
      [42, 46],
      [26, 51],
      [-27, -1],
      [-80, -52],
      [-71, -28],
      [-62, -4],
      [-15, 53],
      [34, 108],
      [-32, 39],
      [-156, -54],
      [-19, -19],
      [0, -14],
      [17, -11],
      [-2, -10],
      [-21, -10],
      [-35, 10],
      [0, 23],
      [17, 37],
      [-19, 17],
      [-48, -12],
      [-65, -23],
      [-140, -74],
      [2, 16],
      [37, 46],
      [149, 59],
      [73, 27],
      [46, -4],
      [28, 21],
      [26, 0],
      [22, -21],
      [56, 9],
      [89, 38],
      [38, -13],
      [-13, -64],
      [23, -5],
      [59, 52],
      [39, 18],
      [17, -16],
      [-125, -86],
      [-35, -42],
      [-7, -24],
      [21, -32],
      [50, 5],
      [205, 171],
      [99, 56],
      [99, 33],
      [99, 8],
      [82, 30],
      [67, 51],
      [284, 113],
      [116, 70],
      [110, 125],
      [27, 8],
      [18, -17]
    ],
    [
      [15434, 4329],
      [48, 70],
      [-3, 28],
      [-118, -43],
      [-39, -27],
      [-10, -24],
      [23, -15],
      [56, -6],
      [43, 17]
    ],
    [[15847, 4126], [-1, -34], [14, 4], [31, 44], [1, 17], [-30, -10], [-15, -21]],
    [
      [15885, 4480],
      [10, 34],
      [-52, -12],
      [-50, -12],
      [-101, -20],
      [-21, -38],
      [-32, -38],
      [-9, -26],
      [15, -15],
      [3, -17],
      [-11, -22],
      [16, -5],
      [46, -5],
      [58, 55],
      [12, 30],
      [24, 51],
      [92, 40]
    ],
    [[16042, 4578], [5, 17], [-44, 0], [-26, -21], [-5, -18], [8, -19], [31, 13], [31, 28]],
    [
      [9011, 1404],
      [41, 57],
      [47, 39],
      [91, 45],
      [24, 29],
      [18, 76],
      [14, 121],
      [29, 78],
      [45, 33],
      [12, 38],
      [-20, 41],
      [17, 40],
      [54, 40],
      [46, 12],
      [37, -17],
      [17, -40],
      [3, -105],
      [-9, -32],
      [-6, -26],
      [11, -25],
      [21, -28],
      [132, -39],
      [28, -7],
      [21, 4],
      [14, -40],
      [5, -97],
      [-15, -70],
      [-34, -44],
      [-53, -27],
      [-74, -12],
      [-278, 33],
      [-104, -13],
      [-43, -24],
      [-33, -23],
      [-30, -34]
    ],
    [
      [7351, 6945],
      [-9, 3],
      [-34, -31],
      [-60, -10],
      [-86, 11],
      [-30, -13],
      [26, -38],
      [-22, -58],
      [-70, -79],
      [-40, -81],
      [-11, -83],
      [11, -73],
      [33, -63],
      [48, -42],
      [64, -18],
      [13, -26],
      [-37, -31],
      [-66, -18],
      [-95, -5],
      [-95, -35],
      [-97, -66],
      [-56, -59],
      [-16, -52],
      [-3, -76],
      [13, -186],
      [-8, -72],
      [-41, -132],
      [-54, -172],
      [-65, -108],
      [-12, 15],
      [10, 52],
      [54, 173],
      [8, 73],
      [-5, 39],
      [-18, 3],
      [-22, -33],
      [-23, -71],
      [-6, -38],
      [14, -6],
      [-10, -48],
      [-34, -90],
      [-29, -19],
      [-25, 54],
      [2, 69],
      [41, 140],
      [-15, 27],
      [24, 87],
      [10, 63],
      [-12, 65],
      [10, 60],
      [31, 55],
      [-10, 17],
      [-50, -22],
      [-73, -102],
      [-95, -183],
      [-74, -80],
      [-53, 23],
      [-42, -15],
      [-32, -54],
      [-35, -27],
      [-38, 1],
      [-23, -64],
      [-8, -129],
      [-37, -73],
      [-111, -44],
      [-22, -33],
      [-2, -48],
      [18, -62],
      [0, -124],
      [-16, -186],
      [-37, -163],
      [-57, -139],
      [-63, -115],
      [-104, -152],
      [-2, -30],
      [52, -118],
      [21, -88],
      [10, -176],
      [-2, -29],
      [-100, -217],
      [-2, -31],
      [107, -364],
      [95, -261],
      [11, -28],
      [29, -105],
      [52, -204],
      [27, -203],
      [2, -201],
      [-13, -195],
      [-27, -187],
      [-40, -180],
      [-53, -171],
      [-55, -135],
      [-113, -224],
      [-55, -149],
      [-57, -121],
      [-59, -92],
      [-107, -102],
      [-275, -188],
      [-86, -37],
      [-95, -13],
      [-103, 8],
      [-56, 21],
      [-12, 35],
      [-23, 4],
      [-29, 23],
      [-38, 40],
      [-42, 71],
      [-46, 148],
      [-49, 223],
      [-43, 133],
      [-38, 42],
      [-35, 67],
      [-30, 91],
      [-15, 205],
      [0, 319],
      [10, 217],
      [21, 116],
      [-6, 88],
      [-32, 61],
      [-19, 76],
      [-7, 90],
      [-17, 62],
      [-26, 32],
      [-5, 46],
      [14, 58],
      [0, 40],
      [-14, 21],
      [-3, 33],
      [7, 43],
      [-13, 127],
      [4, 77],
      [17, 84],
      [57, 152],
      [19, 74],
      [4, 76],
      [27, 96],
      [75, 197],
      [1, 46],
      [-25, 89],
      [-6, 67],
      [5, 88],
      [24, 110],
      [42, 130],
      [55, 104],
      [69, 78],
      [22, 89],
      [-24, 99],
      [9, 146],
      [41, 195],
      [46, 145],
      [50, 97],
      [35, 92],
      [18, 88],
      [-8, 39],
      [-38, 42],
      [-36, 61],
      [-61, 13],
      [-75, -40],
      [-11, 3],
      [-38, -6],
      [-72, -89],
      [-65, -155],
      [-46, -64],
      [-64, -44],
      [-47, -48],
      [-31, -54],
      [-40, -12],
      [-50, 31],
      [-9, 73],
      [30, 118],
      [46, 113],
      [95, 184],
      [4, 44],
      [23, 54],
      [147, 80],
      [56, 62],
      [20, 80],
      [5, 117],
      [10, 33],
      [21, 14],
      [65, 99],
      [108, 183],
      [97, 188],
      [85, 195],
      [63, 110],
      [40, 25],
      [29, 59],
      [46, 181],
      [34, 83],
      [23, 32],
      [11, -19],
      [-22, -100],
      [0, -83],
      [19, -105],
      [36, -36],
      [54, 33],
      [47, 69],
      [39, 106],
      [45, 44],
      [50, -19],
      [43, 7],
      [33, 32],
      [64, 16],
      [12, -40],
      [-8, -80],
      [-20, -39],
      [-33, 3],
      [-37, -43],
      [-42, -87],
      [2, -70],
      [45, -52],
      [22, -3],
      [0, 47],
      [21, 47],
      [67, 85],
      [7, 28],
      [128, 75],
      [48, 48],
      [13, 54],
      [67, 113],
      [63, 35],
      [86, 16],
      [83, -10],
      [81, -35],
      [45, 1],
      [11, 36],
      [69, 25],
      [126, 12],
      [93, 52],
      [60, 91],
      [77, 47],
      [93, 1],
      [126, -32],
      [159, -66],
      [129, -91],
      [170, -181],
      [45, -18]
    ],
    [[7357, 7056], [-6, -111]],
    [
      [6621, 6894],
      [-46, 24],
      [-35, -60],
      [-34, -129],
      [-5, -74],
      [26, -18],
      [35, 8],
      [42, 34],
      [17, 61],
      [-9, 138],
      [9, 16]
    ],
    [
      [6060, 5939],
      [-7, -15],
      [22, -69],
      [22, -22],
      [28, -2],
      [12, 23],
      [-4, 46],
      [-19, 30],
      [-54, 9]
    ],
    [
      [5186, 6348],
      [-12, 17],
      [-29, -91],
      [20, -35],
      [55, -6],
      [30, 16],
      [5, 67],
      [-5, 22],
      [-16, 9],
      [-48, 1]
    ],
    [
      [5063, 6171],
      [-39, -38],
      [-34, -86],
      [-45, -55],
      [-54, -25],
      [-38, -57],
      [-21, -89],
      [-36, -80],
      [-81, -121],
      [-6, -34],
      [4, -47],
      [17, -61],
      [27, -40],
      [39, -20],
      [35, 13],
      [93, 138],
      [8, 30],
      [-12, 15],
      [13, 52],
      [37, 91],
      [28, 44],
      [17, -3],
      [19, 18],
      [19, 38],
      [-3, 36],
      [-24, 32],
      [2, 14],
      [27, -6],
      [10, 34],
      [-7, 73],
      [12, 29],
      [32, -13],
      [18, 20],
      [6, 54],
      [-19, 33],
      [-44, 11]
    ],
    [[7771, 8112], [-15, -42]],
    [
      [7756, 8070],
      [-29, -4],
      [-37, -18],
      [-43, -42],
      [-61, -44],
      [-63, 11],
      [-44, 62],
      [-62, 13],
      [-82, -35],
      [-72, 2],
      [-63, 38],
      [-121, 19],
      [-14, 29],
      [3, 30],
      [18, 32],
      [11, 75],
      [2, 120],
      [17, 84],
      [30, 49],
      [-47, 20],
      [-124, -9],
      [-142, -40],
      [-161, -73],
      [-175, -25],
      [-188, 22],
      [-143, -4],
      [-158, -36],
      [-27, 12],
      [-100, -48],
      [-171, -108],
      [-129, -101],
      [-106, -112],
      [-17, -2],
      [-25, 26],
      [-54, 34],
      [-57, 1],
      [-61, -32],
      [-63, 20],
      [-63, 74],
      [-57, 21],
      [-52, -32],
      [-85, -15],
      [-188, 17],
      [-19, 24],
      [-17, 80],
      [-20, 35],
      [-37, 24],
      [-58, 89],
      [-128, 226],
      [-36, 27],
      [-104, 70],
      [-100, 31],
      [-252, 44],
      [-55, 31],
      [-32, -6],
      [-167, -122],
      [-7, -26],
      [-10, -26],
      [-28, -43],
      [-22, -20],
      [-13, 20],
      [-14, 36],
      [3, 78],
      [18, 121],
      [-10, 79],
      [-39, 36],
      [-16, 44],
      [5, 51],
      [-28, 41],
      [-61, 31],
      [-28, 41],
      [4, 50],
      [-6, 34],
      [-16, 19],
      [-65, -40],
      [-114, -96],
      [-78, -87],
      [-41, -76],
      [-61, -43],
      [-81, -10],
      [-65, -33],
      [-115, -111],
      [-85, -58],
      [-142, -41],
      [-199, -27],
      [-165, -72],
      [-133, -119],
      [-161, -100],
      [-188, -82],
      [-127, -29],
      [-68, 24],
      [-80, 50],
      [-120, 91],
      [-18, -8],
      [9, -19],
      [-2, -15],
      [-23, -26],
      [-153, -67],
      [-20, 25],
      [57, 92],
      [22, 69],
      [-1, 63],
      [25, 78],
      [53, 92],
      [19, 63],
      [-14, 32],
      [-58, 54],
      [-46, 3],
      [-60, -17],
      [-181, -136],
      [-18, 10],
      [-58, 22],
      [-15, -4],
      [-139, -118],
      [-142, -71],
      [-337, -101],
      [-86, -12],
      [-66, 7],
      [-86, 59],
      [-34, 41],
      [3, 43],
      [40, 47],
      [419, 365],
      [259, 260],
      [430, 489],
      [111, 97],
      [249, 162],
      [401, 156],
      [197, 86],
      [171, 95],
      [110, 72],
      [47, 51],
      [39, 23],
      [29, -3],
      [57, 34],
      [98, 125],
      [11, 11],
      [93, 35],
      [31, 60],
      [18, 71],
      [45, 119],
      [30, 49],
      [8, 55],
      [-16, 59],
      [5, 55],
      [24, 53],
      [115, 74],
      [312, 135],
      [10, -17],
      [-54, -104],
      [-58, -198],
      [-21, -5],
      [-19, -24],
      [-17, -44],
      [45, 3],
      [105, 49],
      [49, 25],
      [12, 49],
      [17, 97],
      [22, 54],
      [29, 10],
      [19, 46],
      [10, 84],
      [23, 49],
      [37, 16],
      [16, 29],
      [-3, 40],
      [10, 110],
      [7, 51],
      [33, 75],
      [45, 23],
      [62, -4],
      [43, -38],
      [24, -72],
      [6, -59],
      [-11, -46],
      [-60, -84],
      [-173, -202],
      [-18, -56],
      [15, -30],
      [10, -42],
      [21, -6],
      [38, 37],
      [16, 34],
      [-5, 29],
      [14, 38],
      [32, 47],
      [25, 14],
      [19, -19],
      [31, 8],
      [43, 35],
      [78, 18],
      [15, 2],
      [94, 136],
      [30, 79],
      [0, 68],
      [-38, 70],
      [-78, 74],
      [-22, 79],
      [-2, 39],
      [-15, 50],
      [16, 19],
      [47, -13],
      [40, 7],
      [33, 26],
      [58, 5],
      [47, -4],
      [57, -64],
      [205, -61],
      [76, -33],
      [94, -33],
      [24, -11],
      [15, -36],
      [35, -15],
      [54, 5],
      [56, -11],
      [96, -56],
      [88, -38],
      [139, 5],
      [76, 16],
      [121, -26],
      [68, 0],
      [69, 45],
      [59, 11],
      [49, -21],
      [22, -27],
      [40, -39],
      [122, 19],
      [27, -10],
      [36, -43],
      [6, -27],
      [20, -42],
      [46, -81],
      [10, -27],
      [32, -34],
      [26, -78],
      [19, -123],
      [22, -75],
      [24, -28],
      [11, -32],
      [22, -128],
      [55, -148],
      [81, -145],
      [109, -141],
      [136, -94],
      [164, -46],
      [145, -11],
      [127, 26],
      [169, 14],
      [210, 2],
      [109, -17],
      [8, -34],
      [-20, -50],
      [-48, -66],
      [-5, -43],
      [-3, -55],
      [-38, -74],
      [-20, -45],
      [6, -31],
      [-11, -41],
      [-28, -52],
      [16, -70],
      [62, -87],
      [76, -70],
      [89, -53],
      [55, -56],
      [23, -57],
      [25, -31],
      [27, -6],
      [25, -26],
      [22, -48],
      [-24, -84],
      [-70, -120],
      [-31, -77],
      [7, -36],
      [-60, -131],
      [-3, -52],
      [25, -36],
      [95, -51],
      [30, -42],
      [25, -6],
      [18, 30],
      [38, 13],
      [56, -2],
      [54, -21],
      [49, -40],
      [18, -41],
      [-15, -42],
      [-25, -15],
      [-77, 6],
      [-43, -21],
      [-26, -45],
      [-8, -120],
      [21, -51],
      [13, 0],
      [17, 24],
      [22, 49],
      [24, 5],
      [26, -38],
      [-28, -77],
      [-81, -117],
      [-21, -37],
      [0, -53],
      [68, -95],
      [25, 9],
      [16, 24],
      [27, 29],
      [89, 14]
    ],
    [[1353, 8797], [-11, -22], [5, -22], [39, -26], [14, 11], [-5, 22], [-25, 32], [-17, 5]],
    [
      [1455, 8581],
      [65, 44],
      [-8, 21],
      [-23, 2],
      [-133, -90],
      [-30, -35],
      [-9, -30],
      [24, -4],
      [57, 22],
      [30, 22],
      [3, 25],
      [24, 23]
    ],
    [[1562, 8725], [12, 56], [-8, 13], [-33, -9], [-75, -32], [-1, -23], [72, -14], [33, 9]],
    [
      [1699, 8922],
      [0, 32],
      [-12, 11],
      [-26, -7],
      [-16, -29],
      [-8, -48],
      [12, -12],
      [33, 25],
      [17, 28]
    ],
    [
      [4054, 11756],
      [-30, -46],
      [11, -106],
      [19, -38],
      [26, 1],
      [33, 41],
      [57, 21],
      [82, 3],
      [43, 15],
      [4, 28],
      [30, 15],
      [12, 25],
      [5, 48],
      [-48, 28],
      [-102, 9],
      [-81, -12],
      [-61, -32]
    ],
    [
      [3199, 10333],
      [287, 167],
      [31, 47],
      [-3, 25],
      [15, 36],
      [97, 90],
      [-15, 13],
      [-89, -27],
      [-167, -87],
      [-392, -238],
      [-46, -38],
      [-7, -19],
      [31, 0],
      [0, -23],
      [-30, -47],
      [4, -33],
      [38, -17],
      [75, 15],
      [112, 49],
      [40, 23],
      [6, 14],
      [-36, 7],
      [5, 26],
      [44, 17]
    ],
    [
      [4321, 9503],
      [49, 35],
      [-52, 54],
      [-85, 26],
      [-116, 0],
      [-116, -18],
      [-211, -81],
      [-74, -52],
      [-168, -184],
      [-47, -87],
      [-6, -68],
      [42, -48],
      [139, -34],
      [8, 23],
      [-7, 11],
      [-20, 13],
      [37, 63],
      [11, -26],
      [-4, -69],
      [-29, -51],
      [-30, -35],
      [0, -47],
      [29, -58],
      [37, 14],
      [44, 85],
      [122, 135],
      [-4, 31],
      [19, 51],
      [206, 174],
      [41, 54],
      [-8, 20],
      [-44, 30],
      [6, 9],
      [231, 30]
    ],
    [
      [4421, 11619],
      [69, 11],
      [-5, 37],
      [-41, 65],
      [-35, 31],
      [-29, -1],
      [-17, -10],
      [-4, -16],
      [21, -48],
      [1, -24],
      [-13, -23],
      [13, -15],
      [40, -7]
    ],
    [
      [6507, 10026],
      [10, 23],
      [-45, 40],
      [-124, 39],
      [-62, -4],
      [-56, -21],
      [-48, -39],
      [-20, -30],
      [9, -21],
      [83, -24],
      [193, 23],
      [60, 14]
    ],
    [
      [7771, 8112],
      [108, 41],
      [70, 9],
      [47, -37],
      [32, -51],
      [3, -55],
      [-22, -66],
      [-3, -59],
      [17, -51],
      [58, -27],
      [100, -3],
      [52, -14],
      [19, -42],
      [20, -7],
      [182, -19],
      [68, -19],
      [21, -28],
      [113, -36],
      [434, -89],
      [131, 11],
      [176, 30],
      [55, -49],
      [26, -6],
      [40, 8],
      [248, 25],
      [33, -48],
      [122, -46],
      [80, -3],
      [104, -5],
      [195, -22],
      [48, -16],
      [3, -52],
      [1, -22],
      [15, 0],
      [23, 24],
      [38, 9],
      [27, 22],
      [16, 36],
      [31, 6],
      [71, -45],
      [3, -25],
      [-99, -55],
      [11, -23],
      [5, -14],
      [27, -7],
      [2, -25],
      [-32, -29],
      [30, -18],
      [95, 42],
      [53, 2],
      [146, -32],
      [67, 0],
      [40, 26],
      [5, 24],
      [18, 16],
      [12, -6],
      [4, -32],
      [-3, -59],
      [64, -21],
      [129, 17],
      [77, -2],
      [73, 0],
      [13, -13],
      [18, -37],
      [55, -66],
      [19, -30],
      [17, -67],
      [40, -41],
      [13, -36],
      [2, -32],
      [12, -32],
      [36, -71],
      [48, -50],
      [25, -44],
      [9, -52],
      [26, -20],
      [41, 11],
      [17, 22],
      [-2, 43],
      [15, -8],
      [26, -73],
      [-3, -116],
      [6, -66],
      [36, -94],
      [47, -25],
      [32, -37],
      [93, -25],
      [10, 7],
      [-18, 60],
      [9, 16],
      [32, 4],
      [53, -9],
      [35, -30],
      [15, -51],
      [-8, -65],
      [-62, -66],
      [-10, -57],
      [14, -35],
      [15, -7],
      [18, 23],
      [25, -11],
      [34, -45],
      [9, -35],
      [-16, -28],
      [-61, -15],
      [-2, -27],
      [64, -59],
      [50, -79],
      [26, -8],
      [14, -25],
      [3, -42],
      [26, -28],
      [74, -32],
      [29, -82],
      [15, -48],
      [13, -55],
      [22, -6],
      [17, 26],
      [0, 49],
      [37, 19],
      [4, -26],
      [-21, -59],
      [4, -48],
      [31, -35],
      [-39, -12],
      [-110, 11],
      [-50, 48],
      [-15, 17],
      [-33, 6],
      [-11, 16],
      [15, 49],
      [-9, 14],
      [-67, -12],
      [-68, -35],
      [-44, -42],
      [1, -49],
      [29, -52],
      [57, -53],
      [24, -90],
      [-10, -128],
      [-29, -85],
      [-48, -45],
      [-142, 32],
      [-234, 108],
      [-133, 96],
      [-56, 131],
      [-21, 9],
      [-120, -42],
      [-64, -54],
      [-56, -77],
      [-16, 24],
      [25, 125],
      [10, 78],
      [-5, 31],
      [-24, 25],
      [-44, 19],
      [-58, -10],
      [-89, -59],
      [-20, 7],
      [3, 20],
      [91, 96],
      [30, 53],
      [1, 48],
      [13, 44],
      [26, 38],
      [-6, 17],
      [-39, -3],
      [-20, -17],
      [-42, -68],
      [-60, 7],
      [-8, 15],
      [17, 39],
      [-13, 25],
      [-42, 10],
      [-14, 44],
      [-14, 10],
      [-52, 31],
      [-5, 84],
      [-17, 67],
      [-31, 51],
      [-3, 51],
      [30, 84],
      [-11, 11],
      [-104, 11],
      [-195, 10],
      [-102, -17],
      [-9, -43],
      [23, -33],
      [54, -22],
      [39, -37],
      [13, -32],
      [36, -32],
      [42, -69],
      [48, -105],
      [36, -44],
      [38, -5],
      [22, -10],
      [-3, -105],
      [17, -65],
      [41, -60],
      [22, -67],
      [-2, -116],
      [3, -115],
      [-20, -74],
      [-55, -86],
      [-34, -72],
      [-16, -57],
      [-29, -38],
      [-43, -19],
      [-45, -41],
      [-46, -62],
      [-30, -77],
      [-14, -92],
      [-38, -95],
      [-61, -97],
      [-24, -112],
      [14, -124],
      [1, -537],
      [12, -133],
      [-8, -132],
      [-26, -131],
      [-46, -102],
      [-63, -74],
      [-64, -54],
      [-62, -35],
      [-59, -60],
      [-55, -84],
      [-94, -75],
      [-131, -66],
      [-105, -41],
      [-15, 34],
      [-29, 67],
      [-55, 314],
      [-81, 560],
      [-66, 348],
      [-53, 134],
      [-72, 102],
      [-90, 69],
      [-69, 31],
      [-48, -10],
      [-75, -56],
      [-140, -50],
      [-57, -33],
      [-33, -44],
      [-39, -25],
      [-44, -8],
      [-5, -12],
      [34, -16],
      [0, -43],
      [-94, -191],
      [-22, -43],
      [-57, -16],
      [-47, -38],
      [-48, -96],
      [-47, -24],
      [-67, 9],
      [-113, 58],
      [-68, 59],
      [-22, 58],
      [-2, 70],
      [20, 82],
      [-2, 66],
      [43, 138],
      [58, 49],
      [83, 11],
      [47, 22],
      [8, 32],
      [26, 23],
      [43, 14],
      [37, 87],
      [30, 161],
      [35, 88],
      [25, 15],
      [23, -5],
      [2, -19],
      [96, 103],
      [36, 71],
      [2, 73],
      [9, 213],
      [12, 138],
      [13, 63],
      [-5, 161],
      [-28, 133],
      [-109, 145],
      [-15, 49],
      [3, 51],
      [22, 54],
      [49, 13],
      [84, -42],
      [15, 10],
      [-24, 25],
      [-14, 49],
      [-4, 72],
      [-11, 40],
      [-18, 7],
      [-47, 83],
      [-6, 43],
      [5, 40],
      [-18, 40],
      [-63, 63],
      [-115, 64],
      [-276, 122],
      [-49, 47],
      [-69, 28],
      [-86, 8],
      [-55, 37],
      [-23, 65],
      [-39, 59],
      [-53, 53],
      [-84, 34],
      [-115, 16],
      [-122, 57],
      [-128, 97],
      [-72, 30]
    ],
    [
      [7357, 7056],
      [30, -2],
      [14, 14],
      [0, 18],
      [-12, 21],
      [8, 76],
      [-8, 23],
      [50, 119],
      [44, 31],
      [50, -11],
      [31, -26],
      [13, -40],
      [16, -15],
      [31, -1],
      [17, -6],
      [10, 3],
      [16, 22],
      [107, -10],
      [362, -60],
      [50, 9],
      [28, 21],
      [7, 33],
      [-24, 30],
      [-94, 65],
      [-29, 42],
      [-8, 40],
      [14, 38],
      [-20, 34],
      [-54, 30],
      [-118, 9],
      [-21, 31],
      [9, 30],
      [35, 39],
      [1, 57],
      [-11, 33],
      [-15, 98],
      [-26, 66],
      [-52, 97],
      [-21, 35],
      [-31, 21]
    ],
    [
      [8770, 7261],
      [-47, -7],
      [-45, -39],
      [-29, -43],
      [-14, -47],
      [27, -35],
      [67, -25],
      [54, 19],
      [63, 104],
      [-2, 24],
      [-25, 24],
      [-49, 25]
    ],
    [
      [8586, 7310],
      [-51, 89],
      [-28, 26],
      [-76, 16],
      [-15, -17],
      [19, -40],
      [-3, -21],
      [-32, -31],
      [-60, -42],
      [-45, -9],
      [-31, 22],
      [-20, -3],
      [-4, -41],
      [9, -19],
      [43, -37],
      [92, -20],
      [143, -3],
      [70, 26],
      [25, 52],
      [-10, 35],
      [-26, 17]
    ],
    [
      [8286, 7506],
      [21, 30],
      [13, 75],
      [-27, 62],
      [-69, 49],
      [-80, 29],
      [-90, 11],
      [-36, -35],
      [19, -81],
      [39, -94],
      [60, -106],
      [49, -28],
      [23, 8],
      [20, 35],
      [11, 45],
      [19, 16],
      [28, -16]
    ],
    [
      [7946, 7785],
      [-20, -1],
      [-7, -32],
      [11, -38],
      [50, -61],
      [14, 8],
      [5, 24],
      [-4, 38],
      [-18, 34],
      [-31, 28]
    ],
    [
      [7982, 7811],
      [0, 34],
      [-31, 76],
      [-1, 74],
      [27, 72],
      [-5, 45],
      [-38, 17],
      [-35, -9],
      [-59, -77],
      [40, -68],
      [25, -95],
      [21, -38],
      [25, -23],
      [31, -8]
    ],
    [
      [7685, 6967],
      [-128, 22],
      [-6, -26],
      [76, -80],
      [46, -30],
      [29, 3],
      [22, 16],
      [16, 30],
      [-6, 35],
      [-36, 42],
      [-13, -12]
    ],
    [
      [9154, 7009],
      [25, -35],
      [201, -23],
      [110, -38],
      [116, -71],
      [101, -47],
      [86, -22],
      [67, -37],
      [48, -52],
      [118, -58],
      [15, -33],
      [20, -14],
      [27, 8],
      [208, 198],
      [61, 10],
      [0, -20],
      [-105, -118],
      [-38, -12],
      [-20, 1],
      [-54, -53],
      [0, -13],
      [31, 7],
      [32, -14],
      [33, -34],
      [53, 17],
      [73, 70],
      [69, 104],
      [93, 213],
      [-8, 15],
      [-84, 1],
      [-14, 21],
      [26, 47],
      [8, 43],
      [-10, 39],
      [-23, -4],
      [-37, -48],
      [-28, -67],
      [-33, -128],
      [-8, 5],
      [-35, 71],
      [-3, 55],
      [14, 62],
      [-13, 37],
      [-38, 11],
      [-18, 28],
      [1, 45],
      [-14, 36],
      [-29, 28],
      [-49, -4],
      [-67, -34],
      [-52, -57],
      [-36, -79],
      [-33, 3],
      [-29, 84],
      [-31, 34],
      [-33, -16],
      [-20, 15],
      [-6, 45],
      [-38, 13],
      [-68, -17],
      [-44, -26],
      [-18, -33],
      [-18, -6],
      [-18, 22],
      [-29, -8],
      [-41, -39],
      [-12, -29],
      [24, -41],
      [-4, -16],
      [40, -17],
      [5, -25],
      [-16, -43],
      [-23, 2],
      [-29, 47],
      [-25, 17],
      [-19, -14],
      [-18, 6],
      [-15, 27],
      [-31, 5],
      [-46, -17],
      [-28, 10],
      [-9, 35],
      [-19, 27],
      [-30, 19],
      [-9, 42],
      [11, 64],
      [-24, 22],
      [-60, -20],
      [-40, -24],
      [-21, -28],
      [-30, -3],
      [-40, 22],
      [-34, -4],
      [-28, -31],
      [-27, 5],
      [-26, 40],
      [-26, 9],
      [-24, -23],
      [-15, -35],
      [-6, -48],
      [62, -43],
      [198, -66]
    ],
    [
      [10404, 6524],
      [-22, 7],
      [-94, -52],
      [-9, -22],
      [5, -18],
      [18, -15],
      [27, 4],
      [61, 55],
      [14, 41]
    ],
    [
      [11945, 5370],
      [-1, 20],
      [17, 62],
      [-14, 27],
      [-44, 11],
      [-21, -12],
      [4, -33],
      [17, -33],
      [42, -42]
    ],
    [
      [3786, 4219],
      [-9, -205],
      [-50, -114],
      [-48, -91],
      [-33, 8],
      [-27, 199],
      [-27, 142],
      [34, 129],
      [67, 132],
      [60, 3],
      [37, -25],
      [-4, -178]
    ],
    [
      [11671, 7719],
      [38, 12],
      [63, -7],
      [32, -35],
      [62, -22],
      [106, 13],
      [29, 31],
      [-25, 38],
      [17, 35],
      [72, 35],
      [198, -10],
      [77, 13],
      [126, 17],
      [159, -75],
      [59, -119],
      [-51, -81],
      [-87, -29],
      [-199, 74],
      [-105, -31],
      [-110, -51],
      [-19, 8],
      [14, 52],
      [-23, 53],
      [-48, 18],
      [-231, -30],
      [-154, 91]
    ],
    [
      [12671, 4601],
      [-31, -63],
      [-1, -53],
      [-19, -15],
      [-20, 75],
      [-18, 125],
      [-37, 63],
      [-10, 19],
      [17, 42],
      [64, 49],
      [58, 75],
      [32, 82],
      [21, 64],
      [6, 60],
      [14, 91],
      [18, -13],
      [13, -84],
      [76, -117],
      [54, -43],
      [35, -70],
      [14, -66],
      [-9, -48],
      [-30, -55],
      [-141, -65],
      [-106, -53]
    ]
  ],
  "transform": {
    "scale": [0.0009990424944179362, 0.0006327663375915449],
    "translate": [-92.10703124999999, 41.398095703124994]
  },
  "objects": {
    "ne_50m_lakes": {
      "type": "GeometryCollection",
      "geometries": [
        { "arcs": [[0, 1], [2], [3]], "type": "Polygon" },
        { "arcs": [[4], [5], [6], [7], [8]], "type": "Polygon" },
        { "arcs": [[-1, 9]], "type": "Polygon" },
        { "arcs": [[10, 11], [12], [13], [14], [15]], "type": "Polygon" },
        {
          "arcs": [[16, 17], [18], [19], [20], [21], [22], [23], [24], [25], [26]],
          "type": "Polygon"
        },
        {
          "arcs": [[27, -12, 28, -17], [29], [30], [31], [32], [33], [34], [35], [36], [37]],
          "type": "Polygon"
        },
        { "arcs": [[38]], "type": "Polygon" },
        { "arcs": [[39]], "type": "Polygon" },
        { "arcs": [[40]], "type": "Polygon" }
      ]
    }
  }
}

},{}],17:[function(_dereq_,module,exports){
module.exports={"type":"Topology","arcs":[[[107575,102770],[-61,-55],[-178,-2],[-57,22],[-60,3],[-39,-21],[-76,78],[-64,16],[-89,68],[-100,-14],[-77,21],[-273,-5],[-59,29],[122,13],[129,25],[59,86],[43,-49],[62,-15],[52,41],[79,-34],[50,17],[181,-68],[11,83],[84,-43],[90,2],[-9,51],[51,24],[29,-64],[-21,-129],[33,-33],[79,-19],[9,-28]],[[114511,104543],[32,64],[81,-38],[38,4],[-2,39],[49,14],[3,54],[62,135],[106,-28],[11,26],[71,25],[-8,60],[96,119],[55,39],[32,-30],[70,-26],[42,-88],[94,19],[29,-34],[33,21],[72,-71],[52,-107],[71,-20],[101,-64],[17,-21],[269,-135],[98,-58],[22,-36],[17,-76],[101,-105],[25,-56],[2,-65],[18,-41],[36,-18],[24,-47],[-61,-136],[-35,-8],[-34,30],[-48,15],[-26,-24],[47,-111],[-48,-35],[10,-103],[-17,-34],[-48,-30],[-81,-26],[-60,5],[-37,39],[-81,29],[-126,20],[-84,-12],[-59,-61],[5,-34],[-35,-57],[12,-99],[-34,-44],[-63,-8],[-74,25],[-99,-9],[-37,21],[60,52],[-30,73],[9,58],[-53,43],[-36,67],[-63,167],[-43,5],[-1,41],[-30,65],[-35,29],[-76,35],[-17,68],[17,58],[-19,74],[-35,28],[-19,50],[-68,73],[-98,77],[-157,94],[-12,35]],[[132844,90738],[54,-10],[9,-36],[-7,-55],[-26,-74],[-27,-36],[-28,-1],[-205,-119],[-30,-48],[-40,-41],[-52,-33],[-40,-42],[-27,-51],[-23,-22],[-282,66],[-51,-4],[-17,-111],[-22,-30],[-27,-95],[-1,-41],[21,-37],[112,7],[23,-8],[21,-53],[19,-98],[5,-99],[-8,-102],[-5,-337],[-11,-116],[-22,-43],[-2,-73],[20,-101],[3,-72],[-12,-43],[-23,-24],[-32,-6],[-20,-31],[25,-26],[-4,-32],[-65,-49],[-54,-18],[-104,-14],[-6,20],[21,68],[48,73],[0,19],[-198,229],[-41,65],[-20,50],[-2,34],[35,105],[70,172],[38,116],[5,58],[-8,71],[-21,85],[-15,165],[-30,16],[-3,38],[-60,71],[-35,77],[-35,-4],[-15,-58],[-28,0],[1,40],[-23,12],[-17,-101],[-30,-59],[-51,-52],[-36,-56],[-22,-59],[-7,-67],[10,-74],[16,-48],[20,-22],[-7,-33],[-37,-43],[-25,-45],[-13,-46],[7,-33],[82,59],[-5,-73],[-26,-28],[-4,-99],[14,-83],[-34,-2],[-19,28],[-28,-19],[-1,-48],[40,-24],[-68,-122],[-66,-80],[-57,-30],[-31,-3],[-36,52],[-6,40],[9,49],[22,58],[-1,65],[-25,73],[-11,60],[3,49],[20,73],[62,184],[15,70],[27,65],[38,59],[26,67],[17,74],[23,60],[16,108],[38,123],[30,44],[67,41],[24,-6],[32,40],[55,-18],[38,1],[69,41],[26,1],[19,-26],[-7,-37],[-30,-49],[-18,-50],[-5,-50],[15,-42],[36,-33],[18,1],[14,53],[44,-39],[29,39],[73,39],[18,112],[30,101],[-17,108],[-62,113],[-29,64],[104,4],[46,-28],[22,-112],[-3,-52],[-13,-35],[44,-13],[97,45],[44,6],[78,-28],[50,4],[192,-30],[45,28],[49,63],[6,70],[-42,103],[-13,45],[-105,-74],[-58,40],[-32,38],[-54,7],[-76,-22],[-55,11],[-34,45],[-14,49],[8,54],[37,54],[68,55],[47,21],[27,-11],[25,-49],[22,-86],[19,-13],[18,59],[24,18],[61,-31],[71,45],[1,46],[-84,73],[22,25],[71,-15],[31,-46],[-8,-78],[9,-59],[28,-42],[54,-32],[80,-22],[57,-9]],[[115468,49544],[-51,-116],[9,-341],[-7,-37],[-26,-21],[-23,-86],[-21,-51],[-60,-76],[-21,-17],[-52,-22],[-42,-35],[-45,-74],[-63,-129],[-34,-81],[-5,-30],[-73,-140],[-31,-19],[-23,27],[2,67],[-7,48],[2,53],[-16,26],[-40,-1],[-7,38],[31,79],[23,79],[29,62],[35,44],[49,90],[24,25],[90,146],[39,53],[27,16],[116,131],[47,75],[50,129],[50,91],[24,-3]],[[48437,86840],[17,-10]],[[48454,86830],[-17,-57],[-7,-88],[5,-46],[100,-61],[60,-8],[85,36],[50,-10],[21,-24],[33,-63],[29,118],[29,42],[59,51],[89,60],[67,33],[92,15],[28,21],[18,68],[47,73],[137,135],[53,36],[89,22],[125,9],[124,-18],[122,-45],[85,-20],[44,3],[177,-29],[-6,23],[-117,17],[-37,18],[-49,5],[33,58],[18,7],[84,85],[71,26],[103,19],[104,34],[89,-5],[122,26],[123,-8],[31,-14],[38,4],[71,43],[-7,41],[-53,65],[-16,79],[92,-33],[5,-34],[-15,-34],[4,-47],[19,-58],[5,-80],[-24,-34],[-44,-30],[-40,-48],[-37,-64],[-49,-53],[-62,-41],[-59,-51],[-56,-61],[-113,-85],[-217,-144],[-152,-99],[-312,-137],[-210,-115],[-128,-93],[-112,-114],[-50,-36],[-55,-9],[-87,6],[-91,-13],[-94,-34],[-116,-55],[-26,8],[-68,56],[-19,-12],[-87,10],[-40,-20],[-40,-3],[-7,33],[159,43],[-91,53],[-33,-37],[-61,23],[-31,35],[-102,57],[-64,49],[-53,10],[-8,17],[56,128],[26,41],[32,21],[4,37],[46,59],[-7,33],[8,58],[19,76],[22,55]],[[50844,87499],[18,-3],[27,50],[-44,21],[-15,-11],[14,-57]],[[48710,86336],[-3,64],[-18,-9],[-6,-66],[27,11]],[[105680,59539],[-47,-90],[-24,9],[-60,-11],[-65,-75],[-35,-21],[-42,24],[-9,95],[14,9],[27,96],[8,43],[42,43],[48,28],[34,0],[21,-23],[46,-18],[42,-109]],[[117685,33970],[13,-82],[-28,-34],[-62,105],[-48,105],[-20,18],[-32,9],[-10,72],[-21,46],[-29,7],[-4,-92],[-15,-35],[-2,-62],[-30,-36],[-23,-13],[-30,51],[-24,18],[-22,34],[-31,103],[-5,59],[5,44],[28,86],[15,89],[-9,76],[-34,65],[-29,30],[-23,-4],[-28,40],[-52,140],[-5,31],[15,76],[3,100],[-8,200],[-18,81],[-14,31],[-7,88],[-29,43],[-10,39],[0,151],[-10,52],[-18,30],[-24,8],[-20,30],[-21,96],[-3,38],[16,150],[24,68],[37,47],[17,42],[-1,38],[18,48],[56,102],[3,24],[-20,214],[-16,106],[-18,78],[-4,57],[10,34],[-4,54],[-18,429],[-8,23],[23,68],[-1,30],[19,34],[-7,45],[-26,75],[-28,60],[-58,92],[-24,55],[-18,69],[-25,140],[-14,30],[-8,86],[7,36],[22,50],[-3,87],[28,65],[54,-39],[136,-181],[4,-23],[86,-180],[29,-77],[0,-39],[27,-158],[5,-73],[-5,-101],[11,-64],[-4,-27],[32,-94],[14,-79],[3,-77],[-8,-75],[-28,-110],[-4,-44],[6,-40],[17,-44],[30,-47],[22,-57],[15,-66],[13,-30],[45,-10],[23,-39],[28,-65],[12,-107],[-1,-147],[-13,-127],[-26,-107],[-20,-58],[-40,-44],[-36,-58],[-15,-41],[6,-23],[0,-131],[-10,-54],[57,-206],[10,-115],[-21,-97],[-2,-86],[19,-77],[5,-91],[-6,-109],[7,-69],[20,-32],[9,-88],[-2,-145],[6,-80],[12,-15],[82,-12],[16,-12],[24,-77],[49,-223],[-2,-23],[49,-176]],[[117311,36036],[11,49],[-29,40],[-22,-10],[-7,-34],[22,-44],[25,-1]],[[117334,36004],[13,-16],[36,28],[6,54],[-23,26],[-21,-6],[-11,-86]],[[117367,47229],[56,-172],[-9,-19],[-56,-27],[-74,-4],[-57,-23],[-12,-14],[10,-89],[-11,-29],[-19,-6],[-28,40],[-28,-1],[-65,17],[-25,-26],[-56,-103],[-11,-61],[5,-75],[22,-68],[33,-76],[-24,-89],[-30,-21],[-4,-34],[-77,-144],[-17,-56],[3,-29],[26,-40],[-6,-17],[-57,-3],[-9,-19],[-4,-58],[54,-66],[-15,-16],[-57,25],[-36,6],[-37,-25],[-12,-61],[17,-63],[-43,-49],[-5,-65],[-27,-37],[-34,-20],[-51,9],[-41,-49],[-20,-54],[90,7],[15,-12],[28,-64],[-66,-1],[-22,-33],[-44,9],[-39,-46],[-3,-22],[29,-27],[141,-35],[39,-5],[98,52],[24,-11],[10,-64],[-25,-43],[-101,-112],[-49,-67],[-46,-28],[6,-40],[-16,-27],[-120,39],[-84,74],[-40,8],[-31,-18],[-22,-27],[-9,-32],[-5,-168],[46,-173],[-41,-12],[-66,-103],[-9,98],[49,65],[8,37],[-28,104],[6,96],[-7,27],[-35,2],[-39,20],[-24,49],[-56,-64],[-16,51],[-48,8],[-48,47],[-22,84],[-22,20],[-24,0],[-39,-78],[-35,-10],[-11,-16],[26,-86],[-6,-48],[-58,20],[-26,-29],[-10,-57],[5,-64],[15,-75],[-19,5],[-31,73],[-19,-25],[-22,-5],[-17,-86],[-20,45],[-2,52],[20,68],[-27,22],[-20,64],[30,60],[-12,61],[-16,40],[-21,27],[-16,48],[0,184],[-7,59],[1,79],[57,412],[28,68],[9,55],[23,250],[-10,71],[-38,6],[-16,114],[-19,91],[56,208],[28,81],[75,149],[3,62],[-35,70],[6,46],[48,38],[13,32],[61,16],[30,45],[28,-13],[16,24],[39,-5],[19,48],[31,25],[7,-35],[21,-3],[28,16],[18,28],[14,88],[16,39],[18,16],[20,-59],[-3,-49],[9,-27],[52,10],[43,47],[15,0],[10,-53],[11,-16],[25,11],[37,37],[29,44],[4,42],[36,12],[37,127],[-8,38],[50,32],[55,-9],[-21,-79],[13,-17],[63,-15],[14,-38],[-26,-21],[25,-50],[29,-10],[28,14],[34,74],[24,14],[19,-7],[9,-73],[10,-21],[67,-21],[50,30],[17,-4],[0,-116],[4,-60],[28,-82],[18,-17],[32,2],[-4,-51],[6,-28],[63,-90],[35,-85],[15,-5],[25,40],[10,65],[28,44],[70,50],[105,48]],[[115809,45138],[3,60],[-9,20],[-85,-38],[-47,86],[-16,-40],[-5,-81],[35,-56],[3,-69],[30,11],[11,66],[70,20],[10,21]],[[115952,46949],[12,50],[-29,44],[-45,-38],[-8,19],[1,71],[-56,4],[-12,-22],[23,-63],[61,-38],[-2,-33],[-23,-41],[-20,-79],[82,50],[16,76]],[[116464,45297],[5,40],[-25,44],[-25,97],[-74,-8],[-13,67],[-38,-2],[-31,-31],[9,-39],[3,-76],[28,-23],[53,-18],[79,-50],[29,-1]],[[116580,47516],[9,53],[-41,24],[-55,-19],[-4,-82],[-19,-37],[47,-6],[5,-21],[39,-39],[-14,98],[33,29]],[[161330,99227],[-4,-87],[-39,-108],[-27,-178],[12,-94],[-48,-177],[-22,-40],[-23,-122],[-36,-44],[-12,-66],[-20,-167],[3,-110],[34,-83],[-9,-54],[-22,-51],[-44,-59],[-106,-121],[-27,-71],[1,-124],[-12,-21],[-87,1],[-12,27],[30,153],[-39,35],[-61,-25],[-31,-33],[-128,-209],[-50,-65],[117,53],[64,15],[46,-19],[28,-54],[12,-54],[-4,-53],[-30,-46],[-96,-59],[-61,19],[-158,-129],[-53,-50],[19,-68],[-59,-86],[-133,-101],[-34,-53],[-108,-62],[-108,-12],[-112,-50],[-169,-99],[-94,-63],[-20,-26],[-35,-84],[-54,-38],[-89,51],[-59,11],[-108,-56],[-26,-38],[-7,-48],[13,-57],[-11,-62],[-43,-137],[-72,-88],[-57,-83],[-99,-68],[-212,-81],[-217,-54],[-94,-32],[-46,-30],[-100,-7],[-154,16],[-137,33],[-199,92],[-61,57],[46,34],[200,36],[123,13],[113,1],[74,20],[34,40],[54,20],[72,0],[68,11],[66,22],[46,30],[24,39],[78,47],[89,73],[71,101],[98,160],[118,117],[135,75],[93,60],[180,199],[15,34],[57,63],[-5,53],[-54,14],[42,55],[128,129],[94,70],[177,114],[24,103],[30,58],[95,94],[160,129],[94,119],[31,56],[50,128],[65,97],[27,68],[16,74],[27,59],[76,101],[37,67],[56,169],[30,42],[-9,50],[21,30],[58,41],[70,155],[15,64],[9,135],[-13,54],[2,67],[26,37],[44,31],[153,153],[108,2],[51,-18],[63,-77],[-55,-91],[2,-34]],[[160103,97104],[17,117],[-30,6],[-94,-108],[-135,-69],[-209,-96],[-21,-27],[2,-51],[205,41],[73,33],[117,110],[75,44]],[[40675,94732],[8,-36],[37,-63],[-56,5],[-36,25],[-52,67],[-27,-60],[-12,-104],[18,-78],[-11,-25],[-51,-41],[-45,0],[-61,19],[-39,40],[-18,102],[-12,106],[17,150],[-8,78],[24,145],[50,43],[39,86],[48,32],[38,65],[-2,17],[-56,22],[-107,-150],[-41,-1],[-2,43],[40,75],[101,209],[3,36],[-35,30],[-18,80],[-31,13],[-108,-56],[-37,-54],[2,-66],[-9,-32],[-47,-48],[-33,38],[-8,95],[8,86],[22,76],[-6,93],[-34,108],[-29,55],[-52,-12],[-29,-25],[-32,-71],[-29,59],[43,69],[18,43],[-21,27],[-124,-1],[4,-84],[-33,-123],[-18,-23],[-81,36],[-44,50],[-4,66],[-16,45],[-39,58],[-3,45],[-26,65],[-141,96],[-15,31],[-8,66],[-50,66],[-150,157],[-6,54],[-44,69],[0,55],[72,46],[159,50],[-14,44],[-101,16],[-57,-5],[-115,-27],[-35,23],[-111,128],[23,134],[68,119],[25,71],[35,152],[-11,32],[42,107],[70,28],[115,-51],[92,9],[91,-17],[171,-116],[58,-9],[5,36],[-102,78],[-75,48],[-34,36],[-22,49],[1,49],[25,49],[47,139],[89,179],[52,-33],[-29,-55],[81,1],[-12,-48],[-128,-60],[-76,-164],[12,-30],[77,-58],[16,43],[55,48],[75,12],[29,-27],[-122,-32],[19,-83],[62,-144],[2,-89],[23,-55],[48,-63],[32,-107],[42,-102],[62,-84],[22,-77],[50,-86],[-14,-43],[-28,-14],[15,-81],[38,-89],[79,-104],[50,-182],[24,-61],[41,-134],[-10,-61],[36,-119],[57,-78],[15,-48],[-6,-71],[43,-10],[37,-38],[0,-38],[-29,-10],[20,-56],[71,-100],[30,-61],[2,-51],[30,-59],[88,-98],[32,-54],[1,-40],[-29,-27],[-33,-85],[-5,-64],[21,-151],[5,-91],[-8,-29]],[[39759,96217],[30,37],[-33,162],[-30,10],[-30,-77],[2,-58],[38,-77],[23,3]],[[40523,95123],[20,38],[-65,10],[-16,-58],[-27,-42],[-35,-24],[-1,-37],[48,-16],[48,29],[28,100]],[[33274,106021],[61,-18],[-44,-46],[-7,-46],[-144,-87],[-7,57],[-22,56],[-55,20],[-173,36],[-32,19],[-116,16],[-199,12],[-221,1],[-55,-61],[-124,-34],[-176,-62],[-2,-41],[255,46],[176,-6],[73,24],[320,-11],[75,-22],[-59,-40],[-71,-29],[-134,-21],[-83,15],[-74,-25],[-57,-5],[-123,8],[68,-60],[-34,-81],[-26,-8],[-117,15],[-65,-21],[-137,-93],[-7,-45],[-48,-38],[-116,-58],[-30,-66],[80,-9],[-70,-93],[-130,-87],[-144,-113],[-67,-85],[-2,-19],[-71,-45],[-203,-107],[-172,-67],[-22,54],[-49,33],[-33,-16],[-124,-21],[-23,-26],[-68,-36],[-95,-72],[-16,-51],[13,-30],[-25,-48],[45,-5],[11,-48],[-58,-73],[-61,-29],[-82,-17],[-66,3],[-51,22],[-96,-16],[-143,-55],[-84,-21],[-28,15],[-82,-13],[-246,-60],[-82,-7],[-157,19],[-78,-9],[-106,19],[-215,79],[-46,38],[-15,33],[-48,-11],[-63,20],[-217,12],[-80,21],[-31,22],[-132,133],[19,27],[237,-111],[159,52],[41,33],[-42,53],[20,30],[99,-5],[84,-37],[147,-17],[45,-36],[101,-42],[-18,58],[64,36],[-20,40],[-78,51],[70,39],[54,-11],[63,27],[-16,67],[52,121],[50,43],[9,50],[28,15],[100,21],[72,-6],[30,-34],[185,30],[74,41],[49,53],[-8,51],[-65,30],[-11,36],[-40,35],[-55,22],[-65,60],[5,51],[-71,66],[-72,-7],[-10,16],[30,70],[-31,39],[11,104],[-241,3],[-97,139],[-73,44],[-27,42],[20,37],[53,-1],[88,-41],[5,-50],[48,-3],[224,-53],[52,-60],[54,-42],[120,-30],[78,-50],[33,-34],[52,-16],[109,13],[50,33],[40,4],[33,-25],[-11,-56],[53,-44],[78,-10],[76,-58],[80,-96],[106,-79],[155,-13],[115,22],[45,-9],[66,51],[19,-34],[211,34],[134,57],[128,89],[178,130],[46,48],[33,72],[76,92],[21,61],[49,57],[79,54],[104,52],[239,86],[87,17],[274,6],[97,-26],[108,-12],[39,-15],[169,-26],[177,-36],[14,-20]],[[31682,105343],[22,46],[-71,-17],[-92,-77],[-50,-53],[-138,-54],[-127,-38],[-101,-7],[-160,-34],[-24,-21],[53,-44],[97,2],[14,-107],[24,-42],[55,-20],[82,20],[73,42],[48,99],[65,76],[123,74],[82,80],[25,75]],[[52718,88895],[-2,-51],[-25,-31],[-237,-200],[-67,-72],[42,-56],[49,35],[25,-36],[-12,-48],[43,9],[2,-34],[-58,-86],[-51,-5],[39,-105],[5,-81],[9,-24],[-5,-66],[-51,-35],[-50,-2],[-57,-26],[-63,-50],[-43,-45],[-23,-42],[-112,-43],[-54,-34],[-60,9],[-112,6],[-87,-6],[-64,-19],[-56,11],[-49,41],[-76,32],[-197,35],[-165,-7],[-136,-24],[-213,-71],[-93,-41],[-73,-41],[-37,10],[-74,1],[-117,40],[-28,63],[69,124],[52,118],[72,64],[28,-2],[56,42],[79,84],[62,52],[42,19],[115,19],[72,25],[52,8],[31,-6],[47,10],[65,27],[87,20],[110,12],[69,14],[31,16],[59,-1],[27,33],[57,-6],[18,-45],[30,-26],[48,-8],[55,16],[21,-45],[58,23],[-17,-45],[40,-31],[56,16],[71,43],[14,33],[-78,-12],[11,33],[80,75],[-1,30],[-88,-48],[-37,-2],[-8,31],[20,64],[-19,24],[-91,-33],[-15,-38],[-20,6],[3,69],[44,10],[27,-12],[85,28],[22,-8],[7,-41],[-31,-58],[41,-16],[120,101],[58,33],[58,20],[57,5],[49,17],[39,30],[166,67],[68,42],[64,74],[26,-5]],[[52191,88571],[26,57],[-69,-26],[-15,-38],[58,7]],[[52455,88659],[6,20],[-120,-25],[-30,-45],[8,-51],[27,-2],[34,32],[20,48],[55,23]],[[28114,109401],[176,56],[16,-22],[-93,-116],[-3,-77],[-113,-97],[-67,-41],[47,-36],[-10,-22],[-87,-7],[-74,-122],[-6,-46],[30,-23],[165,14],[-16,-78],[-24,-38],[-53,-27],[-202,-28],[-43,18],[-68,-36],[-85,46],[-70,57],[-64,28],[-194,11],[-134,-4],[-67,-9],[-50,-29],[4,-35],[-161,-90],[-9,-30],[132,-87],[130,-79],[19,-38],[-28,-22],[-94,-1],[-98,-17],[-124,-113],[-37,-26],[-180,-96],[-231,-110],[-123,-52],[-93,-68],[-80,6],[-61,43],[-5,50],[48,34],[181,32],[51,33],[6,40],[258,194],[85,114],[22,91],[-9,46],[-48,25],[-116,19],[-184,-49],[-119,-53],[-109,-57],[-105,-67],[-34,-125],[-6,-127],[-90,-53],[-158,-71],[-304,58],[-194,43],[-281,42],[-44,2],[1,64],[150,30],[196,57],[50,61],[-1,61],[-69,74],[59,42],[241,116],[98,90],[25,61],[-214,118],[-11,40],[94,19],[360,-55],[174,8],[106,41],[26,50],[-17,35],[-148,29],[-96,52],[-243,52],[-190,16],[-339,4],[13,-77],[-27,-44],[-67,-58],[-26,37],[-116,30],[-101,3],[-117,-13],[-165,5],[-122,25],[-140,-4],[-53,-40],[0,-88],[-31,-8],[-52,31],[-50,-93],[-47,24],[-72,84],[-41,94],[16,85],[150,29],[164,-3],[71,31],[495,111],[224,84],[52,-41],[142,56],[185,25],[3,37],[105,-15],[212,49],[308,97],[129,86],[317,70],[58,22],[217,44],[202,49],[71,-37],[50,2],[268,-56],[67,-22],[-91,-75],[-298,-124],[-188,-40],[-27,-50],[-117,-58],[3,-63],[-50,-14],[-39,-32],[19,-24],[250,-14],[32,-24],[216,-44],[61,28],[225,-11],[204,5],[89,15],[156,48],[85,109],[110,81],[162,54],[26,-17],[-30,-84],[-77,-85],[13,-35]],[[48437,86840],[52,57],[52,26],[26,62],[7,72],[17,47],[26,19],[6,71],[59,30],[31,-34],[2,-62],[-9,-34],[19,-31],[93,-28],[20,-21],[3,-58],[-8,-41],[-51,-42],[-44,-7],[-162,20],[-61,-8],[-61,-48]],[[105269,102539],[31,-2],[-12,-48],[-46,-40],[-10,-112],[-103,-134],[-36,-18],[-129,-37],[-93,-43],[-69,-96],[-40,-3],[-21,28],[7,52],[-63,6],[-52,-30],[-48,-64],[-83,-48],[-60,-65],[-69,14],[-101,-22],[-3,56],[93,75],[69,97],[3,22],[-50,19],[-13,31],[-2,65],[-44,35],[28,36],[-5,75],[19,0],[18,-70],[42,0],[64,178],[22,27],[72,7],[100,-54],[28,-64],[63,-45],[28,30],[-39,156],[-32,94],[-11,72],[23,64],[44,-30],[35,43],[107,12],[98,19],[53,-11],[19,-44],[46,9],[77,-6],[49,-20],[1,-67],[-31,-79],[-37,-62],[9,-22],[54,14]],[[49885,72446],[-47,-135],[-19,-42],[-25,-6],[-43,30],[-22,27],[-21,55],[-27,42],[-28,19],[-26,38],[-6,36],[15,34],[67,80],[41,66],[39,32],[26,0],[41,-67],[20,-53],[20,-98],[-5,-58]],[[47437,57840],[13,-51],[13,-83],[-19,-49],[-55,-28],[-26,16],[-38,5],[-61,22],[-88,40],[-88,47],[-62,39],[-72,63],[-46,65],[-46,122],[-22,75],[-17,94],[-21,36],[-1,86],[-8,148],[2,47],[24,89],[43,8],[71,-40],[51,-44],[86,-107],[47,-44],[60,-68],[74,-112],[72,-127],[54,-106],[60,-143]],[[118869,58339],[-19,24],[-49,-3],[-33,36],[-20,-24],[-22,21],[-23,51],[4,114],[12,127],[11,43],[28,51],[41,-5],[17,-27],[44,-6],[37,45],[74,29],[41,-42],[15,-33],[6,-55],[16,-45],[-9,-92],[-16,-71],[-22,-30],[-10,-70],[-32,-76],[-31,-114],[-36,4],[-15,30],[12,38],[-21,80]],[[56475,32151],[-71,57],[22,47],[1,77],[-105,58],[-71,75],[-38,55],[-24,-12],[12,-60],[-15,-14],[-89,86],[-3,89],[35,81],[21,12],[60,-54],[0,49],[-30,99],[-35,99],[-5,27],[49,50],[12,32],[42,22],[12,17],[32,4],[36,-74],[49,-27],[55,-95],[51,-34],[36,-103],[17,-6],[47,-48],[51,-81],[37,-28],[15,-50],[0,-57],[11,-15],[53,-16],[43,-4],[20,-23],[-18,-43],[-36,-62],[10,-48],[24,-51],[80,-28],[34,-1],[25,-17],[2,-59],[-61,-22],[-26,-53],[-29,-23],[-69,-16],[5,-37],[17,-34],[7,-44],[-9,-39],[-41,6],[-63,53],[-29,39],[27,62],[8,50],[-3,78],[57,53],[36,-15],[-3,32],[-33,28],[-40,61],[-37,27],[-14,-5],[-3,-42],[-16,-93],[-51,-21],[-60,10],[-24,19]],[[38933,97040],[12,-36],[-35,-28],[-3,-43],[-38,-29],[-34,3],[-51,99],[-47,-21],[-66,-3],[-45,-32],[-57,-5],[-134,51],[-43,32],[-80,79],[-77,35],[-141,94],[-34,17],[-29,41],[22,29],[125,-6],[28,22],[-24,47],[143,13],[39,9],[6,45],[-41,110],[-29,30],[11,30],[49,39],[-10,31],[-75,-2],[-5,50],[89,9],[17,14],[-2,117],[-31,96],[14,20],[53,-11],[47,39],[24,-6],[75,-50],[73,-7],[-34,-29],[-147,-105],[-28,-25],[51,-20],[27,-49],[-33,-75],[18,-22],[91,36],[69,48],[109,100],[33,-22],[-74,-132],[-27,-32],[-78,-8],[-10,-62],[-87,-61],[-56,-15],[-29,-44],[118,-93],[24,-39],[28,-78],[-11,-49],[32,-94],[55,-2],[63,43],[60,20],[89,16],[15,-24],[4,-51],[32,-54]],[[118362,104699],[4,-41],[-31,-73],[13,-61],[-18,-43],[-50,-43],[-105,-48],[-128,-107],[18,-42],[-48,-34],[-64,36],[-99,31],[-59,79],[-29,16],[-126,-33],[-54,38],[-88,23],[-43,-40],[-23,11],[-20,54],[25,91],[80,17],[59,-5],[77,-27],[0,-18],[-76,-7],[2,-40],[36,-29],[175,-14],[50,-37],[34,11],[28,69],[-12,83],[-44,89],[-46,60],[-90,68],[-88,52],[-223,108],[-83,51],[-58,46],[-84,86],[-38,55],[100,-45],[40,16],[19,47],[-5,59],[-28,58],[-166,158],[1,87],[-93,109],[-70,28],[-20,86],[64,-5],[49,-75],[112,-141],[8,-61],[44,-38],[102,-57],[36,11],[-12,69],[15,59],[-29,87],[51,24],[-87,125],[55,-9],[124,-182],[23,-100],[80,-102],[50,-22],[30,25],[-13,44],[-54,79],[26,5],[144,-64],[66,29],[58,52],[8,48],[-16,63],[-63,106],[-47,51],[-48,21],[-72,13],[-119,78],[-61,-48],[-59,28],[-69,86],[-82,150],[30,39],[43,-39],[57,-27],[122,-28],[101,-69],[142,-70],[80,-57],[24,18],[47,-45],[91,-52],[37,-56],[11,-101],[-13,-66],[-41,-133],[-32,-56],[8,-37],[76,-83],[30,-90],[68,-93],[44,-94],[109,-151],[72,-57],[28,-37]],[[117723,105206],[12,58],[-23,74],[-33,3],[-33,-86],[-20,-120],[24,-19],[42,29],[31,61]],[[31548,85640],[8,-57],[-80,-57],[-55,-104],[-45,-27],[-38,-39],[-22,2],[-38,108],[-6,55],[-19,64],[-18,-23],[-22,-74],[-20,8],[-30,74],[-57,85],[-13,37],[-21,108],[-16,47],[-86,143],[-21,105],[0,81],[-10,20],[62,37],[64,28],[33,-7],[20,-52],[-17,-63],[17,-37],[61,-61],[45,-10],[22,-105],[30,-84],[28,4],[14,24],[-21,172],[17,23],[76,-60],[89,14],[16,-7],[-8,-50],[-22,-11],[-74,9],[-26,-26],[-10,-32],[129,-204],[64,-58]],[[31432,85521],[9,75],[-24,94],[-39,2],[8,-88],[46,-83]],[[47467,90117],[-60,-22],[-50,6],[-15,-64],[-41,-47],[-24,-48],[0,-92],[19,-38],[66,-35],[-15,-34],[-93,-13],[-56,-21],[-56,-39],[-34,-34],[-10,-77],[7,-109],[-4,-43],[-56,-180],[-39,-23],[32,101],[-31,-19],[-26,-82],[-31,20],[1,42],[24,82],[-9,17],[20,88],[-2,75],[-17,29],[-42,-61],[-56,-108],[-43,-48],[-55,6],[-19,-32],[-42,-16],[-14,-38],[-5,-76],[-22,-43],[-65,-26],[-13,-48],[10,-36],[1,-74],[-11,-110],[-21,-96],[-33,-83],[-37,-67],[-61,-90],[42,-140],[4,-121],[-58,-129],[-1,-18],[63,-215],[61,-171],[48,-183],[15,-120],[2,-119],[-8,-115],[-16,-111],[-23,-106],[-31,-102],[-32,-79],[-66,-133],[-32,-88],[-34,-71],[-34,-54],[-63,-61],[-160,-112],[-51,-21],[-55,-9],[-60,6],[-33,12],[-60,60],[-24,42],[-27,88],[-29,131],[-25,79],[-23,26],[-37,93],[-9,121],[0,189],[6,128],[12,69],[-3,52],[-30,80],[-4,54],[-28,83],[8,58],[-13,132],[12,95],[33,90],[14,89],[15,57],[44,116],[-14,80],[0,92],[38,142],[32,62],[40,45],[14,53],[-15,59],[6,87],[24,114],[26,86],[50,112],[6,75],[-43,61],[-36,8],[-44,-24],[-28,-1],[-43,-53],[-37,-92],[-27,-38],[-107,-93],[-29,18],[-5,44],[17,69],[27,67],[56,109],[15,58],[87,47],[32,37],[12,47],[3,70],[56,85],[63,108],[56,112],[50,115],[37,65],[41,50],[26,107],[28,-2],[10,-112],[21,-21],[32,19],[28,41],[22,63],[26,26],[30,-11],[81,33],[3,-72],[-53,-46],[-24,-52],[1,-41],[26,-31],[25,54],[44,66],[74,45],[28,28],[8,33],[39,66],[37,21],[50,9],[49,-5],[47,-21],[73,36],[74,7],[55,32],[34,53],[45,27],[55,2],[73,-19],[93,-40],[76,-53],[99,-107],[26,-11]],[[47470,90183],[-3,-66]],[[47040,90087],[-27,14],[-20,-35],[-20,-76],[-3,-44],[36,-6],[24,20],[10,36],[0,91]],[[46713,89523],[8,-50],[36,-1],[-13,45],[-31,6]],[[46202,89764],[-24,-44],[11,-20],[50,5],[3,40],[-40,19]],[[46129,89660],[-22,-22],[-20,-52],[-26,-32],[-32,-15],[-22,-34],[-12,-52],[-21,-47],[-48,-72],[-1,-48],[26,-60],[43,-4],[54,82],[5,57],[39,80],[32,31],[10,126],[29,4],[4,32],[-38,26]],[[47712,90808],[-8,-25]],[[47704,90783],[-39,-13],[-61,-51],[-37,7],[-25,36],[-37,8],[-47,-21],[-43,1],[-36,23],[-71,11],[11,169],[10,50],[-10,41],[-72,-5],[-83,-24],[-94,-44],[-102,-15],[-110,14],[-83,-2],[-93,-22],[-16,7],[-58,-28],[-100,-64],[-75,-60],[-63,-66],[-55,34],[-34,1],[-36,-20],[-36,13],[-37,44],[-34,12],[-30,-19],[-50,-8],[-109,9],[-22,62],[-67,87],[-75,134],[-82,57],[-58,19],[-148,26],[-50,15],[-98,-72],[-39,-69],[-15,33],[11,118],[-5,47],[-33,47],[-13,54],[-52,43],[-10,61],[-39,-23],[-66,-58],[-45,-51],[-24,-45],[-36,-25],[-47,-7],[-39,-19],[-67,-65],[-50,-34],[-83,-25],[-115,-16],[-97,-43],[-78,-69],[-94,-60],[-110,-49],[-74,-17],[-40,15],[-117,83],[-19,-40],[-90,-40],[-12,15],[34,55],[12,77],[56,138],[-41,51],[-63,-8],[-105,-81],[-53,17],[-82,-70],[-83,-41],[-197,-60],[-89,-4],[-50,36],[-18,49],[23,27],[246,217],[151,154],[251,289],[65,57],[145,96],[235,92],[115,51],[100,56],[64,43],[50,44],[51,17],[64,81],[54,21],[18,35],[59,174],[-9,35],[17,64],[67,44],[182,80],[-60,-189],[-7,-41],[90,44],[17,86],[42,65],[5,49],[35,39],[18,136],[19,44],[63,11],[25,-22],[17,-77],[-6,-27],[-35,-50],[-101,-120],[-11,-33],[27,-46],[32,42],[24,67],[43,2],[25,21],[55,12],[55,80],[17,47],[1,40],[-23,42],[-45,43],[-23,100],[79,23],[62,0],[33,-38],[120,-35],[113,-46],[29,-31],[65,-3],[107,-55],[126,12],[71,-16],[39,1],[41,26],[63,-6],[37,-38],[70,11],[37,-32],[82,-171],[24,-117],[20,-35],[13,-76],[32,-88],[48,-85],[64,-84],[79,-55],[96,-28],[85,-6],[74,15],[99,8],[122,2],[64,-10],[-7,-50],[-28,-39],[-5,-59],[-34,-69],[-19,-74],[10,-41],[36,-52],[44,-41],[52,-32],[33,-32],[27,-53],[44,-47],[-15,-49],[-41,-71],[-13,-67],[-35,-78],[12,-52],[56,-30],[17,-24],[48,21],[33,-1],[60,-36],[10,-24],[-23,-34],[-45,4],[-40,-40],[-5,-70],[42,12],[30,-19],[-17,-46],[-59,-91],[0,-31],[39,-56],[41,36],[51,9]],[[44021,91085],[-58,-13],[-9,-41],[51,26],[16,28]],[[44083,91170],[3,41],[-64,-24],[61,-17]],[[45540,92963],[-17,-27],[6,-63],[26,-21],[52,36],[74,10],[30,69],[-28,17],[-60,5],[-48,-7],[-35,-19]],[[45040,92121],[168,99],[25,64],[57,53],[-61,-8],[-98,-52],[-229,-141],[-27,-22],[14,-25],[-18,-28],[26,-29],[43,9],[66,29],[8,42],[26,9]],[[45696,91631],[29,20],[-31,32],[-50,15],[-67,0],[-69,-10],[-122,-48],[-44,-31],[-98,-109],[-27,-51],[-4,-40],[25,-29],[81,-20],[-20,-42],[17,-62],[22,8],[26,50],[70,81],[10,48],[120,103],[24,32],[-27,35],[135,18]],[[45754,92882],[38,28],[-44,57],[-30,-16],[5,-56],[31,-13]],[[46974,91940],[-21,37],[-73,23],[-36,-2],[-60,-36],[-7,-30],[49,-14],[148,22]],[[47712,90808],[63,24],[41,5],[47,-52],[-12,-72],[8,-64],[35,-17],[89,-9],[22,-30],[107,-11],[52,-27],[66,-22],[254,-52],[76,6],[103,18],[32,-30],[183,17],[19,-28],[72,-28],[108,-4],[114,-13],[27,-10],[3,-44],[45,20],[42,38],[43,-42],[-57,-32],[26,-41],[0,-28],[55,25],[30,1],[86,-19],[39,0],[44,36],[1,-54],[37,-13],[75,10],[88,-1],[61,-86],[10,-40],[24,-24],[8,-40],[29,-61],[42,-56],[6,-30],[71,-15],[2,-108],[21,-55],[46,-37],[55,-15],[-5,40],[55,6],[30,-48],[-6,-38],[-36,-39],[3,-54],[33,2],[21,-26],[-5,-37],[-37,-25],[38,-35],[52,-66],[18,-42],[43,-19],[33,-110],[34,2],[-2,-56],[-65,7],[-63,51],[4,38],[-39,-7],[-66,-46],[0,-30],[51,-61],[14,-53],[-6,-76],[-17,-51],[-28,-26],[-83,19],[-137,64],[-77,57],[-33,77],[-13,5],[-70,-24],[-79,-64],[20,120],[-16,34],[-61,5],[-51,-35],[-10,16],[53,57],[41,108],[-27,8],[-36,-50],[-35,4],[-2,47],[-72,56],[-2,50],[-29,69],[10,87],[-175,12],[-60,-9],[9,-46],[54,-35],[53,-78],[49,-88],[35,-10],[-1,-62],[10,-38],[36,-75],[1,-136],[-12,-45],[-32,-50],[-29,-77],[-68,-57],[-27,-37],[-18,-46],[-8,-54],[-22,-56],[-36,-58],[-14,-66],[8,-74],[1,-317],[7,-79],[-5,-77],[-15,-78],[-27,-61],[-37,-43],[-73,-53],[-35,-35],[-32,-50],[-55,-44],[-76,-39],[-62,-24],[-25,59],[-33,186],[-47,331],[-39,206],[-31,79],[-41,60],[-54,41],[-40,19],[-27,-6],[-44,-33],[-115,-49],[-43,-41],[-8,-47],[-68,-139],[-61,-31],[-28,-57],[-28,-15],[-39,6],[-66,34],[-39,35],[-15,76],[11,87],[25,82],[34,29],[76,19],[20,33],[25,8],[22,52],[17,95],[20,52],[30,-6],[56,61],[21,42],[6,170],[15,118],[-3,95],[-17,79],[-63,86],[-7,60],[13,32],[28,7],[50,-25],[-16,93],[-44,77],[-1,48],[-47,61],[-68,38],[-161,72],[-29,28],[-40,17],[-51,4],[-31,22],[-36,74],[-31,31],[-50,20],[-67,10],[-71,34],[-75,57],[-42,17]],[[47470,90183],[19,30],[0,58],[29,71],[26,19],[48,-23],[17,-32],[42,11],[63,-6],[212,-36],[29,6],[20,31],[-69,56],[-21,49],[-4,43],[-31,18],[-69,5],[-13,18],[26,41],[-15,111],[-57,118],[-18,12]],[[48296,90305],[-27,-5],[-44,-48],[8,-49],[39,-14],[32,11],[37,62],[-45,43]],[[48189,90333],[-46,69],[-45,9],[0,-46],[-53,-43],[-53,-30],[24,-22],[55,-12],[83,-2],[41,16],[14,31],[-20,30]],[[48013,90449],[20,62],[-16,37],[-40,29],[-47,18],[-53,6],[-20,-21],[33,-103],[36,-64],[42,-11],[18,47],[27,0]],[[47836,90630],[-18,64],[-2,44],[14,69],[-43,5],[-34,-45],[22,-40],[28,-80],[33,-17]],[[47662,90131],[-75,13],[-3,-16],[44,-47],[44,-16],[22,27],[-32,39]],[[48521,90155],[14,-21],[118,-13],[64,-22],[68,-43],[59,-27],[50,-13],[40,-22],[28,-31],[104,-57],[122,117],[36,-6],[-62,-70],[-33,-6],[-32,-32],[56,-32],[31,10],[43,42],[40,61],[54,127],[-61,21],[19,54],[-19,20],[-38,-68],[-19,-75],[-25,45],[6,69],[-41,45],[-7,47],[-45,15],[-40,-21],[-30,-34],[-21,-46],[-20,1],[-16,51],[-49,19],[-26,34],[-66,-25],[-10,-20],[-38,5],[-31,-41],[-2,-54],[-45,21],[-61,-1],[-39,73],[-8,50],[-88,-44],[-23,13],[-52,-18],[-31,29],[-26,-62],[36,-26],[116,-39]],[[49251,89868],[-67,-26],[23,-30],[44,56]],[[50152,89186],[9,48],[-46,16],[12,-39],[25,-25]],[[40058,71013],[-17,-75],[-47,-179],[-45,-281],[-43,-285],[-20,-208],[-22,-317],[-4,2],[-6,154],[10,72],[10,217],[-16,43],[-21,7],[-17,68],[20,27],[44,-42],[21,152],[-6,43],[-2,126],[35,59],[-7,28],[-44,52],[18,155],[-4,62],[-29,30],[44,83],[21,16],[20,-35],[39,-13],[8,15],[7,100],[26,5],[21,-32],[6,-49]],[[36893,101836],[-65,-23],[-31,-58],[36,-42],[-1,-23],[-67,-46],[-27,9],[-75,-50],[-30,-90],[8,-107],[33,-39],[-33,-84],[-23,37],[-4,47],[-26,31],[-65,14],[-53,-7],[-98,-68],[-39,-18],[-22,16],[103,117],[-9,18],[-62,-6],[-9,30],[69,69],[-18,84],[29,42],[21,57],[-24,47],[29,117],[36,43],[32,7],[65,-29],[85,37],[10,50],[47,8],[111,0],[18,-52],[69,-62],[7,-28],[-27,-48]],[[33151,102581],[-27,0],[-102,-41],[-253,-72],[-110,-54],[-118,-82],[-75,-87],[-21,-48],[41,-38],[-24,-42],[-153,0],[-44,-16],[-23,-64],[-30,-31],[-22,109],[-39,-6],[-28,-27],[-27,-66],[-41,23],[5,64],[-96,82],[-23,49],[12,31],[77,-13],[14,-37],[62,27],[86,69],[152,190],[76,65],[76,43],[186,131],[94,56],[57,52],[-14,30],[56,55],[93,62],[87,26],[80,2],[52,-12],[80,-83],[121,12],[19,-7],[-39,-90],[67,-19],[123,45],[186,36],[156,-36],[231,-37],[187,-44],[247,-20],[121,-7],[344,-27],[-12,-36],[-391,10],[-289,20],[-145,19],[-80,-11],[-21,36],[-63,5],[-61,-113],[-71,-11],[-148,-40],[-315,-38],[-112,-5],[-95,-16],[-76,27]],[[54407,95121],[-51,-39],[-74,-100],[-202,-190],[38,92],[78,80],[42,86],[35,34],[-26,49],[-70,-77],[-60,-77],[-104,-148],[-33,-11],[-99,-214],[-11,-89],[18,-102],[-32,14],[3,-60],[22,-98],[-11,-23],[-50,-47],[-61,-38],[-40,-73],[-44,-4],[-28,35],[59,15],[50,81],[79,37],[30,29],[-31,44],[-11,53],[12,26],[23,242],[-28,-47],[-23,-70],[-3,-53],[-42,19],[-60,-74],[1,46],[88,256],[63,155],[-39,9],[12,78],[146,91],[154,151],[103,75],[80,66],[48,-12],[-33,-64],[88,42],[46,6],[6,-24],[-29,-29],[-127,-91],[-14,-27],[29,-22],[135,84],[51,10],[4,-35],[-107,-67]],[[100488,96446],[2,-73],[-44,-54],[-175,-164],[-105,37],[-162,49],[-14,19],[31,81],[-2,161],[43,3],[61,41],[3,56],[-48,52],[-28,47],[-29,82],[117,97],[64,35],[12,-22],[13,-156],[74,-27],[65,5],[5,-32],[-24,-33],[-7,-91],[24,-51],[106,-19],[18,-43]],[[100449,96355],[20,34],[-22,58],[-81,22],[-69,-28],[-95,-63],[-122,-105],[43,-41],[106,-34],[44,4],[71,69],[70,49],[35,35]],[[119685,101629],[-174,13],[-15,-20],[47,-61],[-13,-19],[-61,24],[-30,28],[-5,120],[-18,35],[-42,11],[-119,49],[-28,-9],[-228,228],[-53,73],[-58,41],[-95,35],[8,29],[-11,115],[92,-83],[67,-31],[-7,-53],[29,-34],[46,-20],[15,-55],[60,18],[-26,50],[85,-8],[133,-63],[47,-8],[12,87],[-63,35],[-25,45],[46,24],[-17,108],[-45,-3],[-28,36],[-54,27],[-79,-6],[-36,-19],[-47,-51],[-20,28],[19,65],[102,76],[13,37],[-30,86],[116,-2],[85,20],[23,-45],[52,-24],[37,10],[2,84],[65,-39],[66,-8],[-44,-56],[-54,-28],[-40,2],[-153,35],[-30,-25],[153,-113],[180,-194],[190,-187],[139,-48],[-11,-51],[21,-78],[-79,-25],[-34,-75],[3,-108],[-23,0],[-51,67],[-42,-5],[46,-64],[-11,-23]],[[142865,87123],[-6,-23],[-78,-13],[-168,-114],[-10,-58],[-81,-118],[-70,-37],[-47,-17],[-222,-13],[-66,4],[-60,16],[-200,95],[-86,26],[-99,8],[-21,57],[0,67],[132,26],[172,72],[128,44],[218,53],[115,-17],[87,47],[223,24],[108,6],[-76,-86],[11,-32],[96,-17]],[[42033,93464],[2,-98],[-36,-5],[-42,63],[-24,59],[-33,13],[14,60],[-65,11],[-164,-36],[-58,-33],[-19,-65],[80,13],[37,-52],[67,-62],[64,32],[30,-12],[37,-68],[82,-13],[-11,-38],[-46,-28],[-40,17],[-35,-41],[-124,-117],[-103,-76],[-62,4],[-57,54],[-50,36],[-21,-5],[-83,-58],[-55,-10],[-23,39],[-2,40],[28,54],[-5,72],[9,29],[37,28],[126,27],[28,42],[-4,75],[-39,20],[-4,25],[43,29],[70,104],[76,47],[27,39],[-25,21],[-63,-24],[-67,-69],[-115,-59],[-18,8],[13,67],[-21,74],[31,16],[81,-17],[20,-42],[50,34],[-17,47],[27,26],[138,31],[17,37],[-10,41],[-60,136],[53,33],[29,-53],[5,-74],[46,-107],[61,-13],[21,-36],[-13,-43],[11,-37],[-32,-34],[39,-45],[76,-17],[23,-45],[-6,-20],[54,-51]],[[41653,93294],[-51,2],[33,-95],[44,-1],[37,-51],[38,5],[13,72],[35,57],[-40,33],[-63,-22],[-46,0]],[[45555,94014],[-67,-46],[22,-41],[-13,-78],[10,-68],[-12,-72],[-19,-51],[5,-68],[-13,-39],[-73,-6],[-37,-24],[-22,-42],[-20,4],[-32,65],[-40,-1],[-35,-32],[-88,29],[-7,28],[44,92],[70,19],[-11,37],[-79,-15],[-24,-59],[-49,-44],[-53,-18],[-24,26],[-10,50],[57,-22],[8,55],[-13,71],[-1,57],[-45,-38],[-22,11],[-28,46],[14,45],[38,-35],[37,32],[-3,56],[33,49],[-9,85],[7,31],[72,146],[76,60],[52,-36],[50,-1],[64,-70],[45,-3],[-36,62],[-43,37],[18,24],[89,-73],[35,-73],[-20,-22],[-10,-49],[6,-61],[46,39],[60,-69]],[[45249,94071],[-45,20],[-14,-88],[-54,-66],[18,-27],[74,-13],[-27,57],[6,50],[42,67]],[[45348,93721],[-5,48],[-56,-36],[-1,-46],[54,5],[8,29]],[[37544,101444],[105,9],[1,-97],[-51,-55],[-61,-104],[-55,-63],[-128,-89],[-34,-53],[7,-56],[121,-58],[30,-32],[-17,-19],[-54,-11],[-45,-44],[-9,-72],[51,-84],[-97,-103],[-29,-47],[-79,-24],[-16,-27],[-34,-105],[-40,-83],[-245,-126],[-54,-97],[-39,9],[-25,38],[-29,-6],[-38,-53],[-23,-74],[-28,24],[18,161],[-24,19],[3,44],[30,29],[29,-17],[23,60],[14,135],[14,82],[21,-3],[-6,-161],[17,-101],[43,-7],[68,114],[75,31],[33,53],[30,12],[32,41],[16,93],[-7,65],[-17,44],[-93,51],[29,102],[-56,69],[-76,11],[-8,65],[25,80],[114,106],[13,32],[-11,84],[99,40],[49,63],[137,124],[33,58],[-1,42],[57,65],[36,-18],[-13,-80],[25,-36],[34,25],[26,59],[38,-31],[-25,-115],[18,-51],[48,28],[5,35]],[[52740,87718],[-29,-42],[-94,-4],[-79,58],[20,23],[117,-14],[65,-21]],[[6978,103182],[-81,-49],[-71,-28],[-82,-21],[-96,-35],[-17,-33],[5,-41],[-16,-70],[-18,-21],[-172,-18],[-489,-91],[-81,40],[-18,74],[39,76],[343,116],[156,41],[389,80],[209,-20]],[[5814,101406],[-35,-30],[16,-69],[-27,-55],[-54,-61],[10,102],[-27,45],[-260,77],[-42,8],[-91,41],[-30,63],[40,63],[104,18],[210,-74],[177,-91],[9,-37]],[[26886,83791],[-14,-63],[-39,-11],[-28,24],[-19,84],[-15,28],[4,49],[43,73],[72,31],[10,-44],[-14,-171]],[[117722,33650],[-32,-106],[-12,-6],[-30,54],[0,95],[6,45],[34,14],[34,-96]],[[146500,92115],[-52,-68],[-18,-38],[-23,-79],[-35,-31],[-32,7],[-144,128],[-134,101],[-83,57],[-106,34],[-108,43],[-103,78],[-33,32],[7,35],[84,-1],[55,12],[93,8],[37,33],[20,65],[-82,117],[-3,50],[18,41],[3,136],[28,108],[70,97],[92,72],[73,22],[96,45],[49,45],[-12,74],[-147,167],[-110,25],[-71,82],[-35,6],[-67,82],[15,25],[128,-43],[76,23],[88,84],[42,-17],[-84,-64],[-52,-92],[327,-279],[42,-71],[-18,-34],[-105,-36],[-167,-86],[-51,-18],[-98,-61],[-56,-90],[-15,-89],[17,-95],[56,-168],[94,-25],[85,-12],[-7,-22],[-132,-19],[-50,-47],[-9,-34],[93,-65],[83,-42],[100,-64],[165,-57],[78,-50],[28,-37]],[[45383,88505],[-5,-121],[-57,-122],[-20,5],[-15,118],[-17,84],[21,76],[39,79],[56,-14],[-2,-105]],[[55262,109563],[134,248],[9,77],[-19,52],[31,27],[101,-8],[102,-24],[52,-67],[63,-51],[36,-69],[62,-38],[64,-16],[38,11],[201,-57],[122,-56],[144,-44],[14,-31],[-36,-42],[-69,-25],[0,-32],[130,-11],[100,20],[19,-42],[-70,-66],[-77,-43],[-66,-81],[-58,38],[-50,-38],[-71,48],[-50,4],[-122,-36],[-113,-23],[-120,88],[-21,-7],[-54,-73],[-82,-45],[-148,-154],[-71,-35],[-56,31],[58,35],[25,55],[-9,63],[29,60],[10,109],[-60,81],[-44,91],[-78,76]],[[118517,49687],[-30,-106],[-36,-18],[-21,40],[-65,267],[-54,131],[-52,86],[-25,7],[-27,61],[-24,120],[-86,188],[-17,107],[-38,108],[-14,166],[5,94],[19,151],[18,218],[4,106],[14,158],[26,11],[29,-24],[17,15],[14,63],[20,28],[25,-20],[13,-188],[11,-264],[-1,-88],[7,-163],[-3,-380],[53,-249],[68,-182],[65,-82],[66,-62],[26,-123],[-7,-176]],[[37421,106408],[74,43],[35,65],[35,31],[62,26],[20,55],[54,-13],[53,-43],[42,17],[-23,48],[33,22],[82,-56],[75,31],[62,63],[115,15],[35,-73],[-16,-89],[-87,-71],[-95,6],[-33,-11],[-43,-112],[6,-66],[26,-39],[103,-49],[0,-48],[-85,-56],[-18,-36],[22,-51],[-4,-39],[-64,-45],[-89,-44],[-53,16],[-7,35],[41,101],[-4,58],[-50,19],[-11,-62],[-38,-98],[-53,34],[-97,23],[-40,75],[-97,5],[-27,36],[-8,74],[-28,68],[-45,42],[8,40],[47,27],[85,26]],[[41697,107285],[221,-34],[13,-44],[-50,-59],[-84,-15],[-99,32],[-138,32],[-258,13],[-118,-20],[-165,-1],[-95,66],[-164,62],[-27,101],[41,21],[557,-47],[125,-33],[143,-20],[83,-19],[15,-35]],[[108618,98195],[-27,45]],[[108591,98240],[94,99],[37,54],[28,60],[56,-35],[72,10],[115,-2],[-176,-93],[-102,-61],[-97,-77]],[[108618,98195],[-14,-21],[-151,-104],[-79,-31],[-56,40],[17,26],[131,40],[66,42],[59,53]],[[54130,88580],[-19,-88],[21,-167],[-19,-9],[-28,204],[4,102],[52,83],[-13,104],[-13,35],[-26,194],[26,128],[9,151],[22,76],[13,11],[16,-108],[24,27],[30,73],[34,49],[34,-25],[-13,-35],[-48,-49],[12,-35],[-13,-199],[-11,-112],[-29,-42],[20,-45],[0,-45],[-25,-46],[-10,-71],[-50,-161]],[[29419,78417],[-2,-49],[-48,-79],[-53,8],[-72,146],[-88,121],[5,66],[25,18],[36,-3],[38,-63],[131,-122],[28,-43]],[[28421,106396],[-24,-58],[-114,-31],[-108,0],[-27,18],[-188,85],[-56,55],[-155,26],[-94,25],[-1,36],[38,49],[129,80],[142,53],[59,-24],[50,-73],[88,-58],[111,-30],[82,-69],[68,-84]],[[54367,8725],[-12,54],[-3,73],[12,58],[36,63],[37,43],[43,-14],[42,-34],[38,-69],[23,-79],[-18,-28],[-57,13],[-87,-88],[-54,8]],[[166444,78321],[-19,-29],[-48,-134],[-35,-23],[-97,-45],[-22,-27],[-33,9],[12,33],[-29,42],[-32,-26],[-2,44],[27,21],[56,-35],[5,52],[-30,87],[-120,3],[-37,-15],[-17,29],[15,22],[47,6],[71,25],[36,-3],[30,-64],[29,-91],[38,-48],[39,15],[37,111],[-9,30],[-77,20],[-43,37],[13,43],[27,1],[43,65],[32,64],[72,-32],[22,4],[-35,56],[-35,97],[3,77],[30,21],[16,-24],[21,-98],[33,-48],[65,-19],[50,-81],[-48,-168],[-37,-81],[-22,32],[-15,47],[-27,-2]],[[156019,95356],[-9,-97],[-23,-145],[-28,-111],[3,-34],[-106,-163],[-36,-140],[16,-31],[-17,-54],[-43,-28],[-104,-42],[-8,36],[8,58],[20,56],[12,83],[29,342],[-7,133],[42,94],[26,38],[-13,70],[106,134],[112,-153],[20,-46]],[[151064,94332],[-47,104],[65,58],[259,133],[60,38],[43,-18],[15,-29],[-8,-64],[58,-104],[88,-60],[15,-34],[42,-17],[21,-26],[-35,-24],[-50,-4],[-70,-49],[-107,-133],[-112,-43],[-70,-5],[-60,45],[-68,83],[-27,80],[-12,69]],[[114097,39218],[-12,-73],[-51,-170],[-75,-158],[-72,-169],[-43,-116],[-61,-10],[-41,3],[-10,-35],[3,-41],[-21,-19],[-56,3],[-26,49],[-5,81],[0,118],[25,127],[161,312],[131,238],[51,-7],[72,-75],[30,-58]],[[174754,89425],[-39,-102],[-76,-216],[-90,-93],[-55,-77],[-52,36],[-20,58],[-101,105],[23,211],[-55,99],[-6,130],[39,69],[93,38],[100,1],[104,-31],[85,-62],[49,-75],[1,-91]],[[158233,58930],[-62,74],[-37,3],[-84,51],[-68,101],[-21,109],[-57,15],[-66,59],[-52,74],[-31,86],[-32,157],[19,41],[65,-15],[81,-52],[69,-134],[59,-178],[21,-76],[19,-23],[59,8],[42,-22],[16,-79],[2,-75],[49,-42],[9,-82]],[[156000,81722],[-18,-62],[0,-67],[-11,-78],[-45,-31],[-61,23],[-54,56],[-119,4],[-157,31],[-77,57],[-74,106],[-29,55],[17,33],[80,8],[10,31],[-17,73],[-10,77],[28,32],[70,32],[87,-15],[100,-48],[101,-117],[81,-114],[98,-86]],[[138733,110865],[-102,6],[-55,15],[-24,49],[20,32],[76,4],[168,-73],[-22,-35],[-61,2]],[[122365,78840],[-38,60],[-53,116],[-18,184],[-5,136],[10,53],[19,13],[72,-175],[73,-92],[61,-87],[15,-58],[-10,-75],[1,-73],[-34,-22],[-93,20]],[[122727,77868],[-10,-71],[-41,-71],[-31,-4],[-42,20],[-108,84],[-27,68],[-8,89],[5,55],[25,47],[-4,85],[24,46],[24,19],[19,-38],[81,-216],[93,-113]],[[117874,76967],[-23,-341],[-14,-30],[-58,-1],[-11,21],[-1,135],[12,118],[24,86],[34,42],[22,9],[15,-39]],[[113143,101700],[41,28]],[[113184,101728],[14,-37],[71,-9],[89,8],[63,-38],[73,-83],[30,-65],[12,-56],[-40,-51],[-64,7],[-47,27],[-17,38],[-47,50],[-72,26],[-32,-4],[-74,159]],[[54902,92604],[-65,27],[-45,87],[-37,42],[-15,45],[24,61],[77,21],[93,-30],[90,-60],[41,-46],[-3,-86],[-24,-40],[-60,-30],[-76,9]],[[51975,87175],[-29,-116],[-56,-78],[32,123],[20,0],[21,75],[12,-4]],[[52079,86962],[-18,39],[-4,124],[-22,130],[0,130],[43,-174],[1,-249]],[[52289,87041],[-86,95],[-62,163],[13,139],[30,-157],[20,-84],[85,-156]],[[99797,57004],[-51,-159],[1,-77],[10,-103],[25,-111],[-39,-17],[-16,34],[-17,117],[-40,66],[25,45],[-3,25],[-25,21],[13,31],[-11,71],[-39,66],[0,26],[22,13],[13,40],[6,155],[13,17],[44,-39],[22,13],[30,-8],[14,-22],[-18,-63],[10,-100],[11,-41]],[[123772,82593],[47,-39],[5,-115],[60,-142],[-15,-80],[-99,-49],[-29,-72],[-45,-71],[-51,13],[-49,73],[-43,197],[-20,166],[-2,65],[-11,74],[-87,113],[-4,96],[45,54],[18,58],[-9,34],[-91,1],[6,45],[46,67],[147,0],[59,-86],[18,-93],[-24,-131],[-7,-118],[28,-50],[107,-10]],[[34891,101244],[87,-39],[-2,-55],[16,-57],[-5,-66],[-26,-47],[-67,-37],[-70,-5],[-98,-40],[-39,10],[-53,-9],[-62,26],[-47,-32],[-50,-63],[-38,-23],[-32,6],[10,84],[-6,33],[32,61],[-24,46],[85,-12],[54,26],[103,101],[137,47],[61,53],[34,-8]],[[38783,106163],[151,200],[27,-17],[-18,-52],[17,-48],[85,-54],[10,-57],[-58,-52],[-85,-53],[-110,-9],[-91,11],[-20,51],[92,80]],[[37061,103699],[38,-8],[62,50],[37,53],[-30,17],[-41,54],[1,79],[49,45],[53,15],[27,-32],[33,-96],[29,-52],[65,-77],[65,-127],[22,-118],[-21,-121],[-57,25],[-24,-75],[-36,83],[-111,180],[-62,37],[-56,-13],[-65,9],[-7,59],[29,13]],[[130909,86361],[-28,-69],[-78,-60],[-92,-38],[-106,-21],[-58,75],[-42,76],[-66,-9],[-5,46],[24,182],[52,203],[52,110],[71,29],[92,-69],[45,-104],[-9,-41],[11,-37],[30,-52],[10,-91],[85,-70],[12,-60]],[[141955,98567],[72,47],[116,34],[39,4],[17,-34],[48,-29],[28,32],[19,73],[61,49],[82,30],[74,-5],[16,-21],[90,-15],[90,-32],[10,-44],[-65,-34],[-94,-30],[-78,-53],[-4,-31],[33,-26],[28,-75],[-28,-26],[-55,-18],[-15,30],[-81,35],[-47,9],[-64,28],[-113,8],[-133,-29],[-53,15],[11,65],[-4,43]],[[145106,90091],[-224,101],[-114,87],[-75,84],[-13,90],[3,78],[-6,65],[26,81],[51,57],[70,3],[65,-23],[46,-33],[42,-42],[38,-67],[28,-88],[7,-98],[16,-127],[40,-168]],[[44154,75612],[-24,-71],[-41,-18],[-22,17],[-3,36],[13,45],[37,60],[36,-20],[4,-49]],[[53770,99682],[-48,-47],[-66,-17],[-89,63],[-133,40],[-99,22],[-39,43],[-15,72],[4,56],[27,32],[63,26],[109,-17],[97,-52],[37,-94],[52,-17],[62,-6],[40,-37],[-2,-67]],[[49992,90575],[59,3],[55,-34],[62,8],[16,19],[-4,43],[42,20],[115,-5],[119,17],[93,-44],[35,-71],[-30,-48],[-50,-17],[-117,45],[-61,-19],[-65,-30],[-16,66],[-28,11],[-136,-17],[-89,53]],[[165639,92728],[-67,-23],[-15,17],[10,47],[-6,70],[-48,76],[-15,52],[41,61],[98,90],[95,127],[82,123],[62,67],[42,-19],[-12,-94],[18,-38],[35,-18],[21,-34],[51,-13],[63,-5],[16,-32],[-41,-34],[-52,-1],[-30,24],[-31,-1],[-22,-27],[-10,-40],[-118,-131],[-109,-143],[-11,-54],[-47,-47]],[[114604,36948],[-34,-117],[-11,-79],[-48,-74],[-49,-34],[-49,1],[-47,20],[-22,136],[-23,16],[-9,34],[40,36],[104,211],[18,-2],[43,43],[48,8],[42,-56],[11,-71],[-14,-72]],[[113698,112114],[109,36],[-59,-101],[20,-32],[-69,-65],[-28,-66],[-88,-17],[-93,-52],[1,-26],[100,-55],[55,-46],[-37,-19],[-178,26],[-58,-17],[59,-69],[-71,-42],[-53,46],[-73,-56],[-112,27],[-119,-23],[98,61],[-7,37],[-52,40],[-91,20],[-5,17],[188,75],[45,55],[108,62],[114,46],[84,11],[-55,60],[65,42],[22,-58],[55,11],[65,69],[60,3]],[[112723,107502],[30,45],[-3,62],[24,22],[143,-6],[60,-25],[2,-77],[26,-42],[60,-36],[62,8],[63,38],[87,16],[116,-14],[2,-70],[-75,-3],[-36,-21],[15,-34],[230,-80],[4,-36],[-53,0],[-85,22],[-161,75],[-53,46],[-33,-4],[-50,-38],[-57,-18],[-24,-34],[3,-40],[-202,101],[-12,36],[-75,79],[-8,28]],[[113143,101700],[-42,111],[-60,92],[-51,47],[-47,183],[-111,148],[5,73],[78,66],[111,38],[132,11],[93,-25],[39,-64],[38,-318],[4,-81],[-139,-167],[-9,-86]],[[114137,45778],[31,-86],[23,-116],[40,-146],[-11,-39],[-33,-34],[-24,-63],[-54,-108],[-43,-63],[-48,-14],[-48,-68],[-35,-14],[3,71],[35,79],[-8,44],[-16,12],[2,70],[22,65],[-4,70],[45,159],[27,178],[23,39],[73,-36]],[[114267,46726],[-31,-9],[-18,58],[-5,45],[19,72],[80,254],[86,57],[82,-1],[61,-49],[-11,-121],[-132,-192],[-78,-97],[-33,5],[-20,-22]],[[167612,76372],[-121,-109],[-93,-32],[-78,35],[-73,110],[-23,151],[28,130],[71,93],[64,38],[22,-21],[10,-60],[63,-15],[-4,-41],[-24,-34],[-3,-28],[38,-86],[23,-25],[5,-44],[-26,-49],[39,-19],[46,43],[36,-37]],[[52274,99811],[75,133],[12,37],[-16,41],[20,103],[29,29],[21,-59],[42,-65],[-25,-18],[42,-53],[155,-122],[-8,-43],[-113,-17],[-176,11],[-58,23]],[[60649,18548],[-104,-60],[-95,-22],[-91,-50],[-79,108],[-41,48],[-28,55],[29,40],[61,24],[61,37],[45,59],[127,-36],[32,10],[50,-38],[31,-77],[2,-98]],[[116459,93133],[-2,-24],[-50,10],[-20,55],[-36,19],[-69,4],[-55,27],[-87,77],[-124,96],[-67,42],[-72,54],[-158,102],[-8,39],[40,14],[93,-27],[138,-67],[65,-40],[99,-42],[30,18],[-13,60],[3,68],[28,20],[24,-69],[62,-73],[91,-118],[105,-152],[-17,-93]],[[52542,97477],[-7,66]],[[52535,97543],[37,0],[62,27],[11,44],[36,39],[175,35],[152,4],[30,-50],[33,-21],[225,68],[39,2],[5,-36],[-56,-17],[-41,-41],[-93,-39],[-129,-41],[2,-23],[189,-12],[100,-17],[104,13],[81,32],[-3,-36],[-128,-51],[-75,-4],[-75,-25],[-53,-1],[16,43],[-20,15],[-51,-15],[-68,-4],[-232,18],[-14,-16],[33,-36],[-19,-41],[-48,-38],[-128,-30],[-98,9],[-13,45],[35,43],[70,51],[22,27],[-22,31],[-84,-18]],[[52542,97477],[-159,-34],[-14,-32],[-76,-36],[-190,-54],[-204,-6],[-21,20],[-102,6],[-53,37],[-3,43],[17,30],[-92,19],[-37,19],[43,40],[112,71],[46,12],[31,-19],[159,29],[67,24],[79,58],[21,63],[-92,75],[28,29],[106,-58],[33,27],[55,-7],[21,33],[-59,54],[9,16],[111,28],[61,-48],[67,6],[-20,-36],[-42,-28],[-29,-104],[-75,-61],[-32,-43],[-63,-59],[-90,-61],[31,-22],[83,12],[266,23]],[[5470,101309],[-32,-66],[-102,-121],[5,-55],[-38,-60],[-108,65],[7,50],[67,13],[30,29],[2,40],[90,95],[55,28],[24,-18]],[[55221,64579],[-16,-34],[-41,-12],[-109,68],[9,35],[55,19],[58,-35],[44,-41]],[[3647,109448],[1,-38],[-81,-3],[-84,11],[-111,73],[12,22],[112,-7],[114,-22],[37,-36]],[[29741,103705],[-55,-126],[-76,-32],[-5,-48],[-107,27],[-75,33],[-160,-10],[-10,36],[54,50],[160,7],[106,30],[75,40],[65,12],[28,-19]],[[52221,97097],[-32,-45],[2,-39],[-56,-22],[-35,-48],[3,-57],[23,-50],[-32,-12],[-27,25],[-22,116],[33,91],[-6,120],[28,47],[72,69],[10,-46],[25,-2],[56,63],[30,-16],[10,-50],[-15,-32],[-67,-18],[-21,-19],[21,-75]],[[54944,370],[-10,-19],[-66,-19],[-95,-35],[-97,-5],[-90,14],[-45,14],[-72,-12],[-104,-62],[-19,-58],[25,-48],[98,9],[7,-21],[-40,-52],[-41,-68],[-20,-8],[-32,98],[-58,140],[-35,31],[-62,-12],[-6,13],[41,58],[109,-19],[33,26],[-6,33],[-106,73],[-80,-34],[13,59],[-5,48],[58,82],[30,-45],[43,-96],[68,-44],[85,-15],[41,9],[170,13],[190,38],[44,-28],[34,-58]],[[149911,75846],[-19,46],[-12,116],[-2,139],[37,32],[132,-48],[43,10],[36,51],[89,-8],[46,36],[19,-2],[21,-32],[1,-46],[-48,-42],[-52,-71],[-76,-71],[-68,-43],[-147,-67]],[[116630,83632],[39,-39],[43,-13],[44,-81],[-33,-125],[-49,-52],[-64,34],[-110,101],[-46,98],[36,98],[44,85],[9,84],[-2,65],[26,-6],[43,-60],[18,-60],[-9,-50],[11,-79]],[[35067,103395],[-35,-4],[-32,-32],[-5,-63],[15,-44],[-28,-33],[-37,-1],[-23,48],[-42,-43],[-57,16],[29,97],[-4,82],[-39,20],[15,41],[39,12],[46,-61],[28,0],[34,102],[61,45],[85,34],[98,-42],[-38,-50],[-151,-18],[-24,-26],[14,-31],[62,9],[35,-20],[2,-30],[-48,-8]],[[113165,104274],[-113,59],[15,30],[47,-37],[77,-9],[65,26],[19,62],[50,-32],[24,-56],[-62,-32],[-122,-11]],[[114298,106079],[53,-33],[12,51],[-27,82],[34,5],[89,-99],[44,-75],[-6,-59],[-85,-47],[-25,54],[-49,32],[-48,55],[8,34]],[[114156,105926],[-6,63],[43,30],[57,-79],[44,13],[34,-27],[-41,-61],[-67,17],[-70,-8],[-19,32],[25,20]],[[190942,112043],[-28,-20],[-103,9],[-170,41],[-92,61],[20,29],[91,-20],[191,-57],[91,-43]],[[138804,110753],[-36,19],[41,53],[107,54],[91,31],[46,-42],[-6,-40],[-73,-51],[-170,-24]],[[138294,112830],[14,-34],[-83,-38],[-81,-4],[-61,61],[27,56],[113,-18],[71,-23]],[[138506,112715],[-60,61],[3,31],[49,64],[61,3],[1,-79],[67,-62],[-18,-32],[-103,14]],[[138338,112859],[-112,28],[-16,44],[249,15],[29,-24],[-72,-51],[-78,-12]],[[117540,106674],[77,-22],[134,9],[4,-54],[-165,-48],[-133,76],[-25,-28],[33,-131],[40,-105],[-4,-139],[-28,20],[-37,207],[-18,47],[-60,68],[-88,56],[-58,154],[-56,68],[-83,-19],[77,164],[47,-3],[127,-55],[-17,-21],[141,-43],[3,-30],[75,-97],[-22,-42],[36,-32]],[[116621,106616],[84,41],[30,-11],[35,-62],[88,27],[50,-19],[34,-113],[44,-64],[-29,-26],[-70,-12],[-114,2],[-119,24],[-113,64],[-19,30],[6,61],[-15,39],[57,1],[51,18]],[[114960,106086],[-9,57],[51,33],[47,71],[44,-21],[10,-67],[-62,-85],[-31,-16],[-32,-67],[-78,4],[60,91]],[[114603,106475],[-20,-62],[21,-81],[84,-156],[112,-135],[-59,-22],[-51,-78],[-30,7],[74,130],[-47,57],[-61,28],[-104,120],[-24,54],[-129,12],[-38,20],[-79,3],[-56,68],[-35,64],[-114,141],[29,92],[55,-55],[199,-129],[151,-77],[77,20],[45,-21]],[[115292,109105],[62,-142],[10,-47],[-82,9],[-119,47],[-175,82],[42,126],[38,40],[64,31],[80,6],[41,-23],[41,-100],[-2,-29]],[[116079,108568],[86,-51],[29,-58],[-46,2],[-89,48],[-122,34],[-162,59],[-31,50],[-123,45],[-94,60],[-38,57],[108,57],[170,31],[40,-12],[39,-40],[13,-46],[110,-79],[19,-31],[10,-83],[28,-34],[53,-9]],[[54441,123662],[290,44],[1617,178],[62,-29],[-315,-44],[-1290,-133],[-289,-26],[-75,10]],[[44682,96745],[52,32],[7,29],[82,53],[93,-37],[100,20],[-54,-109],[-62,-16],[-191,14],[-27,14]],[[148173,86433],[-68,40],[-91,7],[-122,-29],[-83,-42],[-14,23],[4,50],[35,67],[19,91],[36,33],[250,-85],[34,-102],[0,-53]],[[115115,105242],[-11,102],[71,-17],[48,-89],[-19,-38],[-48,-7],[-41,49]],[[115503,101743],[-50,-23],[-30,-59],[-140,27],[-164,40],[-11,23],[95,75],[89,37],[74,110],[45,-44],[80,-38],[73,-1],[6,-74],[25,-99],[-92,26]],[[165333,74549],[1,-32],[-43,-1],[-67,39],[-15,-19],[-10,-75],[15,-84],[16,-54],[3,-46],[-14,-39],[-25,1],[-35,94],[-36,0],[-11,52],[-135,135],[-29,79],[-49,113],[36,133],[38,96],[15,91],[-20,46],[18,22],[49,24],[22,-80],[-16,-20],[-5,-38],[-23,-80],[-9,-61],[-20,-25],[19,-38],[32,-23],[-35,-49],[9,-52],[126,-29],[-4,-64],[27,21],[19,56],[17,13],[100,-12],[14,-10],[6,-63],[19,-21]],[[165189,75407],[35,-62],[-13,-22],[-64,1],[-58,26],[-18,51],[13,29],[27,7],[35,-16],[32,7],[11,-21]],[[165034,75391],[27,-5],[11,-52],[-4,-39],[-25,-25],[-97,-25],[-39,66],[4,60],[78,-15],[19,57],[26,-22]],[[163954,75539],[26,38],[7,29],[52,20],[4,20],[-26,51],[26,1],[59,-62],[12,-39],[11,-117],[-11,-34],[-47,86],[-20,5],[-20,-33],[-17,-69],[-21,21],[-35,83]],[[149310,77033],[-18,-6],[-45,-102],[-43,-72],[-46,-6],[-66,27],[1,71],[-12,26],[-45,-2],[-50,-13],[-47,17],[-69,106],[25,34],[55,-12],[33,-31],[32,7],[30,43],[60,29],[89,-5],[86,-38],[33,-45],[-3,-28]],[[151339,92347],[-96,-75],[-105,-63],[-52,-61],[-43,-107],[-47,-77],[-68,45],[-55,109],[-8,72],[75,42],[85,70],[0,125],[22,83],[80,7],[35,-19],[-48,-96],[2,-75],[44,-17],[37,0],[58,51],[-2,128],[55,-18],[40,-54],[-9,-70]],[[200060,11037],[-37,-50],[-44,-34],[-77,-43],[-24,4],[26,92],[-10,46],[-19,22],[8,43],[27,46],[79,0],[66,-41],[15,-34],[-10,-51]],[[149873,85340],[39,-187],[18,-116],[38,-119],[-37,-38],[-30,11],[-31,112],[0,111],[-63,201],[-46,68],[31,53],[43,39],[18,-13],[20,-122]],[[101078,90715],[-17,-44],[-106,16],[-115,-13],[-77,-25],[-105,-146],[-14,0],[15,91],[108,145],[143,66],[135,-58],[33,-32]],[[93260,98464],[100,-13],[27,-45],[-6,-53],[-42,-95],[-43,-2],[-80,14],[-13,43],[57,151]],[[105792,102260],[5,-26],[-29,-153],[42,-49],[-10,-26],[-152,-109],[-38,-90],[-23,-91],[-48,-108],[-94,-143],[-56,-115],[-50,-13],[-22,120],[54,211],[104,240],[68,144],[73,111],[176,163],[21,-5],[-21,-61]],[[122206,83156],[-39,42],[-24,81],[-15,7],[-44,-24],[-57,-2],[-84,31],[-109,10],[7,45],[55,88],[-1,46],[40,50],[115,35],[72,11],[52,0],[85,14],[49,35],[-3,21],[-42,11],[-7,33],[45,29],[166,38],[73,-10],[29,-28],[-46,-46],[-87,-58],[-53,-22],[-48,-90],[-14,-106],[59,-51],[19,-56],[-49,-75],[-52,-67],[-35,-10],[-57,18]],[[32586,88755],[-43,-15],[-40,59],[-31,9],[2,31],[-21,20],[-61,-30],[-15,30],[16,23],[46,9],[40,71],[32,9],[33,-21],[4,-91],[36,-52],[2,-52]],[[7581,113276],[-196,-47],[-222,-36],[-81,6],[-132,52],[-48,53],[-91,58],[129,22],[131,38],[117,11],[350,-125],[43,-32]],[[30271,92118],[-22,127],[35,21],[51,-71],[13,-76],[12,-140],[-31,-68],[-54,2],[-8,25],[22,41],[-13,26],[-60,23],[-2,20],[35,30],[22,40]],[[30575,80806],[-20,11],[-72,156],[-59,65],[-55,80],[-29,-1],[-43,-123],[-57,-25],[-50,23],[-15,44],[-31,22],[-45,-20],[-32,25],[-79,0],[-43,-77],[-32,-3],[-25,51],[-3,44],[95,30],[83,37],[47,41],[26,229],[30,35],[1,-71],[-14,-102],[-8,-117],[45,-93],[40,-57],[39,26],[76,157],[19,-1],[128,-183],[46,-84],[39,-92],[-2,-27]],[[31639,102051],[6,-40],[-21,-86],[32,-88],[-12,-26],[-52,3],[-107,38],[-63,62],[-28,43],[17,37],[-2,47],[-26,25],[-35,-12],[-46,14],[-32,37],[10,61],[55,89],[108,43],[65,-5],[3,-30],[-45,-8],[-3,-39],[78,-21],[83,7],[74,-58],[-59,-93]],[[60096,97910],[-16,-52],[58,-66],[-56,-45],[-61,33],[-108,20],[-46,-29],[22,-48],[-71,-4],[-70,-22],[-33,-41],[-107,3],[-67,33],[-8,60],[-22,48],[-187,63],[-40,42],[-51,-10],[-22,-30],[28,-47],[-31,-66],[-34,-25],[-99,-16],[-54,-26],[-96,-66],[-50,-21],[-20,-34],[1,-33],[61,-46],[1,-34],[24,-69],[-9,-23],[21,-37],[90,-43],[-27,-31],[-125,32],[-172,32],[-32,49],[114,16],[42,42],[-38,103],[31,100],[-6,54],[-67,33],[-31,33],[-6,44],[-39,60],[-85,24],[-57,46],[-168,197],[-83,74],[-13,49],[-51,59],[-45,19],[7,-50],[46,-79],[-33,-80],[11,-71],[33,-78],[-29,-14],[-65,103],[-68,200],[-8,126],[51,5],[16,24],[-13,108],[12,31],[-27,76],[-26,36],[19,34],[-18,66],[27,23],[11,43],[110,-113],[-4,-61],[-40,-12],[12,-45],[-17,-37],[41,-92],[3,-26],[90,-103],[5,-94],[98,-122],[17,-5],[33,-65],[39,9],[47,-73],[46,-36],[18,11],[104,-112],[61,-46],[32,34],[89,-46],[52,17],[-43,59],[28,8],[156,-28],[28,21],[-74,71],[-41,67],[-36,109],[-53,89],[-82,23],[15,59],[4,79],[28,40],[51,-12],[-48,-76],[46,-42],[47,-108],[58,-87],[-6,-63],[29,-60],[93,-26],[48,-1],[56,40],[89,1],[58,-19],[67,11],[-38,73],[-24,16],[-82,-19],[-9,20],[32,54],[15,167],[20,98],[-19,120],[45,-44],[15,-36],[-20,-86],[5,-45],[107,-157],[112,-100],[45,-69],[71,-68],[97,-18],[68,19],[55,40],[115,-10]],[[156158,116353],[-43,31],[90,92],[-24,96],[-96,150],[-77,33],[-156,16],[-52,27],[-113,9],[-153,-14],[-145,7],[-73,19],[-45,38],[58,18],[214,-10],[116,29],[36,39],[153,90],[4,49],[-39,86],[-67,73],[-130,56],[21,33],[78,38],[89,11],[26,55],[-104,67],[-142,23],[39,53],[251,-81],[43,-34],[2,-51],[-41,-47],[16,-102],[78,-61],[79,-35],[128,-37],[437,-37],[396,1],[451,32],[140,32],[157,48],[276,49],[-14,41],[-107,25],[-87,40],[20,54],[120,65],[17,54],[69,8],[118,-27],[95,-75],[138,-91],[259,-56],[15,-38],[-59,-28],[-423,-41],[-345,-101],[-104,0],[-101,-37],[-646,-72],[-248,-13],[-32,-82],[51,-86],[68,-77],[-123,-3],[-168,50],[-155,61],[-67,45],[31,71],[-89,11],[-209,-26],[16,-36],[128,-35],[-2,-54],[33,-79],[71,-95],[95,-94],[53,-84],[2,-105],[-239,8],[-59,11]],[[45060,62001],[-54,-71],[-31,-137],[-33,-43],[-60,-25],[-33,-6],[-38,-27],[-50,100],[41,31],[124,75],[79,80],[64,72],[-9,-49]],[[184412,114648],[-4,-27],[-73,-81],[-104,4],[-55,101],[42,112],[54,8],[85,-33],[55,-84]],[[149481,82362],[-56,-9],[-87,17],[-130,32],[-36,45],[184,43],[52,-4],[47,-14],[46,-35],[0,-61],[-20,-14]],[[105296,106443],[-12,62],[91,-19],[101,-40],[69,3],[59,-56],[-25,-40],[-67,-8],[-19,-60],[28,-128],[-5,-109],[-20,-9],[-47,92],[49,51],[-42,101],[-88,60],[-88,39],[16,61]],[[43616,94956],[41,93],[36,-2],[139,-52],[135,102],[83,88],[49,-21],[52,-49],[2,-19],[-59,-48],[-128,-50],[-111,-56],[-107,-8],[-116,-26],[-29,12],[13,36]],[[116078,39801],[4,-48],[-205,119],[-70,-8],[-54,-14],[-22,26],[-28,70],[-30,94],[-16,77],[3,44],[20,32],[37,32],[75,-45],[247,-333],[39,-46]],[[157243,64632],[-31,-1],[-62,-27],[-43,-44],[-18,25],[-26,78],[-26,102],[-12,72],[32,-11],[13,-36],[45,-31],[51,-12],[4,-48],[22,-28],[51,-39]],[[53453,92612],[-9,-21],[-54,-24],[-39,6],[-29,27],[-65,104],[-47,15],[-13,-53],[15,-30],[-10,-115],[-71,-22],[-34,-48],[23,176],[34,34],[20,47],[-6,45],[-30,-15],[-56,-66],[-10,-43],[5,-61],[-56,-6],[-69,33],[-1,38],[63,106],[59,4],[41,69],[40,34],[42,-17],[92,60],[53,25],[96,16],[21,-43],[-34,-108],[15,-30],[40,-26],[44,8],[3,-39],[22,-58],[31,-33],[54,-22],[72,-50],[-80,-24],[-44,28],[-49,59],[-39,29],[-40,-9]],[[166913,77967],[-9,-39],[-54,-4],[-27,-11],[-31,-41],[-16,-6],[-28,105],[-26,53],[29,22],[58,4],[14,57],[-5,66],[41,38],[-7,29],[-39,64],[-33,-5],[-15,-79],[-32,-35],[-78,11],[-49,37],[103,-19],[26,11],[18,70],[21,39],[50,14],[29,-34],[26,-81],[-2,-51],[15,-22],[15,-79],[6,-114]],[[148633,112275],[-9,-14],[-122,-5],[-122,33],[-20,26],[16,78],[-22,76],[-44,97],[-80,57],[-72,18],[36,107],[138,63],[96,24],[66,-7],[26,-27],[-17,-33],[-165,-61],[14,-60],[146,-131],[38,-64],[25,-116],[72,-61]],[[112874,104951],[-59,7],[-44,24],[-54,-27],[64,-38],[6,-54],[-65,7],[-12,34],[-119,136],[-87,48],[-76,53],[79,103],[18,-13],[-38,-86],[24,-6],[63,33],[30,36],[32,-48],[36,2],[48,-25],[20,38],[-65,27],[-31,35],[-3,40],[105,-59],[32,-44],[15,-138],[16,-22],[59,-24],[6,-39]],[[106653,102751],[-55,-32],[-228,-75],[-61,-43],[-67,2],[-45,16],[-21,53],[-41,29],[-97,26],[152,31],[67,-17],[76,5],[78,-16],[115,1],[127,20]],[[112635,106307],[-22,-22],[6,-51],[27,-49],[-11,-29],[-98,59],[-19,139],[79,30],[12,36],[-79,64],[-18,26],[23,32],[96,-41],[65,-60],[4,-28],[-45,-59],[-20,-47]],[[54242,98815],[16,26],[135,50],[162,32],[76,7],[124,47],[14,-42],[-5,-56],[-40,-38],[-67,-40],[-91,-31],[-45,-37],[-199,-59],[-38,0],[-60,19],[-117,-43],[-83,-17],[26,36],[48,26],[-3,35],[-47,26],[13,19],[74,-7],[107,47]],[[118838,103551],[31,64],[67,52],[83,28],[114,10],[95,-27],[48,-57],[14,-74],[-4,-53],[-17,-16],[-148,-8],[-144,9],[-111,41],[-28,31]],[[64941,51785],[-27,-53],[10,-49],[34,-17],[-5,-51],[-34,-18],[-30,26],[-38,-37],[-23,-75],[-21,-25],[-56,28],[-24,31],[-44,2],[3,43],[30,44],[-4,73],[66,168],[13,80],[14,32],[19,0],[44,-26],[30,-49],[18,-14],[57,8],[25,-29],[-24,-86],[-33,-6]],[[117968,43916],[10,78],[17,68],[33,60],[20,-33],[-8,-117],[-18,-117],[-37,-52],[-19,46],[2,67]],[[119320,53463],[48,-1],[3,-16],[-22,-104],[-95,-230],[-59,-120],[-61,-56],[-22,1],[5,68],[21,45],[70,79],[13,40],[-14,179],[3,34],[39,53],[48,39],[23,-11]],[[107694,110387],[-202,101],[-154,32],[-58,22],[-135,28],[-53,58],[-14,53],[178,-74],[398,-142],[164,-93],[19,-40],[-143,55]],[[107344,91027],[-64,-39],[-117,-3],[-3,37],[122,37],[109,49],[93,55],[148,134],[50,-41],[-12,-45],[-275,-161],[-51,-23]],[[51855,87334],[-27,-43],[-15,15],[25,48],[17,-20]],[[1422,103784],[-36,-20],[-28,-47],[12,-67],[-42,-53],[-83,-44],[-15,41],[13,71],[-36,9],[-121,-20],[-6,18],[118,53],[-52,35],[51,42],[33,-15],[44,21],[130,-14],[18,-10]],[[31949,86603],[-22,-134],[-28,-6],[-20,98],[8,134],[47,16],[17,-36],[-2,-72]],[[31684,84881],[-50,-30],[-13,-59],[-24,20],[-15,56],[37,71],[-5,66],[-12,50],[52,-7],[17,-30],[43,-122],[-30,-15]],[[23099,109882],[34,68],[81,-1],[52,-65],[-31,-54],[-4,-81],[-50,-39],[-55,28],[-23,69],[-4,75]],[[41920,91425],[-34,-88],[-26,-21],[-50,2],[-51,62],[-27,-22],[-32,43],[18,33],[123,37],[57,27],[-15,-81],[37,8]],[[44232,93054],[3,-17],[-31,-88],[-21,-34],[-132,24],[-19,20],[-9,78],[-23,26],[37,28],[20,-13],[26,-61],[39,-10],[70,8],[40,39]],[[52773,91149],[-3,-67],[-31,-100],[-25,-3],[-41,38],[-99,4],[36,56],[10,52],[73,25],[54,89],[26,-94]],[[54779,9683],[-5,-39],[-108,9],[-52,19],[-4,67],[15,27],[60,34],[55,7],[34,-36],[5,-88]],[[26117,103815],[-89,-44],[-28,40],[-64,19],[48,51],[98,66],[9,33],[-20,138],[93,-46],[40,-54],[10,-68],[-34,-73],[-63,-62]],[[54755,803],[-184,112],[-52,2],[-26,-17],[-161,-15],[32,54],[3,26],[-29,74],[51,64],[31,0],[38,-36],[171,-69],[151,-55],[76,-53],[18,-38],[-7,-48],[-41,-22],[-71,21]],[[107607,46421],[-55,97],[-14,34],[-9,75],[-34,21],[-19,39],[101,18],[36,-4],[26,-22],[18,-47],[-8,-64],[13,-170],[-16,-7],[-39,30]],[[115056,106980],[17,83],[44,-18],[14,-33],[74,-30],[85,-19],[-39,-23],[-8,-33],[30,-43],[-15,-26],[10,-64],[-127,89],[2,55],[-20,21],[-78,14],[11,27]],[[165168,75457],[-63,89],[-2,25],[26,14],[50,-36],[56,0],[-8,-110],[-59,18]],[[165845,91887],[-31,49],[1,35],[23,56],[40,61],[41,44],[54,2],[40,-80],[-7,-36],[-30,-50],[-45,-51],[-71,-49],[-15,19]],[[145748,76724],[-17,1],[-123,47],[-46,2],[-32,23],[45,66],[23,14],[74,23],[24,-51],[87,-83],[30,-44],[-26,-14],[-39,16]],[[147823,76446],[-7,-29],[-35,-12],[-45,-35],[-37,-132],[-41,-124],[-22,-15],[-22,11],[-13,60],[70,174],[53,153],[27,137],[52,-24],[13,-79],[7,-85]],[[143724,96824],[-23,-64],[11,-37],[-34,-28],[-102,27],[-68,80],[-8,28],[7,63],[43,56],[71,40],[60,-24],[24,-25],[42,-90],[-23,-26]],[[139990,96919],[-111,127],[-64,60],[-44,1],[-32,-19],[-21,29],[38,81],[71,26],[55,2],[118,-34],[11,-23],[-3,-171],[-18,-79]],[[144050,99228],[-11,-43],[-122,-76],[-66,-21],[-82,7],[-14,51],[34,25],[120,67],[58,13],[70,-1],[13,-22]],[[119668,104654],[33,71],[52,30],[57,13],[14,-21],[-35,-113],[22,-64],[-34,-45],[-82,-8],[-20,17],[-7,120]],[[113945,109296],[129,-11],[-48,-53],[-126,-76],[-97,-92],[-153,-27],[-172,1],[-16,13],[-7,106],[51,-35],[34,-49],[193,60],[32,46],[159,44],[1,28],[-63,64],[47,4],[36,-23]],[[117060,110754],[51,15],[32,-27],[17,-50],[109,-150],[19,-48],[-28,-38],[-127,-6],[1,145],[-48,66],[6,21],[-32,72]],[[115961,109474],[24,-71],[-137,42],[-76,1],[-5,26],[73,19],[7,41],[-128,166],[-190,16],[-13,23],[39,37],[51,1],[34,33],[-50,24],[38,26],[68,-24],[51,-43],[81,4],[76,-9],[12,-18],[-34,-51],[77,-52],[17,-30],[-70,-67],[55,-94]],[[196896,6507],[-23,-73],[-25,-107],[0,-79],[-38,-63],[0,140],[11,35],[-17,205],[18,17],[7,-90],[33,-35],[34,50]],[[180244,115139],[-137,-73],[-158,-19],[-61,99],[222,36],[95,-13],[39,-30]],[[150254,74349],[-72,15],[-24,-23],[-10,-37],[-20,7],[-50,-41],[-24,-39],[-32,-9],[-48,68],[63,14],[19,26],[46,101],[-11,108],[-37,38],[-32,-1],[-87,-55],[68,76],[35,22],[41,-18],[56,-58],[1,-74],[34,-56],[33,20],[43,8],[43,-40],[1,-80],[-36,28]],[[139180,99847],[-8,-55],[-24,-32],[-100,-20],[-16,-22],[-43,3],[-30,103],[56,19],[58,-8],[71,21],[36,-9]],[[136650,99472],[-49,-47],[-53,-35],[-63,5],[-58,25],[0,30],[30,69],[20,17],[59,-59],[65,14],[49,-19]],[[116309,108982],[38,-33],[66,-28],[-9,-72],[-74,-36],[-86,-30],[-33,7],[-22,57],[81,43],[-3,34],[-54,33],[3,25],[50,64],[43,-64]],[[190037,112830],[-30,-17],[-37,-69],[-39,-6],[-105,68],[-9,54],[89,13],[126,-18],[5,-25]],[[190676,112236],[43,27],[-26,49],[97,-17],[61,-28],[70,-55],[-139,-19],[24,-54],[-112,39],[-18,58]],[[166546,76673],[-19,47],[-15,77],[85,3],[43,8],[14,-47],[-3,-36],[-64,-67],[-41,15]],[[163429,75300],[1,-24],[-50,-121],[-33,-33],[-21,28],[-23,64],[7,53],[49,-10],[17,52],[53,-9]],[[112350,101833],[-1,-54],[-33,-117],[-20,2],[-27,69],[-44,80],[22,71],[31,14],[60,-21],[12,-44]],[[166055,76777],[-17,1],[-48,45],[-27,9],[-45,-51],[-49,-65],[-35,-9],[-56,20],[-88,116],[-12,55],[11,45],[30,10],[28,-14],[20,-78],[22,-16],[34,-2],[29,17],[34,42],[48,5],[44,-36],[64,-66],[13,-28]],[[104602,106811],[18,31],[48,25],[116,12],[-73,-70],[20,-90],[55,-63],[85,-34],[-14,-21],[-91,40],[-99,69],[-43,52],[-22,49]],[[59049,88678],[-21,-66],[-57,-31],[-45,6],[-14,25],[19,69],[37,15],[15,-41],[41,2],[25,21]],[[58382,90191],[-26,-14],[-23,23],[66,67],[29,90],[50,47],[7,-61],[-34,-86],[-69,-66]],[[26662,86560],[-41,-182],[-42,11],[-10,61],[-3,111],[8,72],[38,40],[51,-63],[-1,-50]],[[27184,84612],[-32,-23],[-35,45],[-52,96],[-35,88],[32,63],[38,17],[25,-26],[11,-92],[19,-101],[29,-67]],[[123376,85185],[-10,25],[25,67],[47,-22],[86,-75],[77,-29],[134,-150],[31,-71],[-3,-45],[-63,-37],[-135,-6],[-38,45],[-37,176],[-52,71],[-62,51]],[[164240,75624],[-56,-22],[11,-49],[-48,16],[-3,68],[23,27],[54,15],[17,-17],[2,-38]],[[164101,75974],[23,-1],[4,51],[58,22],[21,-35],[3,-36],[-53,-46],[-44,-5],[-12,50]],[[165695,79644],[-44,-80],[-23,-3],[-22,42],[-26,111],[6,59],[23,-6],[52,-50],[30,-42],[4,-31]],[[164101,75677],[-69,82],[30,4],[39,-34],[18,20],[-14,59],[34,6],[21,21],[32,-33],[-2,-45],[-41,-30],[-16,-48],[-32,-2]],[[163926,75379],[-29,-49],[-20,-66],[-27,18],[-16,62],[1,89],[42,55],[18,0],[-6,-54],[8,-36],[29,-19]],[[163842,75236],[-26,11],[-16,47],[-56,-52],[-16,13],[24,47],[-37,36],[46,19],[32,24],[16,-4],[33,-141]],[[27085,81952],[-23,-46],[-18,19],[41,27]],[[26925,111397],[41,57],[-12,39],[62,79],[34,-2],[29,-43],[50,-37],[34,-47],[-74,-31],[-32,-97],[35,-70],[46,-33],[-9,-27],[-121,-72],[-67,-1],[7,216],[-23,69]],[[56057,93526],[-34,48],[0,70],[-34,49],[-51,22],[-40,-48],[13,-62],[-15,-19],[-66,26],[-21,-66],[-40,11],[-39,-6],[-36,-36],[-21,13],[-31,88],[-38,-10],[-39,-26],[51,120],[52,-62],[17,-32],[28,6],[48,38],[67,93],[17,11],[60,-12],[21,28],[-2,44],[-30,79],[41,-17],[31,-108],[38,-14],[41,6],[16,-44],[71,-82],[31,-9],[32,-68],[33,-95],[2,-30],[-42,4],[-19,60],[-8,69],[-26,34],[-21,0],[-47,-40],[-10,-33]],[[63689,93363],[26,62],[16,5],[20,-68],[-24,-35],[-49,-11],[-45,-46],[-64,-102],[-115,-139],[-45,-46],[-38,-68],[-106,-106],[-58,-27],[-91,19],[137,103],[67,79],[60,33],[112,109],[59,82],[5,77],[24,8],[109,71]],[[23133,110193],[-80,-7],[-31,167],[-39,75],[-46,12],[-23,36],[69,5],[70,-25],[56,-60],[117,-79],[-5,-93],[-88,-31]],[[23487,110233],[13,-15],[-58,-162],[-52,-42],[-84,-19],[-93,24],[17,79],[-26,43],[283,92]],[[27632,108070],[-6,76],[-87,128],[97,42],[126,27],[105,-86],[11,-59],[-19,-43],[25,-98],[-52,-77],[4,-28],[46,-26],[49,6],[6,-53],[-109,-14],[-144,98],[-34,41],[-18,66]],[[44194,97570],[21,65],[40,32],[48,-30],[41,31],[58,15],[82,-50],[9,-17],[139,-41],[-41,-36],[-79,-40],[-73,-10],[-56,5],[-78,58],[-111,18]],[[42552,95051],[-42,-35],[-49,38],[-64,75],[-31,50],[15,53],[50,35],[51,-21],[40,18],[43,-24],[40,-107],[-37,-34],[-16,-48]],[[41925,97541],[-67,-4],[-98,-20],[-96,23],[-102,48],[-45,39],[-9,30],[9,69],[105,53],[70,-3],[8,-38],[-38,-28],[18,-37],[38,-16],[66,5],[202,26],[49,-9],[27,-33],[-7,-45],[-52,-38],[-78,-22]],[[55806,96850],[65,-21],[-46,-67],[5,-146],[27,-99],[-38,-9],[-21,50],[0,95],[-32,17],[-37,-63],[-81,69],[-2,84],[-38,-19],[-55,-58],[-24,52],[34,61],[50,5],[52,18],[1,31],[-47,18],[5,42],[41,-21],[34,-37],[55,-10],[52,8]],[[35981,98873],[-65,-25],[-49,19],[-58,-26],[-27,-33],[-50,-90],[-180,-65],[-63,30],[-12,110],[24,76],[-12,59],[19,28],[72,50],[74,-21],[71,25],[93,-14],[31,14],[60,-47],[70,-8],[18,-61],[-16,-21]],[[42335,90544],[13,-68],[-16,-24],[-75,-13],[-33,20],[-30,40],[-11,51],[21,49],[47,34],[41,15],[33,-18],[10,-86]],[[55348,88101],[44,-45],[-9,-44],[-26,-7],[-55,60],[-43,29],[8,48],[33,40],[25,-21],[23,-60]],[[50576,88732],[-18,-38],[-24,4],[-10,73],[-27,49],[9,25],[38,29],[34,45],[31,86],[30,32],[45,-69],[31,-26],[28,-80],[-22,-61],[-145,-69]],[[42108,98569],[136,-39],[-7,-32],[-78,0],[-45,-14],[-44,-68],[-99,-33],[-79,-16],[-72,-30],[-29,-61],[-51,-39],[-124,-63],[-2,39],[18,43],[50,38],[84,120],[18,123],[120,40],[95,6],[109,-14]],[[150479,80586],[-29,-63],[-131,34],[-54,65],[25,18],[71,-15],[91,-7],[27,-32]],[[124908,76001],[-20,-25],[-48,0],[-50,13],[-177,59],[-37,2],[-127,-16],[-19,72],[-42,32],[-4,28],[21,33],[66,13],[72,-55],[46,-27],[24,-2],[109,-51],[40,10],[41,28],[42,-3],[17,-38],[46,-73]],[[27624,87778],[-24,1],[-58,-47],[-25,21],[-42,85],[4,56],[80,11],[75,-19],[30,-100],[-40,-8]],[[64012,92461],[-60,-63],[-61,-13],[75,124],[22,11],[33,-32],[-9,-27]],[[41127,99298],[-67,-129],[-32,-32],[-87,-22],[13,63],[46,37],[13,40],[114,43]],[[302,108262],[-12,-32],[-113,-69],[-66,-15],[-49,10],[-56,43],[-6,26],[46,24],[97,-3],[159,16]],[[1229,104466],[10,-35],[-45,-27],[-46,64],[3,22],[61,45],[164,13],[-147,-82]],[[1984,104140],[-38,-33],[-104,39],[-3,29],[58,10],[87,-45]],[[4116,102651],[-46,-20],[-10,62],[-23,43],[31,21],[39,-59],[9,-47]],[[4491,103344],[-20,-8],[-163,25],[-170,-47],[-85,-14],[-71,3],[-23,48],[77,-3],[141,30],[127,34],[70,-30],[67,-4],[50,-34]],[[3943,103580],[16,23],[124,49],[83,-40],[-2,-37],[-177,20],[-44,-15]],[[4207,102990],[-69,-20],[13,-49],[82,-24],[18,-36],[-102,17],[-37,34],[-72,32],[-31,30],[-50,97],[22,2],[87,-63],[89,16],[85,-7],[-35,-29]],[[151931,93236],[-11,-24],[-59,-32],[-84,-33],[-81,-41],[-91,100],[-36,14],[-184,106],[-19,26],[18,50],[115,15],[117,-28],[123,-82],[108,9],[77,-52],[7,-28]],[[147303,76224],[-21,-40],[-4,-56],[-25,-63],[-26,76],[-38,28],[-31,-1],[-42,-25],[-32,4],[-44,23],[-36,33],[-14,36],[5,44],[24,40],[43,10],[50,-16],[96,0],[88,-47],[7,-46]],[[121889,58230],[-4,-35],[-24,5],[0,45],[28,-15]],[[42696,92586],[-45,16],[-32,27],[-4,31],[112,-48],[45,-6],[53,58],[-140,44],[-7,21],[-76,4],[-22,-35],[-95,50],[-35,-7],[-25,15],[6,39],[-19,45],[23,40],[-31,35],[-47,26],[37,31],[-17,52],[82,30],[36,-102],[35,-32],[1,-42],[-38,-56],[25,-15],[61,-2],[81,19],[49,-5],[60,12],[8,-17],[-53,-40],[2,-18],[79,-25],[45,10],[30,-32],[21,-91],[45,-11],[33,-61],[-13,-37],[-37,44],[-47,24],[-186,9]],[[25389,99736],[-60,6],[-126,-27],[-62,-3],[-51,-16],[-74,-62],[-170,10],[-167,28],[-39,-18],[-5,-50],[18,-74],[26,-14],[33,-90],[123,-120],[126,-185],[37,-30],[6,-148],[-8,-34],[-89,220],[-52,94],[-54,74],[-78,66],[-100,103],[-65,32],[-75,88],[10,73],[-40,12],[-63,0],[-103,47],[18,16],[81,-26],[34,32],[-110,114],[-80,142],[-206,328],[-1,42],[45,-14],[175,-253],[59,-118],[64,-77],[63,-24],[40,17],[12,-30],[-12,-65],[21,-50],[40,-37],[106,11],[193,-22],[129,-19],[70,18],[16,26],[70,37],[125,22],[70,19],[42,-1],[84,-28],[93,-41],[24,-25],[0,-34],[-72,7],[-91,51]],[[40864,99785],[-21,-53],[-93,13],[-121,-6],[-59,-22],[-12,22],[40,60],[132,25],[42,40],[89,50],[87,27],[24,-35],[-141,-58],[36,-44],[-3,-19]],[[34915,91793],[-69,186],[-28,60],[-41,29],[-39,-48],[-21,-73],[-25,2],[-126,-34],[-28,-79],[-25,38],[-51,-26],[-22,20],[-89,18],[-88,0],[-56,-47],[-108,-10],[-62,-57],[-120,-145],[9,59],[-12,67],[-67,35],[-74,4],[-2,27],[75,1],[67,-25],[24,-46],[29,-12],[75,82],[15,31],[-16,39],[151,-24],[44,50],[152,9],[77,-9],[35,9],[-2,36],[95,23],[36,32],[66,34],[17,45],[-13,48],[91,64],[39,-12],[24,-28],[-19,-55],[8,-44],[53,-60],[-8,-19],[16,-51],[14,-104],[-1,-40]],[[37297,91777],[-53,20],[-57,-3],[-44,-21],[-6,24],[71,63],[-22,93],[-132,27],[-32,31],[-22,92],[24,70],[-26,40],[-61,46],[-81,34],[-93,5],[-74,-27],[-30,-31],[-59,-32],[-71,-13],[-36,11],[-36,65],[-98,-74],[-5,53],[54,63],[53,11],[30,-43],[43,-38],[57,8],[83,65],[55,24],[51,9],[66,-9],[74,-24],[68,-37],[45,-49],[24,-79],[-13,-54],[2,-62],[28,-20],[70,-16],[14,13],[-19,111],[61,25],[28,-4],[-18,-46],[33,-136],[41,-34],[-25,-49],[30,-50],[82,-36],[80,-22],[97,5],[54,12],[21,17],[183,28],[43,16],[37,-37],[-152,-93],[-111,8],[-192,-48],[-67,62],[-67,36]],[[19776,103030],[-37,17],[-32,59],[-60,137],[-49,81],[-123,120],[-131,79],[-72,61],[-32,47],[-103,123],[-34,63],[63,-21],[93,-121],[106,-98],[56,-36],[105,-52],[93,-105],[78,-128],[16,-74],[67,-118],[-4,-34]],[[28077,93434],[-30,-5],[-61,71],[-45,103],[-10,75],[16,140],[26,21],[-15,-178],[13,-64],[23,-55],[80,-76],[3,-32]],[[29066,92330],[-113,-50],[-19,-29],[-8,-85],[-58,-13],[33,135],[24,50],[49,37],[-16,37],[-40,15],[-33,-20],[-26,7],[41,48],[39,1],[42,-17],[85,-116]],[[26846,92042],[-63,27],[-32,31],[-18,45],[-191,188],[-3,38],[218,-229],[55,-70],[34,-30]],[[35397,102582],[-93,-28],[-63,-41],[-72,-123],[-49,16],[36,85],[150,112],[196,136],[118,-31],[38,-24],[-38,-60],[-223,-42]],[[11364,105577],[-24,-47],[-29,10],[-82,110],[12,88],[43,-84],[80,-77]],[[41495,92211],[110,-59],[-6,-56],[-61,-16],[-125,-3],[-81,17],[-13,33],[12,43],[38,42],[58,29],[44,-5],[24,-25]],[[41466,92297],[27,25],[78,31],[88,17],[46,-3],[29,-20],[2,-35],[-24,-32],[-88,-37],[-84,-2],[-67,8],[-7,48]],[[39741,105779],[-95,108],[-43,24],[-64,-66],[-21,-3],[-4,86],[-36,51],[-24,62],[61,-16],[86,-2],[9,39],[55,-8],[68,21],[128,-1],[29,25],[20,53],[36,-22],[-26,-56],[29,-188],[56,-56],[25,-65],[92,-67],[-74,-37],[-106,65],[-90,16],[-108,-32],[-49,49],[46,20]],[[55725,107777],[6,23],[-98,59],[-66,18],[-79,-23],[-26,-58],[-75,-59],[-130,-16],[-30,-16],[-92,22],[-83,45],[-59,45],[1,46],[-50,53],[-34,84],[111,-6],[-13,-69],[32,-19],[40,34],[97,22],[43,20],[-14,39],[-61,35],[49,62],[-9,73],[-71,106],[15,61],[-58,84],[17,10],[113,-91],[35,-20],[107,-22],[40,-51],[81,-45],[91,-36],[98,-51],[45,-43],[35,-109],[128,-152],[22,-35],[-61,-40],[38,-77],[-183,53],[48,44]],[[38016,74891],[-11,-16],[-60,-24],[-65,7],[-23,23],[-57,13],[15,53],[-29,37],[20,43],[56,-64],[63,-65],[69,51],[2,-33],[20,-25]],[[52313,91376],[26,-16],[45,-3],[16,-33],[-36,-26],[-30,8],[-64,53],[-29,112],[-13,27],[-48,-26],[-4,65],[-17,40],[-27,29],[13,60],[-31,64],[174,162],[3,-17],[-106,-96],[-10,-43],[64,44],[-16,-57],[102,-32],[30,-28],[-8,-43],[-32,-64],[-21,-150],[19,-30]],[[36967,66190],[-1,-30],[-43,-35],[-87,41],[-43,14],[-95,-4],[-50,19],[-50,30],[-18,27],[216,64],[64,-9],[72,-113],[35,-4]],[[46683,58823],[-1,-98],[-41,-21],[-50,-12],[-32,5],[10,60],[-15,36],[-85,-6],[-35,28],[-25,37],[-11,58],[71,21],[40,33],[36,16],[26,-6],[35,-85],[43,-63],[34,-3]],[[65358,24219],[-65,-39],[-9,-64],[51,-94],[-54,27],[-26,0],[-21,-51],[4,-33],[29,-49],[11,-42],[-35,-95],[52,-34],[52,-22],[-82,-47],[-54,-18],[19,-58],[-14,-16],[-60,-7],[-32,-41],[-31,22],[0,107],[16,33],[35,30],[17,42],[-6,33],[-27,39],[-86,23],[6,42],[107,-9],[13,82],[-12,28],[-157,52],[15,11],[92,-16],[48,7],[33,27],[5,41],[-35,50],[21,37],[52,19],[8,66],[-19,39],[-1,34],[19,31],[-7,52],[11,35],[29,11],[12,34],[-1,43],[-32,103],[10,84],[16,28],[33,24],[-10,-73],[8,-76],[11,-57],[-4,-333],[45,-62]],[[50446,55862],[-12,-38],[-63,-28],[-18,-26],[-15,9],[26,47],[-53,45],[-14,-71],[-20,-52],[-31,35],[-10,-32],[-22,15],[31,52],[36,89],[17,70],[67,45],[6,-58],[-26,-29],[21,-47],[53,-43],[27,17]],[[114014,31789],[-26,-45],[-6,-28],[-38,12],[-22,-89],[-39,-35],[-38,32],[-5,60],[-20,9],[-85,-55],[-21,-32],[-11,-67],[-43,65],[-24,6],[-111,-47],[-20,-20],[-19,-81],[11,-43],[-16,-50],[-17,-17],[-52,24],[-69,-99],[-52,-27],[-27,-31],[-21,-43],[-9,-38],[-45,-147],[-44,-53],[-27,-23],[-1,-33],[-59,-32],[-39,-66],[-99,-224],[-23,-41],[-33,-2],[-1,30],[44,73],[7,72],[38,77],[49,73],[48,57],[20,43],[-16,41],[16,66],[35,84],[104,211],[46,7],[-1,23],[18,77],[1,29],[48,10],[18,15],[-24,61],[2,16],[62,46],[44,-46],[40,1],[37,59],[36,45],[4,48],[42,24],[22,27],[8,31],[31,7],[25,-24],[60,42],[33,4],[15,29],[26,10],[81,-35],[55,-37],[51,-4],[26,-32]],[[116647,48850],[-89,-18],[-63,-36],[-39,-40],[-21,-36],[9,-20],[69,-11],[-29,-38],[-7,-28],[27,-13],[44,32],[5,-60],[-39,-124],[-18,29],[-62,55],[-61,64],[-70,46],[-18,-24],[14,-23],[1,-64],[-20,19],[-21,52],[-35,6],[-13,-51],[-95,67],[-15,24],[-10,56],[5,23],[38,-7],[73,29],[39,33],[6,-29],[74,13],[35,-15],[14,12],[46,69],[16,-17],[63,51],[44,-13],[27,20],[33,50],[19,-3],[6,-41],[18,-39]],[[120028,103194],[48,-16],[61,-51],[34,-3],[42,-40],[73,-29],[42,3],[-33,-55],[44,-50],[-37,-42],[-107,62],[-88,73],[-122,115],[-82,85],[24,29],[58,-56],[43,-25]],[[102749,91739],[-27,-20],[-136,7],[-118,79],[-40,68],[-40,95],[11,9],[57,-46],[166,-101],[127,-91]],[[110652,104722],[-113,-29],[-117,-58],[-4,15],[125,82],[121,30],[-12,-40]],[[117838,76456],[-14,-68],[-48,0],[-23,69],[13,49],[27,1],[45,-51]],[[117902,77956],[-32,-57],[-26,77],[4,49],[15,34],[20,7],[21,-61],[-2,-49]],[[172544,32225],[-52,-120],[8,-28],[-17,-63],[-30,-62],[-29,13],[-15,-26],[-6,-55],[4,-43],[-8,-130],[-16,-79],[-13,-10],[-23,67],[-66,96],[0,68],[-27,6],[-4,26],[22,59],[31,122],[9,60],[24,24],[-18,60],[7,15],[47,10],[21,43],[10,-79],[12,0],[78,51],[27,0],[24,-25]],[[154323,79878],[-15,-14],[-97,8],[-7,-38],[-47,6],[-33,73],[4,25],[42,48],[15,69],[65,57],[19,-51],[4,-95],[50,-88]],[[154078,80007],[-63,-103],[-42,-29],[-34,-8],[-86,76],[-10,84],[96,24],[70,-3],[69,-41]],[[139934,86348],[-52,-13],[-56,11],[-79,-2],[-78,-15],[-93,-71],[-26,22],[32,38],[11,37],[-29,100],[15,24],[77,-2],[73,-19],[85,-46],[89,-25],[31,-39]],[[144624,75950],[-31,-26],[-26,3],[-50,48],[35,54],[0,78],[8,45],[20,10],[36,-42],[1,-26],[-26,-34],[33,-86],[0,-24]],[[144785,76010],[-18,-42],[-34,-42],[-36,-18],[-24,18],[-11,46],[1,92],[20,19],[67,-19],[29,-23],[6,-31]],[[97215,53554],[2,-25],[-26,-30],[-45,-18],[-12,-138],[3,-78],[-21,-26],[-13,13],[-4,49],[-19,1],[-30,-24],[-29,67],[-23,36],[-91,89],[-101,77],[-121,80],[-23,45],[8,32],[47,-47],[17,24],[19,-16],[16,-53],[14,-13],[28,7],[7,80],[8,11],[34,-104],[38,5],[44,18],[11,-16],[19,-71],[28,-18],[23,1],[19,110],[-26,39],[1,33],[14,21],[32,-2],[38,23],[15,37],[-26,4],[-10,31],[48,34],[-21,37],[8,28],[-10,88],[-43,37],[14,13],[38,-7],[49,68],[-6,36],[-22,9],[-19,27],[-31,84],[-2,42],[-28,22],[-36,6],[-1,-59],[-37,26],[-17,35],[-50,7],[21,28],[42,6],[86,-10],[-17,46],[-28,46],[-68,5],[-36,-7],[-111,-42],[-65,6],[71,45],[39,15],[36,2],[19,16],[4,34],[54,1],[36,32],[-33,22],[3,58],[-23,107],[-72,79],[-21,-36],[-26,-2],[9,45],[-39,41],[-31,-32],[-16,5],[0,61],[-10,24],[-28,-50],[-30,5],[-66,-4],[4,31],[55,7],[32,29],[-72,100],[-42,34],[6,47],[-75,176],[-51,71],[-24,14],[-51,-12],[-14,-76],[-39,100],[-20,9],[-57,52],[-32,14],[-22,-7],[-60,-63],[-69,-33],[0,42],[42,20],[126,103],[63,-74],[83,-20],[38,16],[-4,27],[-38,58],[-7,57],[-17,67],[0,120],[26,2],[5,-137],[27,-71],[29,-9],[17,-64],[38,-114],[34,-63],[33,13],[-4,100],[32,-48],[36,-29],[-57,-97],[14,-49],[35,-62],[56,-68],[42,-73],[40,-39],[53,12],[23,19],[41,10],[23,-7],[26,-30],[5,-26],[-75,-8],[-6,-18],[119,-137],[15,-34],[12,-92],[18,-61],[16,-17],[46,-1],[9,27],[-7,39],[0,53],[12,79],[-10,35],[5,37],[17,12],[41,-7],[-10,70],[-27,75],[18,12],[47,102],[-73,37],[-4,136],[25,24],[-8,-148],[75,-34],[24,-49],[-41,-12],[-18,-101],[28,-75],[21,-29],[-40,-14],[-22,-19],[-30,-113],[34,-70],[53,39],[-6,-96],[17,-37],[3,-176],[-25,3],[-6,-39],[12,-43],[41,49],[-25,-177],[10,-70],[-13,-42],[5,-31],[-5,-58],[16,-62],[-11,-23],[-10,-93],[-14,-48],[-2,-120],[11,-89],[10,2]],[[97132,53558],[30,43],[-13,39],[-26,26],[-31,12],[-26,-27],[-3,-58],[45,-54],[24,19]]],"transform":{"scale":[0.0017093077481475233,0.0010698559416761953],"translate":[-165.898486328125,-50.620019531249994]},"objects":{"ne_50m_lakes":{"type":"GeometryCollection","geometries":[{"arcs":[[0]],"type":"Polygon"},{"arcs":[[1]],"type":"Polygon"},{"arcs":[[2]],"type":"Polygon"},{"arcs":[[3]],"type":"Polygon"},{"arcs":[[4,5],[6],[7]],"type":"Polygon"},{"arcs":[[8]],"type":"Polygon"},{"arcs":[[9],[10],[11]],"type":"Polygon"},{"arcs":[[12],[13],[14],[15],[16]],"type":"Polygon"},{"arcs":[[17],[18]],"type":"Polygon"},{"arcs":[[19],[20],[21]],"type":"Polygon"},{"arcs":[[22],[23]],"type":"Polygon"},{"arcs":[[24],[25],[26]],"type":"Polygon"},{"arcs":[[27]],"type":"Polygon"},{"arcs":[[-5,28]],"type":"Polygon"},{"arcs":[[29]],"type":"Polygon"},{"arcs":[[30]],"type":"Polygon"},{"arcs":[[31]],"type":"Polygon"},{"arcs":[[32]],"type":"Polygon"},{"arcs":[[33]],"type":"Polygon"},{"arcs":[[34]],"type":"Polygon"},{"arcs":[[35],[36]],"type":"Polygon"},{"arcs":[[37],[38]],"type":"Polygon"},{"arcs":[[39,40],[41],[42],[43],[44]],"type":"Polygon"},{"arcs":[[45,46],[47],[48],[49],[50],[51],[52],[53]],"type":"Polygon"},{"arcs":[[54,-41,55,-46],[56],[57],[58],[59],[60],[61],[62],[63]],"type":"Polygon"},{"arcs":[[64]],"type":"Polygon"},{"arcs":[[65]],"type":"Polygon"},{"arcs":[[66]],"type":"Polygon"},{"arcs":[[67]],"type":"Polygon"},{"arcs":[[68],[69]],"type":"Polygon"},{"arcs":[[70]],"type":"Polygon"},{"arcs":[[71]],"type":"Polygon"},{"arcs":[[72],[73]],"type":"Polygon"},{"arcs":[[74],[75],[76]],"type":"Polygon"},{"arcs":[[77]],"type":"Polygon"},{"arcs":[[78]],"type":"Polygon"},{"arcs":[[79]],"type":"Polygon"},{"arcs":[[80]],"type":"Polygon"},{"arcs":[[81]],"type":"Polygon"},{"arcs":[[82]],"type":"Polygon"},{"arcs":[[83]],"type":"Polygon"},{"arcs":[[84]],"type":"Polygon"},{"arcs":[[85]],"type":"Polygon"},{"arcs":[[86]],"type":"Polygon"},{"arcs":[[87]],"type":"Polygon"},{"arcs":[[88]],"type":"Polygon"},{"arcs":[[89,90]],"type":"Polygon"},{"arcs":[[-90,91]],"type":"Polygon"},{"arcs":[[92]],"type":"Polygon"},{"arcs":[[93]],"type":"Polygon"},{"arcs":[[94]],"type":"Polygon"},{"arcs":[[95]],"type":"Polygon"},{"arcs":[[96]],"type":"Polygon"},{"arcs":[[97]],"type":"Polygon"},{"arcs":[[98]],"type":"Polygon"},{"arcs":[[99]],"type":"Polygon"},{"arcs":[[100]],"type":"Polygon"},{"arcs":[[101]],"type":"Polygon"},{"arcs":[[102]],"type":"Polygon"},{"arcs":[[103]],"type":"Polygon"},{"arcs":[[104]],"type":"Polygon"},{"arcs":[[105]],"type":"Polygon"},{"arcs":[[106]],"type":"Polygon"},{"arcs":[[107,108]],"type":"Polygon"},{"arcs":[[109]],"type":"Polygon"},{"arcs":[[110]],"type":"Polygon"},{"arcs":[[111]],"type":"Polygon"},{"arcs":[[112]],"type":"Polygon"},{"arcs":[[113]],"type":"Polygon"},{"arcs":[[114]],"type":"Polygon"},{"arcs":[[115]],"type":"Polygon"},{"arcs":[[116]],"type":"Polygon"},{"arcs":[[117]],"type":"Polygon"},{"arcs":[[118]],"type":"Polygon"},{"arcs":[[119]],"type":"Polygon"},{"arcs":[[120]],"type":"Polygon"},{"arcs":[[121]],"type":"Polygon"},{"arcs":[[122]],"type":"Polygon"},{"arcs":[[123]],"type":"Polygon"},{"arcs":[[124]],"type":"Polygon"},{"arcs":[[125]],"type":"Polygon"},{"arcs":[[126]],"type":"Polygon"},{"arcs":[[127]],"type":"Polygon"},{"arcs":[[-108,128]],"type":"Polygon"},{"arcs":[[129]],"type":"Polygon"},{"arcs":[[130]],"type":"Polygon"},{"arcs":[[131]],"type":"Polygon"},{"arcs":[[132]],"type":"Polygon"},{"arcs":[[133]],"type":"Polygon"},{"arcs":[[134]],"type":"Polygon"},{"arcs":[[135,136]],"type":"Polygon"},{"arcs":[[-136,137]],"type":"Polygon"},{"arcs":[[138]],"type":"Polygon"},{"arcs":[[139]],"type":"Polygon"},{"arcs":[[140]],"type":"Polygon"},{"arcs":[[141]],"type":"Polygon"},{"arcs":[[142]],"type":"Polygon"},{"arcs":[[143]],"type":"Polygon"},{"arcs":[[144]],"type":"Polygon"},{"arcs":[[145]],"type":"Polygon"},{"arcs":[[146]],"type":"Polygon"},{"arcs":[[147]],"type":"Polygon"},{"arcs":[[148]],"type":"Polygon"},{"arcs":[[149]],"type":"Polygon"},{"arcs":[[150]],"type":"Polygon"},{"arcs":[[151]],"type":"Polygon"},{"arcs":[[152]],"type":"Polygon"},{"arcs":[[153]],"type":"Polygon"},{"arcs":[[154]],"type":"Polygon"},{"arcs":[[155]],"type":"Polygon"},{"arcs":[[156]],"type":"Polygon"},{"arcs":[[157]],"type":"Polygon"},{"arcs":[[158]],"type":"Polygon"},{"arcs":[[159]],"type":"Polygon"},{"arcs":[[160]],"type":"Polygon"},{"arcs":[[161]],"type":"Polygon"},{"arcs":[[162]],"type":"Polygon"},{"arcs":[[163]],"type":"Polygon"},{"arcs":[[164]],"type":"Polygon"},{"arcs":[[165]],"type":"Polygon"},{"arcs":[[166]],"type":"Polygon"},{"arcs":[[167]],"type":"Polygon"},{"arcs":[[168]],"type":"Polygon"},{"arcs":[[169]],"type":"Polygon"},{"arcs":[[170]],"type":"Polygon"},{"arcs":[[171]],"type":"Polygon"},{"arcs":[[172]],"type":"Polygon"},{"arcs":[[173]],"type":"Polygon"},{"arcs":[[174]],"type":"Polygon"},{"arcs":[[175]],"type":"Polygon"},{"arcs":[[176]],"type":"Polygon"},{"arcs":[[177]],"type":"Polygon"},{"arcs":[[178]],"type":"Polygon"},{"arcs":[[179]],"type":"Polygon"},{"arcs":[[180]],"type":"Polygon"},{"arcs":[[181]],"type":"Polygon"},{"arcs":[[182]],"type":"Polygon"},{"arcs":[[183]],"type":"Polygon"},{"arcs":[[184]],"type":"Polygon"},{"arcs":[[185]],"type":"Polygon"},{"arcs":[[186]],"type":"Polygon"},{"arcs":[[187]],"type":"Polygon"},{"arcs":[[188]],"type":"Polygon"},{"arcs":[[189]],"type":"Polygon"},{"arcs":[[190]],"type":"Polygon"},{"arcs":[[191]],"type":"Polygon"},{"arcs":[[192]],"type":"Polygon"},{"arcs":[[193]],"type":"Polygon"},{"arcs":[[194]],"type":"Polygon"},{"arcs":[[195]],"type":"Polygon"},{"arcs":[[196]],"type":"Polygon"},{"arcs":[[197]],"type":"Polygon"},{"arcs":[[198]],"type":"Polygon"},{"arcs":[[199]],"type":"Polygon"},{"arcs":[[200]],"type":"Polygon"},{"arcs":[[201]],"type":"Polygon"},{"arcs":[[202]],"type":"Polygon"},{"arcs":[[203]],"type":"Polygon"},{"arcs":[[204]],"type":"Polygon"},{"arcs":[[205]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[206]],"type":"Polygon"},{"arcs":[[207]],"type":"Polygon"},{"arcs":[[208]],"type":"Polygon"},{"arcs":[[209]],"type":"Polygon"},{"arcs":[[210]],"type":"Polygon"},{"arcs":[[211]],"type":"Polygon"},{"arcs":[[212]],"type":"Polygon"},{"arcs":[[213]],"type":"Polygon"},{"arcs":[[214]],"type":"Polygon"},{"arcs":[[215]],"type":"Polygon"},{"arcs":[[216]],"type":"Polygon"},{"arcs":[[217]],"type":"Polygon"},{"arcs":[[218]],"type":"Polygon"},{"arcs":[[219]],"type":"Polygon"},{"arcs":[[220]],"type":"Polygon"},{"arcs":[[221]],"type":"Polygon"},{"arcs":[[222]],"type":"Polygon"},{"arcs":[[223]],"type":"Polygon"},{"arcs":[[224]],"type":"Polygon"},{"arcs":[[225]],"type":"Polygon"},{"arcs":[[226]],"type":"Polygon"},{"arcs":[[227]],"type":"Polygon"},{"arcs":[[228]],"type":"Polygon"},{"arcs":[[229]],"type":"Polygon"},{"arcs":[[230]],"type":"Polygon"},{"arcs":[[231]],"type":"Polygon"},{"arcs":[[232]],"type":"Polygon"},{"arcs":[[233]],"type":"Polygon"},{"arcs":[[234]],"type":"Polygon"},{"arcs":[[235]],"type":"Polygon"},{"arcs":[[236]],"type":"Polygon"},{"arcs":[[237]],"type":"Polygon"},{"arcs":[[238]],"type":"Polygon"},{"arcs":[[239]],"type":"Polygon"},{"arcs":[[240]],"type":"Polygon"},{"arcs":[[241]],"type":"Polygon"},{"arcs":[[242]],"type":"Polygon"},{"arcs":[[243]],"type":"Polygon"},{"arcs":[[244]],"type":"Polygon"},{"arcs":[[245]],"type":"Polygon"},{"arcs":[[246]],"type":"Polygon"},{"arcs":[[247]],"type":"Polygon"},{"arcs":[[248]],"type":"Polygon"},{"arcs":[[249]],"type":"Polygon"},{"arcs":[[250]],"type":"Polygon"},{"arcs":[[251]],"type":"Polygon"},{"arcs":[[252]],"type":"Polygon"},{"arcs":[[253]],"type":"Polygon"},{"arcs":[[254]],"type":"Polygon"},{"arcs":[[255]],"type":"Polygon"},{"arcs":[[256]],"type":"Polygon"},{"arcs":[[257]],"type":"Polygon"},{"arcs":[[258]],"type":"Polygon"},{"arcs":[[259]],"type":"Polygon"},{"arcs":[[260]],"type":"Polygon"},{"arcs":[[261]],"type":"Polygon"},{"arcs":[[262]],"type":"Polygon"},{"arcs":[[263]],"type":"Polygon"},{"arcs":[[264]],"type":"Polygon"},{"arcs":[[265]],"type":"Polygon"},{"arcs":[[266]],"type":"Polygon"},{"arcs":[[267]],"type":"Polygon"},{"arcs":[[268]],"type":"Polygon"},{"arcs":[[269]],"type":"Polygon"},{"arcs":[[270]],"type":"Polygon"},{"arcs":[[271]],"type":"Polygon"},{"arcs":[[272]],"type":"Polygon"},{"arcs":[[273]],"type":"Polygon"},{"arcs":[[274]],"type":"Polygon"},{"arcs":[[275]],"type":"Polygon"},{"arcs":[[276]],"type":"Polygon"},{"arcs":[[277]],"type":"Polygon"},{"arcs":[[278]],"type":"Polygon"},{"arcs":[[279]],"type":"Polygon"},{"arcs":[[280]],"type":"Polygon"},{"arcs":[[281]],"type":"Polygon"},{"arcs":[[282]],"type":"Polygon"},{"arcs":[[283]],"type":"Polygon"},{"arcs":[[284]],"type":"Polygon"},{"arcs":[[285]],"type":"Polygon"},{"arcs":[[286]],"type":"Polygon"},{"arcs":[[287]],"type":"Polygon"},{"arcs":[[288]],"type":"Polygon"},{"arcs":[[289]],"type":"Polygon"},{"arcs":[[290]],"type":"Polygon"},{"arcs":[[291]],"type":"Polygon"},{"arcs":[[292]],"type":"Polygon"},{"arcs":[[293]],"type":"Polygon"},{"arcs":[[294]],"type":"Polygon"},{"arcs":[[295]],"type":"Polygon"},{"arcs":[[296]],"type":"Polygon"},{"arcs":[[297]],"type":"Polygon"},{"arcs":[[298]],"type":"Polygon"},{"arcs":[[299]],"type":"Polygon"},{"arcs":[[300]],"type":"Polygon"},{"arcs":[[301]],"type":"Polygon"},{"arcs":[[302]],"type":"Polygon"},{"arcs":[[303]],"type":"Polygon"},{"arcs":[[304]],"type":"Polygon"},{"arcs":[[305]],"type":"Polygon"},{"arcs":[[306]],"type":"Polygon"},{"arcs":[[307]],"type":"Polygon"},{"arcs":[[308]],"type":"Polygon"},{"arcs":[[309]],"type":"Polygon"},{"arcs":[[310]],"type":"Polygon"},{"arcs":[[311]],"type":"Polygon"},{"arcs":[[312]],"type":"Polygon"},{"arcs":[[313]],"type":"Polygon"},{"arcs":[[314]],"type":"Polygon"},{"arcs":[[315]],"type":"Polygon"},{"arcs":[[316]],"type":"Polygon"},{"arcs":[[317]],"type":"Polygon"},{"arcs":[[318],[319]],"type":"Polygon"}]}}}
},{}],18:[function(_dereq_,module,exports){
module.exports={"type":"Topology","arcs":[[[31325,5416],[22,-45],[-56,-2],[1,36],[33,11]],[[30948,5749],[17,13],[53,-13],[69,-32],[2,-71],[-27,-29],[-32,-66],[-38,-43],[-45,-12],[-76,0],[-28,19],[-72,0],[-35,-20],[-67,9],[-15,14],[6,118],[-25,72],[26,37],[3,32],[28,12],[68,-15],[172,-9],[16,-16]],[[12551,25254],[-51,-1],[-1,42],[26,97],[51,-57],[24,-55],[-49,-26]],[[12542,25111],[-31,7],[-40,41],[37,58],[54,-3],[21,-30],[-11,-52],[-30,-21]],[[12119,26126],[-40,-11],[-38,19],[-15,35],[35,10],[58,-53]],[[8333,27067],[-69,-6],[11,31],[42,19],[16,-44]],[[3330,27230],[26,-26],[56,7],[30,-17],[11,-46],[-6,-79],[23,-21],[10,-59],[-98,-12],[-42,-21],[-24,-38],[-43,31],[-79,15],[-98,61],[-43,13],[-82,101],[51,14],[115,-8],[14,40],[85,44],[31,-9],[63,10]],[[6986,25919],[13,-25],[43,29],[51,10],[25,-35],[-18,-25],[20,-45],[54,-36],[-54,-62],[-102,-1],[28,-55],[-55,-13],[-47,-42],[-62,-6],[-63,-45],[-40,-56],[1,-35],[-31,-53],[-60,-43],[-27,24],[77,86],[-57,31],[-18,43],[-72,-1],[39,-64],[-21,-48],[-22,4],[-44,60],[-20,86],[-37,67],[9,56],[38,58],[70,41],[69,-4],[62,-140],[-3,122],[34,35],[-51,48],[-7,32],[37,31],[30,-12],[38,-65],[46,38],[44,8],[11,78],[85,-39],[-13,-37]],[[13303,24927],[-24,-55],[-10,-64],[7,-67],[24,-116],[-38,-82],[-11,-57]],[[13251,24486],[-100,-132]],[[13151,24354],[-76,20],[-55,179],[18,44],[64,28],[-29,29],[-6,132],[-11,63],[-58,115],[-44,27],[-96,-41],[-41,-29],[-24,-117],[-31,-58],[-38,18],[-28,78],[60,107],[44,118],[-49,112],[-44,21],[-34,51],[-14,54],[-33,12],[-3,60],[-42,16],[-45,50],[-4,75],[-27,3],[-153,64],[9,84],[-29,112],[-30,45],[26,27],[121,-66],[-110,107],[5,76],[-44,-42],[-24,7],[-42,71],[-38,40],[-9,71],[-40,11],[-76,69],[-49,9],[-31,35],[-46,99],[-6,49],[-47,52],[-23,119],[-39,128],[-35,-9],[-5,-54],[32,-59],[15,-91],[35,-122],[43,-187],[-25,-37],[-45,11],[-40,62],[-34,18],[-83,-9],[3,94],[-51,108],[-50,-12],[-62,42],[-32,48],[-73,22],[-33,-9],[-30,-47],[50,4],[86,-49],[43,-44],[-20,-43],[72,-7],[33,-50],[13,-62],[-20,-17],[-92,-12],[-33,-42],[-32,11],[-103,65],[-131,96],[-5,23],[-52,41],[-76,120],[-90,75],[-57,21],[-40,30],[-176,97],[-127,87],[52,33],[30,59],[-22,63],[18,86],[-92,-95],[-101,-54],[-56,-8],[-117,20],[-156,79],[31,39],[-9,53],[-56,-45],[-38,-13],[-243,59],[-109,6],[-155,-21],[-83,-22],[-94,2],[16,35],[-68,54],[-85,17],[-72,57],[22,102],[-44,1],[-73,-80],[-87,31],[-93,19],[62,89],[-131,5],[-28,45],[-70,-25],[31,57],[-68,44],[15,118],[-80,-63],[-106,-14],[-18,-27],[-75,28],[-17,-36],[-65,-10],[-32,41],[45,111],[-101,-80],[-37,-90],[-59,-26],[80,-53],[-8,-86],[68,18],[24,-42],[-61,-47],[-9,-36],[5,-80],[-37,-23],[-27,-68],[-114,-20],[-77,43],[-40,-18],[-46,1],[-32,-55],[-68,-19],[3,-54],[-64,13],[14,-44],[-64,-64],[-39,8],[-41,-24],[-26,15],[-27,-70],[-48,-44],[-30,-50],[-112,1],[-70,-30],[-33,1],[-62,49],[31,62],[44,29],[81,28],[58,62],[27,78],[-86,-66],[-30,-5],[-70,25],[-24,42],[32,105],[79,110],[38,135],[-12,99],[9,43],[102,50],[141,94],[66,-35],[57,-11],[39,14],[62,-5],[127,-37],[8,29],[-124,29],[-104,66],[39,95],[109,99],[-72,-15],[-36,-30],[-34,-61],[-45,-13],[-100,-4],[-39,21],[-92,-53],[-56,-57],[-123,-55],[-39,-35],[-13,-40],[9,-40],[-68,-37],[-76,-79],[3,-61],[-30,-37],[-79,-50],[-66,1],[63,-59],[10,-44],[-36,-61],[-95,-24],[16,-74],[-88,-36],[-14,41],[-65,-47],[12,-17],[-57,-72],[-76,-57],[6,-13],[-31,-93],[14,-18],[131,-42],[65,-40],[23,-53],[-28,-52],[-50,-52],[-67,-34],[-44,-48],[-18,-66],[-34,-39],[-11,-66],[-89,-21],[-4,-32],[-116,-21],[-26,-52],[-59,-55],[-76,-37],[-52,-85],[-62,-12],[-15,-50],[-52,1],[-65,-62],[23,-45],[-22,-69],[-42,-49],[-42,-2],[-99,-92],[-70,-9],[-27,-25],[-25,-64],[-80,5],[-44,-27],[14,-25],[-58,-33],[-62,-22],[-34,-52],[41,-17],[31,-54],[-69,-63],[-28,42],[-45,-86],[-177,-74],[-26,-17],[-12,-58],[-48,76],[-61,-25],[-78,-69],[-35,-13],[-25,-40],[-61,-9],[-24,-29],[-35,16],[-56,-56],[-78,-17],[-27,13],[1,35],[29,46],[-24,38],[-39,-19],[-23,-32],[-16,-71],[-92,-129],[-26,2],[-46,-48],[-34,14],[8,35],[-31,50],[-29,-13],[9,-75],[-17,-38],[-53,-22],[-36,47],[-34,10],[-3,-76],[-44,-41],[3,146],[46,37],[42,-5],[124,136],[45,74],[53,64],[61,54],[66,42],[133,59],[19,-36],[49,9],[-10,-42],[43,-59],[25,42],[63,5],[51,-33],[10,33],[-47,37],[-16,36],[44,120],[21,38],[143,126],[139,65],[81,88],[26,-22],[61,-11],[4,125],[20,38],[73,94],[76,58],[42,52],[79,46],[-20,19],[-1,54],[18,77],[2,80],[16,46],[31,16],[-22,128],[17,42],[64,70],[70,49],[-18,19],[21,111],[-49,-57],[-145,-67],[-98,-55],[-47,-13],[-31,14],[-55,107],[21,76],[57,19],[-56,27],[-25,-8],[-45,-74],[-23,11],[-28,-119],[25,-101],[-5,-41],[-45,-19],[-36,34],[-76,128],[-86,98],[-69,-47],[-64,44],[-59,74],[-81,-49],[-44,-42],[-30,0],[-80,-37],[-31,-29],[-9,-38],[-109,-29],[-108,16],[80,37],[36,40],[-15,53],[-4,108],[-26,-17],[-25,36],[-12,71],[28,41],[25,77],[1,38],[-22,64],[-64,136],[-29,102],[-49,54],[78,169],[-35,-10],[-56,-102],[-36,-50],[27,-88],[-19,-70],[-45,2],[-40,-36],[-94,-40],[-128,-22],[-62,2],[-65,47],[3,49],[-94,80],[-53,79],[-38,2],[-72,54],[10,45],[-91,21],[93,101],[33,69],[26,9],[120,-49],[29,-36],[-12,-61],[85,82],[29,-11],[46,-79],[55,38],[30,47],[-65,62],[-26,47],[-67,-54],[-124,3],[-86,31],[-86,-5],[-31,23],[83,62],[-55,4],[-87,60],[3,-54],[-51,-2],[-35,101],[-48,19],[-11,34],[30,45],[-46,31],[-52,-3],[14,51],[73,18],[-67,63],[129,34],[-37,73],[11,45],[71,104],[70,87],[97,55],[14,22],[-14,91],[22,87],[24,25],[83,18],[-41,37],[32,45],[83,24],[103,-35],[33,-35],[71,-41],[82,18],[70,65],[52,30],[93,135],[47,3],[39,-42],[129,8],[66,15],[46,31],[75,88],[14,45],[-35,109],[-23,113],[-64,74],[-46,23],[-8,44],[61,-5],[72,32],[27,52],[-14,59],[-49,55],[-34,10],[-77,-66],[-81,11],[-30,-37],[-84,-34],[-45,-33],[-131,-122],[-19,73],[-90,69],[-27,-23],[69,-45],[-26,-49],[-94,65],[-63,19],[-163,-2],[-164,-63],[-66,2],[-85,25],[-193,36],[-50,22],[-43,52],[20,51],[-2,51],[-37,13],[-77,74],[18,19],[64,11],[22,47],[70,30],[-114,24],[-15,-7],[-204,43],[-161,74],[-28,46],[140,52],[49,52],[91,9],[90,89],[97,48],[50,13],[57,-25],[79,-5],[45,28],[-77,40],[18,38],[90,47],[106,14],[108,60],[59,17],[111,11],[67,-14],[11,-26],[-35,-79],[2,-46],[-38,-37],[93,-67],[145,-5],[79,12],[82,-24],[103,10],[110,-9],[72,101],[28,16],[70,-32],[22,59],[-117,37],[-78,-18],[-25,20],[8,42],[-83,103],[-77,21],[-38,82],[67,27],[30,-15],[34,-60],[31,-9],[-9,-60],[39,-55],[88,-51],[71,19],[79,-12],[73,-46],[37,-5],[116,24],[-9,78],[-27,20],[-78,-4],[-61,34],[-52,-9],[-95,-52],[-48,21],[-79,55],[-6,53],[71,90],[-95,36],[-118,-15],[-168,4],[-147,38],[-52,48],[-61,144],[-51,67],[-348,224],[-158,57],[-76,62],[-94,22],[54,54],[54,190],[-7,45],[193,-9],[129,8],[205,28],[120,51],[91,68],[104,112],[18,114],[39,75],[167,173],[129,121],[67,-50],[178,36],[100,73],[142,74],[43,-12],[-40,-49],[28,-12],[-25,-59],[67,-5],[11,89],[19,17],[-60,53],[78,78],[101,47],[88,-39],[104,-1],[37,21],[134,2],[109,49],[113,78],[114,116],[87,-14],[182,-53],[63,-34],[-62,-65],[-94,-34],[76,-48],[81,30],[113,105],[62,-8],[75,-49],[-31,-48],[51,-23],[112,-24],[76,39],[116,7],[74,21],[123,-29],[81,3],[72,-37],[-57,-39],[-9,-41],[82,-50],[112,2],[-49,-55],[199,-17],[27,17],[128,30],[72,-34],[69,0],[78,33],[166,-4],[161,-43],[57,-51],[64,20],[105,-27],[45,-44],[177,-24],[88,11],[128,-3],[127,-14],[106,-56],[67,-21],[159,-14],[56,29],[98,8],[87,25],[108,-7],[38,14],[141,-42],[114,-84],[165,-51],[80,-60],[112,-2]],[[10271,31976],[0,-299],[0,-449],[0,-449],[0,-449],[0,-449],[0,-299],[0,-299],[0,-300],[0,-449],[0,-299],[0,-299],[0,-449],[0,-300],[132,-42],[20,42],[132,-60],[82,75],[165,7],[0,-32],[-29,-101],[39,-46],[93,-47],[21,-63],[287,-283],[42,-165],[45,43],[128,83],[65,1],[31,66],[0,91],[33,-1],[27,39],[-20,37],[107,30],[126,67],[118,-110],[-6,-71],[35,-78],[74,-46],[65,-58],[40,-96],[22,-25],[109,-74],[116,-151],[-6,-37],[41,-59],[99,-177],[101,-193],[88,-154],[-29,-61],[85,-28],[-20,-87],[65,-32],[12,-104],[69,4],[136,-99],[94,-33],[25,-39],[48,-17],[18,-56],[54,-20],[33,14],[20,-49],[0,-65]],[[4065,24463],[32,-126],[-62,-56],[-135,-3],[-45,-25],[-63,-74],[-35,-12],[-64,8],[-23,64],[5,33],[50,43],[49,96],[29,17],[41,-7],[112,71],[72,6],[37,-35]],[[12397,24751],[6,-24],[-40,-34],[-10,-28],[-51,-55],[4,74],[-28,43],[28,22],[50,-9],[41,11]],[[12940,24513],[29,-90],[-48,-5],[-7,53],[-45,48],[5,56],[18,32],[30,-25],[18,-69]],[[12726,25041],[-5,-85],[-43,-8],[-37,26],[-8,40],[-58,11],[-11,61],[28,22],[60,132],[18,-6],[69,-124],[-13,-69]],[[13039,24723],[-28,-114],[-29,-31],[-36,32],[-29,4],[-15,50],[-92,-77],[-10,99],[55,86],[6,126],[98,64],[41,-51],[43,-96],[-4,-92]],[[12380,25499],[82,-15],[32,-48],[-2,-105],[-20,-28],[-34,-5],[3,-43],[-19,-32],[-107,-1],[-13,27],[-9,144],[-39,65],[-43,44],[31,30],[44,-3],[94,-30]],[[12519,24418],[25,17],[43,-18],[-8,-69],[-17,-37],[-28,12],[-55,75],[-32,59],[-35,105],[-46,18],[-3,48],[36,11],[55,-57],[28,-53],[4,-40],[33,-71]],[[8782,27264],[59,-11],[15,-19],[-136,-60],[-23,69],[39,37],[46,-16]],[[8412,26938],[-38,8],[16,37],[57,69],[39,31],[96,136],[37,-55],[-92,-92],[-50,-95],[-65,-39]],[[8433,27264],[0,-50],[-20,-67],[-39,4],[21,118],[38,-5]],[[6956,25561],[-41,-16],[-39,-48],[-22,26],[25,69],[97,-9],[-20,-22]],[[7100,26257],[-42,29],[48,40],[29,-12],[-35,-57]],[[6249,24893],[-47,4],[32,43],[15,-47]],[[6493,25208],[-13,34],[59,50],[20,-14],[-66,-70]],[[4836,24634],[28,-4],[26,-63],[-34,-13],[-59,8],[-6,70],[45,2]],[[4390,24394],[-37,43],[38,27],[17,-26],[-18,-44]],[[5060,24538],[-46,-43],[-31,7],[-6,33],[37,36],[39,15],[7,-48]],[[3411,23996],[-53,-12],[-15,63],[54,19],[52,-44],[-38,-26]],[[3489,24030],[-26,60],[29,16],[23,-53],[-26,-23]],[[4771,26305],[-42,-4],[-17,50],[40,35],[75,30],[-56,-111]],[[2348,23370],[-8,-28],[-72,30],[63,16],[17,-18]],[[2061,23234],[-26,10],[40,49],[26,-28],[-40,-31]],[[1582,23075],[-42,1],[41,59],[43,-30],[-42,-30]],[[2331,25310],[77,-21],[-31,-27],[-46,48]],[[600,22937],[-33,1],[-4,49],[34,-9],[3,-41]],[[289,22791],[-23,-12],[-96,0],[67,42],[35,34],[12,34],[28,-3],[-19,-49],[-4,-46]],[[11937,25677],[41,-100],[56,-224],[2,-59],[-12,-41],[7,-113],[-33,-32],[-27,42],[-37,100],[-10,73],[-39,33],[-1,51],[-46,-1],[4,56],[32,49],[-39,28],[-10,54],[-35,29],[-29,-89],[-44,15],[-15,67],[9,44],[57,38],[26,64],[39,9],[104,-93]],[[12017,26092],[70,-12],[52,3],[47,-77],[45,-105],[-35,16],[-60,117],[-16,-7],[11,-77],[91,-155],[-10,-67],[14,-59],[-49,-18],[-44,-79],[-48,-46],[-43,18],[4,72],[25,128],[-48,80],[-26,182],[-19,78],[-24,60],[20,29],[43,-81]],[[11727,26134],[40,-50],[-8,-79],[74,68],[95,-38],[21,-50],[-12,-70],[-37,-12],[-35,12],[-5,-43],[75,-4],[29,-69],[-16,-56],[-43,16],[-114,77],[-30,-6],[-3,-87],[-20,-31],[-61,14],[-46,116],[-105,102],[-30,51],[15,63],[40,24],[62,6],[98,59],[16,-13]],[[12325,25158],[100,-10],[27,-43],[2,-74],[94,-49],[44,-51],[47,-107],[39,-66],[-62,19],[-22,-67],[48,10],[56,-51],[16,-43],[-14,-39],[64,-8],[-7,-89],[6,-117],[-12,-41],[-31,-7],[-36,46],[-21,54],[-57,23],[-11,51],[-54,-2],[-88,160],[40,18],[-35,56],[-2,55],[-57,-3],[-20,42],[-51,1],[-35,49],[27,26],[38,-20],[34,24],[22,39],[-13,50],[-22,9],[-87,-48],[-19,29],[62,75],[-18,36],[8,63]],[[12208,25417],[44,-32],[25,-67],[-32,-35],[-8,-50],[0,-99],[-18,-84],[-33,3],[-33,-29],[-16,65],[14,107],[31,22],[-57,64],[-26,73],[3,59],[27,41],[36,7],[43,-45]],[[7119,26194],[28,27],[65,-52],[27,-35],[-35,-42],[-44,46],[-11,-60],[-63,-17],[-68,-44],[-55,-9],[-110,46],[111,107],[57,-9],[-19,70],[46,18],[36,-12],[35,-34]],[[2825,23625],[-84,-55],[-65,-100],[-155,-107],[75,108],[10,67],[33,41],[56,1],[11,77],[30,47],[83,30],[50,-38],[-16,-50],[-28,-21]],[[3198,23910],[32,-9],[35,59],[39,-34],[-71,-76],[-17,-43],[29,-25],[-71,-59],[-19,-31],[-106,-35],[-63,-34],[-39,-33],[-40,-9],[-39,33],[27,24],[50,12],[109,68],[12,56],[21,32],[-28,89],[39,46],[67,20],[29,-3],[4,-48]],[[1282,23006],[88,-35],[-63,-14],[-58,11],[4,42],[29,-4]],[[971,22954],[-148,-21],[27,27],[124,45],[100,42],[13,27],[-49,24],[74,53],[41,-45],[-10,-44],[-34,-23],[16,-34],[-61,-29],[-93,-22]],[[442,22868],[43,-24],[0,-34],[-92,-64],[-18,24],[-10,59],[38,25],[10,75],[41,-22],[-12,-39]],[[87,22756],[-50,12],[20,47],[-57,61],[22,17],[44,1],[43,-40],[8,-32],[-30,-66]],[[28579,19353],[192,0],[230,0],[230,1],[229,0],[28,99],[25,46],[53,-14],[37,38],[29,-38],[27,73],[30,12],[-3,63],[29,47],[48,49],[4,32],[31,54],[-3,77],[16,99],[50,98],[16,137],[212,386],[53,-19],[0,-78],[31,-32],[30,-4],[96,42],[47,36],[59,-41],[80,-99],[3,-299],[3,-408],[53,-42],[45,-9],[2,-37],[-19,-29],[18,-49],[-12,-52],[29,-52],[32,13],[35,-16]],[[30674,19437],[12,-92],[-7,-53],[33,-30],[-57,-78],[-47,11],[-54,-26],[-26,-43],[-52,8],[-18,-63],[-42,-36],[-17,55],[-44,8],[-26,-32],[-22,32],[-19,-65],[-25,-36],[-50,37],[23,31],[-38,22],[-32,-26],[2,-44],[-30,-90],[-1,-39],[-44,-57],[-32,8],[-38,-49],[-69,-27],[-17,-24],[-48,25],[-27,-10],[-32,-35],[-25,-49],[19,-23],[-88,-142],[-34,-110],[-25,-33],[-26,-125],[25,-80],[21,-27],[-47,-33],[-27,-62],[-32,-51],[13,-16],[72,-37],[33,-96],[-11,-28],[30,-24],[9,-70],[61,-38],[81,50],[-24,79],[31,-10],[10,-78],[-2,-67],[-31,0],[-95,-26],[-21,-23],[-49,-25],[-3,90],[-9,3],[-75,-85],[-60,-17],[-4,64],[-35,62],[-26,-66],[-5,-92],[-22,-38],[-68,-25],[-44,5],[-93,-25],[-29,11],[-30,-19],[-102,-5],[-21,10],[-28,-35],[-43,-21],[-111,-79],[-74,-98],[-27,4],[-21,-42],[-56,-76],[-10,-41],[6,-37],[53,-13],[26,-52],[-20,-131],[-25,-69],[-16,-108],[-22,-58],[-37,-65],[-11,-51],[-48,-49],[-12,-43],[-41,-106],[-44,-26],[16,100],[-42,33],[-24,-1],[-26,39],[-34,28],[-47,77],[-14,-6],[45,-100],[5,-97],[23,-65],[35,-75],[26,-22],[-11,-74],[26,-66],[-1,-39],[-26,-22],[-7,-66],[-19,-7],[-35,-90],[-67,-223],[3,-37],[-31,-22],[-32,-46],[-11,-66],[-22,-74],[-18,57],[6,69],[24,113],[27,70],[20,34],[17,68],[-53,9],[15,59],[-37,42],[18,47],[-52,9],[-44,42],[-23,69],[9,54],[57,10],[-67,101],[36,18],[-40,80],[35,-9],[-14,103],[23,63],[55,74],[-1,64],[-29,-12],[-9,-66],[-50,-56],[-35,-10],[-5,-40],[-42,15],[43,-93],[-15,-22],[-13,-68],[-5,-79],[10,-108],[24,-58],[-7,-30],[30,-113],[-72,45],[-82,33],[-27,52],[-25,24],[-43,-52],[51,9],[39,-82],[72,-32],[27,-20],[21,-43],[33,-25],[25,-35],[-12,-88],[-10,-24],[-41,3],[-83,131],[21,-65],[26,-28],[38,-65],[49,-29],[10,-39],[2,-71],[-38,15],[0,-89],[20,-33],[12,-49],[-32,-31],[-56,77],[-8,41],[-20,-2],[-84,51],[23,-48],[70,-26],[10,-65],[36,-44],[5,-33],[24,-3],[43,32],[68,-21],[30,-131],[36,-219],[55,-179],[-6,-4],[-40,119],[-23,86],[-23,151],[-18,2],[-9,-50],[34,-105],[-4,-47],[-73,53],[-1,-69],[-21,-15],[-43,9],[-27,-54],[-22,-6],[-33,30],[-13,-60],[32,-7],[70,5],[42,19],[37,-10],[3,-47],[-7,-97],[23,17],[6,89],[35,32],[22,-30],[8,-69],[-8,-61],[-53,-71],[-38,-66],[-19,-13],[-60,24],[-28,-2],[-11,57],[-22,11],[-7,-39],[-29,-11],[63,-83],[-26,-60],[-6,-41],[-42,-42],[10,-26],[79,25],[27,-27],[-43,-81],[-53,-13],[-7,-23],[-87,-5],[-24,5],[-32,-47],[-30,2],[-5,-48],[-39,-38],[-64,-86],[-48,-108],[-24,-83],[-109,3],[-47,-22],[-73,-77],[-22,-33],[-60,-130],[-15,-83],[-23,-55],[-40,-48],[-46,-21],[-8,-47],[-52,-63],[-36,12],[11,-42],[-35,-55],[-58,-23],[-36,-38],[-58,16],[27,-56],[-3,-37],[-34,-30],[-11,53],[-32,-50],[19,-42],[-17,-38],[-26,-15],[-20,-86],[-53,-34],[14,-33],[-9,-31],[-28,-25],[6,-89],[-33,-82],[-23,-9],[25,-46],[-21,-47],[-21,14],[-3,-57],[-19,-109],[5,-74],[12,-46],[20,-190],[14,-66],[24,-178],[40,-172],[56,-209],[93,-253],[11,-36],[-27,-85],[-30,73],[12,60],[-22,27],[-1,51],[-19,13],[25,-193],[27,-102],[117,-499],[28,-63],[10,-46],[11,-95],[2,-123],[-19,-224],[-4,-152],[-26,-47],[-22,-63],[-8,-99],[-10,-49],[-33,-52],[-20,2],[-50,-39],[-34,10],[-42,-22],[-27,2],[-15,47],[19,46],[-36,136],[-33,76],[-5,50],[-56,31],[-41,47],[-26,84],[-16,148],[-47,60],[-11,78],[3,96],[-8,36],[-23,-13],[0,-50],[-31,16],[-41,97],[-50,175],[-26,50],[22,13],[65,159],[-13,35],[-20,-13],[-16,41],[-9,-93],[-14,-30],[-23,-5],[-27,70],[26,201],[24,127],[3,206],[-33,84],[-144,205],[-111,243],[-97,91],[-73,-20],[-13,-18],[-7,-63],[-47,-5],[-69,-64],[-24,3],[-39,-29],[-80,-21],[-16,8],[19,52],[-13,40],[-41,51],[-48,75],[20,85],[-38,-22],[-4,-40],[-115,85],[-78,34],[60,15],[-5,33],[-53,2],[-64,-48],[-144,-33],[22,30],[38,17],[5,37],[-29,-1],[-27,20],[-4,-44],[-27,-59],[-54,-23],[-10,38],[-30,-53],[-47,14],[-18,59],[-19,22],[-6,91],[-23,28],[-29,-162],[-59,4],[-95,-9],[-58,30],[-88,-42],[-27,7],[-34,-63],[-40,-29],[-101,52],[-25,43],[-50,14],[-29,-52],[-23,-70],[36,-39],[30,-18],[50,15],[28,34],[22,-1],[21,25],[20,-28],[-42,-57],[20,-39],[43,-8],[7,45],[19,28],[22,-24],[16,-46],[1,-51],[-49,-25],[-53,-78],[13,-42],[45,-60],[74,-45],[17,1],[18,-44],[28,-24],[-1,-30],[-25,-23],[-13,-42],[-22,33],[-26,-42],[-6,36],[-25,64],[-49,55],[-47,16],[-8,43],[-16,21],[-78,41],[5,-30],[25,-25],[0,-48],[-14,-79],[-31,-40],[-24,78],[-21,20],[-35,2],[-23,-14],[-25,-62],[-20,-10],[-70,32],[-79,48],[37,32],[-24,55],[-24,29],[-51,21],[-44,98],[-42,2],[-19,44],[-34,-19],[-33,-51],[15,-55],[-49,-18],[-113,20],[-34,20],[-44,40],[-62,33],[-143,-4],[-36,-23],[-16,43],[30,54],[-19,33],[-29,-84],[16,-64],[-58,-10],[-180,-129],[62,66],[-22,10],[-47,-10],[15,55],[-5,49],[-25,1],[-16,-39],[-37,13],[8,-88],[29,-82],[-69,-104],[-4,-45],[-33,-59],[-106,-112],[-54,-54],[-46,-27],[-71,35],[-48,-31],[-35,62],[-17,-4],[31,-112],[29,-16],[-39,-47],[-32,-13],[-26,41],[-10,-103],[-23,-32],[-21,16],[-52,-23],[4,-43],[30,17],[-11,-55],[-27,-54],[-22,-13],[-34,8],[-16,-17],[39,-85],[-25,-129],[-16,-47],[-24,-7],[-43,41],[-3,-55],[57,-25],[3,-62],[-22,-77],[24,-141],[9,-105],[9,-46],[52,-168],[18,-2],[1,-53]],[[22383,9599],[-38,-10],[-21,-36],[-22,10],[-41,48],[-59,29],[-78,12],[-53,24],[-28,36],[-30,22],[-31,7],[-46,52],[-30,21],[-39,9],[-26,25],[-18,60],[-16,101],[-20,63],[-39,78],[1,68],[-19,87],[7,65],[-6,42],[-25,45],[-44,48],[-37,70],[-31,93],[-30,64],[-31,36],[-20,43],[-11,52],[1,38],[-60,163],[-24,78],[-6,49],[-27,58],[-72,112],[-6,30],[-73,89],[-22,56],[-16,17],[-152,10],[-48,15],[-29,24],[-20,-3],[-12,-30],[-25,-20],[-39,-8],[-32,-56],[-28,-102],[-16,-117],[-36,-43],[-19,-46],[-21,-22],[-25,1],[-46,36],[-66,69],[-52,44],[-38,16],[-34,32],[-29,48],[-51,48],[-28,54],[-33,90],[-16,70],[0,74],[-43,160],[-39,102],[-82,82],[-66,89],[-83,133],[-58,81],[-34,27],[-30,48],[-25,69],[-30,46],[-303,3],[-182,2],[-1,-230],[-293,-1],[-195,-1],[-293,-1],[-262,152],[-262,151],[-262,152],[-262,151],[14,29],[17,77],[-189,-27],[-238,-33],[-237,-33]],[[16864,12965],[-2,60],[-29,7],[-8,73],[2,68],[-15,82],[-41,101],[-89,124],[-45,42],[-36,52],[-50,19],[-9,-24],[-32,16],[5,59],[-31,81],[-26,9],[-64,-5],[-86,45],[-26,27],[-8,47],[-41,42],[-53,41],[-29,-10],[-39,6],[-55,30],[-33,3],[-62,-8],[-24,6],[-45,56],[6,118],[-11,71],[8,65],[-20,41],[-41,27],[-8,33],[7,47],[-11,30],[-34,29],[-32,65],[-40,35],[-17,59],[-88,186],[-59,90],[-9,53],[-3,71],[36,82],[-5,61],[-20,45],[-78,26],[-64,111],[-4,85],[-25,87],[-4,117],[36,9],[4,-68],[39,-47],[29,-10],[-27,96],[-21,30],[-24,87],[46,41],[36,5],[102,-8],[-9,23],[-36,-2],[-42,24],[-33,-29],[-67,40],[-25,-18],[-3,-80],[8,-59],[-81,55],[-57,78],[-5,92],[-17,14],[-21,74],[-46,45],[-38,71],[-76,119],[-5,104],[-28,132],[12,75],[-2,53],[-13,81],[-14,43],[-62,120],[-60,81],[-13,122],[24,113],[17,33],[25,99],[-2,96],[20,117],[-14,121],[-12,50],[-23,35],[9,101],[-39,71],[-18,132],[3,104],[-36,117],[23,102],[30,169],[16,35],[29,124],[10,20],[5,187],[8,142],[15,47],[-5,49],[5,65],[-4,67],[31,319],[5,90],[-9,136],[4,153],[10,21],[66,0],[42,21],[-46,39],[-58,-16],[-17,17],[-31,-11],[6,108],[29,-30],[15,117],[-50,43],[-11,61],[74,51],[-40,11],[-15,23],[-34,-7],[-9,99],[-31,100],[-18,130],[-23,65],[-45,61],[-11,36],[-11,91],[6,69],[13,45],[56,-38],[70,-30],[55,-38],[188,-25],[47,16],[28,-35],[45,4],[23,25],[11,-64],[23,-68],[-34,-73],[-12,28],[-63,-123],[19,3],[62,73],[47,79],[17,-44],[-29,-39],[21,-148],[-7,-57],[-64,18],[-30,-20],[-19,-60],[21,-21],[28,8],[30,-18],[27,28],[17,56],[34,19],[18,30],[-4,119],[-11,25],[4,86],[42,98],[-31,52],[-38,9],[-8,50],[31,34],[-22,41],[-50,46],[48,21],[-5,84],[-13,56],[-25,-8],[-38,118]],[[15301,21396],[382,0],[235,0],[236,0],[235,0],[353,0],[354,0],[235,0],[236,0],[235,0],[353,0],[354,0],[235,0],[236,0],[235,0],[236,0],[235,0],[353,0],[354,0],[235,0],[354,0],[235,0],[353,0],[236,0],[353,0],[353,0],[419,-1],[2,194],[59,-11],[24,-23],[3,-95],[11,-59],[50,-134],[57,-19],[99,-23],[57,-27],[39,-42],[68,19],[23,28],[34,6],[72,-9],[73,-40],[64,-50],[23,-81],[32,27],[35,4],[46,-18],[41,-54],[94,-71],[36,0],[46,24],[49,45],[56,4],[27,-49],[117,-3],[63,10],[27,-53],[115,-7]],[[24473,20888],[-69,-52],[-102,-58],[-111,-51],[-99,-84],[-119,-158],[-71,-84],[-116,-119],[-2,-42],[23,-19],[42,1],[93,33],[39,23],[39,38],[25,-9],[50,44],[29,5],[20,-28],[-27,-76],[-16,-80],[43,21],[9,23],[74,-54],[87,36],[81,71],[45,23],[94,22],[56,55],[57,28],[32,53],[50,44],[5,-34],[32,-36],[17,-68],[-6,-65],[72,59],[121,-32],[39,-32],[77,-155],[52,-5],[37,15],[51,-37],[33,10],[26,-19],[29,36],[83,68],[79,23],[91,-6],[48,9],[45,23],[73,16],[5,-22],[-10,-120],[70,-19],[23,11],[46,-28],[29,28],[27,0],[27,-64],[7,-61],[45,-81],[40,-45],[-23,-20],[-100,20],[-50,-3],[-30,30],[-47,-92],[-47,59],[-80,51],[-60,10],[-38,-45],[-95,-41],[-45,15],[-41,-17],[-35,-69],[-36,-25],[-32,-66],[-13,40],[36,54],[-1,39],[-65,-26],[-10,-35],[-38,-21],[-31,3],[-36,-62],[-50,-125],[-57,-106],[-6,-64],[-57,-46],[-54,-166],[16,-34],[50,51],[51,100],[51,10],[23,-46],[-42,-137],[-11,-63],[-2,-108],[-34,-59],[-18,-78],[7,-94],[-21,-64],[-35,-181],[2,-105],[29,-119],[-6,-136],[4,-170],[28,-65],[25,-115],[25,-71],[28,-33],[44,-10],[50,17],[75,61],[46,63],[62,160],[30,99],[19,119],[3,128],[-8,66],[-22,100],[-29,94],[-30,118],[28,80],[-2,66],[-19,77],[28,49],[33,83],[15,113],[1,86],[31,14],[19,87],[29,25],[26,-2],[47,85],[28,17],[-5,-98],[-12,-68],[45,2],[27,98],[3,125],[43,41],[70,19],[-33,58],[0,51],[30,51],[7,36],[50,10],[91,-61],[32,-5],[38,-28],[17,-41],[179,-99],[23,-33],[21,-69],[-35,-59],[3,-32],[30,-47],[9,-96],[-10,-158],[-36,-56],[-14,3],[-28,-109],[-22,-22],[-51,-27],[-12,-44],[-5,-71],[25,-38],[31,-19],[32,5],[13,31],[29,18],[36,101],[36,33],[72,38],[44,-33],[34,-76],[19,-113],[37,-283],[12,-33]],[[26450,18335],[-19,-142],[-16,-59],[-28,-34]],[[26387,18100],[4,53],[-15,18],[-28,-17],[-2,-38],[-21,-36],[-15,-74],[-50,-45],[-43,-174],[-28,-35],[-26,-70],[151,-105],[22,22],[36,-30],[49,-59],[111,44],[110,9],[76,83],[61,50],[99,63],[148,76],[174,133],[53,46],[26,34],[57,50],[41,64],[51,61],[-12,77],[3,62],[-43,19]],[[27376,18381],[1,90]],[[27377,18471],[103,40],[64,14],[78,3],[93,-19],[59,-40],[26,-6],[72,14],[81,-9],[78,43],[61,75],[51,15],[24,19],[-9,115],[10,39],[27,47],[-27,60],[-23,-19],[-20,31],[32,39],[112,110],[9,48]],[[28278,19090],[110,144],[61,65],[89,54],[41,0]],[[29186,17295],[-1,-37],[63,56],[37,14],[14,-19],[-65,-61],[-60,-35],[-57,-25],[-120,-63],[-19,4],[-98,-32],[-41,-4],[-15,33],[-53,-39],[-5,30],[42,78],[85,66],[55,12],[26,-10],[46,17],[135,18],[31,-3]],[[30380,19009],[-62,-20],[1,36],[30,47],[17,-9],[14,-54]],[[29739,17495],[-76,-26],[30,63],[16,4],[30,-41]],[[23861,11412],[-10,-7],[-46,44],[33,35],[30,-29],[-7,-43]],[[22419,10593],[10,57],[22,53],[16,-14],[-48,-96]],[[22376,9701],[-27,87],[-37,251],[14,-10],[41,-257],[9,-71]],[[16527,13116],[-34,6],[-10,62],[44,-68]],[[16059,13675],[-34,-1],[-23,50],[49,6],[25,-27],[-17,-28]],[[16104,13757],[88,-37],[-68,-20],[-21,14],[1,43]],[[16528,13402],[13,-38],[-41,3],[-12,49],[40,-14]],[[15303,21232],[-31,-40],[-25,16],[24,43],[32,-19]],[[15361,20967],[57,-87],[-13,-35],[-40,38],[-13,70],[-40,57],[22,57],[27,4],[8,-33],[-41,-27],[33,-44]],[[27066,10274],[-21,49],[-31,138],[6,18],[46,-205]],[[25171,19506],[-3,-103],[-33,-61],[-2,-31],[-26,-45],[-28,15],[-4,46],[23,39],[15,55],[38,45],[20,40]],[[26143,19875],[3,-34],[-58,-7],[-38,18],[71,60],[22,-37]],[[25977,20037],[-16,10],[-23,65],[36,23],[-6,-62],[9,-36]],[[24966,20586],[-63,-10],[1,-37],[-57,-56],[-4,-27],[-56,-76],[2,58],[-50,26],[14,50],[47,60],[78,43],[64,6],[24,-37]],[[24657,20854],[-16,-28],[-52,-20],[3,57],[154,105],[-10,-59],[-79,-55]],[[26740,6212],[23,7],[43,-8],[-48,-38],[-28,3],[10,36]],[[29421,7450],[-47,13],[7,31],[38,-9],[2,-35]],[[29360,7488],[-38,40],[24,18],[14,-58]],[[29235,7494],[32,-28],[-32,-20],[0,48]],[[32317,3041],[-29,27],[3,40],[24,35],[16,-32],[-3,-47],[-11,-23]],[[32394,3381],[-16,-53],[-30,34],[-1,67],[33,81],[16,-32],[-2,-97]],[[31915,5132],[-58,51],[12,32],[23,-25],[23,-58]],[[27345,917],[-73,-2],[-36,-53],[-32,-27],[-35,-108],[18,-23],[-104,-125],[-20,-19],[-47,-13],[-25,-38],[0,-70],[26,-41],[29,-75],[51,-94],[18,-86],[-28,-34],[-48,-4],[-42,-77],[-63,-25],[-50,-3],[-15,29],[-4,82],[-33,141],[-8,96],[-26,-23],[-10,-95],[-20,-22],[-65,49],[-48,151],[-14,62],[-37,14],[-31,26],[-34,4],[-18,-15],[-17,18],[-3,42],[-36,-19],[-46,7],[-41,17],[-28,-9],[-24,-29],[4,-76],[-7,-14]],[[26323,436],[-1,30],[-18,65],[-21,30],[7,27],[38,44],[1,93],[-17,54],[10,33],[39,48],[0,27],[-39,53],[-16,2],[0,200],[38,73],[44,-44],[22,36]],[[26410,1207],[53,-76],[-1,-47],[10,-65],[29,-35],[-3,-56],[46,-50],[82,12],[19,31],[46,-98],[53,-24],[41,3],[39,14],[62,38],[45,69],[36,31],[116,65],[41,69],[34,17],[36,52],[20,41],[21,20],[61,-15],[40,-19],[44,-6],[-35,-261]],[[26675,58],[-46,19],[-23,42],[11,72],[23,14],[16,-50],[-6,-50],[22,-23],[3,-24]],[[26246,3981],[-16,-31],[-18,-61],[-18,46],[-19,-40],[11,-30],[21,-9],[31,-209],[-7,-38],[-37,-107],[-18,-31],[-23,-132],[-20,-215],[14,-192],[-7,-178],[8,-95],[-23,-8],[-9,33],[10,61],[-9,17],[-36,-57],[29,-110],[-5,-30],[-4,-104],[-23,18],[-1,-66],[-17,-36],[35,-19],[14,-93],[-25,-39],[-24,-71],[-10,-65],[10,-87],[32,-101],[20,-8]],[[26112,1894],[-5,-42],[-42,-47],[-29,-5],[-69,24],[-24,61],[-26,41],[-39,6],[-58,31],[-58,-55],[-174,113],[-23,10],[-33,-63]],[[25532,1968],[-24,70],[-36,68],[-141,208],[-51,125],[-28,89],[-26,47],[-76,96],[-17,38],[-75,127],[-58,75],[0,32],[23,40],[12,-2],[33,-61],[24,30]],[[25092,2950],[76,6],[28,32],[11,87],[13,22],[40,9],[5,15],[-17,197],[10,33],[34,6],[65,-10],[12,8],[50,115],[31,28],[15,-43],[55,-62],[29,71],[28,23],[28,41],[29,59],[54,61],[-4,53],[9,41],[34,42],[0,40],[14,34],[35,29],[19,-10],[17,-39],[23,-27],[30,-14],[23,5],[17,25],[35,3],[42,17],[13,28],[26,-7],[29,12],[32,31],[45,26],[14,36],[34,16],[71,-8]],[[22383,9599],[-5,-106],[-17,-86],[-55,-181],[-23,-112],[-44,-320],[-14,-209],[-3,-128],[-10,-218],[5,-186],[-19,-85],[-12,-78],[4,-58],[17,-119],[5,-89],[50,-152],[27,-53],[47,-72],[-6,-65],[-24,14],[11,52],[-20,24],[-37,76],[-45,136],[32,-217],[11,-35],[23,-29],[4,-41],[39,-151],[45,-155],[3,-43],[18,-52],[114,-219],[69,-163],[25,-155],[14,-48],[7,-66],[46,-74],[14,-48],[24,-27],[20,-80],[29,-48],[-3,-25],[53,-33],[18,15],[105,-10],[46,-66],[60,-29],[70,-174],[42,-5],[56,15],[88,56],[30,28],[58,37],[90,7],[28,-11],[67,23],[32,29],[16,45],[74,32],[63,5],[66,16],[26,-40],[-9,-56],[30,-30],[56,-12],[19,5],[25,44],[46,42],[-1,49],[-24,44],[6,49],[79,101],[29,27],[59,103],[13,193],[12,34],[39,60],[7,58],[0,272],[7,80],[29,152],[47,57],[100,79],[260,71],[36,17],[45,47],[33,16],[59,-2],[33,25],[35,-6],[88,-35],[57,-32],[63,-13],[18,15],[-8,54],[24,20],[26,-15],[34,-66],[24,-21],[2,-96],[13,-43],[-26,-136],[-17,-50],[-37,-80],[-44,-63],[-56,-142],[-12,-66],[0,-53],[9,-51],[-9,-38],[-33,-23],[-28,-73],[12,-43],[60,15],[-2,-41],[-14,-30],[-48,-47],[10,-52],[20,36],[13,-17],[-25,-124],[-17,-127],[-22,-73],[-8,-107],[-25,-91],[-8,3],[-21,85],[-27,43],[13,104],[-18,55],[-34,-59],[1,-39],[-28,-88]],[[24827,5763],[-46,3],[-17,-16],[-17,-80],[-61,-166],[-25,-26],[-43,43],[-22,-15],[-8,-79]],[[24588,5427],[-282,0],[-223,0],[-1,-289],[-115,2],[51,-73],[31,-70],[38,-56],[44,-41],[29,-40],[14,-40],[7,-61],[45,-37],[15,-25],[-12,-117],[4,-46],[-147,-1],[-209,0],[-125,-384],[-5,-42],[36,-84],[-19,-37],[-12,-94],[8,-65],[-21,-75]],[[23739,3752],[-82,151],[-76,153],[-60,88],[-39,70],[-104,155],[-53,71],[-50,84],[-45,47],[-45,31],[-37,41],[-14,-30],[34,-29],[30,-12],[46,-43],[7,-22],[-129,86],[-53,7],[20,59],[-17,18],[-28,-31],[-11,38],[-30,27],[-34,-57],[20,-49],[-51,-17],[-91,-104],[-85,-44],[-122,-100],[-54,-5],[-28,-16],[-82,38],[-104,94],[-157,29],[-106,123],[-106,50],[-67,118],[-40,5],[-26,19],[-95,42],[-95,29],[-92,102],[-61,33],[-52,41],[-115,70],[-42,39],[-40,60],[-66,62],[-28,51],[-31,19],[-45,97],[-23,42],[-41,25],[-62,-7],[-91,43],[-42,11],[-88,63],[-117,71],[-71,158],[-59,100],[-37,42],[-64,51],[-35,41],[-55,32],[-92,80],[-30,69],[-17,61],[-49,73],[-54,139],[-14,50],[-11,78],[-27,81],[7,26],[28,32],[46,7],[32,34],[2,46],[-21,43],[-45,13],[10,34],[60,141],[1,134],[7,56],[-61,66],[-27,102],[-33,87],[1,174],[-41,154],[-42,76],[-21,27],[-59,119],[-46,69],[-46,130],[-45,82],[-57,139],[-41,69],[-188,233],[11,0],[55,-57],[10,33],[-6,33],[-25,1],[-59,25],[-37,38],[-16,41],[-2,46],[-53,98],[4,28],[31,-6],[-3,39],[-78,51],[-25,37],[-79,81],[-10,58],[-67,-29],[-11,36],[27,19],[-10,28],[-52,-26],[-30,21],[-22,48],[-11,157],[42,107],[20,32],[-21,117],[-57,91],[-75,-4],[-21,34],[-17,58],[-15,101],[-92,42],[-55,82],[-23,66],[-15,104],[10,73],[14,35],[-64,26],[-44,-14],[-56,40],[-44,76],[-52,138],[-58,44],[-42,91],[-30,90],[-29,37],[-32,63],[-17,129],[-43,40],[-10,96],[-44,94],[-22,77],[-23,59],[-6,72],[-18,89],[-35,108],[-29,72],[-15,73],[6,74],[-5,71],[18,5],[-1,48],[-51,39],[-69,19],[-40,27],[-2,62],[-35,47],[-52,36],[-15,-53],[-41,-9],[-32,24],[-78,92],[-81,33],[14,-78],[-16,-56],[-9,-195],[33,-101],[16,-99],[5,-74],[14,-59],[-4,-137],[5,-42],[22,-69],[41,-64],[8,-34],[54,-49],[32,-64],[65,-87],[20,-37],[58,-136],[12,-90],[33,10],[17,-68],[32,-8],[31,-147],[16,-20],[45,-24],[1,-68],[20,-44],[-4,-58],[17,-50],[4,-87],[50,-86],[62,-69],[13,-90],[25,-82],[52,-55],[-2,-58],[35,-65],[6,-83],[30,-54],[15,4],[-33,90],[6,63],[62,-88],[7,-65],[21,-38],[2,-50],[36,-145],[1,-99],[10,-73],[39,-111],[32,-23],[6,-55],[33,-141],[38,-77],[19,-65],[3,-41],[-14,-60],[-2,-41],[21,-128],[31,-65],[35,-16],[14,-42],[20,48],[-8,57],[17,20],[67,-87],[11,-33],[47,-88],[28,-69],[9,-73],[46,-32],[25,-61],[2,-38],[-12,-98],[-11,-28],[-37,-41],[-27,-50],[-54,-49],[-23,4],[-22,57],[-26,172],[-17,36],[-13,54],[-20,45],[-74,68],[-37,72],[-36,47],[-39,69],[-106,115],[-44,58],[-29,58],[-33,-7],[-12,66],[-62,87],[-13,48],[-3,56],[20,228],[-16,73],[-1,69],[-15,75],[-57,154],[-77,53],[-100,137],[-51,140],[-19,-34],[-42,3],[-48,-39],[-28,36],[-39,90],[-43,11],[-31,60],[-25,20],[-37,7],[-31,31],[-14,81],[-12,28],[-48,55],[-40,62],[-37,39],[-12,54],[59,-6],[70,-23],[33,6],[41,39],[0,-50],[19,-29],[27,4],[-11,54],[5,45],[-30,8],[22,40],[26,107],[12,105],[-27,91],[-45,64],[-98,187],[-91,148],[-48,22],[-40,54],[-70,76],[-30,39],[-21,92],[-16,12],[5,64],[-7,112],[-12,29],[-38,28],[-18,200],[-65,85],[-10,119],[-35,81],[-56,104],[-2,69],[13,30],[1,48],[-63,75],[-18,103],[-42,74],[-18,97]],[[25201,6701],[-22,41],[26,86],[47,0],[-51,-127]],[[19799,7371],[-8,-42],[-21,17],[-9,69],[28,-11],[10,-33]],[[18676,9109],[-25,-18],[-9,61],[34,-43]],[[17961,11182],[-30,23],[-64,108],[-23,54],[-2,82],[22,-7],[26,-38],[11,-75],[48,-19],[12,-128]],[[17405,10679],[-4,-17],[-46,35],[25,59],[-4,63],[22,-8],[13,-80],[-6,-52]],[[18225,11158],[-21,-121],[-21,2],[-44,38],[-5,24],[17,140],[51,38],[23,-121]],[[16557,11126],[-12,-20],[-32,107],[24,9],[18,-44],[2,-52]],[[18364,8796],[-4,-24],[-83,96],[43,2],[44,-74]],[[18265,8874],[-29,54],[-4,40],[-33,31],[20,83],[26,-114],[1,-38],[19,-56]],[[27874,5756],[34,-19],[64,-16],[57,-68],[97,-54],[38,-122],[-25,-17],[-61,-7],[-40,31],[-51,24],[-25,-64],[-26,3],[-9,24],[-24,-29],[-13,-67],[-43,61],[-29,12],[-57,2],[-27,9],[-31,72],[-23,15],[-30,73],[-61,14],[-13,35],[4,33],[30,50],[34,-2],[32,12],[29,28],[117,-28],[52,-5]],[[25092,2950],[1,54],[-22,22],[-24,75],[3,40],[-32,17],[-29,-13],[-29,20]],[[24960,3165],[21,27],[-12,35],[18,150],[-24,39],[-52,-5],[-25,41],[-19,14],[-35,-23],[-47,-47],[-18,66],[-42,19],[-50,56],[-6,19],[-37,46],[-15,39],[-85,45]],[[24532,3686],[21,58],[32,39],[3,32],[-17,101],[14,89],[8,17],[46,36],[81,111],[126,189]],[[24846,4358],[26,-14],[71,91],[71,16],[36,-61],[56,22],[104,-37],[42,17],[76,3],[34,-9],[49,52],[68,35],[-5,25],[47,0],[82,-52],[89,9],[52,37],[91,-55],[24,-41],[11,36],[71,-26],[134,-197],[-54,42],[-39,-14],[5,-47],[46,-3],[35,-54],[11,-35],[58,37],[50,-27],[22,-83],[37,-44]],[[29388,6401],[6,-15],[14,-136],[-11,-70],[28,-46],[-3,-33],[-42,-74],[21,-67],[-3,-63],[-34,-61],[-37,-8],[35,-93],[31,-39],[7,-36],[-8,-34],[-1,-85]],[[29391,5541],[-49,76],[-31,21],[-123,-4],[-36,-22],[-67,-13],[-140,51],[-72,-11],[-28,-20],[-38,-76],[-56,89],[-30,27],[-62,40],[-16,53],[25,89],[44,20],[35,-11],[66,-34],[48,-5],[27,-22],[186,-34],[35,-11],[26,17],[22,42],[66,13],[8,51],[-32,35],[-51,77],[-45,91],[20,31],[-8,56],[18,102],[-97,87],[-72,14],[-22,11],[-12,32],[10,44],[24,24],[55,26],[66,12],[66,-14],[57,-45],[59,-35],[73,-12],[33,-13],[15,11]],[[29105,5920],[-5,-36],[-70,42],[-58,55],[3,29],[29,7],[69,-54],[32,-43]],[[29144,6565],[7,-26],[-57,25],[15,29],[35,-28]],[[24588,5427],[-6,-269],[-8,-391],[-7,-324],[95,-2]],[[24662,4441],[26,-14],[54,-51],[18,44],[-15,52],[101,-114]],[[24532,3686],[-49,-3],[-10,-33],[8,-54],[-34,-30],[-22,-54],[-34,-16],[-48,-73],[-16,-35],[3,-50]],[[24330,3338],[-106,84],[-35,14],[-149,-1],[-64,33],[-73,63],[-49,59],[-115,162]],[[32413,3726],[-66,-14],[-7,32],[22,37],[-36,26],[-22,78],[26,36],[28,-25],[29,-47],[-4,-25],[18,-37],[12,-61]],[[32275,4615],[-54,-1],[-5,36],[11,31],[-8,38],[16,38],[18,-20],[14,-54],[50,-55],[-42,-13]],[[32202,4501],[-22,-23],[-25,51],[-9,122],[12,28],[30,-15],[24,-28],[-6,-23],[3,-92],[-7,-20]],[[32301,4441],[-22,2],[-2,31],[19,26],[5,-59]],[[33672,20292],[-18,40],[-15,78],[28,1],[-8,-60],[13,-59]],[[31810,5557],[-31,0]],[[31779,5557],[16,23],[15,-23]],[[24960,3165],[-2,-59],[-30,-53],[-69,-9],[-66,26],[-26,-15],[-98,51],[-113,99],[-68,16],[-78,26],[-46,63],[-34,28]],[[29388,6401],[12,9],[19,58],[30,24],[32,-1],[57,-23],[78,34],[33,-14],[55,-57],[55,-2],[36,-49],[49,-21],[47,19],[19,-42],[3,-60],[15,-54],[24,-35],[114,14],[26,-28],[-25,-36],[-79,2],[-5,-45],[64,-16],[64,-30],[72,-20],[60,-43],[84,-120],[11,-30],[-5,-38],[-24,-62],[-33,-32],[-26,-68],[-17,-2],[-27,62],[-32,35],[-38,-4],[-56,21],[-68,-13],[-69,15],[-69,-36],[-45,-63],[-82,-17],[-45,61],[-32,4],[-74,-37],[-15,-25],[-7,-79],[-44,-113],[-36,-104],[-12,-5],[-22,46],[-31,24],[-12,93],[-26,43]],[[32287,4112],[-26,-11],[-11,89],[-18,64],[7,55],[38,-25],[11,-29],[8,-79],[-9,-64]],[[26610,8166],[73,-24],[84,6],[25,-35],[48,15],[98,7],[28,-44],[51,-38],[80,0],[32,-34],[78,-68],[35,-85],[76,-87],[26,-10],[77,4],[60,-15],[150,-132],[76,-107],[64,-53],[25,-11],[13,49],[55,-81],[43,-35],[-11,-24],[-51,8],[44,-69],[30,57],[73,-97],[40,-31],[12,-28],[53,-6],[54,-24],[51,-48],[97,-11],[25,-26],[10,-34],[-18,-49],[18,-31],[-40,-13],[5,-50],[22,9],[84,-16],[35,6],[91,-32],[42,-40],[60,-96],[36,-28],[41,-2],[19,-18],[3,-63],[-27,-45],[-105,-11],[-60,-29],[-42,-38],[-31,-14],[-17,30],[-31,-34],[-73,-1],[-59,35],[-108,16],[-172,-26],[-60,-24],[-59,0],[-70,-16],[-69,-3],[44,116],[94,112],[31,54],[-1,44],[-23,35],[-11,42],[-66,24],[-146,13],[-33,24],[-54,85],[-26,24],[-23,40],[-13,125],[-28,112],[-25,40],[-26,13],[-101,-34],[-24,5],[-175,92],[-89,67],[-22,31],[-48,97],[-4,-30],[-128,-6],[-33,23],[-28,87],[-11,-64],[-36,-20],[-24,41],[-104,9],[-43,46],[-29,50],[89,40],[20,35],[-17,50],[-44,24],[-230,5],[-34,-48],[-41,-42],[-37,-65],[-41,-43],[-24,-41],[-30,-18],[-32,11],[-23,-15],[-59,-5],[-18,-40],[-18,-76],[-58,-23],[-58,-55],[-14,71],[-50,-16],[-43,-36],[-11,34],[95,70],[44,3],[13,19],[-16,93],[6,63],[23,49],[44,74],[21,24],[217,155],[164,39],[87,56],[69,19],[72,-14]],[[27762,7546],[-24,7],[-40,36],[-5,26],[40,-2],[38,-15],[-9,-52]],[[27663,7716],[-42,19],[-26,52],[34,4],[34,-75]],[[27496,7853],[38,-11],[39,4],[19,-43],[-29,-5],[-44,10],[-35,23],[12,22]],[[27704,7635],[-37,20],[1,42],[32,-17],[4,-45]],[[26410,7351],[-25,-27],[-55,-39],[-29,-1],[-30,14],[-21,32],[47,31],[-30,102],[25,78],[65,-17],[11,-10],[24,-63],[18,-100]],[[26323,436],[-19,57],[-49,87],[-1,78],[-9,43],[-64,66],[-21,-5],[13,-45],[34,-58],[2,-52],[-70,21],[-17,17],[-35,69],[25,59],[8,38],[-1,80],[-6,39],[-27,59],[-44,64],[-61,53],[-29,43],[-72,32],[-27,22],[-22,40],[5,73],[-20,56],[-86,111],[-48,41],[-10,-45],[21,-46],[54,-43],[21,-58],[-31,-62],[-16,-16],[-16,-44],[-55,117],[-85,47],[-16,29],[-32,89],[-15,81],[6,55],[46,121],[-1,56],[-13,24],[-54,55],[43,44],[3,40]],[[26112,1894],[15,-53],[4,-41],[35,-138],[28,-77],[61,-140],[27,-26],[44,-112],[25,-52],[46,-28],[13,-20]],[[12576,24026],[26,-31],[60,20],[23,-16],[26,-62],[3,-41],[-16,-31],[-97,-65],[9,-18],[96,17],[21,110],[-7,69],[54,4],[75,51],[-10,-80],[-32,-74],[-28,-130],[-2,-107],[-23,-58],[-93,-39],[-47,2],[-63,60],[15,30],[72,-9],[-66,55],[-48,25],[-19,62],[-38,77],[-13,72],[10,115],[112,-8]],[[12825,23548],[28,-47],[5,-93],[-74,-28],[49,-57],[75,-28],[-38,-63],[58,-83],[34,-13],[-19,-37],[51,-14],[7,-37],[-29,-33],[-55,43],[-39,83],[-68,72],[-78,108],[-20,16],[-26,77],[32,34],[-90,37],[-15,37],[49,-4],[85,22],[40,27],[38,-19]],[[14084,22240],[137,-64],[137,-32],[101,-37],[62,-12],[36,-21],[69,-156],[83,-144],[1,-45],[27,-58],[52,-52],[126,-66],[53,-39],[48,-72],[19,-68],[35,-65],[36,-124],[29,45],[29,-110],[-14,-25],[-30,11],[-41,-48],[-89,27],[-214,108],[-49,29],[-70,55],[-4,31],[55,69],[3,30],[-67,-12],[-54,3],[-27,-29],[-23,6],[-32,39],[-46,32],[35,25],[-37,70],[-32,-15],[-19,62],[-22,27],[-40,11],[-17,-31],[-31,34],[-36,-15],[2,88],[53,36],[-67,53],[-25,57],[-61,35],[-21,-32],[-37,-1],[-26,58],[4,55],[-43,-25],[-33,57],[-57,0],[-52,-18],[9,43],[-6,52],[-31,17],[16,51],[73,17],[42,-27],[-10,87],[-62,6],[5,-37],[-37,-19],[-53,0],[-58,57],[-22,69],[29,43],[39,15],[51,1],[56,-20],[143,-92]],[[10271,31976],[165,-24],[118,10],[220,-55],[136,-102],[110,-50],[45,-34],[240,-96],[150,-39],[164,-3],[70,-26],[140,-69],[-13,68],[-63,32],[-66,13],[1,43],[78,3],[-89,53],[61,103],[53,13],[67,-7],[32,74],[79,9],[100,-9],[11,101],[58,2],[46,-57],[54,-26],[-79,-116],[60,10],[131,54],[86,14],[89,111],[140,52],[47,-24],[64,25],[136,84],[84,-5],[58,74],[118,47],[89,-30],[64,3],[59,52],[53,-45],[-30,-55],[-224,-111],[-69,-49],[-70,-27],[-249,-42],[-35,-17],[-38,-59],[-56,-48],[-98,-25],[-36,-27],[-76,-91],[-80,-73],[21,-20],[165,-17],[-17,101],[14,35],[221,108],[36,72],[62,15],[152,-16],[8,-90],[89,117],[85,89],[65,33],[151,55],[85,14],[35,-54],[46,-19],[73,50],[88,77],[30,76],[142,58],[-58,36],[-4,36],[-50,29],[12,54],[37,25],[66,-29],[104,-76],[66,-66],[58,-91],[47,-107],[40,-63],[152,-135],[92,-44],[95,-16],[59,40],[1,27],[-52,75],[38,54],[125,132],[-51,27],[109,56],[31,-21],[-8,-98],[18,-78],[78,-39],[-4,-19],[-94,-117],[15,-24],[226,-1],[109,84],[28,101],[24,37],[264,3],[148,-21],[158,-59],[185,-122],[121,-41],[272,-43],[213,-110],[266,-69],[30,4],[265,-31],[-52,53],[172,7],[176,-63],[103,-53],[57,-44],[124,-134],[-36,-67],[-50,-10],[-135,11],[-24,-38],[-76,-32],[-11,-59],[-72,-48],[81,-54],[130,-12],[65,-22],[148,-23],[187,-3],[93,-11],[122,2],[55,24],[183,13],[106,34],[56,-16],[170,83],[91,12],[114,-125],[152,-11],[24,-38],[12,-73],[37,-29],[28,74],[34,4],[67,-96],[99,-75],[22,-48],[-49,-48],[-101,-3],[97,-104],[94,-90],[14,29],[-10,113],[34,21],[129,-62],[-53,89],[18,25],[-55,49],[-44,95],[-2,65],[-89,113],[4,44],[54,45],[0,64],[113,15],[63,24],[80,11],[46,45],[67,-8],[6,61],[37,32],[64,13],[71,61],[2,43],[-44,13],[-145,-62],[-41,-81],[-94,9],[-55,-36],[-131,14],[-33,-24],[2,-57],[-146,-12],[-117,65],[12,52],[91,119],[151,19],[91,21],[167,62],[184,56],[101,-28],[72,-75],[30,-133],[65,-66],[78,-41],[76,-18],[-7,-34],[48,-44],[81,-16],[149,19],[50,24],[216,-156],[133,-40],[90,7],[91,-27],[193,54],[111,21],[188,-13],[83,-16],[90,-31],[63,1],[139,52],[-79,71],[20,55],[59,-43],[98,-112],[36,-30],[131,-48],[83,56],[-14,66],[-94,54],[-53,9],[-103,-39],[-90,83],[16,35],[-70,102],[46,26],[69,-33],[117,35],[-40,63],[82,-14],[64,7],[40,-22],[80,-101],[137,-7],[-68,-93],[37,-6],[42,52],[100,44],[11,-40],[-45,-226],[-47,-84],[0,-22],[51,-69],[5,-44],[43,-9],[81,23],[-18,-63],[70,8],[32,-28],[-1,-72],[-97,-24],[-69,26],[-94,-15],[131,-127],[-44,93],[20,15],[96,-18],[75,33],[26,144],[-57,178],[-51,65],[52,146],[93,32],[73,-22],[32,15],[135,116],[91,91],[76,35],[-8,176],[-24,43],[-36,6],[-26,-83],[-62,-31],[-73,-9],[-33,31],[10,56],[132,108],[-47,15],[0,83],[27,14],[151,31],[-52,-62],[64,8],[10,94],[-33,20],[-101,-37],[-70,4],[-100,100],[-53,-37],[-129,46],[-115,57],[-71,12],[-43,36],[-103,129],[-14,104],[88,110],[49,13],[-86,57],[-52,67],[39,238],[74,63],[59,-2],[55,-33],[45,4],[43,79],[-117,18],[-2,37],[61,45],[29,46],[86,65],[129,41],[49,-2],[21,-67],[47,-43],[138,1],[13,-66],[140,-90],[85,-102],[18,-99],[-6,-78],[-21,-33],[94,-69],[70,-35],[93,-133],[64,-25],[70,-83],[-14,-16],[-105,-2],[-81,49],[-61,-69],[106,-10],[41,-24],[-85,-75],[-128,-91],[121,-21],[28,-36],[82,-37],[52,8],[94,53],[59,-15],[-32,-40],[202,-19],[-74,-110],[-114,-3],[136,-92],[76,-124],[-20,-22],[8,-67],[-17,-110],[90,-106],[67,67],[35,63],[-2,53],[28,57],[32,143],[87,108],[82,20],[115,-106],[89,-53],[76,-73],[64,-209],[-12,-70],[-38,-13],[-98,24],[9,-173],[32,-95],[55,-72],[138,-138],[22,-73],[42,-15],[125,112],[57,25],[20,33],[17,101],[29,57],[115,127],[61,191],[-1,95],[60,71],[67,-14],[113,16],[-66,36],[6,33],[46,28],[15,57],[-69,47],[-76,34],[-12,62],[4,100],[-24,53],[16,65],[-17,34],[135,-21],[103,23],[112,-7],[159,-69],[186,-8],[103,3],[67,-25],[-74,-94],[92,-25],[44,-82],[55,14],[149,-40],[23,-40],[-78,-60],[-94,-48],[158,-31],[32,-35],[-5,-63],[-68,-52],[-132,-48],[-56,22],[-95,-30],[44,-56],[0,-27],[47,-72],[58,15],[-24,-94],[34,-66],[74,-71],[82,-65],[39,-70],[-8,-52],[-46,-147],[-45,-35],[-68,-7],[-89,-115],[-123,-91],[-84,-18],[-97,-82],[-62,-12],[-110,159],[-36,31],[-75,22],[14,29],[-73,69],[-56,-20],[116,-114],[27,-15],[35,-80],[74,-121],[-22,-40],[-107,41],[-51,-58],[-170,75],[-96,109],[-45,16],[-127,-25],[-157,6],[-31,-59],[44,-49],[79,-26],[83,-43],[11,-34],[-23,-50],[-251,-253],[-36,-45],[-58,-44],[-107,-10],[-80,7],[-194,135],[-20,34],[-95,31],[-92,62],[-91,39],[-27,-33],[-74,5],[-55,23],[-55,-3],[-137,23],[-81,-1],[5,-36],[84,13],[106,-17],[152,-41],[72,-35],[151,-149],[74,-50],[214,-35],[223,-10],[100,-32],[0,-69],[-70,-122],[-167,-218],[-31,-80],[-30,-31],[-151,-89],[-45,-9],[-67,22],[-43,-26],[-52,29],[-59,-5],[-41,35],[-13,42],[-74,-27],[35,-25],[26,-52],[-17,-18],[-61,20],[36,-91],[-64,-83],[-35,-9],[-62,19],[-31,-33],[-66,-5],[-164,71],[-107,16],[-30,28],[-44,-19],[-98,28],[-116,55],[-127,47],[-46,6],[-19,-51],[29,-39],[40,-14],[10,57],[235,-96],[103,-36],[-84,-62],[47,-4],[124,69],[43,-19],[98,-78],[99,-33],[62,-47],[15,-123],[-10,-48],[-39,-37],[-68,-12],[-64,-53],[-28,-7],[-116,15],[-46,15],[-90,-23],[5,-43],[28,-26],[79,-20],[-14,-53],[-55,23],[-95,-20],[-12,-39],[-75,-84],[-11,-47],[-44,-33],[13,-30],[-98,-84],[-58,7],[42,-66],[32,-17],[-29,-38],[-80,-53],[-64,-81],[14,-51],[-47,-14],[-2,-91],[-67,-137],[-72,-169],[-53,-22],[32,-42],[-38,-237],[2,-244],[-3,-108],[-23,-92],[26,-34],[18,-60],[37,-18],[44,-63],[36,-15],[64,22],[76,7],[81,-15],[58,6],[34,-31],[15,-105],[71,-212],[39,-153],[59,-158],[13,-99],[-48,-121],[86,3],[15,-19],[64,29],[250,90],[59,8],[85,-16],[142,-89],[79,-36],[233,-66],[33,-19],[105,-106],[103,-72],[54,-64],[88,-146],[76,-42],[101,-31],[150,-72],[108,-40],[85,-49],[75,-87],[52,-42],[1,-64],[42,32],[98,-14],[220,29],[67,-26],[125,-26],[63,9],[33,-36],[80,-6],[51,-42],[23,-35],[25,-95],[-55,-234],[7,-90],[42,-70],[28,-116],[-14,-106],[12,-126],[-3,-52],[-28,-53],[-9,-48],[9,-54],[66,-76],[44,-82],[72,-113],[8,-33],[-25,-38],[54,-45],[50,-59],[87,-60],[86,-110],[45,-119],[14,-70],[-63,-64],[84,34],[73,-24],[74,-57],[27,7],[14,82],[53,115],[44,29],[20,-39],[48,-25],[24,-39],[29,-79],[20,81],[26,34],[-12,35],[-56,87],[42,117],[34,36],[58,112],[-18,123],[-46,69],[-11,68],[12,50],[-49,95],[-26,188],[-33,157],[30,62],[-20,49],[22,36],[-19,25],[-49,24],[29,30],[-44,30],[-37,62],[-65,185],[209,94],[167,96],[114,86],[157,163],[155,226],[31,57],[29,128],[7,179],[-2,94],[-13,149],[-22,102],[-65,193],[-74,134],[-146,139],[-91,55],[-93,93],[-31,12],[-14,41],[23,112],[123,176],[40,30],[22,39],[-27,49],[35,54],[2,48],[66,5],[43,57],[-12,66],[17,60],[-65,10],[26,43],[-63,144],[46,70],[-55,7],[-21,32],[36,77],[-58,-2],[-80,34],[91,119],[29,75],[-3,106],[61,61],[-50,35],[-69,70],[-36,83],[-17,94],[1,89],[18,38],[129,90],[64,21],[136,-24],[72,-31],[221,-76],[39,-34],[92,32],[180,-66],[13,14],[121,56],[122,76],[37,-8],[120,-89],[87,-79],[58,-6],[11,-44],[-27,-97],[62,43],[77,-46],[40,-53],[16,-60],[43,39],[63,-36],[-5,-34],[-51,-44],[27,-66],[89,-91],[68,-6],[125,-47],[50,-7],[72,14],[26,-15],[25,-82],[52,-10],[55,28],[23,69],[42,-20],[21,-65],[-10,-36],[-57,-63],[-23,-63],[-9,-92],[34,-86],[-11,-74],[-34,-23],[-46,-6],[-151,-2],[-31,-16],[225,-20],[49,-63],[14,-75],[-31,-96],[2,-75],[78,-2],[14,-31],[-28,-50],[10,-47],[-32,-112],[-40,-19],[-21,55],[-32,-42],[-81,-35],[78,-47],[63,16],[74,63],[58,23],[64,-4],[67,8],[78,-62],[25,-96],[20,-20],[5,-80],[34,37],[36,-6],[35,-96],[32,57],[14,-72],[33,-16],[52,45],[100,68],[82,30],[31,29],[26,76],[42,48],[54,-33],[34,-62],[32,3],[-32,107],[52,14],[-18,35],[62,60],[25,40],[0,52],[26,71],[6,50],[39,15],[-35,71],[12,82],[75,27],[-3,40],[66,176],[32,41],[31,3],[57,-35],[22,-50],[-30,-39],[67,-16],[28,-46],[9,-92],[26,15],[35,-91],[50,-68],[-54,-68],[112,-24],[34,-71],[-42,-76],[67,18],[39,-115],[37,11],[49,-106],[-64,-64],[-32,-14],[40,-54],[66,34],[67,-3],[-22,-79],[-38,-61],[-70,-38],[23,-17],[108,33],[29,-2],[49,-93],[52,-4],[61,-53],[-27,-47],[18,-52],[-10,-29],[-103,-65],[10,-22],[60,6],[65,-37],[-35,-68],[32,-26],[64,-7],[82,-89],[-16,-80],[5,-89],[-44,-13],[-146,23],[-84,35],[-8,-19],[77,-33],[34,-39],[38,-3],[-24,-82],[118,-17],[4,-33],[-57,-19],[-2,-30],[43,-12],[53,4],[18,-86],[-41,-27],[72,-20],[27,-46],[26,-2],[45,49],[67,-73],[69,-16],[-19,-69],[60,-86],[-46,-59],[-25,-85],[143,30],[48,26],[19,-58],[47,1],[9,-35],[-70,-96],[86,58],[47,74],[72,-26],[31,-101],[29,-58],[106,-33],[48,20],[8,27],[64,5],[9,-53],[57,-51],[66,-17],[22,-36],[-44,-66],[-37,-28],[-53,-1],[-72,-17],[-19,-33],[-60,-30],[-54,-91],[-58,-3],[-54,-41],[-45,-18],[-82,-49],[-90,0],[-53,-37],[-36,-85],[13,-56],[-53,-48],[5,-39],[34,-5],[88,102],[58,16],[82,55],[112,118],[73,53],[90,38],[5,77],[32,13],[160,-19],[55,-14],[74,-162],[-3,-46],[-91,-64],[28,-37],[15,-47],[88,105],[48,34],[103,13],[31,-72],[67,-19],[40,-60],[30,-65],[17,-67],[-3,-40],[-23,-69],[17,-62],[-12,-73],[7,-58],[-88,-25],[90,-34],[37,-41],[-35,-60],[39,-35],[-1,-54],[-89,-106],[-73,-68],[-74,-59],[-118,-115],[-33,-8],[-56,19],[-130,-27],[-70,-53],[-68,-14],[-66,0],[-29,-30],[-7,-33],[-115,-150],[-30,-51],[-59,-54],[-65,-93],[-56,-38],[-19,-52],[-54,-32],[-146,-17],[-55,14],[-41,-22],[-92,-2],[-120,-50],[-54,66],[-68,3],[-54,20],[-98,12],[-116,-4],[-28,-26],[-96,8],[-41,24],[-78,-1],[-43,-17],[-93,20],[-99,-17],[-87,11],[-24,11],[-136,-31],[-53,18],[-47,-48],[-33,10],[-70,-4],[-34,-29],[-33,-46],[-56,-37],[-81,-200],[-7,-78],[-31,-52],[-49,-8],[-140,-39],[-62,-31],[17,-24],[-54,-26],[-59,-47],[-11,-35],[-72,-56],[-83,-131],[-40,-96],[-48,-69],[-34,-27],[-49,4],[-42,33],[-30,4],[-75,45],[-176,45],[50,-45],[93,-6],[98,-57],[78,-37],[25,-38],[-36,-136],[-25,-47],[-84,-121],[-41,-41],[-71,-145],[-114,-111],[-42,-66],[-98,-50],[-37,-12],[-34,7],[-41,-41],[-48,-25],[-15,-38],[-117,-101],[-44,-13],[-39,-27],[-11,-46],[-34,-27],[-38,-86],[-54,-82],[-65,-14],[-51,-75],[-39,-32],[18,-81],[-34,-10],[-72,-61],[-97,-103]],[[28278,19090],[-71,-64],[-120,-63],[-27,-2],[-55,-29],[-56,-56],[-19,67],[-16,-64],[26,-47],[-9,-54],[-26,-8],[-39,29],[-9,25],[-49,-5],[-23,39],[-115,-30],[-51,-7],[-95,-31],[-39,-2],[-108,-34],[-67,-75],[-39,-21],[-34,-36],[-57,-132],[13,-34],[55,-23],[53,-6],[76,44]],[[27376,18381],[7,-44],[29,-58],[-34,-23],[-91,9],[-57,-14],[-42,3],[-131,-43],[-64,-83],[96,-21],[3,-13],[-144,25],[-59,24],[-58,10],[-59,-4],[-67,-32],[-65,-74],[-43,-89],[-44,-8],[-73,-51],[-42,-51],[-14,-65],[-25,48],[-24,5],[-40,-19],[-28,4],[-47,33],[8,105],[30,26],[28,5],[77,-11],[44,27],[3,54],[-54,27],[-9,17]],[[26450,18335],[66,35],[40,51],[69,73],[20,75],[-2,260],[-3,40],[37,129],[21,45],[33,31],[34,94],[0,75],[-38,101],[-48,92],[-26,18],[24,31],[83,-7],[-5,-47],[15,-66],[51,-59],[15,-57],[36,-21],[-10,-65],[38,34],[39,11],[15,-43],[37,-31],[64,-35],[39,-10],[22,42],[3,41],[-31,80],[31,25],[47,-45],[15,26],[-16,60],[-71,89],[18,13],[-25,63],[20,43],[-14,26],[-50,-17],[-31,50],[-1,59],[-34,9],[-40,102],[-45,83],[-77,-5],[-39,16],[-99,12],[-5,63],[-43,-8],[-118,15],[-42,31],[-87,-9],[-15,16],[-85,-13],[-119,28],[-56,27],[-50,6],[-11,16],[-58,15],[1,74],[-21,29],[-74,-21],[-19,-20],[-19,48],[36,75],[-34,4],[2,38],[52,38],[-19,26],[-52,-4],[-35,30],[-6,29],[49,145],[-33,54],[-61,58],[-22,51],[44,159],[-30,6],[-180,-10],[-45,15],[-37,30],[-53,93],[-31,109],[-11,64],[-38,94],[-52,11],[-17,21],[-49,-11],[-52,8],[-59,-7],[-51,31],[-154,63],[-16,21],[-66,-13],[11,-55],[32,-46],[-9,-48],[-26,-44],[-92,-35],[64,92],[2,34],[-19,36],[-29,-6],[-18,-99],[-38,-84],[-8,-47],[-43,-24],[32,126],[-86,-44],[-32,-24],[-6,-72],[-34,-97],[-26,-11],[-38,-50]],[[15301,21396],[-37,42],[-35,-42],[-40,79],[-10,58],[45,16],[-29,29],[-29,-3],[15,127],[-28,-7],[-3,-60],[-54,-32],[-90,44],[-47,61],[-7,30],[36,29],[29,-67],[15,37],[-31,40],[-27,79],[-46,-51],[-39,-10],[-56,18],[-82,109],[-1,27],[-75,177],[53,112],[-3,96],[-30,-119],[-21,-40],[-41,-20],[-74,11],[-61,-17],[-21,23],[-53,-11],[-121,22],[9,70],[-41,30],[34,42],[-65,55],[-91,-11],[-27,-13],[-58,25],[-48,37],[-73,83],[-2,60],[22,82],[59,92],[39,24],[86,8],[-18,24],[-84,-5],[-29,-15],[-56,-102],[-19,6],[-29,58],[-14,88],[11,53],[28,58],[-36,-1],[17,103],[99,84],[54,20],[15,73],[53,62],[-10,12],[-134,-165],[-76,-36],[-30,-71],[-26,-102],[-29,-83],[-25,107],[-46,82],[89,82],[-7,56],[-59,4],[46,216],[1,25],[-72,-41],[-24,26],[-20,135],[-35,53],[-60,43],[-59,20],[-25,85],[16,55],[50,38],[23,-25],[5,-60],[80,-65],[20,25],[-57,43],[-20,55],[10,76],[50,25],[-48,31],[-51,-46],[-19,6],[-69,-102],[-21,-128],[-49,-23],[-28,-49],[-71,82],[-25,71],[-45,13],[-52,40],[-20,36],[28,73],[41,56],[12,81],[69,16],[-29,31],[-52,-23],[-57,46],[-38,77],[16,102],[42,56],[91,223],[30,59],[61,93],[-39,-12],[-30,48],[-8,-109],[-40,-106],[-30,-2],[27,80],[3,48],[-30,140],[2,61],[17,70]],[[17584,33483],[50,-1],[106,31],[126,22],[14,87],[67,60],[125,3],[195,-50],[215,-90],[-23,-72],[-71,-70],[-78,-41],[60,-29],[101,80],[47,-46],[99,86],[95,51],[64,14],[2,50],[-98,53],[-37,42],[7,33],[181,-12],[110,-56],[39,-34],[95,-42],[117,-117],[64,-184],[53,-122],[18,-77],[79,-57],[160,92],[11,39],[-65,67],[-43,63],[-35,150],[-30,40],[-85,300],[10,40],[62,12],[-30,49],[102,13],[167,-70],[45,43],[129,-41],[185,-95],[110,-114],[50,-191],[98,-223],[31,-57],[105,-149],[10,-73],[-60,-169],[111,-91],[59,-84],[102,-68],[132,-61],[98,6],[210,-122],[71,1],[29,-19],[18,-59],[114,7],[49,-61],[18,-82],[-8,-78],[-77,-19],[-51,36],[-68,-26],[-59,20],[-89,55],[-114,-65],[17,-50],[-24,-35],[-83,6],[-121,62],[-29,-21],[89,-75],[28,-67],[43,-13],[121,69],[81,6],[49,-41],[-24,-65],[66,-33],[10,-43],[-19,-56],[-175,-69],[-112,-33],[-158,-8],[-166,29],[-78,32],[-61,-28],[-147,24],[-18,18],[41,65],[-141,27],[-168,14],[-56,32],[-5,80],[-82,30],[-55,-47],[-49,-85],[-112,-92],[-117,-24],[-139,-10],[-51,-19],[-109,-71],[-145,-42],[-135,-24],[-245,-27],[-77,5],[-51,-23],[-275,-14],[-154,-20],[-73,9],[-58,54],[-77,122],[7,63],[-28,120],[-174,38],[-231,-2],[-127,9],[-133,28],[-114,45],[-55,77],[-108,117],[-25,128],[16,24],[161,38],[283,42],[259,28],[117,-2],[139,-20],[125,-3],[158,-19],[124,25],[112,-1],[16,41],[-107,48],[-285,87],[-169,38],[-158,-8],[-141,-27],[-318,-18],[-92,19],[-184,-10],[-164,13],[-187,132],[-31,41],[29,34],[126,64],[138,27],[300,75],[68,57],[68,8],[41,31],[-330,-52],[-153,-5],[-107,-22],[-128,23],[11,46],[74,18],[38,50],[-174,-10],[-137,24],[-33,40],[12,113],[98,93],[104,49],[-74,85],[29,54],[227,152],[82,43],[189,72],[282,81],[252,82],[93,-22],[48,-51],[21,-82],[-16,-107],[-114,-138]],[[16144,34262],[73,45],[74,-1],[49,-79],[78,75],[55,27],[160,11],[125,-17],[87,-31],[131,-74],[212,-143],[123,-66],[33,-60],[-15,-32],[-151,-60],[-135,-35],[-161,-75],[-110,-36],[-144,-69],[-270,-112],[-47,-38],[-76,-128],[-99,-59],[-114,-16],[-4,-44],[-48,-123],[-21,-131],[-49,-64],[-86,-31],[-149,-29],[-77,28],[-112,-91],[-109,-37],[-107,-54],[-75,18],[-51,49],[-56,104],[-114,131],[-207,80],[-148,71],[-130,-6],[1,90],[105,147],[109,84],[-4,115],[20,26],[96,32],[-2,41],[-68,36],[67,86],[46,89],[86,56],[23,60],[65,63],[-81,46],[-47,49],[-121,202],[339,45],[234,15],[242,39],[119,-5],[210,-91],[169,-51],[105,-11],[-48,-61]],[[21594,34176],[233,-109],[278,44],[88,31],[69,5],[99,-23],[62,-91],[-43,-38],[-66,-14],[34,-91],[-145,-52],[-177,-146],[20,-34],[118,52],[59,1],[109,-51],[-13,-37],[81,-53],[9,-64],[105,39],[36,-7],[28,-109],[-22,-48],[-77,-41],[57,-56],[-48,-81],[39,-40],[3,-68],[-92,-22],[-76,-61],[-100,-22],[-165,17],[-14,-53],[9,-61],[-59,-47],[-69,-24],[-89,35],[-50,-1],[-66,97],[-91,102],[-163,127],[-75,76],[-107,29],[-62,55],[-80,-20],[-62,19],[-69,88],[-119,55],[-86,97],[58,114],[83,36],[78,-10],[46,-58],[58,-31],[54,-70],[72,-17],[167,31],[24,104],[83,-7],[-51,81],[-69,9],[25,59],[120,-7],[-70,51],[-174,-33],[-155,80],[-20,28],[55,44],[96,14],[83,-54],[19,42],[-95,59],[10,58],[119,20],[133,52]],[[25759,29728],[22,-3],[63,98],[31,6],[65,-47],[50,-76],[100,-25],[87,-33],[51,-70],[64,-33],[89,-63],[110,-30],[61,-40],[72,-112],[31,-109],[-11,-98],[105,28],[91,-22],[49,29],[61,-60],[8,-32],[88,-55],[-12,-21],[-205,-154],[-253,104],[-115,22],[-24,113],[-29,17],[-128,32],[-9,70],[-66,-8],[-52,-23],[-34,-44],[-12,-72],[-19,-30],[-81,-79],[-79,-38],[-55,-100],[-35,-42],[-90,-57],[-119,-40],[-29,10],[-19,68],[-41,198],[-25,25],[-137,-26],[-76,3],[-75,-44],[-91,10],[-3,39],[43,81],[40,48],[128,66],[47,43],[-28,123],[-5,97],[51,228],[20,208],[42,147],[41,65],[72,45],[31,-38],[55,-25],[38,-89],[-37,-57],[51,-38],[37,-90]],[[23481,34287],[187,-51],[74,-46],[164,28],[150,-9],[127,-30],[75,-42],[-58,-94],[-102,-81],[-136,-170],[-53,-87],[-79,-90],[-59,-49],[-75,-18],[-262,42],[-241,-23],[122,-45],[55,-56],[4,-71],[-87,-86],[-46,-115],[-29,-14],[-98,22],[-141,-16],[-51,8],[7,78],[-23,165],[-82,143],[-18,82],[7,119],[-15,196],[3,71],[68,31],[157,-48],[-67,73],[-25,52],[30,60],[154,46],[193,3],[103,31],[67,-9]],[[22302,31972],[56,16],[100,-84],[158,-84],[97,-165],[55,-64],[46,-32],[58,29],[29,-34],[-54,-40],[-61,-6],[-58,-56],[-35,-10],[-105,-70],[-54,-5],[-113,40],[-128,2],[-114,66],[-97,35],[-24,53],[-60,-23],[-80,20],[-21,49],[-96,-36],[-52,28],[-34,60],[30,49],[150,19],[52,26],[74,60],[-28,64],[3,43],[71,43],[6,65],[52,40],[80,14],[104,-63],[-7,-49]],[[32336,19834],[47,21],[18,-19],[-46,-34],[18,-22],[-33,-46],[9,-23],[35,23],[35,-5],[29,18],[63,88],[-64,-7],[-8,13],[44,53],[-4,28],[23,45],[57,55],[20,-59],[37,5],[62,-24],[5,-23],[-24,-48],[26,-40],[-48,-31],[-52,-71],[-50,-45],[-80,-33],[-55,10],[-113,-19],[-46,73],[-13,116],[4,60],[20,57],[47,68],[102,253],[70,92],[57,14],[-4,-41],[26,-80],[-15,-79],[-27,-102],[-3,-74],[-23,-50],[-46,-41],[-47,-24],[-53,-52]],[[28579,19353],[39,19],[83,75],[60,27],[79,79],[57,15],[11,18],[15,88],[26,65],[32,54],[26,75],[47,48],[71,40],[65,87],[71,46],[15,36],[79,61],[63,12],[113,52],[30,32],[44,17],[131,93],[36,43],[47,88],[41,45],[15,48],[59,78],[61,103],[31,73],[45,41],[88,117],[48,47],[72,46],[34,44],[53,43],[97,53],[90,65],[123,55],[143,83],[116,44],[82,7],[99,21],[35,-2],[155,-36],[74,-45],[84,-94],[15,-59],[-45,17],[-12,-18],[47,-58],[-3,-72],[-26,-65],[-78,-32],[-37,-68],[-54,-35],[-21,-27],[-61,-44],[-60,5],[-77,41],[-48,40],[-43,-44],[-45,7],[-31,-29],[41,-37],[100,-39],[26,-27],[24,-84],[51,-4],[70,63],[98,-7],[47,-37],[-41,-79],[-16,-103],[-48,-69],[-65,-68],[35,-27],[42,20],[58,-14],[-21,-88],[25,-97],[25,-14],[22,-162],[25,-23],[4,-36],[91,-10],[18,-14],[64,-14],[23,-30],[-62,-44],[50,-32],[48,-52],[53,9],[42,-33],[18,-30],[51,15],[54,-3],[58,-18],[-13,-47],[45,6],[28,-20],[17,22],[57,34],[72,70],[22,-79],[28,-30],[33,-7],[45,23],[39,-58],[19,-68],[-49,-38],[98,-9],[20,-29],[-19,-30],[-50,1],[-29,-26],[-50,-16],[-127,-81],[-65,-30],[-69,-48],[-143,-66],[-34,-1],[-42,-36],[-66,7],[0,-57],[-18,-33],[-41,4],[-36,30],[-30,51],[-28,-81],[-19,51],[-32,-19],[-15,-54],[17,-57],[-28,-15],[-25,-61],[-30,-22],[-31,-62],[-37,-47],[-11,-31],[-62,-72],[-41,-1],[-26,-30],[-4,-60],[-38,-16],[-23,18],[-47,3],[-42,120],[-25,11],[-16,-38],[-25,37],[-18,136],[0,33],[26,115],[64,103],[-14,23],[65,19],[26,39],[-19,17],[208,184],[42,31],[84,39],[24,-9],[2,-52],[33,-12],[39,53],[95,48],[80,5],[-43,37],[-81,-8],[-50,17],[-68,-11],[-73,11],[-41,-44],[-46,26],[24,51],[74,77],[68,108],[-47,-15],[-15,30],[-66,-116],[-33,-7],[-44,-42],[-63,-36],[-72,-70],[-94,-58],[-61,17],[-92,-82],[-20,24],[-54,-31],[-46,-8],[-13,40],[-56,12]],[[33895,22699],[-34,-76],[-124,-31],[9,-66],[51,-4],[4,-81],[-20,-65],[-57,-64],[-16,-67],[-68,-120],[-20,-19],[-8,-56],[-58,-108],[-16,-90],[-27,-61],[26,-58],[71,112],[34,33],[36,80],[34,-14],[-11,-54],[65,31],[45,-30],[65,19],[0,-31],[-52,-55],[-89,-64],[47,-17],[-6,-46],[-42,-64],[63,30],[38,-41],[83,28],[6,-56],[43,15],[-17,-73],[4,-34],[61,54],[32,5],[53,37],[38,53],[27,-11],[5,-66],[41,43],[99,9],[26,-8],[67,-54],[14,-29],[-1,-63],[-51,-55],[-24,-56],[-70,-72],[68,14],[-31,-54],[74,-14],[-27,-78],[67,-23],[39,41],[53,8],[45,42],[10,-45],[-32,-87],[-55,-24],[-54,-63],[-11,-74],[-23,-32],[-12,-53],[-46,-63],[27,-76],[45,6],[27,42],[61,130],[100,76],[14,-27],[-23,-33],[-43,-111],[-18,-82],[1,-72],[14,-29],[49,50],[19,34],[25,77],[20,-12],[16,-101],[-8,-63],[-63,-165],[6,-66],[-20,-80],[-30,-70],[-27,-18],[-34,36],[-25,-3],[-43,-40],[-22,24],[10,142],[-15,92],[-62,-101],[-60,-61],[-27,31],[23,105],[61,154],[-33,189],[-56,53],[-11,-45],[-48,-110],[-8,-48],[-30,-19],[-86,-21],[-65,-108],[-14,-66],[-48,-71],[-60,4],[-71,-24],[-44,31],[-2,23],[32,51],[96,45],[36,50],[47,98],[60,34],[29,28],[-57,45],[-80,2],[-19,-82],[-39,-19],[-79,33],[5,132],[-49,-12],[-36,6],[-19,-53],[-66,-28],[-73,-13],[-14,-14],[-49,5],[-144,29],[-52,-3],[-73,26],[-112,1],[-168,-50],[-49,-5],[-39,33],[-29,118],[25,67],[86,84],[69,85],[61,96],[-64,14],[-123,-10],[0,19],[72,47],[18,-23],[35,-4],[61,208],[25,42],[60,-12],[24,8],[14,52],[-39,41],[-6,65],[29,58],[25,22],[44,-35],[16,18],[-27,47],[-8,51],[59,166],[68,225],[37,72],[9,48],[46,46],[35,78],[8,45],[47,60],[16,67],[18,29],[46,34],[86,46],[50,41],[34,-2],[10,-29],[48,-19],[7,55],[40,9],[17,-30]],[[25298,32673],[128,20],[133,53],[152,-1],[8,-37],[76,7],[32,177],[4,76],[-50,37],[-113,8],[-40,27],[-56,71],[-87,62],[101,59],[39,58],[107,1],[104,-14],[-55,61],[60,35],[-62,15],[-85,-12],[-51,20],[-82,108],[-2,72],[52,78],[54,14],[277,-80],[-5,20],[-226,85],[-80,22],[-20,31],[121,118],[111,27],[55,35],[91,1],[99,60],[295,79],[197,0],[95,-18],[55,-31],[46,-79],[24,-85],[151,-99],[3,-99],[88,-93],[-111,-96],[-23,-52],[43,-16],[-93,-111],[69,-54],[-63,-16],[-8,-80],[40,6],[115,113],[75,33],[20,93],[73,22],[34,-54],[63,2],[29,28],[89,-61],[50,-3],[33,44],[75,-37],[86,0],[108,-51],[52,20],[-207,72],[-46,39],[-7,37],[36,47],[96,45],[120,25],[186,-8],[54,-13],[140,-63],[134,-1],[143,-70],[37,-106],[-94,-96],[135,31],[106,-7],[70,-25],[15,-51],[-28,-49],[-160,-60],[53,-25],[-46,-88],[30,-85],[44,126],[54,45],[73,5],[41,34],[77,-15],[0,-68],[64,-55],[160,102],[89,-37],[60,2],[135,-25],[114,-50],[57,-52],[10,-63],[-76,-69],[-100,0],[-46,-20],[-76,-65],[-93,-55],[89,1],[84,92],[73,18],[103,-37],[51,1],[39,37],[70,27],[43,-49],[-5,-70],[-102,-100],[-155,-56],[-84,-69],[57,-7],[32,45],[80,-2],[-10,-149],[92,190],[116,80],[168,55],[71,-30],[145,-11],[77,-40],[109,-39],[36,-66],[-32,-34],[-86,-47],[-140,-24],[-110,-45],[-100,-60],[73,-13],[86,58],[156,24],[54,-31],[-79,-85],[73,-20],[71,64],[51,11],[-4,38],[32,62],[40,18],[56,-18],[136,-127],[52,-120],[-24,-40],[-151,28],[-85,-12],[-71,-56],[-82,0],[-126,-36],[14,-28],[80,23],[75,7],[126,-52],[192,-2],[85,-26],[78,-43],[26,-75],[-35,-20],[-146,17],[-42,-9],[-125,41],[-101,-42],[152,-65],[1,-54],[-37,-77],[-122,31],[-110,-17],[-41,-26],[149,-14],[50,-26],[31,-67],[45,-14],[61,14],[55,-27],[95,-19],[71,6],[26,-25],[-65,-25],[-10,-25],[84,-101],[55,62],[88,16],[-18,-123],[14,-42],[66,56],[3,-57],[31,-18],[45,33],[58,-26],[52,2],[90,53],[39,-27],[-53,-50],[2,-54],[52,-2],[53,29],[31,-14],[102,-97],[38,16],[47,-45],[-63,-36],[19,-76],[-81,3],[-15,-50],[144,4],[3,26],[64,32],[76,-21],[77,-52],[-34,-31],[-27,-77],[-122,-103],[204,65],[161,-23],[71,73],[43,-14],[184,-190],[-63,-47],[-104,62],[-51,-27],[121,-63],[22,-66],[-80,-31],[-82,13],[-60,37],[-8,-49],[-35,-45],[80,-41],[69,-57],[-139,-32],[31,-72],[-63,-56],[-13,-43],[-86,-9],[-125,44],[-63,-17],[87,-30],[-8,-37],[0,-161],[-10,-57],[-57,-97],[-80,92],[-45,7],[-25,-29],[-88,96],[-24,-70],[-30,27],[-123,151],[-18,53],[-63,101],[7,29],[53,61],[92,39],[50,90],[-92,-59],[-83,-35],[-66,-11],[-78,5],[-8,46],[-57,24],[-40,37],[-75,38],[-65,92],[-20,49],[-147,-18],[35,-45],[1,-57],[-34,-9],[-69,50],[-89,34],[49,-102],[144,-120],[-24,-40],[-77,-18],[-77,23],[-88,84],[-86,61],[-2,-39],[69,-48],[-11,-72],[39,-57],[69,-27],[-25,-77],[14,-31],[93,41],[47,-28],[26,-54],[51,-22],[-57,-43],[72,-58],[26,-81],[31,3],[35,-118],[27,72],[113,-75],[-21,-39],[41,-22],[59,77],[47,-17],[71,-75],[22,10],[43,-48],[-65,-54],[92,-9],[29,-37],[-38,-68],[-82,7],[137,-158],[92,10],[-2,-35],[30,-32],[46,-98],[-42,-14],[18,-111],[-5,-102],[-50,3],[-52,143],[-29,49],[-27,-20],[23,-139],[-18,-36],[110,-175],[-69,-17],[-67,22],[33,-118],[-37,-27],[-24,34],[-131,111],[-34,-5],[-32,40],[-120,14],[-52,98],[-10,-71],[-93,90],[-12,31],[-45,-6],[-87,91],[-59,78],[-31,0],[27,-108],[-31,13],[-112,96],[-69,46],[-101,13],[-15,-25],[67,-97],[82,-83],[64,-92],[62,-30],[66,-10],[-16,-43],[74,-32],[70,-53],[65,-73],[92,-51],[101,-153],[83,-48],[-43,-55],[21,-96],[-55,-24],[-81,30],[-70,44],[-104,31],[-72,40],[-259,43],[-194,87],[-78,86],[-50,88],[-70,24],[-75,-17],[-46,3],[-92,54],[-64,21],[-81,56],[-44,12],[-100,74],[-98,124],[76,0],[66,42],[-15,38],[-59,62],[-31,9],[-106,-8],[-19,25],[18,69],[-62,-38],[-38,42],[-14,46],[-76,65],[-72,85],[-74,72],[45,70],[-96,22],[-78,-13],[-11,-57],[-41,-12],[-9,85],[-79,13],[-27,18],[-33,79],[-79,-30],[77,-106],[-15,-31],[-103,-21],[-72,17],[-44,26],[-63,-8],[-14,-68],[-97,-8],[-80,-38],[-43,0],[-44,-31],[-155,20],[-95,42],[-45,2],[-70,68],[-36,61],[-4,68],[26,96],[60,69],[119,45],[23,30],[-27,55],[37,64],[151,-20],[83,-22],[180,-73],[50,-44],[40,-66],[-24,-55],[38,-33],[26,86],[-41,65],[-74,62],[18,30],[108,-21],[42,5],[51,50],[87,-8],[76,16],[42,44],[100,18],[56,-18],[31,30],[-51,115],[-79,57],[-67,69],[-44,65],[16,35],[122,77],[97,76],[84,86],[68,28],[24,79],[43,75],[117,53],[40,62],[-36,45],[-62,162],[-90,147],[-44,83],[-74,82],[6,51],[-93,-37],[-65,69],[27,76],[-19,63],[-74,0],[42,-63],[-84,-19],[-42,20],[-62,65],[-9,43],[-64,17],[44,35],[-45,46],[-82,-37],[-86,22],[-27,-33],[-216,-100],[-60,14],[18,148],[48,22],[118,-11],[79,65],[6,28],[-39,54],[-71,35],[-100,28],[-31,39],[27,87],[-46,-16],[-42,-42],[-107,40],[72,48],[-26,41],[-75,15],[-72,-10],[-39,33],[-24,139],[-33,39],[-122,-9],[-93,63],[-123,128],[-92,-85],[14,-33],[91,-24],[58,-75],[8,-67],[-30,-36],[-115,-43],[-86,0],[-151,52],[-258,49],[-128,10],[9,-32],[91,-33],[75,-68],[1,-61],[-177,108],[-159,-52],[-96,13],[-167,70],[-121,-20],[-91,-1],[-183,22],[-107,37],[-144,17],[-96,-38],[-150,56],[-84,131],[-137,5],[16,-44],[-121,0],[-110,-38],[-138,98],[-124,41],[-99,121],[-45,122],[1,36],[235,-19],[57,-24],[128,-20],[184,22],[11,12],[-119,50],[-83,51],[-198,17],[-136,24],[-169,69],[-42,36],[-50,233],[26,52],[74,56],[-44,25],[-11,105],[45,80],[95,121],[26,135],[79,94],[59,31],[16,46],[148,99],[124,65],[263,57],[100,11],[125,-3],[233,-21],[37,-58],[-63,-47],[-132,-73],[-112,-103],[-109,-153],[-41,-46],[-11,-56],[21,-52],[75,-102],[-5,-173],[34,-115],[50,-66],[138,-110],[126,-80],[16,-24],[-106,-64],[-149,-33],[-178,-78]],[[32144,21447],[-116,-7],[-92,31],[-68,16],[-67,28],[-145,89],[-30,69],[-58,63],[-153,87],[-12,31],[31,20],[66,8],[103,-34],[129,-30],[120,-57],[62,-42],[138,-111],[25,-9],[62,-54],[34,-74],[-29,-24]],[[31589,20103],[20,6],[16,42],[86,-28],[23,-25],[43,-19],[124,19],[143,14],[38,-33],[-81,-74],[-57,-38],[-8,-19],[6,-97],[-58,-5],[-38,17],[-4,63],[-62,66],[-42,-27],[-2,-24],[-100,40],[-61,91],[-69,9],[1,70],[-36,35],[-42,5],[10,66],[36,67],[63,82],[4,-76],[-30,-70],[51,-70],[11,-34],[-11,-33],[26,-20]],[[26565,28547],[15,-36],[-22,-79],[-24,-40],[-76,-68],[-49,-60],[-124,-99],[-100,15],[-89,-40],[-17,73],[-39,51],[1,38],[44,47],[100,172],[24,9],[90,-25],[71,37],[69,-5],[62,21],[64,-11]],[[27243,28269],[72,-84],[-10,-114],[-108,-212],[-28,-9],[-77,78],[-31,16],[-19,41],[0,101],[11,50],[59,97],[42,32],[89,4]],[[28312,31296],[144,-45],[25,-48],[-17,-90],[11,-109],[-11,-110],[-21,-40],[-55,-47],[-105,-43],[-252,-24],[-86,15],[-42,73],[-41,152],[21,74],[79,124],[56,71],[104,45],[53,-5],[84,14],[53,-7]],[[27246,34027],[47,-7],[298,13],[62,-9],[236,-76],[56,-74],[68,-23],[52,-78],[78,-40],[-6,-42],[61,-60],[-86,-31],[-169,12],[-227,27],[-133,-8],[-226,-56],[-101,-8],[-121,55],[-31,58],[-37,127],[-120,24],[-57,56],[-7,174],[121,23],[81,-30],[161,-27]],[[26916,23319],[-77,6],[-95,47],[-134,55],[-55,47],[53,69],[141,20],[55,-10],[103,-144],[15,-47],[-6,-43]],[[27442,25059],[-8,45],[15,77],[37,28],[3,-62],[-15,-53],[-32,-35]],[[27412,25121],[-41,-102],[-25,-93],[-28,19],[37,110],[-22,20],[-51,-124],[-25,-30],[-30,0],[-44,-35],[74,146],[-19,15],[-60,-93],[-34,-35],[-26,31],[58,93],[53,67],[38,-17],[19,33],[24,96],[9,70],[41,-41],[13,-34],[31,-16],[16,-52],[-8,-28]],[[27124,25091],[-31,3],[23,53],[90,45],[-2,-40],[-80,-61]],[[20134,33651],[-40,9],[-176,132],[-30,56],[-96,56],[-109,34],[25,61],[60,49],[304,36],[110,-11],[129,-66],[26,-38],[-8,-96],[-58,-95],[-78,-88],[-59,-39]],[[13054,24206],[-8,32],[13,46],[36,-7],[-8,-38],[-33,-33]],[[13025,22916],[-20,27],[-4,45],[30,-3],[-6,-69]],[[13244,23939],[-40,-59],[-43,38],[-43,17],[55,63],[27,6],[44,-65]],[[13351,23534],[-23,5],[-60,86],[-43,32],[-58,70],[18,44],[71,-41],[69,-58],[52,-99],[-26,-39]],[[13251,24486],[2,-41],[-40,-68],[-62,-23]],[[13709,23417],[13,-163],[-1,-52],[-31,-92],[-30,9],[-10,123],[-37,39],[-51,76],[-25,100],[26,63],[14,71],[81,-51],[29,-34],[22,-89]],[[14389,21581],[8,-47],[-46,17],[6,43],[32,-13]],[[14924,21671],[-57,30],[-75,75],[58,-4],[31,-31],[43,-70]],[[14237,21710],[-28,4],[-54,53],[4,61],[27,18],[24,-15],[25,-50],[2,-71]],[[14640,21961],[-21,17],[-28,93],[16,53],[29,-13],[34,-86],[-30,-64]],[[13603,23197],[-9,-23],[-50,72],[-27,60],[17,35],[53,-66],[16,-78]],[[13883,22666],[-46,19],[-16,60],[32,42],[32,-63],[-2,-58]],[[13760,23141],[-18,82],[-1,70],[20,43],[32,-20],[-20,-145],[-13,-30]],[[13499,23444],[-26,16],[-25,54],[13,24],[38,-94]],[[13539,23508],[-43,13],[20,73],[21,-13],[2,-73]],[[32113,20521],[-3,72],[82,86],[42,27],[0,-42],[-30,-2],[-84,-111],[-7,-30]],[[31396,20829],[-7,-37],[-24,-32],[-12,57],[43,12]],[[30774,19188],[-38,-28],[15,70],[27,14],[-4,-56]],[[28895,19591],[-58,-1],[-21,25],[45,23],[46,72],[13,2],[-25,-121]],[[34272,21781],[32,-43],[-72,-36],[-14,36],[54,43]],[[31309,28344],[85,-4],[27,-33],[-17,-35],[-49,-18],[-67,19],[21,71]],[[30367,27157],[-36,7],[7,54],[55,103],[63,-11],[26,-46],[-7,-30],[-47,-45],[-61,-32]],[[29786,28339],[-56,2],[-62,23],[-24,55],[-37,42],[-41,15],[0,31],[127,-36],[64,-37],[42,-61],[-13,-34]],[[31307,27733],[-155,101],[-11,44],[84,19],[93,-12],[34,-36],[-8,-62],[-37,-54]],[[31252,27996],[51,-60],[-23,-15],[-65,33],[-13,67],[50,-25]],[[30112,26542],[-40,-33],[-7,74],[37,13],[10,-54]],[[31424,27221],[-10,-35],[-101,58],[-7,46],[52,7],[66,-76]],[[32367,25005],[-55,4],[29,48],[33,-22],[-7,-30]],[[32160,25781],[23,-15],[6,-56],[-44,3],[-49,38],[-10,27],[45,16],[29,-13]],[[30455,31920],[-79,20],[59,51],[65,-24],[-45,-47]],[[31901,30648],[-40,8],[55,53],[63,1],[-78,-62]],[[27288,22911],[-38,1],[-34,29],[86,45],[17,-14],[-31,-61]],[[28879,31020],[-207,3],[-55,20],[-48,80],[11,42],[90,14],[138,-37],[82,-3],[41,-23],[8,-91],[-60,-5]],[[27704,28811],[48,-20],[46,-74],[-16,-50],[-96,-38],[-86,64],[-59,60],[-11,48],[75,20],[99,-10]],[[22483,33659],[-86,31],[7,69],[57,30],[73,-59],[-51,-71]],[[22072,34137],[-194,-29],[-114,34],[11,24],[122,42],[233,50],[100,2],[-15,-53],[-143,-70]],[[24301,31858],[6,-32],[-51,-48],[-36,55],[81,25]],[[24220,31756],[-22,-6],[-47,76],[40,6],[29,-76]],[[28532,31310],[-119,64],[-25,26],[31,83],[84,-20],[51,-97],[-22,-56]],[[27399,31230],[-53,21],[6,52],[55,9],[34,-44],[-42,-38]],[[28035,28772],[-29,-5],[-75,33],[-85,71],[63,48],[103,-55],[30,-36],[-7,-56]],[[27275,32047],[18,-39],[-52,-42],[-90,-11],[-21,-44],[-57,7],[-79,58],[-91,14],[35,38],[63,24],[58,2],[67,-32],[104,36],[45,-11]],[[27662,32009],[17,-39],[-139,-69],[-88,10],[58,59],[49,6],[71,46],[32,-13]],[[27769,29078],[-20,-24],[-59,9],[-1,27],[80,-12]],[[27336,31564],[-54,40],[45,64],[95,44],[26,52],[39,14],[18,55],[73,8],[29,-42],[-84,-81],[-88,-118],[-99,-36]],[[27948,31717],[-61,-3],[-45,69],[16,73],[37,13],[122,-19],[16,-39],[-85,-94]],[[25209,32213],[83,-4],[27,-45],[-61,-21],[-126,22],[-36,43],[59,22],[54,-17]],[[26089,30002],[136,-51],[-41,-34],[-70,32],[-49,49],[-84,-8],[-7,74],[-85,76],[12,34],[126,-53],[69,-55],[-7,-64]],[[25296,30996],[-73,38],[-13,29],[17,73],[-31,46],[71,105],[37,-9],[41,-53],[10,-131],[-18,-57],[-41,-41]],[[25827,29889],[-30,-3],[-41,61],[-39,25],[-18,43],[-3,99],[58,-3],[68,-110],[25,-70],[-20,-42]],[[21535,31544],[-41,-40],[-70,28],[2,106],[54,20],[65,-64],[-10,-50]],[[21509,32409],[-86,26],[-16,51],[68,2],[34,-79]],[[22834,31937],[36,-34],[-32,-81],[-87,50],[-50,12],[33,70],[46,9],[54,-26]],[[21271,31847],[-33,22],[53,51],[27,-40],[-47,-33]],[[21085,31431],[-116,40],[-1,22],[71,61],[51,-14],[29,-38],[-34,-71]],[[20341,31339],[-44,6],[-74,38],[16,46],[85,-10],[36,-30],[-19,-50]],[[19413,30825],[-70,14],[7,102],[38,-4],[0,-55],[25,-57]],[[19360,30621],[79,-3],[-8,-40],[-36,-18],[-35,61]],[[26060,19938],[-41,-19],[-33,91],[35,8],[41,-25],[-2,-55]],[[25570,20755],[-70,-12],[9,44],[51,-12],[10,-20]],[[24893,21316],[39,14],[41,-12],[-14,-38],[-71,-14],[5,50]],[[26194,19859],[21,-24],[-32,-40],[-26,20],[37,44]],[[26300,19777],[-55,22],[10,32],[44,4],[33,23],[22,-68],[51,19],[69,45],[52,-58],[24,44],[32,12],[12,-35],[37,-46],[47,-15],[-45,-103],[-35,-28],[-40,15],[-82,61],[-114,58],[-62,18]],[[28157,18958],[-25,-13],[-26,-44],[-17,29],[15,25],[53,3]],[[24827,5763],[-15,-58],[61,-4],[12,-64],[-3,-53],[-31,-141],[-17,-121],[18,-48],[-17,-64],[-8,-102],[9,-118],[-14,-169],[-25,-74],[-44,-101],[-36,-22],[-51,-118],[-4,-65]],[[24922,5483],[-1,39],[16,79],[13,-8],[-28,-110]],[[27765,8722],[-27,-44],[20,-23],[18,50],[15,-40],[12,-137],[-15,-67],[-55,7],[-22,147],[-17,26],[-24,66],[83,35],[12,-20]],[[27884,9570],[-30,47],[24,82],[0,68],[11,102],[-9,37],[-25,30],[-50,116],[17,-5],[23,-45],[28,-30],[2,-26],[53,-56],[8,-100],[-36,-48],[-6,-146],[-10,-26]],[[29044,7156],[-9,-37],[-29,-72],[-66,-18],[-72,-3],[-5,85],[26,11],[18,34],[27,5],[52,-24],[27,26],[22,55],[13,-7],[-4,-55]],[[27741,8957],[0,-125],[-38,-48],[-45,-42],[-9,39],[-32,52],[-48,40],[-19,43],[53,13],[-15,52],[31,83],[7,54],[-14,87],[13,5],[36,-30],[16,-30],[1,-41],[36,-108],[27,-44]],[[27534,9992],[34,-16],[28,13],[50,-4],[40,15],[5,-41],[-85,-14],[-121,-67],[-35,12],[45,67],[10,71],[29,-36]],[[27850,9114],[-31,-11],[-18,33],[55,13],[-6,-35]],[[28759,7941],[-11,-30],[-49,24],[-11,47],[26,4],[12,-27],[33,-18]],[[28656,8629],[-22,-55],[-11,5],[6,69],[27,-19]],[[28718,7680],[22,56],[36,64],[17,15],[5,46],[-17,34],[7,42],[21,-18],[11,-73],[-38,-91],[-33,-28],[-31,-47]],[[28043,9356],[46,-57],[39,-22],[60,-98],[3,-24],[-7,-107],[-18,-104],[-10,37],[10,93],[14,46],[-2,48],[-57,99],[-36,15],[-34,46],[-44,-3],[19,71],[17,-40]],[[28315,8313],[-32,10],[-48,63],[2,28],[53,-81],[25,-20]],[[28543,8028],[-37,89],[-44,25],[-7,113],[-42,129],[26,-22],[30,-109],[12,-94],[35,-32],[25,-45],[2,-54]],[[28414,8697],[-44,-31],[16,65],[-24,33],[-46,133],[2,46],[38,-130],[58,-116]],[[29074,7317],[-37,-25],[-3,30],[18,24],[22,-29]],[[32167,5028],[-9,-20],[-30,8],[-8,43],[19,37],[37,-36],[-9,-32]],[[32159,5304],[-27,11],[-2,60],[28,-27],[1,-44]],[[31810,5557],[-31,0]]],"transform":{"scale":[0.0036210230499170754,0.0019523957147595781],"translate":[-178.19453125,7.2200683593749915]},"objects":{"ne_50m_admin_0_countries_lakes":{"type":"GeometryCollection","geometries":[{"arcs":[[0]],"type":"Polygon"},{"arcs":[[1]],"type":"Polygon"},{"arcs":[[[2]],[[3]],[[4]],[[5]],[[6]],[[7]],[[8,9,10,11]],[[12]],[[13]],[[14]],[[15]],[[16]],[[17]],[[18]],[[19]],[[20]],[[21]],[[22]],[[23]],[[24]],[[25]],[[26]],[[27]],[[28]],[[29]],[[30]],[[31]],[[32]],[[33]],[[34]],[[35]],[[36]],[[37]],[[38]],[[39]],[[40]],[[41]],[[42]],[[43]],[[44]],[[45]],[[46]],[[47]],[[48]],[[49]],[[50,51,52,53,54,55,56,57,58,59,60]],[[61]],[[62]],[[63]],[[64]],[[65]],[[66]],[[67]],[[68]],[[69]],[[70]],[[71]],[[72]],[[73]],[[74]],[[75]],[[76]],[[77]],[[78]]],"type":"MultiPolygon"},{"type":null},{"arcs":[[79]],"type":"Polygon"},{"type":null},{"arcs":[[[80]],[[81]],[[82]]],"type":"MultiPolygon"},{"type":null},{"arcs":[[83]],"type":"Polygon"},{"arcs":[[84]],"type":"Polygon"},{"arcs":[[85]],"type":"Polygon"},{"arcs":[[[86,87,88]],[[89]]],"type":"MultiPolygon"},{"arcs":[[90,91,92,93]],"type":"Polygon"},{"type":null},{"arcs":[[[-53,94,95,96,97]],[[98]],[[99]],[[100]],[[101]],[[102]],[[103]],[[104]],[[105]],[[106]]],"type":"MultiPolygon"},{"arcs":[[107]],"type":"Polygon"},{"arcs":[[-94,108,109,110,111]],"type":"Polygon"},{"arcs":[[[112,113]],[[114]],[[115]]],"type":"MultiPolygon"},{"arcs":[[-97,116,117,-111,118,119]],"type":"Polygon"},{"arcs":[[[120]],[[121]],[[122]],[[123]]],"type":"MultiPolygon"},{"arcs":[[124]],"type":"Polygon"},{"arcs":[[125,126]],"type":"Polygon"},{"type":null},{"arcs":[[-110,127,-119]],"type":"Polygon"},{"arcs":[[-113,128]],"type":"Polygon"},{"arcs":[[129]],"type":"Polygon"},{"arcs":[[[130]],[[131]],[[132]],[[133]],[[134]],[[135]]],"type":"MultiPolygon"},{"arcs":[[-88,136,-92,137]],"type":"Polygon"},{"arcs":[[[138]],[[139]],[[140]],[[-12,141,-61,142,-59,143,-57,144,-55,145]],[[146]],[[147]],[[148]],[[149]],[[150]],[[151]],[[152]],[[-51,153]],[[154]],[[155]],[[156]],[[157]],[[158]],[[159]],[[160]],[[161]],[[162]],[[163]],[[164]],[[165]],[[166]],[[167]],[[168]],[[169]],[[170]],[[-10,171]],[[172]],[[173]],[[174]],[[175]],[[176]],[[177]],[[178]],[[179]],[[180]],[[181]],[[182]],[[183]],[[184]],[[185]],[[186]],[[187]],[[188]],[[189]],[[190]],[[191]],[[192]],[[193]],[[194]],[[195]],[[196]],[[197]],[[198]],[[199]],[[200]],[[201]],[[202]],[[203]],[[204]],[[205]],[[206]],[[207]],[[208]],[[209]],[[210]],[[211]],[[212]],[[213]],[[214]],[[215]],[[216]],[[217]],[[218]],[[219]],[[220]],[[221]],[[222]],[[223]],[[224]],[[225]],[[226]],[[227]],[[228]],[[229]],[[230]]],"type":"MultiPolygon"},{"arcs":[[[-96,231,-117]],[[232]]],"type":"MultiPolygon"},{"type":null},{"arcs":[[[233]],[[234]],[[235]],[[236]],[[237]],[[238]],[[239]],[[240]],[[241]],[[242]],[[243]],[[244]],[[245]],[[246]]],"type":"MultiPolygon"},{"arcs":[[[247]],[[248]]],"type":"MultiPolygon"},{"arcs":[[-126,249]],"type":"Polygon"}]}}}
},{}],19:[function(_dereq_,module,exports){
module.exports={"type":"Topology","arcs":[[[18577,12344],[-320,0],[-348,0],[-232,0],[-207,0]],[[17470,12344],[2,22],[-36,66],[-45,21],[-64,106],[3,63],[-38,12],[20,82],[2,33],[-13,70],[9,26],[-18,72],[-20,116],[-55,94],[-37,20],[-21,-16],[-29,45],[-9,44],[-79,64],[13,27],[-9,27],[-47,54],[-56,30],[-10,45],[-66,65],[-9,51],[-27,48],[-43,54],[-29,71],[-23,-1],[-19,-37],[-30,10],[-27,76],[-56,81],[-21,55],[-73,-7],[-44,34],[-20,32],[25,56],[-74,52],[-16,-37],[-35,-14],[-10,49],[-34,71],[9,26],[-50,108],[-28,4],[-20,66],[-24,20],[-10,44],[-52,42],[-8,-49],[-65,35],[-38,73],[-88,23],[-26,46],[-20,13],[2,50],[39,0],[-19,42],[-20,5],[-31,69],[0,194],[0,388],[0,389],[0,389],[0,388],[0,389],[0,388],[0,195],[0,388]],[[15853,17901],[170,0],[341,0],[340,0],[171,0],[340,0],[341,0],[170,0],[340,0],[341,0],[170,0]],[[18577,17901],[0,-347],[0,-521],[0,-521],[0,-347],[0,-521],[0,-521],[0,-348],[0,-521],[0,-347],[0,-347],[0,-348],[0,-347],[0,-521]],[[17470,12344],[-142,0],[-232,0],[-167,0]],[[16929,12344],[-269,0]],[[16660,12344],[-145,0],[-232,0],[-349,0],[-232,0],[-232,0],[-233,0],[-144,0]],[[15093,12344],[-37,41],[-21,-7],[-13,-34]],[[15022,12344],[-10,0]],[[15012,12344],[-17,64],[-13,14],[-9,57],[44,15],[-29,29],[-28,-2],[15,124],[-28,-6],[-3,-60],[-27,-4],[-26,-27],[-89,43],[-46,60],[-7,30],[16,38],[19,-10],[29,-66],[15,37],[-30,39],[-8,30],[-1,93],[-18,-5],[0,-40],[-46,-50],[-38,-10],[-55,18],[-82,107],[0,27],[-41,93],[-14,49],[-20,32],[53,111],[-3,95],[-16,-8],[1,-46],[-15,-64],[-20,-40],[-41,-19],[-73,11],[-30,39],[-6,-45],[-24,-11],[-21,22],[-52,-10],[-74,18],[-46,3],[-3,39],[13,31],[-40,6],[0,23],[32,41],[-34,30],[10,48],[-21,2],[-19,-25],[-89,-11],[-27,-13],[-57,24],[-47,37],[-73,82],[-2,59],[23,81],[58,91],[38,23],[85,8],[-18,24],[-83,-5],[-29,-14],[-55,-102],[-18,6],[-29,58],[-14,86],[11,53],[28,57],[-36,-1],[17,101],[98,83],[53,20],[31,-12],[-22,64],[6,20],[52,61],[-10,12],[-132,-163],[-75,-35],[-30,-70],[-25,-101],[-29,-82],[-25,106],[-45,81],[87,80],[-6,56],[-33,-10],[-25,14],[45,213],[1,25],[-71,-41],[-23,25],[-20,134],[-35,52],[-59,42],[-58,20],[-24,84],[15,54],[49,38],[23,-26],[-15,-53],[20,-5],[79,-65],[20,25],[-56,43],[-21,53],[-3,44],[13,32],[50,24],[-47,30],[-51,-45],[-19,6],[-67,-101],[-9,-82],[-12,-43],[-49,-23],[-27,-48],[-71,81],[-24,70],[-44,12],[-52,40],[-19,35],[27,73],[40,54],[12,80],[69,16],[45,33],[-75,-2],[-50,-23],[-56,45],[-29,41],[-9,35],[10,30],[6,70],[17,41],[24,15],[90,219],[30,59],[60,92],[-39,-12],[-30,47],[-7,-108],[-40,-104],[-29,-2],[26,79],[3,47],[-16,58],[-13,80],[2,61],[16,69]],[[13122,15825],[1,64],[-21,47],[-32,-13],[-54,20],[-17,55],[-47,16],[-25,39],[-93,33],[-134,98],[-68,-5],[-12,102],[-64,33],[19,85],[-83,27],[29,60],[-87,152],[-100,191],[-97,174],[-41,58],[6,37],[-114,149],[-108,73],[-22,25],[-39,94],[-64,57],[-73,45],[-35,78],[6,69],[-116,108],[-125,-66],[-105,-29],[20,-36],[-27,-39],[-33,2],[0,-91],[-30,-64],[-64,-1],[-127,-82],[-44,-43],[-42,163],[-283,279],[-20,62],[-95,51]],[[10662,17902],[262,-1],[259,0],[259,0],[259,0],[260,0],[259,0],[259,0],[259,0],[260,0],[259,0],[259,0],[259,0],[260,0],[259,0],[259,0],[260,0]],[[14813,17901],[390,0],[260,0],[390,0]],[[13071,15389],[11,56],[37,82],[-6,57],[-17,57],[-7,66],[9,63],[24,55]],[[13122,15825],[-27,-100],[4,-55],[20,-83],[5,-67],[-53,-131]],[[13071,15389],[2,-39],[-39,-67],[-51,-34],[-11,11]],[[12972,15260],[72,89],[27,40]],[[12876,15113],[-7,32],[13,45],[35,-6],[-8,-38],[-33,-33]],[[12405,14936],[25,-30],[60,19],[23,-15],[26,-61],[3,-41],[-16,-30],[-96,-64],[8,-18],[95,16],[21,109],[-7,68],[53,4],[75,50],[-10,-78],[-32,-74],[-19,-64],[-9,-64],[-1,-105],[-23,-58],[-92,-38],[-47,3],[-62,58],[15,29],[71,-9],[-65,55],[-47,25],[-19,61],[-38,75],[-12,72],[9,113],[111,-8]],[[13064,14851],[-40,-58],[-42,37],[-43,17],[55,62],[27,6],[27,-26],[16,-38]],[[13170,14451],[-23,5],[-59,85],[-42,31],[-58,69],[18,44],[70,-41],[68,-57],[52,-97],[-26,-39]],[[13356,14426],[-43,12],[19,73],[21,-14],[3,-71]],[[13523,14336],[13,-161],[-1,-51],[-18,-34],[-13,-57],[-29,9],[-10,122],[-18,40],[-18,-2],[-51,75],[-3,35],[-22,37],[0,27],[25,61],[14,71],[17,-3],[63,-48],[29,-34],[22,-87]],[[12651,14465],[28,-46],[8,-42],[-3,-50],[-74,-28],[48,-55],[74,-28],[-32,-40],[-4,-22],[56,-82],[34,-13],[-19,-37],[51,-13],[7,-37],[-29,-33],[-55,43],[-38,82],[-67,71],[-77,107],[-20,15],[-25,76],[31,33],[-89,37],[-15,37],[49,-5],[84,22],[39,26],[38,-18]],[[13316,14363],[-26,16],[-25,53],[13,23],[38,-92]],[[13418,14119],[-8,-23],[-50,71],[-27,59],[18,35],[52,-65],[15,-31],[0,-46]],[[13573,14064],[-18,81],[-1,68],[20,43],[32,-20],[-21,-143],[-12,-29]],[[12848,13842],[-20,27],[-4,44],[30,-3],[-6,-68]],[[13694,13596],[-15,-8],[-30,27],[-16,58],[32,41],[31,-62],[-2,-56]],[[13892,13175],[136,-63],[135,-31],[99,-37],[61,-11],[37,-21],[67,-154],[82,-142],[1,-44],[27,-57],[51,-52],[40,-24],[85,-40],[51,-39],[48,-71],[19,-67],[34,-64],[35,-122],[11,47],[19,-2],[28,-109],[-13,-25],[-31,11],[-24,-42],[-16,-5],[-88,26],[-210,107],[-49,28],[-69,54],[-4,31],[55,68],[2,29],[-66,-11],[-53,3],[-26,-29],[-23,6],[-32,39],[-46,31],[35,24],[-7,31],[-30,39],[-31,-15],[-19,61],[-22,27],[-39,11],[-16,-31],[-32,34],[-35,-15],[2,86],[53,36],[-67,52],[-25,57],[-60,34],[-20,-31],[-37,-2],[-25,57],[4,54],[-43,-24],[-32,56],[-57,0],[-26,-21],[-25,3],[9,43],[-6,51],[-30,17],[15,50],[21,13],[51,4],[42,-26],[-10,85],[-62,6],[6,-37],[-36,-18],[-53,0],[-21,11],[-36,45],[-23,44],[1,24],[29,42],[38,15],[50,1],[56,-20],[140,-91]],[[14441,12901],[-21,17],[-27,91],[16,52],[29,-12],[33,-86],[-30,-62]],[[14497,12867],[-13,53],[24,16],[8,-48],[-19,-21]],[[14044,12653],[-28,4],[-53,53],[11,32],[-8,28],[27,17],[24,-14],[24,-49],[3,-71]],[[14721,12615],[-56,29],[-74,74],[3,20],[54,-24],[31,-30],[42,-69]],[[14193,12526],[8,-46],[-45,16],[6,43],[31,-13]],[[14917,12223],[-22,-2],[-18,89],[31,-30],[9,-57]],[[24312,16311],[-281,-321],[-191,-227],[-190,-227],[-250,-316],[-222,-291],[-139,-185],[-209,-233],[-209,-234],[0,-436],[0,-218],[-1,-435],[0,-436],[1,-218]],[[22621,12534],[-2,-191],[-180,1],[-233,0],[-149,0]],[[22057,12344],[-199,0],[-232,0],[-233,0],[-232,0],[-232,0]],[[20929,12344],[-10,215],[-22,431],[-22,431],[-11,215],[-21,431],[-11,215],[-22,431],[-21,430],[-11,215],[-22,431],[0,396],[0,396],[0,396],[0,264],[0,264],[0,396]],[[20756,17901],[245,0],[245,0],[368,0],[368,0],[367,0],[376,1]],[[22725,17902],[-3,-25],[2,-240],[-3,-106],[-22,-91],[25,-34],[17,-59],[37,-18],[43,-62],[37,-15],[1,-28],[-14,-162],[17,20],[17,125],[23,56],[19,12],[75,6],[80,-14],[56,6],[34,-31],[15,-104],[70,-209],[38,-150],[58,-156],[14,-98],[-48,-119],[85,3],[15,-18],[63,27],[247,90],[58,8],[83,-17],[141,-87],[77,-35],[123,-34],[35,-16],[72,-16]],[[31100,10839],[-35,-24],[-43,-66]],[[31022,10749],[-45,-14],[-15,30],[-66,-115],[-32,-6],[-44,-42],[-61,-36],[-72,-68],[-92,-58],[-61,17],[-91,-81],[-19,24],[-53,-30],[-45,-9],[-13,40],[-45,-1],[-11,13]],[[30257,10413],[-34,16],[-32,-14],[-29,52],[12,51],[-18,49],[19,28],[-2,37],[-44,9],[-52,41],[-3,402],[-3,295],[-79,97],[-38,36],[-20,4],[-46,-34],[-52,-17],[-43,-25],[-48,17],[-14,29]],[[29731,11486],[116,64],[36,30],[32,55],[1,172],[72,0],[0,34],[136,0],[9,-29],[58,-43],[16,16],[68,16],[29,-10],[7,33],[60,29]],[[30371,11853],[20,-5],[55,28],[19,-3],[41,-37],[98,-39],[26,-26],[24,-82],[16,-13],[34,8],[38,41],[32,22],[61,18],[-12,-27],[47,2],[47,-37],[-41,-78],[-16,-102],[-48,-68],[-63,-67],[34,-26],[41,20],[58,-14],[-21,-87],[25,-96],[24,-14],[10,-80],[13,-44],[-1,-35],[24,-22],[5,-36],[89,-10],[19,-14],[62,-14],[23,-29],[-54,-48]],[[30970,11785],[-7,-37],[-24,-31],[-12,56],[43,12]],[[30356,10167],[-37,-27],[14,68],[27,14],[-4,-55]],[[32988,13581],[-1,281],[-221,0],[-332,0],[-333,0],[-332,0],[-221,0],[-333,0],[-15,30],[-28,120],[-32,24],[-19,42],[-30,15],[12,34],[57,43],[78,19],[31,20],[-38,44],[0,83],[-33,81],[-35,2],[-32,-19],[-16,-86],[-45,-21],[-11,-20],[-1,-132],[22,-88],[-7,-97],[-36,-44],[-20,-130],[12,-57],[-59,-61],[-25,-4],[-27,61],[19,52],[-37,41],[-46,75],[-56,51],[-71,24],[-21,-3],[-21,-52],[-21,-10],[-41,15],[-28,-25],[-41,28],[-17,50],[-50,47],[-18,-74],[-25,37],[-7,57],[19,39],[-9,122],[16,41],[18,78],[-24,49],[-66,-42],[-15,-92],[-35,-21],[-37,6],[-32,29],[4,83],[23,77],[9,130],[-38,44],[-24,64],[-80,29],[-3,99],[-78,119],[20,55],[47,48],[-30,136],[23,34],[60,-28],[45,50],[-38,105],[-22,47],[-18,60],[37,13],[36,-63],[82,-79],[39,-17],[20,26],[-16,75],[5,25],[-15,67],[6,42],[-20,41],[17,12],[32,-31],[103,-125],[43,-34],[54,-18],[24,-51],[33,-30],[151,-13],[83,7],[129,-4],[118,-46],[43,-11],[25,16],[46,138],[3,35],[-20,50],[11,41],[96,71],[1,104],[-40,14],[-31,27],[-77,134],[0,32],[38,34],[-5,17],[-116,36],[25,67],[-57,35],[-4,48],[49,23],[0,44],[-42,83],[3,24],[53,76],[8,68],[23,54],[9,49],[-11,46],[15,39],[0,93],[35,45],[-62,22],[-22,36],[-33,1],[-24,45],[-26,72],[-43,41],[4,35],[74,96],[69,51],[-10,48],[-42,-5],[-9,49],[30,21],[111,24],[13,20],[-24,30],[-37,16],[-116,-46],[-22,8],[-16,40],[-140,32],[-7,46],[17,14],[103,-25],[17,23],[-40,37],[-11,55],[4,87],[34,39],[-18,26],[-74,-31],[-28,38],[8,139],[34,45],[-46,41],[-13,32],[63,30],[2,26],[-40,18],[-7,36],[55,19]],[[30941,18055],[48,-39],[5,-29],[-30,-39],[67,-15],[27,-46],[9,-90],[26,14],[21,-35],[13,-55],[49,-66],[-53,-67],[45,-16],[65,-8],[34,-70],[-41,-74],[66,17],[38,-113],[37,10],[34,-57],[14,-47],[-62,-64],[-32,-13],[39,-53],[65,33],[66,-3],[-22,-78],[-37,-60],[-69,-37],[23,-17],[107,33],[27,-2],[50,-92],[51,-4],[43,-27],[16,-25],[-26,-46],[18,-51],[-10,-29],[-101,-64],[9,-22],[59,6],[45,-16],[19,-20],[-34,-68],[32,-25],[63,-7],[80,-88],[-15,-79],[5,-87],[-43,-13],[-145,22],[-83,34],[-8,-18],[76,-32],[34,-39],[37,-3],[32,-29],[-44,-11],[-11,-41],[116,-17],[4,-31],[-56,-20],[-2,-29],[42,-12],[53,4],[17,-85],[-40,-26],[71,-20],[27,-45],[25,-2],[45,48],[65,-72],[69,-16],[-18,-68],[58,-85],[-46,-58],[-24,-84],[30,-3],[35,22],[77,11],[46,26],[18,-21],[1,-37],[47,1],[9,-34],[-70,-95],[-2,-28],[87,85],[46,73],[48,-8],[24,-17],[30,-100],[29,-57],[104,-33],[48,20],[7,27],[64,5],[9,-52],[56,-50],[65,-17],[22,-36],[-43,-65],[-37,-27],[-52,-2],[-72,-17],[-18,-32],[-59,-29],[-54,-91],[-56,-2],[-38,-17],[-16,-24],[-45,-17],[-80,-48],[-89,0],[-52,-37],[-36,-83],[-4,-34],[16,-22],[-51,-48],[4,-37],[34,-6],[44,44],[43,56],[57,17],[81,54],[53,50],[57,66],[73,52],[89,38],[65,19],[-25,21],[-38,-9],[2,44],[32,14],[157,-19],[55,-14],[59,-121],[13,-38],[-2,-46],[-90,-63],[28,-36],[15,-47],[87,103],[46,34],[103,13],[30,-71],[66,-19],[40,-59],[29,-64],[16,-66],[-2,-39],[-23,-68],[17,-62],[-12,-72],[19,-47],[-12,-10],[-87,-24],[31,-19],[58,-15],[37,-40],[-35,-60],[39,-34],[-2,-53],[-87,-105],[-73,-67],[-72,-58],[-116,-113],[-34,-7]],[[31723,16666],[23,-15],[6,-55],[-43,3],[-49,37],[-10,27],[44,16],[29,-13]],[[31927,15901],[-54,4],[29,47],[32,-21],[-7,-30]],[[33435,13628],[-34,-75],[-40,-15],[-82,-15],[-4,-34],[12,-32],[51,-4],[8,-30],[-4,-49],[-20,-64],[-56,-64],[-16,-65],[-67,-119],[-19,-18],[-8,-55],[-58,-107],[-15,-47],[-1,-42],[-26,-60],[25,-57],[70,110],[34,32],[35,79],[34,-13],[-11,-54],[64,31],[45,-29],[63,18],[1,-30],[-52,-54],[-87,-63],[46,-18],[-6,-45],[-41,-62],[62,29],[37,-41],[82,28],[6,-55],[42,14],[-16,-71],[4,-34],[60,53],[31,6],[53,36],[37,52],[27,-11],[5,-64],[29,37],[109,13],[26,-8],[66,-53],[13,-29],[-1,-62],[-50,-53],[-24,-56],[-69,-71],[67,14],[-31,-53],[22,-11],[51,-3],[-26,-77],[-74,-54],[73,15],[42,31],[25,-15],[39,41],[52,8],[44,41],[10,-44],[-11,-46],[-20,-40],[-54,-23],[-20,-31],[-34,-32],[-10,-72],[-23,-32],[-12,-53],[-45,-61],[6,-37],[20,-38],[44,6],[27,41],[61,128],[53,36],[45,39],[14,-26],[-23,-33],[-43,-110],[-17,-80],[1,-71],[13,-29],[49,49],[19,34],[25,76],[19,-12],[16,-99],[-8,-62],[-62,-164],[6,-65],[-20,-78],[-29,-70],[-27,-17],[-34,35],[-24,-2],[-42,-40],[-22,24],[9,140],[-8,27],[9,62],[-16,2],[-61,-100],[-40,-51],[-18,-10],[-27,31],[5,44],[17,60],[61,152],[5,27],[-14,34],[-24,125],[-16,25],[-39,27],[-11,-44],[-47,-109],[-8,-47],[-30,-19],[-85,-20],[-64,-107],[-13,-65],[-48,-70],[-59,4],[-70,-24],[-44,31],[-1,23],[31,49],[63,25],[32,20],[36,50],[46,96],[58,34],[29,28],[-20,25],[-35,18],[-80,2],[-19,-81],[-38,-18],[-54,17],[-24,16],[5,130],[-48,-12],[-36,5],[-19,-52],[-65,-27],[-71,-13],[-14,-13],[-49,4],[-142,29],[-51,-3],[-72,25],[-111,1],[-1,27],[-47,-39],[-28,-14],[-90,-23],[-48,-5],[-39,32],[-16,52],[-12,65],[25,66],[85,83],[68,83],[29,44],[73,52],[-43,-1],[-62,14],[-61,-9],[-60,0],[-1,18],[89,95],[-17,-49],[17,-22],[34,-4],[61,205],[25,41],[59,-12],[24,8],[13,52],[-38,40],[-7,64],[29,57],[25,22],[44,-35],[15,18],[-26,46],[-8,51],[58,164],[67,221],[55,94],[-9,25],[45,44],[34,77],[8,45],[47,59],[15,66],[19,29],[44,33],[85,45],[49,41],[34,-2],[10,-29],[48,-18],[6,54],[40,9],[17,-30]],[[33807,12723],[31,-42],[-71,-35],[-13,35],[53,42]],[[31100,10839],[42,-27],[47,-51],[53,8],[41,-32],[18,-30],[50,16],[54,-3],[57,-19],[-13,-46],[44,6],[28,-20],[17,22],[56,34],[71,69],[9,-9],[12,-69],[28,-29],[32,-7],[45,22],[39,-57],[19,-67],[-49,-37],[96,-9],[21,-29],[-19,-30],[-18,12],[-32,-11],[-28,-25],[-50,-16],[-61,-35],[-63,-45],[-65,-30],[-68,-46],[-141,-65],[-34,-2],[-41,-35],[-65,6],[0,-56],[-18,-32],[-41,3],[-35,30],[-30,50],[-27,-79],[-18,50],[-33,-19],[-14,-53],[17,-56],[-28,-15],[-25,-60],[-29,-22],[-31,-61],[-37,-46],[-10,-31],[-61,-71],[-40,0],[-26,-30],[-4,-60],[-23,6],[-15,-22],[-22,18],[-27,-10],[-20,14],[-41,118],[-25,10],[-16,-37],[-24,37],[-18,134],[0,32],[25,113],[63,102],[-13,23],[64,18],[26,39],[-19,17],[205,182],[41,30],[83,38],[24,-9],[2,-51],[32,-12],[39,52],[94,47],[78,6],[25,22],[-67,15],[-79,-9],[-50,17],[-68,-11],[-71,11],[-40,-44],[-46,26],[24,51],[73,75],[66,106]],[[31897,10805],[46,20],[18,-19],[-46,-33],[18,-22],[-33,-46],[9,-22],[35,22],[34,-5],[29,18],[62,87],[-63,-6],[-8,11],[44,53],[-4,28],[22,44],[57,55],[19,-59],[37,5],[62,-23],[4,-24],[-23,-47],[25,-39],[-48,-31],[-51,-69],[-49,-45],[-79,-32],[-54,10],[-57,-14],[-55,-5],[-45,72],[-13,114],[4,60],[20,55],[46,67],[101,250],[69,90],[35,17],[22,-3],[-5,-40],[25,-79],[-14,-78],[-27,-101],[-3,-72],[-22,-50],[-46,-40],[-46,-24],[-52,-50]],[[31936,10575],[-31,24],[13,20],[31,-5],[-13,-39]],[[14813,17901],[-48,20],[-1,44],[-56,114],[-13,61],[-41,10],[-54,85],[-3,20],[29,44],[-21,83],[-61,4],[-21,-51],[-73,-8],[-62,-27],[-90,25],[-42,25],[-38,-41],[-60,24],[-9,-32],[-72,4],[-40,-19],[-45,1],[-16,15],[-15,130],[-40,9],[18,28],[1,56],[-14,67],[-25,48],[-36,24],[-85,11],[-105,90],[-27,65],[-31,12],[-52,86],[-44,42],[-59,-25],[-28,22],[-70,7],[-36,20],[10,27],[-17,38],[15,40],[-13,34],[21,37],[-76,43],[-10,38],[-26,19],[-27,60],[2,50],[23,70],[-54,12],[-21,48],[-42,37],[-12,29],[42,29],[28,56],[-12,34],[-47,41],[-23,64],[-100,80],[-40,15],[-51,50],[-26,46],[3,40],[-48,92],[-68,29],[-62,-40],[-49,0],[-7,23],[20,45],[-48,32],[-43,54],[-73,39],[-64,9],[-13,23],[26,29],[3,38],[39,41],[-45,28],[-17,45],[-50,21],[45,28],[15,30],[54,52],[41,78],[-27,60],[-71,71],[58,49],[-74,38],[-57,-51],[-40,3],[8,47],[-79,-11],[-50,-25],[-48,3],[-24,67],[30,41],[-4,55],[-57,10],[14,68],[23,7],[6,57],[-36,43],[-8,84],[-64,64],[-7,30],[-330,0],[-234,0],[-26,42],[1,47],[31,62],[-20,60],[2,47],[-26,51],[-37,25],[-12,41],[0,294],[0,286]],[[11373,22392],[88,-7],[70,-25],[137,-69],[36,-1],[-48,68],[-63,32],[-65,13],[1,42],[77,4],[-87,52],[59,101],[53,13],[66,-7],[7,51],[24,22],[79,9],[98,-9],[11,100],[57,1],[45,-56],[54,-25],[-19,-40],[-59,-75],[59,11],[129,52],[85,15],[37,37],[31,62],[20,10],[85,-4],[-11,28],[64,27],[46,-23],[63,24],[134,83],[83,-5],[58,73],[71,38],[45,8],[88,-29],[63,2],[58,52],[53,-44],[-30,-55],[-221,-110],[-68,-47],[-69,-27],[-246,-42],[-34,-16],[-38,-58],[-55,-48],[-97,-24],[-36,-27],[-74,-90],[-79,-72],[21,-19],[163,-17],[-18,100],[14,34],[47,31],[51,13],[120,62],[38,50],[-2,22],[61,15],[60,-14],[49,10],[41,-13],[27,-73],[69,100],[84,88],[64,32],[148,54],[84,15],[35,-12],[0,-42],[45,-19],[72,50],[87,75],[30,75],[140,58],[22,19],[-79,16],[-4,36],[-50,28],[12,53],[37,25],[65,-28],[102,-75],[65,-66],[58,-89],[46,-106],[40,-62],[98,-93],[51,-40],[92,-44],[93,-15],[58,40],[1,26],[-51,74],[37,53],[123,130],[-50,27],[56,19],[52,36],[30,-21],[-7,-97],[17,-76],[77,-39],[-3,-19],[-94,-115],[15,-23],[223,-1],[40,22],[27,38],[41,23],[28,99],[23,37],[260,3],[147,-21],[155,-58],[76,-47]],[[15667,22731],[1,-478],[0,-312],[209,-118],[209,-119],[209,-118],[278,-158],[209,-118],[278,-158],[348,-197],[278,-158],[209,-118],[500,-3],[157,-150],[157,-150],[82,-45],[66,-17],[237,-35],[357,-52],[356,-52],[237,-35],[356,-52],[357,-52],[3,-75],[-5,-117],[0,-364],[0,-486],[0,-243],[1,-485],[0,-365]],[[20756,17901],[-409,0],[-408,0],[-272,0],[-409,0],[-409,0],[-272,0]],[[15925,25026],[72,45],[72,-2],[49,-77],[77,74],[55,27],[157,10],[123,-17],[86,-31],[130,-72],[208,-141],[122,-65],[32,-59],[-15,-32],[-148,-59],[-134,-35],[-159,-74],[-108,-35],[-142,-68],[-266,-110],[-46,-38],[-76,-126],[-28,-29],[-70,-30],[-112,-15],[-4,-44],[-47,-120],[-21,-130],[-48,-63],[-84,-30],[-148,-29],[-76,28],[-111,-90],[-107,-37],[-46,-33],[-59,-20],[-75,18],[-49,48],[-55,103],[-113,128],[-205,80],[-146,70],[-128,-6],[1,89],[49,23],[-17,27],[72,94],[107,83],[-11,29],[8,85],[19,25],[95,32],[-2,40],[-68,36],[66,84],[47,89],[84,54],[23,59],[64,63],[-21,30],[-59,15],[-47,48],[-118,200],[334,44],[230,14],[239,39],[118,-5],[207,-90],[166,-50],[104,-10],[-47,-61]],[[18577,24455],[0,-408],[0,-275],[0,-274],[0,-274],[0,-274],[-255,0],[-426,0],[0,-45],[-41,-39],[-61,-1],[-2,85],[-436,0],[-288,0],[-413,1],[-5,-58],[104,-97],[-45,-23]],[[16709,22773],[-67,78],[-25,126],[16,23],[159,38],[279,41],[255,28],[116,-2],[137,-19],[123,-3],[156,-20],[123,26],[110,-2],[16,40],[-106,48],[-281,86],[-166,37],[-157,-8],[-139,-27],[-128,-10],[-185,-7],[-92,18],[-181,-10],[-162,13],[-184,131],[-31,40],[29,34],[124,62],[137,27],[295,75],[68,55],[67,8],[71,-9],[37,22],[-68,18],[-325,-51],[-152,-5],[-105,-22],[-58,0],[-68,23],[10,45],[74,18],[37,49],[-171,-10],[-136,24],[-32,39],[11,112],[97,92],[103,48],[-9,24],[-64,59],[29,54],[66,50],[158,100],[80,42],[78,23],[108,48],[279,80],[249,80],[91,-21],[38,-31],[30,-100],[-15,-106],[-47,-62],[-59,-57],[11,-23],[136,36],[125,22],[13,85],[66,60],[124,2],[192,-49],[162,-62],[50,-26],[-23,-72],[-70,-69],[-77,-40],[60,-28],[99,78],[12,-46],[35,1],[97,85],[93,50],[64,14],[2,49],[-97,52],[-37,42],[8,32],[180,-13]],[[15667,22731],[106,-74],[120,-40],[268,-43],[211,-108],[72,-21],[190,-47],[29,4],[130,-19],[132,-11],[-51,51],[86,0],[14,15],[69,-7],[174,-62],[101,-53],[57,-43],[122,-132],[-36,-67],[-48,-9],[-134,11],[-24,-38],[-75,-32],[-10,-57],[-71,-48],[80,-53],[128,-12],[108,-32],[159,-16],[127,1],[92,-11],[121,3],[54,23],[180,13],[104,33],[55,-16],[169,83],[90,11],[28,-45],[29,-11],[55,-67],[110,-2],[39,-9],[24,-37],[12,-73],[37,-28],[27,73],[33,4],[28,-54],[39,-40],[98,-75],[22,-47],[-49,-47],[-36,-10],[-64,7],[96,-103],[93,-89],[14,29],[-10,111],[33,21],[34,-33],[27,6],[67,-34],[-52,88],[17,24],[-54,48],[-44,94],[-2,64],[-58,65],[-29,46],[3,44],[53,45],[1,63],[112,15],[61,24],[79,10],[46,44],[66,-8],[6,60],[36,32],[63,12],[70,61],[2,42],[-43,13],[-143,-61],[-41,-80],[-47,15],[-45,-6],[-54,-36],[-42,-4],[-87,18],[-34,-23],[2,-57],[-143,-12],[-116,65],[12,51],[90,117],[149,19],[90,21],[164,60],[32,5],[150,51],[100,-28],[30,-18],[40,-55],[30,-132],[64,-64],[77,-41],[75,-18],[-7,-33],[48,-43],[80,-16],[79,4],[67,15],[50,23],[213,-153],[131,-40],[89,6],[90,-26],[190,53],[65,5],[44,16],[186,-13],[82,-15],[89,-31],[61,1],[137,51],[-78,70],[-4,39],[24,15],[59,-42],[96,-111],[36,-29],[88,-41],[42,-7],[81,56],[-14,65],[-92,53],[-53,9],[-101,-38],[-22,10],[-67,71],[15,35],[-68,100],[44,26],[69,-33],[115,35],[-39,62],[30,12],[51,-26],[63,7],[39,-22],[79,-99],[135,-7],[-67,-91],[37,-7],[41,52],[99,43],[11,-40],[-37,-164],[-7,-58],[-47,-82],[0,-23],[50,-68],[5,-43],[42,-9],[80,23],[-17,-62],[69,7],[31,-27],[-1,-71],[-38,-21],[-57,-3],[-68,26],[-93,-14],[129,-126],[-43,92],[19,15],[95,-18],[37,4],[37,28],[9,87],[17,55],[-10,50],[-46,126],[-51,64],[52,144],[21,17],[70,14],[72,-21],[32,15],[133,114],[89,89],[76,35],[55,38],[-53,7],[-9,27],[-2,102],[-22,41],[-36,7],[-26,-82],[-61,-30],[-72,-9],[-33,30],[10,55],[131,107],[-47,14],[0,82],[27,14],[148,31],[-51,-62],[64,9],[9,92],[-32,20],[-99,-37],[-70,4],[-40,32],[-59,67],[-51,-37],[-128,46],[-113,56],[-71,12],[-42,35],[-40,60],[-61,68],[-16,43],[1,59],[87,108],[49,13],[38,66],[-76,-30],[-47,20],[-52,66],[0,41],[22,91],[-15,29],[28,25],[4,49],[73,62],[59,-2],[54,-33],[44,4],[43,79],[-116,17],[-2,36],[61,44],[28,47],[84,64],[128,40],[48,-2],[21,-66],[47,-42],[135,0],[13,-64],[138,-90],[84,-100],[18,-98],[-6,-77],[-21,-32],[54,-28],[39,-40],[69,-34],[47,-80],[45,-52],[-6,-36],[68,12],[70,-82],[-14,-16],[-103,-2],[-81,49],[-60,-69],[105,-10],[41,-23],[-84,-74],[-127,-90],[-14,-31],[84,16],[50,-6],[28,-35],[81,-36],[51,7],[92,53],[59,-15],[-33,-39],[93,-14],[107,-6],[38,-24],[-73,-14],[-38,-70],[-113,-2],[134,-91],[76,-122],[-19,-22],[7,-66],[-17,-109],[89,-104],[66,66],[34,62],[-1,53],[27,55],[32,141],[86,107],[80,20],[114,-105],[88,-52],[75,-72],[63,-206],[-23,-79],[-123,20],[9,-170],[32,-93],[54,-71],[136,-137],[22,-71],[41,-15],[124,110],[55,24],[20,33],[16,100],[29,56],[114,125],[60,188],[0,94],[58,70],[66,-14],[111,16],[-65,36],[7,32],[45,28],[15,56],[-69,46],[-35,-1],[-39,35],[-12,61],[4,98],[-23,53],[15,63],[-17,34],[133,-20],[102,23],[110,-8],[88,-45],[69,-23],[184,-7],[101,3],[67,-25],[-33,-56],[-40,-37],[90,-24],[44,-82],[54,14],[60,-9],[87,-30],[23,-40],[-77,-59],[-93,-47],[155,-30],[32,-35],[-5,-62],[-66,-52],[-131,-47],[-56,22],[-93,-30],[44,-54],[0,-27],[46,-71],[57,14],[-24,-92],[34,-65],[73,-70],[81,-65],[38,-68],[-8,-52],[-45,-145],[-45,-34],[-66,-7],[-36,-34],[-52,-79],[-48,-28],[-73,-61],[-84,-19],[-95,-80],[-61,-12],[-35,59],[-73,98],[-36,30],[-43,-1],[-32,23],[15,29],[-73,67],[-84,28],[-73,-61],[26,-18],[44,35],[33,-4],[114,-112],[26,-14],[35,-79],[73,-120],[-21,-39],[-106,40],[-50,-56],[-118,46],[-51,27],[-94,107],[-44,16],[-51,-18],[-74,-6],[-155,5],[-31,-57],[44,-48],[77,-27],[82,-42],[11,-34],[-23,-49],[-248,-249],[-34,-44],[-58,-44],[-105,-10],[-80,7],[-41,23],[-74,62],[-76,48],[-19,34],[-94,30],[-91,62],[-89,38],[-27,-32],[-73,5],[-55,22],[-53,-3],[-136,23],[-80,-1],[5,-35],[58,-5],[24,17],[106,-17],[183,-53],[88,-67],[98,-102],[73,-49],[211,-35],[220,-9],[99,-32],[0,-68],[-69,-120],[-164,-215],[-31,-79],[-30,-30],[-149,-88],[-45,-9],[-65,22],[-43,-26],[-51,28],[-58,-4],[-40,34],[-13,41],[-73,-26],[34,-24],[25,-51],[-16,-18],[-60,19],[35,-90],[-63,-82],[-34,-8],[-62,18],[-30,-32],[-65,-5],[-162,70],[-105,16],[-31,28],[-43,-19],[-96,27],[-115,55],[-125,46],[-45,6],[-19,-50],[29,-39],[39,-14],[10,56],[231,-94],[102,-35],[-83,-62],[47,-3],[58,38],[64,30],[42,-20],[97,-77],[98,-32],[61,-46],[10,-24],[5,-98],[-10,-47],[-39,-36],[-66,-12],[-64,-52],[-27,-8],[-115,16],[-45,14],[-89,-22],[4,-43],[28,-25],[79,-20],[-15,-53],[-54,23],[-94,-19],[-12,-39],[-73,-83],[-34,-10],[23,-36],[-44,-33],[14,-29],[-42,-30],[-49,-1],[21,-35],[-27,-16],[-57,6],[41,-65],[31,-17],[-28,-37],[-79,-52],[-63,-80],[14,-51],[-47,-13],[-2,-90],[-65,-136],[-71,-166],[-53,-22],[32,-41],[-7,-58],[-28,-151]],[[18577,24455],[106,-54],[39,-34],[94,-41],[42,-36],[-7,-30],[54,-14],[27,-35],[62,-182],[53,-120],[18,-75],[77,-57],[38,-5],[27,44],[93,52],[11,38],[-64,66],[-42,63],[-27,77],[-8,71],[-30,39],[-83,296],[9,39],[61,12],[-29,48],[16,23],[84,-10],[165,-69],[45,43],[127,-41],[110,-63],[73,-30],[51,-42],[57,-71],[49,-188],[97,-220],[31,-57],[67,-85],[36,-61],[10,-72],[-60,-167],[110,-89],[58,-83],[101,-67],[79,-29],[51,-32],[29,16],[68,-10],[104,-55],[103,-65],[70,1],[29,-19],[18,-57],[112,6],[48,-60],[18,-81],[-8,-77],[-29,-24],[-47,6],[-50,35],[-16,43],[-52,-69],[-58,20],[-87,55],[-45,-17],[-68,-48],[17,-49],[-24,-35],[-81,7],[-120,61],[-29,-21],[88,-74],[22,-7],[6,-59],[43,-13],[119,68],[80,6],[48,-41],[-24,-63],[66,-33],[9,-43],[-19,-54],[-172,-68],[-68,-13],[-42,-20],[-73,2],[-84,-10],[-163,29],[-78,31],[-59,-28],[-146,24],[-17,18],[41,64],[-66,6],[-74,20],[-166,14],[-54,32],[-6,79],[-48,31],[-33,-1],[-53,-47],[-49,-84],[-111,-90],[-115,-24],[-137,-10],[-51,-19],[-107,-70],[-143,-42],[-133,-23],[-138,-10],[-104,-16],[-76,5],[-50,-24],[-271,-13],[-152,-19],[-72,8],[-57,53],[-76,121],[6,61],[-27,119],[-104,29],[-68,9],[-227,-3],[-125,9],[-132,28],[-112,44],[-26,44],[-68,70]],[[24953,21806],[-73,38],[-12,29],[16,72],[-30,45],[20,46],[50,58],[36,-9],[40,-53],[11,-129],[-18,-56],[-40,-41]],[[28487,21831],[-205,3],[-54,19],[-47,79],[10,42],[90,13],[135,-36],[81,-4],[41,-22],[7,-90],[-58,-4]],[[19149,21638],[-28,-18],[-41,32],[7,100],[37,-3],[1,-55],[24,-56]],[[31467,21463],[-39,8],[55,53],[62,1],[-78,-62]],[[19097,21438],[78,-4],[-8,-39],[-36,-18],[-34,61]],[[25734,20827],[135,-50],[-41,-34],[-69,31],[-48,49],[-53,-18],[-30,10],[-7,73],[-84,75],[12,33],[124,-52],[68,-54],[-7,-63]],[[25476,20715],[-29,-2],[-41,60],[-38,24],[-18,43],[-10,76],[7,22],[57,-4],[67,-108],[25,-69],[-20,-42]],[[25409,20557],[21,-3],[63,97],[30,6],[64,-46],[50,-76],[98,-25],[86,-32],[51,-69],[63,-32],[88,-63],[108,-30],[60,-38],[72,-111],[30,-107],[-12,-48],[-50,-44],[51,-5],[104,27],[90,-21],[48,29],[60,-60],[8,-32],[86,-54],[-11,-20],[-111,-84],[-66,-58],[-26,-10],[-249,103],[-113,21],[-25,111],[-28,17],[-126,32],[5,52],[-13,17],[-66,-8],[-51,-23],[-34,-43],[-12,-71],[-18,-30],[-80,-77],[-78,-38],[-54,-99],[-35,-41],[-89,-56],[-118,-39],[-28,10],[-19,66],[-40,196],[-25,24],[-135,-25],[-75,3],[-74,-44],[-56,-2],[-34,12],[-3,39],[43,80],[40,47],[126,64],[46,43],[-28,121],[-5,96],[51,224],[20,206],[11,59],[30,86],[41,64],[70,44],[31,-37],[55,-25],[37,-87],[-37,-57],[50,-37],[37,-89]],[[27392,19916],[-20,-23],[-58,8],[-1,27],[86,11],[-7,-23]],[[27655,19614],[-29,-5],[-75,33],[-84,70],[63,48],[101,-55],[30,-35],[-6,-56]],[[27328,19653],[47,-20],[45,-73],[-15,-49],[-54,-30],[-41,-8],[-85,63],[-58,60],[-11,47],[74,20],[98,-10]],[[26205,19392],[14,-35],[-21,-78],[-24,-39],[-75,-67],[-49,-59],[-122,-98],[-31,-2],[-67,17],[-88,-40],[-17,73],[-39,50],[1,37],[44,47],[99,169],[23,9],[89,-24],[70,35],[67,-4],[62,21],[64,-12]],[[29382,19188],[-56,1],[-61,23],[-23,54],[-37,42],[-41,14],[0,32],[126,-36],[63,-37],[41,-60],[-12,-33]],[[28384,19223],[-69,2],[-39,37],[78,-8],[30,-31]],[[30884,19193],[84,-4],[26,-32],[-16,-36],[-49,-17],[-66,19],[21,70]],[[26873,19118],[22,-13],[49,-69],[4,-32],[-14,-80],[-107,-209],[-27,-9],[-76,77],[-30,15],[-20,41],[1,99],[11,50],[58,96],[41,31],[88,3]],[[30827,18849],[19,-3],[31,-56],[-22,-15],[-65,33],[-17,34],[5,32],[39,-1],[10,-24]],[[30881,18590],[-153,100],[-10,44],[82,18],[93,-12],[33,-35],[-8,-61],[-37,-54]],[[29955,18022],[-36,7],[7,54],[54,101],[15,14],[47,-25],[26,-45],[-7,-30],[-46,-43],[-60,-33]],[[30997,18086],[-9,-35],[-81,39],[-19,18],[-8,46],[52,7],[65,-75]],[[26732,17785],[-16,26],[22,24],[39,-8],[-45,-42]],[[29702,17416],[-16,-37],[-22,5],[-8,73],[37,13],[9,-54]],[[27039,16016],[-40,-101],[-25,-92],[-27,19],[36,108],[-22,20],[-50,-122],[-24,-30],[-30,0],[-43,-34],[72,144],[-18,14],[-59,-91],[-34,-35],[-26,31],[58,92],[52,66],[38,-17],[18,33],[24,94],[9,69],[40,-40],[13,-34],[31,-16],[15,-51],[-8,-27]],[[26756,15986],[-31,3],[23,53],[88,43],[-1,-39],[-79,-60]],[[27069,15955],[-8,43],[16,77],[35,28],[-1,-90],[-42,-58]],[[26550,14239],[-19,-7],[-57,14],[-93,46],[-132,53],[-46,27],[-9,20],[24,41],[28,28],[140,19],[54,-10],[64,-85],[37,-57],[15,-46],[-6,-43]],[[26917,13837],[-37,1],[-33,29],[84,44],[17,-13],[-31,-61]],[[23162,25050],[184,-50],[74,-45],[161,28],[148,-9],[125,-29],[75,-42],[-58,-92],[-100,-80],[-36,-57],[-63,-66],[-70,-25],[35,-20],[-53,-86],[-78,-89],[-58,-48],[-74,-18],[-259,42],[-65,-1],[-172,-22],[120,-44],[54,-56],[5,-69],[-86,-85],[-46,-113],[-28,-14],[-97,21],[-139,-15],[-50,7],[7,77],[-23,163],[-81,141],[-18,81],[7,117],[-15,194],[3,69],[67,30],[155,-46],[-66,72],[-25,64],[49,55],[133,36],[148,9],[42,-6],[102,30],[66,-9]],[[21772,24903],[-191,-29],[-112,34],[10,24],[93,19],[28,22],[87,13],[142,36],[99,2],[7,-28],[-51,-43],[-112,-50]],[[21301,24942],[230,-108],[274,43],[87,31],[68,4],[97,-22],[31,-27],[30,-62],[-42,-38],[-65,-14],[33,-90],[-142,-51],[-175,-144],[19,-33],[117,51],[58,1],[108,-50],[-13,-37],[80,-52],[-21,-60],[29,-3],[104,38],[35,-7],[28,-107],[-22,-47],[-75,-40],[55,-55],[-47,-81],[39,-39],[3,-68],[-91,-21],[-75,-60],[-98,-22],[-163,17],[-14,-52],[9,-60],[-58,-46],[-69,-24],[-88,34],[-49,-1],[-30,29],[-34,67],[-91,101],[-161,124],[-73,75],[-106,29],[-61,54],[-79,-19],[-61,18],[-68,87],[-117,55],[-70,63],[-15,32],[23,64],[34,48],[82,36],[77,-10],[45,-58],[58,-30],[52,-69],[72,-17],[164,31],[24,103],[46,-40],[36,33],[-51,79],[-68,9],[26,59],[118,-8],[-69,51],[-172,-33],[-153,79],[-20,28],[55,43],[94,14],[82,-53],[19,41],[-94,58],[10,58],[117,19],[83,-21],[-13,44],[62,29]],[[24954,23460],[11,-11],[115,30],[132,52],[150,0],[8,-37],[53,-39],[22,46],[31,175],[4,74],[-49,37],[-112,8],[-40,26],[-54,71],[-86,60],[100,59],[38,57],[105,1],[103,-14],[-55,60],[60,34],[-62,15],[-83,-11],[-50,19],[-81,107],[-2,71],[52,76],[52,15],[74,-17],[200,-62],[-5,20],[-223,83],[-79,22],[-19,30],[118,117],[110,27],[54,34],[90,1],[97,59],[87,28],[205,50],[194,0],[93,-17],[54,-31],[46,-79],[24,-83],[34,-35],[94,-40],[27,-51],[-3,-69],[87,-92],[-109,-95],[-24,-51],[43,-16],[-59,-83],[-33,-26],[68,-53],[-62,-16],[-8,-79],[40,7],[113,110],[74,33],[20,91],[72,23],[34,-54],[61,3],[29,27],[88,-60],[-2,-116],[39,29],[12,83],[33,44],[74,-36],[85,0],[106,-51],[52,20],[-205,71],[-45,38],[-7,37],[35,46],[95,44],[118,25],[184,-8],[53,-13],[139,-62],[132,0],[111,-46],[61,-80],[5,-48],[-93,-95],[134,31],[104,-7],[69,-25],[15,-50],[-28,-48],[-157,-59],[52,-26],[-46,-86],[30,-84],[44,125],[53,44],[72,5],[40,34],[49,11],[27,-26],[-95,-148],[95,81],[63,-55],[23,-53],[56,-12],[-32,68],[9,27],[102,71],[88,-36],[126,-7],[130,-39],[49,-26],[55,-52],[10,-62],[-75,-68],[-98,0],[-46,-20],[-74,-63],[-92,-55],[88,1],[82,91],[73,17],[101,-36],[50,1],[39,37],[69,27],[42,-49],[-5,-69],[-29,-40],[-71,-59],[-46,-23],[-107,-32],[-83,-68],[56,-7],[32,45],[52,21],[27,-24],[-10,-147],[91,188],[114,79],[166,54],[70,-30],[81,2],[62,-13],[76,-39],[107,-38],[36,-66],[-32,-33],[-85,-46],[-137,-24],[-109,-44],[-99,-60],[72,-12],[44,42],[41,15],[154,24],[53,-31],[-77,-84],[72,-19],[69,62],[50,12],[-3,37],[31,61],[39,18],[56,-18],[134,-125],[52,-119],[-24,-39],[-149,28],[-84,-12],[-70,-55],[-81,0],[-124,-35],[14,-29],[79,24],[74,6],[124,-51],[190,-2],[84,-26],[76,-42],[26,-40],[-1,-34],[-33,-19],[-144,16],[-42,-9],[-77,17],[-47,24],[-71,-23],[-56,15],[-58,-13],[-115,-55],[13,-9],[157,47],[31,-3],[150,-65],[0,-53],[-36,-75],[-120,30],[-59,-1],[-50,-16],[-135,31],[95,-57],[147,-13],[49,-26],[1,-33],[30,-33],[44,-14],[60,14],[55,-27],[94,-18],[70,5],[25,-24],[-64,-25],[-10,-25],[55,-55],[-20,-60],[48,16],[54,61],[86,15],[-17,-121],[14,-41],[65,55],[3,-56],[30,-17],[45,32],[57,-26],[-14,-86],[66,88],[88,53],[39,-27],[-52,-50],[1,-52],[52,-2],[52,28],[30,-14],[102,-96],[37,17],[46,-45],[-62,-36],[19,-75],[-81,3],[-14,-49],[142,4],[3,26],[63,31],[75,-20],[76,-51],[-34,-31],[-26,-76],[-121,-102],[202,64],[158,-22],[70,71],[43,-13],[123,-122],[59,-65],[-2,-21],[-60,-25],[-54,40],[-49,20],[-50,-26],[119,-62],[22,-66],[-80,-30],[-80,13],[-59,36],[-49,-2],[40,-46],[-34,-44],[79,-40],[68,-57],[-69,-15],[-103,5],[48,-36],[18,-56],[-62,-56],[-14,-42],[-84,-9],[-123,44],[-63,-17],[86,-29],[-8,-38],[1,-158],[-11,-56],[-56,-96],[-35,31],[-43,60],[-46,7],[-24,-28],[-53,53],[10,77],[-44,-36],[-23,-69],[-30,26],[-121,149],[-18,53],[-62,99],[7,29],[52,60],[91,38],[49,89],[5,48],[-96,-106],[-82,-35],[-64,-10],[-77,4],[-9,46],[-56,23],[-40,37],[-73,37],[-64,91],[-20,49],[-79,-16],[-66,-2],[35,-45],[0,-56],[-33,-9],[-68,49],[-88,34],[49,-100],[43,-41],[99,-78],[-24,-39],[-76,-18],[-76,22],[-35,25],[-52,59],[-85,60],[-79,-25],[77,-14],[68,-47],[-11,-71],[39,-56],[68,-27],[-24,-76],[13,-30],[52,31],[40,9],[46,-28],[26,-53],[51,-22],[-57,-42],[71,-57],[26,-80],[30,3],[35,-117],[26,71],[32,-14],[80,-60],[-21,-38],[41,-22],[33,57],[25,19],[46,-16],[70,-74],[22,10],[43,-48],[-65,-53],[91,-8],[28,-37],[-37,-68],[-81,8],[44,-43],[54,-74],[-6,-30],[44,-9],[90,9],[-2,-33],[30,-33],[45,-96],[-41,-13],[17,-110],[-4,-100],[-49,2],[-52,141],[-29,48],[-27,-19],[23,-137],[-17,-36],[108,-172],[-22,-18],[-46,1],[-66,22],[32,-116],[-36,-27],[-24,34],[-83,77],[-46,32],[-34,-5],[-31,40],[-67,53],[7,-50],[-59,10],[-50,97],[-10,-70],[-92,88],[-12,31],[-44,-6],[-86,89],[-58,77],[-31,1],[27,-107],[-31,13],[-111,95],[-68,45],[-100,13],[-14,-25],[66,-95],[81,-82],[63,-91],[62,-29],[65,-10],[-17,-43],[73,-31],[70,-53],[63,-72],[92,-50],[50,-70],[49,-81],[82,-47],[-13,-36],[-29,-18],[20,-95],[-54,-23],[-80,29],[-69,44],[-102,30],[-71,40],[-256,42],[-43,11],[-51,32],[-97,43],[-78,84],[-48,87],[-70,24],[-74,-17],[-45,3],[-91,53],[-63,21],[-80,55],[-43,12],[-67,42],[-32,31],[-97,122],[75,1],[66,41],[-16,38],[-57,60],[-31,9],[-105,-8],[-18,25],[17,67],[-61,-37],[-37,42],[-14,45],[-75,64],[-71,84],[-73,71],[45,69],[-95,21],[-28,-17],[-50,5],[-11,-57],[-40,-11],[-9,83],[-77,13],[-27,19],[-33,77],[-77,-29],[26,-25],[50,-80],[-15,-31],[-102,-20],[-71,17],[-43,25],[-62,-8],[-14,-67],[-96,-8],[-79,-37],[-42,0],[-44,-31],[-152,20],[-94,42],[-45,2],[-69,66],[-35,60],[-4,68],[25,94],[60,68],[117,45],[23,29],[-27,55],[37,63],[149,-20],[81,-22],[178,-72],[49,-44],[40,-64],[-24,-55],[37,-32],[26,85],[-40,63],[-73,62],[18,30],[106,-21],[41,4],[51,50],[86,-8],[75,16],[41,43],[98,18],[55,-17],[32,29],[-51,113],[-78,56],[-66,68],[-43,64],[16,35],[120,76],[95,75],[62,65],[43,32],[45,15],[24,78],[43,74],[115,52],[40,61],[-36,44],[-61,160],[-41,78],[-48,67],[-44,82],[-72,80],[6,51],[-92,-37],[-48,35],[-16,34],[26,75],[-19,62],[-72,0],[41,-62],[-6,-15],[-77,-4],[-41,20],[-62,64],[-8,42],[-63,16],[43,35],[-45,46],[-81,-37],[-84,22],[-27,-33],[-213,-99],[-59,15],[17,145],[48,22],[117,-11],[78,65],[5,27],[-38,54],[-71,34],[-98,27],[-31,38],[27,86],[-46,-15],[-41,-42],[-105,40],[71,47],[-26,40],[-74,15],[-71,-9],[-38,32],[-24,137],[-33,38],[-120,-9],[-92,63],[-76,78],[-22,37],[-49,-3],[-51,-47],[0,-56],[90,-23],[57,-74],[8,-66],[-30,-36],[-113,-42],[-85,0],[-149,51],[-255,49],[-125,10],[8,-32],[90,-33],[75,-67],[23,-40],[-23,-20],[-175,107],[-156,-52],[-95,13],[-165,69],[-119,-19],[-90,-1],[-180,21],[-67,15],[-39,22],[-142,17],[-95,-38],[-114,35],[-34,20],[-82,129],[-50,-6],[-86,12],[16,-44],[-47,-8],[-72,8],[-109,-37],[-136,96],[-71,14],[-51,27],[-98,119],[-44,119],[1,36],[93,-11],[90,0],[105,-31],[126,-20],[181,22],[11,12],[-117,49],[-82,51],[-51,9],[-144,7],[-134,24],[-167,68],[-42,35],[-48,230],[25,51],[73,55],[-43,25],[-11,103],[44,79],[93,119],[26,134],[78,93],[59,30],[15,46],[146,97],[123,64],[259,56],[98,11],[124,-3],[229,-21],[37,-57],[-62,-46],[-130,-72],[-111,-101],[-107,-152],[-40,-45],[-11,-55],[20,-51],[74,-101],[-5,-171],[33,-113],[50,-64],[136,-109],[125,-79],[15,-23],[-104,-64],[-147,-32],[-64,-34],[-112,-43]],[[26876,24795],[46,-7],[294,12],[61,-9],[233,-75],[55,-72],[67,-23],[52,-76],[76,-40],[-5,-42],[60,-59],[-85,-30],[-167,11],[-224,27],[-130,-7],[-224,-56],[-99,-8],[-119,54],[-32,58],[-35,125],[-119,23],[-56,55],[7,46],[-16,61],[2,65],[24,19],[95,4],[80,-30],[63,-3],[96,-23]],[[19860,24424],[-39,9],[-174,130],[-29,55],[-94,55],[-108,33],[24,61],[59,48],[301,36],[108,-11],[76,-49],[51,-17],[25,-37],[-7,-94],[-57,-94],[-48,-40],[-29,-46],[-59,-39]],[[22178,24432],[-44,-5],[-41,36],[7,67],[56,29],[72,-58],[-20,-57],[-30,-12]],[[21217,23200],[-85,26],[-16,50],[67,1],[31,-35],[3,-42]],[[24866,23007],[82,-4],[15,-56],[-48,-9],[-124,21],[-36,42],[58,23],[53,-17]],[[21999,22769],[55,16],[98,-83],[50,-20],[107,-63],[32,-44],[63,-118],[55,-64],[45,-32],[57,29],[29,-33],[-54,-40],[-59,-6],[-57,-55],[-36,-10],[-103,-69],[-53,-5],[-112,40],[-69,-6],[-57,8],[-113,65],[-95,34],[-23,52],[-60,-22],[-79,20],[-20,48],[-96,-35],[-50,27],[-34,59],[29,49],[149,18],[51,26],[73,59],[-28,64],[3,42],[43,10],[27,32],[6,64],[51,39],[79,15],[103,-63],[-7,-48]],[[26905,22842],[17,-38],[-51,-41],[-89,-11],[-20,-43],[-56,6],[-79,58],[-90,13],[35,38],[51,-7],[11,30],[57,3],[25,-33],[41,1],[103,36],[45,-12]],[[27286,22806],[17,-39],[-137,-69],[-87,11],[58,58],[48,6],[70,45],[31,-12]],[[30042,22718],[-79,20],[-5,18],[63,31],[64,-24],[-43,-45]],[[22523,22734],[37,-34],[-6,-43],[-26,-36],[-86,49],[-50,12],[33,69],[46,9],[52,-26]],[[20982,22645],[-32,22],[52,51],[27,-40],[-47,-33]],[[27568,22517],[-60,-3],[-45,69],[17,72],[36,12],[120,-18],[16,-39],[-84,-93]],[[23971,22656],[6,-31],[-51,-48],[-35,54],[47,30],[33,-5]],[[26965,22367],[-42,6],[-12,33],[45,64],[93,43],[26,50],[39,14],[-7,32],[25,23],[72,8],[28,-42],[-83,-80],[-17,-33],[-70,-83],[-97,-35]],[[23891,22556],[-22,-6],[-46,75],[39,5],[29,-74]],[[21303,22451],[-40,-22],[-28,36],[57,23],[11,-37]],[[21242,22347],[-40,-39],[-32,29],[-37,-2],[2,105],[53,19],[65,-63],[-11,-49]],[[20798,22236],[-114,39],[-1,22],[70,59],[50,-13],[29,-38],[-3,-36],[-31,-33]],[[28144,22116],[-117,63],[-24,26],[8,56],[22,26],[58,-2],[25,-18],[51,-96],[-23,-55]],[[20064,22145],[-43,6],[-72,37],[-24,34],[39,11],[84,-10],[35,-29],[-19,-49]],[[27027,22037],[-53,21],[6,51],[55,9],[34,-43],[-42,-38]],[[27927,22103],[143,-45],[24,-47],[-17,-89],[11,-108],[-11,-108],[-20,-39],[-54,-47],[-104,-42],[-73,-11],[-176,-13],[-85,16],[-41,71],[-41,150],[21,73],[78,122],[55,71],[103,44],[52,-5],[84,14],[51,-7]],[[24312,16311],[33,-19],[103,-103],[102,-72],[53,-63],[87,-144],[74,-42],[101,-30],[147,-71],[107,-39],[84,-48],[73,-86],[52,-41],[1,-63],[42,31],[97,-14],[51,0],[81,17],[84,11],[66,-25],[124,-25],[62,8],[32,-36],[79,-6],[51,-40],[23,-36],[24,-93],[-54,-231],[6,-88],[42,-69],[27,-115],[-13,-104],[12,-124],[-4,-51],[-27,-53],[-9,-47],[9,-53],[65,-76],[44,-80],[71,-111],[7,-33],[-24,-37],[53,-45],[49,-58],[43,-22],[44,-37],[84,-108],[44,-118],[15,-69],[-14,-22],[-49,-41],[-1,-37],[54,59],[30,11],[72,-23],[73,-56],[50,-83],[36,-45]],[[26880,13316],[1,-456],[0,-366],[1,-366],[0,-245],[0,-267],[-12,-35],[3,-38],[29,-91],[1,-47],[19,-57],[55,-98],[41,-111],[30,-61],[34,-47],[87,-44],[212,-60],[85,-63],[37,-80],[55,-49],[25,-3],[7,38],[25,2],[19,-42],[24,-78],[30,-51],[39,-23],[39,-6],[40,10],[33,-10],[26,-32],[76,15],[125,63],[83,34],[40,4],[38,-9],[49,-33],[-22,-126],[37,-52]],[[28291,10436],[-100,-107]],[[28191,10329],[-40,1],[-38,-18],[-50,-35],[-61,-64],[-108,-143]],[[27894,10070],[-70,-62],[-78,-37],[-40,-25],[-27,-3],[-54,-28],[-56,-55],[-18,66],[-16,-63],[4,-17],[58,11],[-37,-41],[-8,-52],[-26,-9],[-38,29],[-10,24],[-48,-4],[-22,38],[-27,3],[-12,-18],[-28,1],[-47,-16],[-50,-7],[-93,-30],[-39,-2],[-107,-34],[-66,-73],[-39,-22],[-33,-34],[-24,-64],[-32,-67],[13,-34],[54,-21],[52,-6],[75,43]],[[27005,9461],[-1,-90]],[[27004,9371],[7,-42],[25,-36],[3,-22],[-33,-23],[-89,10],[-57,-14],[-41,2],[-130,-42],[-62,-81],[94,-22],[3,-12],[-142,25],[-58,24],[-58,9],[-58,-4],[-41,-12],[-25,-20],[-64,-73],[-42,-87],[-43,-8],[-73,-50],[-41,-51],[-14,-63],[-25,47],[-23,5],[-39,-20],[-28,5],[-47,33],[1,73],[7,30],[29,26],[28,4],[76,-10],[44,26],[3,53],[-53,27],[-9,17]],[[26029,9095],[27,33],[16,58],[19,141]],[[26091,9327],[64,34],[41,51],[67,71],[20,74],[-2,256],[-3,40],[23,67],[14,60],[21,44],[32,31],[33,93],[0,73],[-37,100],[-48,90],[-25,19],[24,30],[81,-7],[-4,-47],[14,-64],[51,-58],[14,-56],[36,-21],[-10,-65],[37,34],[39,11],[15,-42],[36,-31],[64,-34],[38,-10],[21,41],[3,41],[-30,78],[30,25],[47,-45],[15,26],[-16,59],[-70,88],[17,13],[-24,62],[19,42],[-13,26],[-49,-17],[-31,50],[-1,58],[-34,8],[-39,101],[-44,82],[-76,-5],[-38,16],[-98,12],[-12,22],[6,40],[-41,-8],[-117,15],[-42,30],[-85,-9],[-15,16],[-84,-13],[-118,28],[-55,27],[-49,6],[-11,15],[-57,15],[1,73],[-21,28],[-73,-20],[-19,-20],[-18,48],[35,73],[-33,4],[2,38],[51,38],[-19,25],[-30,8],[-21,-12],[-35,30],[-5,28],[22,78],[26,65],[-33,54],[-60,57],[-21,50],[25,77],[18,80],[-30,5],[-177,-10],[-44,15],[-37,30],[-52,91],[-31,108],[-11,63],[-38,92],[-17,17],[-33,-6],[-17,21],[-48,-11],[-52,8],[-58,-6],[-50,30],[-152,62],[-16,20],[-66,-12],[11,-54],[32,-46],[-9,-47],[-25,-43],[-91,-35],[63,91],[2,34],[-19,35],[-29,-6],[-17,-98],[-38,-82],[-8,-47],[-42,-24],[31,125],[-85,-44],[-31,-23],[-6,-71],[-33,-96],[-26,-11],[-38,-49]],[[24140,11843],[-53,7],[-35,-10],[-25,10],[-27,52],[-62,-10],[-78,7],[-37,-4],[-27,48],[-55,-3],[-48,-45],[-46,-23],[-35,0],[-93,70],[-40,53],[-45,18],[-35,-5],[-31,-26],[-24,80],[-63,49],[-72,40],[-71,9],[-33,-7],[-23,-27],[-67,-19],[-26,12],[-13,30],[-55,26],[-98,23],[-56,19],[-23,33],[-2,28],[-25,71],[-11,58],[-3,94],[-23,22],[-59,11]],[[26880,13340],[-53,76],[-6,35],[13,81],[46,100]],[[26880,13632],[0,-292]],[[25943,10748],[-54,21],[9,32],[44,4],[32,23],[22,-67],[67,-5],[-16,24],[67,43],[26,-18],[26,-39],[24,44],[31,12],[12,-35],[36,-45],[47,-15],[-26,-68],[-19,-33],[-34,-28],[-40,15],[-80,60],[-51,22],[-62,35],[-54,7],[-7,11]],[[25838,10829],[21,-24],[-32,-39],[-26,19],[12,29],[25,15]],[[25706,10907],[-21,-26],[-19,7],[-33,90],[35,7],[40,-25],[-2,-53]],[[25222,11712],[-68,-12],[-20,24],[28,19],[51,-11],[9,-20]],[[24555,12265],[38,14],[41,-12],[-14,-37],[-70,-14],[5,49]],[[27775,9940],[-26,-12],[-25,-44],[-17,29],[15,24],[53,3]],[[31160,11069],[20,6],[15,41],[40,-11],[21,-18],[24,2],[23,-24],[42,-19],[45,2],[77,16],[141,14],[38,-33],[-80,-72],[-56,-38],[-7,-19],[5,-95],[-58,-5],[-36,17],[-5,62],[-16,36],[-17,-3],[-27,32],[-42,-27],[-2,-23],[-99,39],[-25,29],[-35,60],[-43,-1],[-25,10],[1,69],[-35,35],[-42,5],[-4,26],[14,39],[35,66],[63,81],[3,-75],[-29,-69],[51,-69],[10,-34],[-11,-32],[26,-20]],[[29731,11486],[1,67],[-52,18],[-106,-193],[-102,-188],[-17,-135],[-30,-50],[-19,-46],[-15,-98],[2,-76],[-30,-53],[-4,-32],[-48,-47],[-29,-47],[4,-62],[-30,-12],[-10,-48],[-17,-24],[-17,36],[-34,-20]],[[29178,10476],[-32,-17],[-34,15],[-25,-45],[-27,-98]],[[29060,10331],[-340,0],[-160,-1]],[[28560,10330],[-180,0],[-189,-1]],[[28191,10329],[39,20],[81,74],[59,26],[78,78],[56,15],[11,17],[15,88],[26,63],[31,54],[26,73],[46,48],[70,39],[65,86],[69,46],[15,35],[78,60],[62,11],[64,34],[48,18],[30,32],[43,17],[129,91],[36,42],[46,88],[41,44],[14,47],[59,77],[60,101],[30,73],[45,40],[87,115],[46,46],[19,5],[53,41],[33,43],[53,42],[95,53],[89,63],[121,55],[141,82],[115,44],[80,6],[98,20],[35,-2],[152,-35],[73,-44],[84,-92],[14,-59],[-44,17],[-12,-18],[46,-57],[-2,-71],[-26,-64],[-77,-32],[-36,-67],[-53,-34],[-21,-27],[-61,-43],[-59,5],[-76,41],[-47,39],[-42,-44],[-45,7],[-21,-9],[-34,11],[-70,-50]],[[32988,13581],[-55,17],[-44,-4],[-84,-22],[-69,-53],[-67,-13],[-65,0],[-29,-30],[-6,-33],[-114,-147],[-30,-50],[-58,-53],[-64,-93],[-55,-37],[-19,-51],[-53,-32],[-98,-7],[-46,-9],[-54,14],[-41,-22],[-61,-7],[-30,5],[-118,-50],[-30,47],[-23,18],[-67,4],[-54,19],[-96,12],[-114,-4],[-28,-26],[-95,8],[-40,24],[-77,-1],[-42,-17],[-93,20],[-97,-17],[-86,11],[-23,12],[-135,-31],[-53,17],[-46,-47],[-32,10],[-69,-4],[-34,-29],[-33,-45],[-54,-36],[-80,-198],[-7,-76],[-30,-52],[-49,-8],[-138,-38],[-61,-30],[16,-24],[-52,-26],[-59,-46],[-11,-35],[-71,-55],[-82,-129],[-39,-95],[-48,-68],[-34,-26],[-48,4],[-41,32],[-30,4],[-74,44],[-173,45],[49,-44],[92,-7],[97,-56],[76,-36],[25,-37],[-36,-135],[-24,-46],[-83,-120],[-40,-40],[-71,-143],[-111,-109],[-42,-65],[-97,-49],[-36,-13],[-33,7],[-41,-40],[-48,-24],[-14,-38],[-115,-100],[-45,-12],[-37,-27],[-12,-45],[-33,-27],[-38,-84],[-52,-82],[-65,-13],[-50,-75],[-38,-31],[17,-79],[-34,-11],[-66,-55]],[[26880,13316],[0,24]],[[26880,13632],[50,42],[20,-39],[30,-13],[31,-31],[39,-97],[19,80],[26,34],[-12,34],[-55,85],[-1,21],[25,35],[17,60],[34,35],[57,111],[-9,30],[-8,91],[-46,69],[-11,66],[12,49],[-36,60],[-12,34],[-26,186],[-33,155],[30,60],[-20,48],[22,36],[-19,25],[-48,23],[28,30],[-43,30],[-36,60],[-46,134],[-32,27],[14,22],[206,92],[165,95],[47,42],[65,43],[154,161],[44,54],[25,47],[85,121],[30,56],[28,127],[8,176],[-2,93],[-13,147],[-22,100],[-65,191],[-72,132],[-144,137],[-89,54],[-92,92],[-31,11],[-14,40],[9,75],[14,36],[51,68],[28,53],[42,53],[39,30],[23,38],[-27,48],[34,53],[2,48],[66,4],[41,57],[-11,65],[17,59],[-64,10],[25,42],[-53,110],[-9,32],[45,69],[-54,7],[-21,32],[36,75],[-58,-2],[-78,34],[90,117],[28,74],[-3,104],[61,60],[-51,36],[-67,68],[-36,82],[-16,93],[1,88],[18,37],[126,89],[63,21],[135,-24],[71,-30],[218,-76],[39,-34],[91,32],[177,-65],[12,14],[120,56],[120,74],[37,-7],[118,-89],[86,-77],[57,-6],[11,-44],[-26,-95],[60,42],[76,-46],[40,-52],[16,-58],[42,38],[62,-36],[-5,-33],[-50,-43],[27,-65],[52,-63],[35,-27],[67,-6],[124,-47],[49,-6],[72,13],[25,-14],[25,-81],[50,-10],[55,28],[23,68],[41,-20],[21,-65],[-10,-34],[-56,-63],[-22,-62],[-10,-90],[13,-52],[21,-33],[-12,-73],[-33,-23],[-45,-6],[-149,-2],[-30,-15],[222,-20],[47,-62],[14,-74],[-31,-95],[3,-74],[77,-2],[14,-30],[-28,-49],[10,-47],[-32,-110],[-39,-19],[-21,54],[-31,-42],[-80,-33],[77,-47],[62,16],[73,61],[57,24],[63,-4],[66,8],[37,-20],[41,-42],[24,-94],[20,-20],[5,-79],[-34,-120],[19,8],[31,113],[17,36],[36,-6],[34,-95],[14,46],[18,11],[13,-72],[33,-15],[51,43],[99,68],[81,29],[31,29],[25,75],[42,47],[52,-32],[34,-61],[32,2],[-26,63],[-5,43],[51,13],[-18,35],[21,31],[38,-3],[2,31],[85,14],[-60,26],[-1,51],[26,70],[6,49],[38,16],[-34,69],[12,81],[20,16],[54,11],[-3,39],[65,174],[31,40],[56,-13]],[[31707,12394],[-114,-7],[-91,31],[-67,15],[-66,28],[-143,88],[-30,68],[-57,62],[-151,86],[-12,30],[31,20],[65,7],[101,-33],[128,-30],[118,-56],[61,-41],[137,-109],[24,-9],[62,-54],[33,-73],[-29,-23]],[[31676,11481],[-25,-25],[22,96],[27,22],[54,63],[42,27],[0,-42],[-29,-2],[-83,-109],[-8,-30]],[[29194,11273],[-19,14],[59,48],[-9,-38],[-31,-24]],[[28502,10564],[-21,-10],[-36,10],[-21,24],[45,23],[45,71],[13,1],[-25,-119]],[[28467,10623],[-44,-6],[36,50],[41,11],[-33,-55]],[[20929,12344],[-349,0],[-232,0],[-146,0]],[[20202,12344],[-202,0],[-233,0],[-232,0],[-232,0],[-233,0],[-232,0],[-261,0]],[[10662,17902],[-35,41],[29,99],[0,32],[-163,-8],[-81,-73],[-131,59],[-19,-41],[-130,41],[0,443],[0,442],[0,443],[0,443],[0,295],[0,295],[0,442],[0,295],[0,443],[0,443],[0,295],[0,442]],[[10132,22773],[162,-24],[117,10],[217,-54],[134,-100],[108,-50],[45,-34],[236,-94],[52,-7],[96,-31],[74,3]],[[12972,15260],[-11,11],[-64,8],[-35,128],[-19,49],[18,43],[63,27],[-28,30],[-7,129],[-10,63],[-57,113],[-4,19],[-40,7],[-95,-40],[-40,-28],[-4,-47],[-19,-69],[-31,-57],[-37,18],[-29,76],[59,106],[45,116],[28,1],[51,23],[-80,11],[-22,22],[-26,54],[-44,21],[-33,50],[-14,53],[-33,12],[-3,59],[-41,15],[-45,50],[-10,33],[7,41],[-27,3],[-150,63],[8,83],[-28,110],[-30,44],[26,27],[119,-65],[-109,105],[5,76],[-17,0],[-26,-42],[-24,7],[-41,70],[-37,39],[-9,71],[-40,11],[-75,67],[-49,9],[-30,35],[-46,97],[-6,49],[-45,51],[-23,118],[-31,82],[-8,43],[-34,-8],[23,-51],[-28,-3],[32,-58],[14,-90],[34,-120],[15,-80],[28,-104],[-25,-37],[-44,11],[-40,61],[-33,18],[-82,-9],[-6,35],[9,57],[-50,107],[1,17],[59,38],[-81,-2],[-7,-53],[-21,-12],[-62,41],[-31,48],[-72,22],[-33,-9],[-29,-47],[49,4],[85,-48],[43,-43],[-20,-42],[70,-8],[33,-48],[13,-62],[-20,-17],[-91,-11],[-32,-42],[-32,10],[-102,65],[-128,94],[-6,23],[-50,41],[-31,40],[-45,78],[-88,74],[-57,21],[-39,30],[-174,95],[-102,64],[-23,22],[51,32],[29,59],[-20,61],[7,43],[29,12],[32,-40],[3,-70],[17,2],[21,49],[-61,81],[-31,8],[-90,-94],[-100,-52],[-56,-8],[-115,19],[-133,63],[-21,15],[31,39],[-9,52],[-55,-45],[-38,-12],[-118,32],[-122,26],[-108,6],[-152,-21],[-82,-22],[-93,2],[16,34],[-26,34],[-41,20],[-84,17],[-57,36],[-13,20],[20,100],[28,65],[-71,-63],[-31,-42],[-41,-37],[-36,4],[-50,26],[-91,19],[24,24],[37,63],[-129,5],[-28,45],[-69,-24],[-18,14],[48,41],[-38,14],[-29,30],[15,116],[-79,-62],[-29,-14],[-34,13],[-41,-12],[-18,-27],[-74,27],[-17,-35],[-64,-10],[-31,40],[44,110],[-56,-53],[-43,-26],[-37,-89],[-58,-26],[78,-52],[-7,-85],[67,18],[23,-41],[-22,-29],[-37,-18],[-9,-35],[5,-79],[-37,-23],[-27,-67],[-58,-16],[-54,-3],[-44,15],[-32,27],[-40,-18],[-40,51],[-5,-50],[-32,-54],[-37,1],[-30,-20],[4,-53],[-64,13],[14,-44],[-63,-62],[-39,7],[-40,-23],[-26,15],[-26,-70],[-48,-43],[-29,-49],[-61,29],[4,-33],[-54,5],[-69,-30],[-32,1],[-61,49],[31,61],[42,28],[80,28],[58,61],[38,68],[-11,8],[-86,-64],[-29,-5],[-69,24],[-24,42],[32,104],[77,108],[15,36],[23,97],[-12,98],[9,42],[101,49],[47,38],[92,54],[25,0],[40,-34],[56,-10],[39,14],[61,-6],[125,-36],[8,28],[-86,15],[-36,14],[-103,65],[-23,25],[35,21],[26,48],[108,97],[-71,-15],[-35,-29],[-34,-60],[-44,-13],[-99,-4],[-38,21],[-91,-52],[-56,-57],[-36,-22],[-85,-32],[-38,-34],[-14,-40],[10,-39],[-67,-37],[-75,-78],[3,-60],[-30,-36],[-77,-50],[-66,1],[63,-58],[10,-43],[-36,-60],[-95,-23],[18,-33],[-2,-41],[-51,-29],[-35,-6],[-14,41],[-64,-47],[11,-17],[-55,-71],[-75,-56],[6,-13],[-31,-92],[14,-18],[129,-41],[64,-39],[23,-52],[-28,-52],[-49,-50],[-66,-34],[-44,-48],[-17,-64],[-34,-39],[-10,-66],[-89,-20],[-3,-32],[-115,-20],[-26,-52],[-58,-53],[-75,-38],[-11,-30],[-40,-53],[-61,-12],[-15,-49],[-51,1],[-64,-61],[22,-44],[-21,-69],[-32,-20],[-10,-28],[-41,-2],[-98,-90],[-69,-9],[-27,-25],[-24,-63],[-79,5],[-43,-27],[13,-25],[-57,-32],[-62,-22],[-33,-50],[41,-18],[30,-53],[-67,-62],[-6,47],[-22,-6],[-22,-56],[-23,-28],[-174,-73],[-26,-17],[-12,-58],[-28,-25],[2,76],[-22,24],[-59,-24],[-77,-68],[-35,-13],[-25,-39],[-60,-9],[-24,-29],[-34,16],[-55,-55],[-77,-17],[-27,13],[1,34],[28,46],[-23,37],[-38,-18],[-23,-32],[-16,-70],[-65,-97],[-26,-30],[-25,2],[-45,-47],[-34,14],[8,34],[-31,49],[-28,-13],[8,-74],[-16,-37],[-52,-22],[-36,47],[-33,9],[-4,-75],[-42,-40],[-14,20],[11,37],[5,87],[45,36],[42,-5],[24,18],[23,39],[35,27],[40,51],[44,73],[53,63],[60,52],[65,42],[131,58],[19,-35],[49,9],[-10,-42],[42,-58],[28,0],[-3,42],[62,4],[50,-32],[10,33],[-46,36],[-16,35],[44,119],[20,37],[141,125],[137,64],[79,86],[26,-22],[60,-10],[-2,75],[6,48],[21,38],[72,92],[75,57],[41,51],[43,19],[55,-33],[-20,60],[-20,18],[-1,54],[17,76],[3,79],[15,45],[31,16],[68,11],[-80,30],[-10,85],[17,41],[63,70],[70,48],[-18,18],[21,109],[-49,-56],[-143,-65],[-97,-55],[-46,-13],[-30,14],[-55,105],[21,75],[56,20],[-55,26],[-25,-8],[-45,-73],[-22,11],[-27,-117],[24,-100],[-5,-40],[-44,-19],[-36,33],[-75,127],[-85,96],[-68,-46],[-63,43],[-57,74],[-80,-49],[-44,-42],[-29,0],[-80,-36],[-30,-29],[-9,-37],[-108,-29],[-106,16],[79,37],[36,39],[-16,52],[-2,60],[39,47],[-40,0],[-27,-17],[-24,36],[-12,69],[28,41],[24,76],[1,37],[-22,63],[-62,135],[-28,100],[-49,53],[36,87],[41,80],[-35,-10],[-55,-101],[-36,-49],[27,-86],[-19,-69],[-44,2],[-40,-36],[-93,-39],[-125,-22],[-62,2],[-64,46],[3,49],[-92,78],[-53,78],[-37,2],[-71,53],[9,44],[-23,13],[-66,8],[92,100],[32,68],[25,9],[119,-48],[28,-36],[-12,-60],[85,80],[28,-10],[45,-78],[55,37],[30,47],[-114,62],[50,-1],[-26,46],[-67,-53],[-122,3],[-85,31],[-85,-5],[-29,22],[81,61],[-55,4],[-85,60],[3,-54],[-51,-2],[-34,100],[-47,18],[-12,34],[30,45],[-45,30],[-36,-23],[-16,20],[15,50],[72,18],[-50,34],[-17,28],[127,33],[-25,30],[-10,42],[10,45],[70,103],[69,86],[53,30],[61,-27],[-18,51],[14,21],[-15,90],[23,86],[23,25],[83,17],[-41,37],[31,44],[82,23],[45,-7],[57,-27],[32,-34],[70,-41],[81,18],[31,17],[39,47],[50,30],[72,94],[21,39],[46,2],[38,-41],[127,8],[66,14],[45,31],[74,87],[13,45],[-34,107],[-23,111],[-63,73],[-45,22],[-8,44],[60,-5],[72,32],[26,52],[-14,57],[-48,55],[-34,10],[-76,-65],[-80,10],[-30,-37],[-82,-32],[-45,-33],[-82,-82],[-20,-37],[-26,-2],[-19,72],[-89,68],[-27,-23],[35,-37],[33,-6],[-25,-49],[-93,64],[-62,19],[-161,-2],[-162,-62],[-65,2],[-84,24],[-190,35],[-50,22],[-42,52],[20,50],[-2,50],[-37,13],[-75,73],[17,19],[63,10],[22,47],[69,29],[-112,24],[-15,-7],[-202,42],[-158,74],[-28,45],[47,12],[92,39],[47,51],[90,9],[89,88],[95,47],[50,13],[57,-25],[77,-4],[45,27],[-77,40],[18,37],[89,46],[105,14],[106,59],[58,17],[110,11],[66,-13],[11,-27],[-35,-77],[3,-46],[-38,-36],[92,-67],[142,-4],[78,12],[82,-24],[101,10],[77,-14],[32,5],[70,99],[28,16],[69,-31],[36,38],[-14,20],[-115,37],[-78,-19],[-24,21],[8,41],[-82,101],[-76,21],[-38,81],[67,26],[30,-14],[33,-59],[31,-9],[-9,-59],[38,-54],[87,-51],[70,19],[78,-12],[72,-45],[36,-6],[114,24],[-8,77],[-27,20],[-77,-4],[-60,34],[-51,-9],[-94,-51],[-48,20],[-77,55],[-6,52],[70,88],[-27,21],[-67,15],[-116,-15],[-101,8],[-65,-4],[-145,38],[-51,47],[-22,38],[-39,104],[-49,65],[-344,222],[-156,55],[-75,62],[-93,21],[2,21],[51,32],[27,75],[27,113],[-7,45],[191,-9],[126,7],[42,10],[160,18],[119,50],[90,68],[102,110],[19,112],[38,74],[164,170],[128,120],[66,-49],[175,35],[95,59],[4,13],[141,73],[41,-12],[-39,-48],[28,-13],[-24,-57],[65,-5],[11,87],[19,17],[-59,53],[76,77],[100,46],[87,-39],[102,-1],[37,21],[133,2],[107,49],[111,76],[113,114],[59,12],[26,-25],[180,-52],[45,-3],[18,-31],[-62,-64],[-93,-34],[75,-47],[80,30],[112,103],[61,-8],[73,-48],[-30,-47],[51,-23],[110,-24],[75,38],[115,7],[72,21],[122,-28],[80,2],[70,-35],[-55,-39],[-9,-40],[54,-20],[27,-29],[110,1],[-48,-54],[196,-17],[26,17],[127,29],[71,-33],[68,0],[77,33],[164,-5],[115,-38],[43,-4],[57,-50],[63,20],[103,-27],[45,-44],[175,-23],[86,10],[127,-2],[61,-16],[63,2],[105,-55],[66,-21],[157,-14],[55,29],[97,8],[86,24],[106,-6],[38,13],[139,-41],[78,-48],[34,-35],[163,-50],[79,-59],[111,-2]],[[1842,19726],[108,-15],[44,4],[55,38],[66,15],[84,-41],[29,-63],[107,-31],[19,-29],[34,-13],[118,0],[76,-19],[-12,-49],[-25,-21],[-70,7],[-70,-7],[-53,-57],[-3,-31],[-29,-21],[-27,70],[-58,40],[-54,11],[-23,45],[-69,48],[-88,33],[-58,0],[-92,-53],[-63,6],[-38,31],[92,102]],[[8244,18318],[-52,9],[6,36],[37,-8],[9,-37]],[[1485,18132],[59,-33],[38,3],[30,-28],[-98,-7],[-88,67],[-31,16],[7,38],[34,19],[17,-51],[32,-24]],[[7130,18083],[17,72],[19,-16],[-36,-56]],[[8318,18128],[0,-49],[-20,-66],[-38,4],[21,116],[37,-5]],[[8663,18128],[58,-11],[14,-18],[-133,-60],[-23,68],[39,36],[45,-15]],[[3285,18094],[25,-25],[55,7],[31,-17],[11,-45],[-7,-78],[23,-21],[10,-58],[-97,-12],[-41,-20],[-24,-38],[-42,30],[-78,15],[-97,61],[-42,12],[-43,44],[-38,56],[50,14],[113,-9],[14,40],[84,44],[31,-9],[62,9]],[[8297,17806],[-37,8],[16,37],[57,68],[38,31],[94,134],[44,-37],[-8,-18],[-90,-90],[-16,-43],[-33,-51],[-65,-39]],[[8219,17934],[-14,-16],[-54,10],[12,30],[41,19],[15,-43]],[[4706,17182],[-41,-4],[-17,50],[40,35],[74,29],[-56,-110]],[[7003,17136],[-41,28],[48,39],[29,-12],[-36,-55]],[[7022,17073],[28,27],[64,-51],[27,-35],[-35,-42],[-43,46],[-18,-22],[6,-37],[-62,-17],[-66,-43],[-55,-10],[-109,46],[73,76],[38,30],[55,-9],[1,33],[-19,36],[45,17],[36,-11],[34,-34]],[[11854,16973],[69,-12],[51,3],[46,-76],[45,-103],[-34,15],[-60,116],[-15,-8],[10,-75],[90,-153],[5,-46],[-15,-20],[14,-58],[-48,-18],[-44,-78],[-47,-45],[-43,18],[-8,52],[12,18],[24,95],[1,31],[-48,79],[-33,222],[-41,127],[26,-4],[43,-80]],[[11954,17007],[-39,-12],[-38,19],[-14,35],[35,10],[56,-52]],[[11568,17014],[39,-49],[-7,-78],[72,67],[94,-37],[20,-50],[-11,-68],[-36,-12],[-35,11],[-5,-42],[74,-4],[28,-68],[-16,-55],[-41,15],[-37,32],[-76,44],[-30,-5],[-3,-86],[-20,-31],[-59,14],[-24,44],[-22,71],[-104,100],[-30,50],[16,62],[39,24],[62,6],[13,22],[40,4],[43,32],[16,-13]],[[6891,16802],[13,-24],[42,28],[50,10],[26,-35],[-18,-24],[20,-45],[47,-15],[6,-20],[-54,-61],[-136,27],[36,-28],[28,-55],[-55,-12],[-47,-42],[-60,-6],[-63,-44],[-39,-55],[1,-35],[-31,-51],[-59,-43],[-26,24],[75,85],[-56,30],[-18,43],[-70,-1],[31,-35],[7,-29],[-21,-47],[-21,4],[-44,59],[-19,85],[-37,65],[9,56],[37,58],[69,39],[45,7],[24,-10],[61,-138],[-4,120],[34,34],[-50,47],[-7,32],[37,31],[29,-12],[37,-64],[46,37],[43,8],[11,77],[59,-18],[25,-20],[-13,-37]],[[6798,16815],[-8,-14],[-68,61],[47,-5],[29,-42]],[[11775,16564],[40,-99],[55,-221],[3,-58],[-13,-41],[8,-111],[-7,-38],[-26,7],[-27,42],[-36,98],[4,39],[-14,33],[-39,32],[-1,50],[-45,-1],[4,56],[31,48],[-38,27],[-9,54],[-35,28],[-29,-87],[-15,-19],[-41,-12],[13,46],[-15,65],[9,44],[29,8],[28,29],[25,64],[38,8],[103,-91]],[[2189,16479],[-29,-24],[-33,34],[74,19],[-12,-29]],[[6861,16449],[-40,-15],[-38,-48],[-22,26],[6,41],[19,27],[95,-9],[-20,-22]],[[12212,16388],[46,0],[34,-15],[32,-47],[0,-85],[-22,-46],[-26,31],[-31,57],[-18,-5],[36,-71],[9,-40],[-18,-51],[-107,0],[-12,26],[-11,98],[2,44],[-39,64],[-42,43],[31,30],[43,-3],[93,-30]],[[12042,16308],[18,-35],[25,3],[25,-66],[-31,-35],[-8,-49],[0,-97],[-18,-83],[-32,2],[-33,-28],[-16,64],[14,106],[30,22],[-56,62],[3,19],[-29,54],[3,58],[27,40],[36,7],[42,-44]],[[12381,16147],[-51,-1],[-1,42],[25,94],[51,-55],[24,-55],[-48,-25]],[[2299,16202],[76,-21],[-30,-26],[-49,33],[3,14]],[[6534,16141],[-34,12],[32,36],[27,-33],[-25,-15]],[[6405,16101],[-13,34],[58,49],[20,-13],[-65,-70]],[[12553,15936],[-5,-83],[-42,-7],[-37,25],[-8,39],[-57,11],[-10,60],[27,22],[60,130],[17,-6],[68,-122],[-13,-69]],[[12372,16006],[-31,6],[-39,41],[5,26],[31,31],[54,-2],[20,-31],[-10,-51],[-30,-20]],[[12157,16052],[52,-10],[47,1],[27,-43],[2,-73],[93,-48],[43,-51],[28,-49],[18,-56],[38,-65],[-2,-25],[-58,43],[-22,-65],[47,9],[55,-50],[16,-42],[-13,-38],[62,-8],[-6,-89],[6,-32],[-1,-83],[-11,-40],[-31,-7],[-36,46],[-20,52],[-57,23],[-10,51],[-53,-2],[-36,70],[-40,57],[-11,30],[39,18],[-34,55],[13,43],[-15,12],[-57,-3],[-19,41],[-50,1],[-35,48],[26,26],[38,-20],[34,24],[22,38],[-13,50],[-23,8],[-45,-18],[-40,-29],[-19,29],[61,74],[-17,35],[7,62]],[[12862,15623],[-28,-112],[-29,-31],[-35,32],[-28,4],[-15,49],[-45,-38],[-33,-67],[-14,29],[-9,97],[54,86],[6,124],[97,63],[41,-50],[42,-95],[-4,-91]],[[6164,15791],[-46,4],[31,42],[15,-46]],[[12228,15651],[7,-23],[-40,-34],[-10,-28],[-51,-54],[5,73],[-28,42],[28,22],[49,-9],[40,11]],[[4770,15535],[28,-3],[26,-62],[-33,-13],[-59,7],[-6,70],[15,36],[29,-35]],[[12349,15323],[25,16],[42,-17],[-8,-69],[-16,-36],[-28,12],[-55,74],[-31,58],[-35,104],[-45,18],[-3,47],[36,11],[54,-57],[28,-52],[3,-39],[33,-70]],[[4991,15441],[-45,-42],[-31,6],[-5,33],[36,36],[38,14],[7,-47]],[[12764,15417],[28,-66],[1,-23],[-47,-5],[-7,52],[-45,47],[5,56],[18,31],[29,-25],[18,-67]],[[4010,15367],[26,-84],[6,-40],[-62,-56],[-133,-2],[-44,-25],[-31,-45],[-31,-28],[-35,-11],[-64,7],[-21,63],[4,32],[49,43],[48,95],[29,16],[41,-6],[35,28],[76,42],[71,6],[36,-35]],[[4330,15299],[-25,13],[-12,30],[39,26],[16,-26],[-18,-43]],[[4261,15074],[-24,-11],[-47,33],[46,9],[25,-31]],[[3441,14940],[-14,2],[-11,57],[28,16],[23,-52],[-26,-23]],[[3365,14907],[-53,-12],[-18,34],[4,28],[53,19],[52,-44],[-38,-25]],[[3154,14821],[32,-8],[35,58],[38,-34],[-70,-74],[-17,-43],[29,-24],[-70,-58],[-19,-31],[-104,-35],[-62,-33],[-39,-33],[-40,-8],[-38,32],[27,24],[49,11],[107,67],[12,56],[21,31],[33,0],[32,18],[-68,18],[-21,22],[-4,30],[39,46],[66,20],[29,-4],[3,-48]],[[2787,14541],[-84,-54],[-64,-99],[-152,-105],[-2,15],[44,44],[31,47],[10,66],[33,41],[55,0],[11,77],[30,46],[58,29],[24,0],[49,-37],[-16,-49],[-27,-21]],[[2316,14290],[-8,-28],[-71,29],[62,17],[17,-18]],[[2033,14155],[-26,10],[39,49],[26,-28],[-39,-31]],[[958,13879],[-146,-20],[26,26],[123,45],[98,41],[13,27],[-48,24],[73,52],[41,-45],[-10,-43],[-34,-23],[16,-33],[-61,-29],[-91,-22]],[[1561,13999],[-42,0],[40,59],[43,-30],[-41,-29]],[[1264,13931],[88,-35],[-63,-13],[-57,10],[-33,21],[37,20],[28,-3]],[[592,13863],[-33,1],[-4,48],[34,-9],[3,-40]],[[436,13794],[43,-23],[0,-33],[-91,-63],[-18,23],[-34,-36],[24,94],[37,25],[11,74],[40,-22],[-12,-39]],[[285,13719],[-22,-12],[-42,6],[-53,-6],[-24,14],[90,28],[34,33],[12,34],[28,-4],[-19,-48],[-4,-45]],[[86,13685],[-13,-16],[-36,27],[20,47],[-57,59],[21,17],[45,2],[42,-40],[7,-31],[-29,-65]],[[595,13767],[-23,-11],[-30,22],[3,26],[50,-37]],[[24709,2945],[-36,-57],[-104,-17],[58,31],[-18,58],[-18,21],[0,57],[-7,33],[-23,27],[-12,-58],[-16,-101],[-72,3]],[[24461,2942],[-6,198],[-12,379],[-5,190],[33,466],[20,280],[33,466],[23,324],[-24,33]],[[24523,5278],[-5,13],[175,-3],[262,-4],[263,-5]],[[25218,5279],[43,-387],[36,-323],[36,-322],[12,-77],[38,-140],[12,-59],[1,-45],[11,-37],[-29,-51],[-24,-134],[-1,-53],[19,-103],[-10,-139],[7,-90],[16,-59]],[[25385,3260],[-264,0],[-176,0],[-264,0],[-2,-71],[9,-33],[40,-63],[-2,-97],[-17,-51]],[[24106,5784],[0,-47],[-16,-10],[-5,-41],[-31,-28],[-13,-36],[4,-45],[-11,-44],[-27,-45],[-11,-51],[4,-55],[-13,-46],[-42,-57]],[[23945,5279],[-7,-51],[-26,-24],[-16,-46],[-18,-17],[-12,-113],[-14,-45],[-22,-11],[-27,-71],[-25,-32],[-9,-30],[-1,-50],[-34,-28],[9,-37],[-4,-37],[-43,-77],[-2,-46],[21,-74],[-7,-72],[16,-39],[-12,-102]],[[23712,4277],[-298,-1],[-196,0],[-294,-1]],[[22924,4275],[-1,274],[-14,11],[-54,6],[-22,-11],[-29,41]],[[22804,4596],[4,275],[5,384],[3,219],[-30,352],[-19,210]],[[22767,6036],[152,0],[303,0],[304,0],[228,0],[227,0],[24,-65],[0,-38],[-20,-40],[-42,-59],[-21,-54],[184,4]],[[18836,3425],[-254,-1],[-289,-1],[-258,149],[-129,75],[-259,150],[-258,149],[-129,75],[13,28],[17,76]],[[17290,4125],[39,10],[18,19],[10,35],[-3,72],[-17,27],[-26,17],[-15,50],[-5,81],[5,44],[27,26],[13,33],[15,76],[-2,110],[24,75],[25,46],[55,73],[-9,34],[-40,39],[-18,28],[-7,52],[-43,103],[-15,57],[0,42]],[[17321,5274],[5,183],[-15,62],[1,78],[-9,50],[-3,95],[-15,48],[7,37],[34,31],[51,0],[34,-44],[27,-15],[17,23],[21,53],[0,416]],[[17476,6291],[171,-1],[170,0],[170,0],[255,0],[170,-1],[255,0],[170,0]],[[18837,6289],[0,-179],[0,-358],[0,-358],[0,-358],[0,-358],[0,-358],[0,-358],[0,-358],[-1,-179]],[[17290,4125],[-265,-38],[-234,-32],[-156,-22],[-2,59],[-29,7],[-8,72],[3,67],[-16,81],[-40,99],[-88,123],[-44,41],[-35,52],[-50,18],[-8,-23],[-32,16],[5,57],[-31,81],[-25,9],[-64,-6],[-85,44],[-25,27],[-9,47],[-39,41],[-53,40],[-29,-9],[-38,6],[-54,29],[-32,3],[-62,-8],[-23,6],[-22,36],[-23,19],[6,117],[-11,69],[8,64],[-20,41],[-41,27],[-7,33],[7,45],[-11,30],[-34,28],[-31,64],[-40,35],[-16,59],[-25,36],[-8,32],[-54,114],[-58,90],[-9,51],[-3,71],[35,80],[-4,60],[-20,45],[-78,26],[-62,109],[-4,84],[-25,85],[-4,116],[35,9],[4,-67],[39,-47],[28,-10],[-26,95],[-20,30],[-11,53],[-14,32],[46,41],[36,4],[100,-7],[-8,23],[-36,-3],[-41,25],[-34,-30],[-65,40],[-25,-18],[-3,-78],[8,-58],[-80,54],[-57,76],[-5,91],[-16,14],[-20,73],[-46,44],[-37,70],[-75,117],[-5,103],[-28,130],[12,74],[-2,53],[-13,79],[-14,43],[-61,118],[-59,79],[-13,121],[24,111],[17,32],[24,99],[-2,94],[20,115],[-14,120],[-12,49],[-22,35],[9,51],[-5,56]],[[14701,8813],[218,0],[287,0],[144,0],[287,0],[216,0]],[[15853,8813],[0,-284],[0,-379],[0,-189],[0,-379],[0,-284],[128,-165],[127,-165],[128,-165],[127,-165],[150,-209],[150,-208],[100,-140],[138,-200],[139,-200],[138,-200],[143,-207]],[[15885,4814],[55,-26],[30,13],[2,-24],[-67,-19],[-21,13],[1,43]],[[15841,4733],[-34,-1],[-23,49],[49,6],[24,-27],[-16,-27]],[[16303,4464],[14,-37],[-41,2],[-12,48],[39,-13]],[[16302,4182],[-33,6],[-9,61],[42,-67]],[[18837,6289],[0,252],[0,253],[0,252],[0,253],[0,252],[0,379],[0,378]],[[18837,8308],[171,0],[171,0],[171,0],[171,0],[171,0],[257,0],[257,0]],[[20206,8308],[272,0],[272,0],[0,-504]],[[20750,7804],[1,-474],[1,-378],[0,-379],[1,-284]],[[20753,6289],[-269,0]],[[20484,6289],[-309,0],[-206,0],[-206,0],[-206,0],[-308,0],[-206,0],[-206,0]],[[28983,8819],[1,-248],[-2,-53],[-10,-41]],[[28972,8477],[-24,3],[-91,-25],[-29,11],[-30,-19],[-100,-5],[-21,10],[-27,-35],[-43,-20],[-109,-78],[-13,-15]],[[28485,8304],[-26,57],[65,58],[-16,38],[17,384]],[[28525,8841],[184,-11],[274,-6],[0,-5]],[[27558,7242],[-25,28]],[[27533,7270],[22,34],[30,-62],[-27,0]],[[28101,7023],[-182,0],[-13,400],[-8,240]],[[27898,7663],[20,41],[21,18],[34,2],[24,-14]],[[27997,7710],[-46,-88],[4,-83],[44,-99],[6,-95],[22,-64],[34,-74],[26,-21],[-10,-73],[15,-27],[9,-63]],[[25385,3260],[33,-122],[1,-20],[199,-21],[174,-19],[260,-29],[87,-9],[13,-87],[22,-11],[14,17],[11,78],[-5,87],[28,44],[72,-40],[46,-4]],[[26340,3124],[13,-46],[19,-187],[13,-65],[24,-176],[40,-170],[55,-205],[92,-250],[11,-35],[-27,-84],[-30,72],[12,59],[-21,27],[-2,50],[-18,13],[25,-190],[26,-101],[116,-492],[27,-62],[10,-45],[11,-94],[2,-121],[-19,-221],[-4,-150],[-9,22],[-16,-69],[-22,-62],[-8,-96],[-10,-49],[-33,-51],[-19,1],[-49,-38],[-34,10],[-41,-22],[-27,3],[-16,45],[19,46],[37,-48],[-5,44],[-36,28],[-31,109],[-32,75],[-5,51],[-56,30],[-40,46],[-26,83],[-15,145],[-25,28],[-3,32],[-19,0],[-11,76],[3,95],[-7,36],[-24,-13],[0,-49],[-29,15],[-42,96],[-48,172],[-26,50],[21,13],[32,77],[24,47],[9,32],[-13,35],[-20,-13],[-16,41],[-26,2],[21,-45],[-4,-49],[-13,-30],[-23,-4],[-27,69],[26,197],[24,126],[2,203],[-32,82],[-142,203],[-110,239],[-95,90],[-73,-20],[-12,-18],[-7,-61],[-46,-5],[-68,-63],[-24,2],[-39,-28],[-42,-7],[-36,-14],[-16,7],[19,52],[-13,40],[-41,50],[-47,74],[4,34],[16,49],[-37,-21],[-4,-39],[-28,24],[-87,59],[-76,34],[59,15],[-5,32],[-52,2],[-63,-47],[-142,-32],[21,29],[37,17],[6,36],[-9,34],[-20,-35],[-27,19],[-4,-43],[-26,-57],[-53,-23],[-4,42]],[[26699,1380],[-21,49],[-31,136],[6,17],[46,-202]],[[26645,302],[9,78],[25,25],[-34,-103]],[[25218,5279],[176,-3],[177,-3]],[[25571,5273],[165,3],[165,3]],[[25901,5279],[-15,-34],[-40,-64],[-11,-50],[52,-54],[31,-45],[56,-38],[8,-35],[62,-176],[65,-91],[26,-45],[13,-43],[55,-71],[19,-38],[1,-29],[30,-59],[12,-39],[22,-35],[33,-29],[24,-68],[15,-104],[15,-61],[24,-26],[32,-89],[11,-53],[-1,-47],[17,-37],[55,-40]],[[26512,3779],[-14,-43],[-52,-33],[13,-33],[-9,-30],[-27,-25],[8,-48],[-2,-40],[-32,-80],[-24,-9],[25,-46],[-20,-46],[-21,14],[-3,-56],[-19,-108],[5,-72]],[[23633,7995],[-29,45],[-37,42],[-13,37],[-206,-3],[-138,-1],[-206,-2],[-275,-3],[-274,-3]],[[22455,8107],[-26,72],[7,73],[-3,97],[-10,104],[-18,41],[4,29],[-16,48],[-18,13],[-3,95],[-23,116],[-39,85],[-16,64],[-8,80],[-26,47]],[[22260,9071],[-16,74],[-22,36],[43,187],[-8,54],[-17,13],[3,71],[-15,63],[39,2]],[[22267,9571],[177,0],[266,0],[266,0],[177,0],[177,0],[356,0]],[[23686,9571],[7,-54],[28,-32],[-12,-116],[12,-111],[13,-47],[30,-34],[46,-22],[28,-34],[10,-49]],[[23848,9072],[17,-35],[28,-31],[38,-85],[49,-57],[4,-60],[-18,-85],[-28,-53],[-7,-39],[-25,-36],[-40,-33],[-46,-23],[-50,-12],[-30,-29],[-9,-44],[4,-38],[18,-30],[12,-66],[-12,-56],[-28,-78],[-31,-51],[-33,-23],[-15,-34],[0,-70],[-13,-9]],[[16929,12344],[0,-499],[51,-82],[35,-69],[18,-105],[-9,-31],[59,-63],[37,-26],[113,-171],[13,-39],[36,-56],[21,8],[21,-45],[22,-7],[46,10],[1,-59],[-38,-198],[2,-62],[17,-27],[3,-54],[-24,-18],[-11,-29],[10,-48],[-9,-56],[63,-44],[97,98],[38,-41],[1,-31],[28,-118],[38,-84],[24,-34],[8,-40],[-8,-32],[36,-61],[41,-5],[25,-42],[26,-119],[24,-31],[26,-6],[27,46],[88,-15],[20,33],[24,14],[43,-7],[50,4],[26,-8],[37,10],[44,-7],[14,84],[30,17],[19,-29],[19,-49],[40,-46]],[[18291,10075],[0,-395],[0,-473],[0,-394]],[[18291,8813],[-255,0],[-153,0],[-203,0],[-204,0]],[[17476,8813],[-304,0],[-304,0],[-202,0]],[[16666,8813],[-1,338],[0,226],[0,351],[17,77],[2,78],[10,40],[-32,39],[-30,7],[-16,27],[-1,47],[23,72],[46,96],[25,69],[5,40],[28,93],[51,145],[20,93],[-12,41],[-26,39],[-41,37],[-36,65]],[[16698,10833],[-14,62],[2,55],[-26,83],[3,16],[-1,324],[0,324],[-2,404],[0,243]],[[23848,9072],[248,-3],[245,-2],[280,-3]],[[24621,9064],[4,-144],[28,-65],[26,-113],[22,-67]],[[24701,8675],[-1,-286],[-1,-302],[-1,-303],[-1,-296],[-21,-43],[2,-29],[-12,-33],[32,-131],[0,-84],[-14,-32],[-15,-69],[-40,-62],[-12,-42],[-41,-44],[8,-36],[-15,-32],[-11,-51],[1,-39],[-15,-30],[15,-44]],[[24559,6687],[-24,-45],[-3,-40],[9,-44],[-23,-35],[-54,-27],[-29,-24],[1,-51],[11,-43],[-5,-45],[-39,11],[-71,38],[-47,-1],[-24,-41],[-5,-55]],[[24256,6285],[-32,44],[-8,-34],[-31,45],[-20,56],[-4,44],[11,68],[-11,36],[-6,70],[-27,52],[-54,69],[-40,42],[-26,17],[-30,34],[-34,53],[-18,41],[-3,28],[16,75],[52,189],[-26,33],[-70,25],[-27,-19],[-23,39],[-20,98],[-50,104],[-79,109],[-46,74],[-24,97],[-10,78],[1,65],[16,78]],[[24701,8675],[29,-36],[44,-9],[49,16],[57,46]],[[24880,8692],[233,0],[330,0],[0,-30]],[[25443,8662],[-2,-328],[-3,-493],[-2,-496]],[[25436,7345],[-9,-91],[18,-22],[-3,-25],[-28,-25],[-73,-40],[-37,24],[-27,-4],[-11,-24],[6,-43],[-10,-34],[-23,-23],[-28,-69],[-38,-32],[-18,-35],[-14,-63],[-15,-39],[-39,-11],[-29,23],[-20,47],[-40,5],[-25,-63],[-7,-50],[-14,-21],[-20,7],[-27,36],[-27,8],[-27,-19],[-47,-67],[-18,25],[-45,35],[-40,-1],[-22,12],[-4,-49],[-27,26],[-52,1],[-2,-43],[-35,-14]],[[20750,7804],[171,0],[286,0],[171,0],[285,0],[114,0],[229,0],[171,0],[229,0],[163,0]],[[22569,7804],[63,-65],[36,9],[14,-15],[10,-40],[-30,-46],[-19,-47],[5,-45],[28,-41],[19,-43],[10,-42],[19,-33],[42,-35],[0,-217],[1,-394],[0,-461]],[[22767,6289],[-314,0],[-252,0],[-252,0],[-126,0],[-252,0],[-252,0],[-125,0],[-252,0],[-189,0]],[[25436,7345],[31,18],[31,-21],[45,2],[38,-59],[17,-55],[31,-34],[45,-12],[33,-21],[37,-42],[42,26],[49,-37],[42,6],[55,38],[33,-2],[10,-42],[23,-43],[40,-47]],[[26038,7020],[6,-14],[3,-75],[-11,-46],[8,-41],[29,-59],[3,-23],[52,-127],[48,-58],[38,-16]],[[26214,6561],[-103,-123],[-92,-80],[-14,-45],[-27,-35],[-17,-37],[-41,-30],[-29,-43],[-72,-41],[-41,-14],[-28,-24]],[[25750,6089],[-17,-9],[-189,10],[-189,10],[-189,9],[-189,10],[-182,-5],[-181,-5],[-8,15],[-62,9],[11,-51],[-1,-46],[-184,0],[-195,0]],[[24175,6036],[14,44],[19,15],[20,-11],[14,15],[13,46],[9,99],[-8,41]],[[24165,6038],[-19,-1]],[[24146,6037],[19,1]],[[23712,4277],[15,-39],[-11,-66],[11,-57],[-15,-35],[25,-102],[31,-47],[-1,-46],[-38,-68],[-1,-30],[-18,-38],[-35,-43],[-24,-56],[-11,-68],[-31,-96],[0,-53],[-28,-64],[9,-54],[-16,-55],[196,0],[197,0],[130,0],[-8,-58],[-19,-79],[4,-60],[37,-75],[26,-121],[19,-15]],[[24156,2852],[-19,-14],[-99,52],[-25,42],[-49,14],[-29,-52],[-22,-69],[35,-38],[30,-18],[49,15],[27,34],[23,-1],[20,24],[20,-28],[-41,-55],[19,-39],[42,-8],[7,44],[19,29],[22,-24],[15,-46],[1,-50],[-48,-25],[-36,-44],[-16,-33],[13,-41],[26,-26],[18,-33],[73,-44],[18,1],[17,-44],[27,-23],[-1,-30],[-24,-23],[-13,-41],[-22,32],[-25,-41],[-7,36],[-24,62],[-48,55],[-47,16],[-7,42],[-16,21],[-77,40],[5,-29],[24,-25],[0,-47],[-13,-79],[-31,-39],[-9,14],[-15,63],[-21,19],[-34,3],[-22,-14],[-25,-62],[-20,-9],[-69,31],[-78,48],[37,32],[-25,53],[-9,54],[-14,-25],[-50,21],[-43,96],[-42,2],[-18,44],[-34,-18],[-32,-51],[20,-42],[-6,-12],[-49,-18],[-111,20],[-33,19],[-44,40],[-61,33],[-141,-5],[-36,-22],[-16,42],[21,21],[10,32],[-6,32]],[[22992,2743],[23,68],[-3,86],[-8,31],[9,65],[-1,43],[26,79],[12,55],[14,113],[-8,67],[-29,72],[0,25],[-20,72],[-22,45],[-2,74],[-21,64],[-38,62],[0,511]],[[23537,2502],[-10,-7],[-45,44],[-3,19],[36,15],[29,-28],[-7,-43]],[[28525,8841],[-7,12],[69,340]],[[28587,9193],[215,-11]],[[28802,9182],[309,-14],[24,13],[28,40],[45,32],[46,3]],[[29254,9256],[-6,-26],[25,-79],[20,-27],[-46,-32],[-27,-61],[-31,-51],[13,-16],[49,-17],[22,-19],[32,-95],[-10,-27],[29,-24],[9,-68],[24,-24],[36,-14],[44,21],[36,28],[-1,23],[-23,55],[31,-10],[10,-77],[-2,-66],[-31,0],[-94,-25],[-21,-23],[-47,-24],[-3,89],[-9,2],[-75,-84],[-58,-16],[-12,96]],[[29138,8665],[-30,46],[-10,53],[-2,57],[-113,-2]],[[29335,8498],[-75,-25],[30,61],[16,5],[29,-41]],[[29480,8442],[-21,-8],[-49,19],[47,21],[23,-32]],[[27533,7270],[-18,13],[-30,42],[-49,30],[0,55],[-15,22],[-52,41]],[[27369,7473],[-20,58],[0,28],[-23,48],[-39,10],[-19,23],[-23,4],[-61,-29],[-25,-47],[-49,8],[-33,32],[-5,-22],[-42,-59],[-29,12],[-59,-83],[-53,-51],[3,258]],[[26892,7663],[251,0],[189,0],[189,0],[188,0],[189,0]],[[28101,7023],[-26,-36],[-7,-65],[-18,-6],[-41,-110]],[[28009,6806],[-66,-13],[-11,-23]],[[27932,6770],[-52,9],[15,58],[-36,42],[17,46],[-7,48],[-21,-33],[-23,-5],[-44,41],[-22,67],[8,53],[57,11],[-67,99],[36,18],[-19,46],[-25,-18],[5,51],[35,-9],[-14,101],[22,62],[54,73],[-1,64],[-28,-12],[-9,-65],[-33,-27],[-16,-29],[-35,-9],[-5,-40],[-41,15],[42,-91],[-15,-22],[-12,-67],[-5,-79],[10,-106],[23,-56],[-6,-30],[8,-47],[20,-65],[-71,45],[-47,17],[-27,64],[-6,-49],[-27,51],[-24,24],[-22,-22],[-21,5],[-2,44],[29,79],[22,28],[10,36],[-4,56]],[[28050,6814],[-7,-2]],[[28043,6812],[7,2]],[[30257,10413],[12,-91],[-7,-53],[32,-28],[-55,-77],[-47,10],[-53,-26],[-26,-41],[-51,7],[-18,-62],[-41,-35],[-16,53],[-44,9],[-26,-32],[-21,32],[-20,-64],[-3,-62],[-21,26],[-50,36],[23,31],[-37,22],[-31,-26],[1,-44],[-29,-88],[-2,-38],[-43,-56],[-32,7],[-37,-48],[-39,-13],[-20,38],[-9,-52],[-17,-24],[-47,25],[-26,-10],[-32,-34],[-25,-48],[18,-23],[-86,-140],[-33,-108],[-25,-33]],[[29274,9353],[-22,47],[-4,38],[-35,76],[-12,352],[-10,280],[-13,330]],[[29967,9991],[-61,-20],[1,36],[30,46],[17,-9],[13,-53]],[[26029,9095],[3,52],[-14,18],[-28,-17],[-2,-38],[-20,-35],[-16,-73],[-48,-45],[-23,-102],[3,-17],[-23,-52],[-27,-34],[-25,-68]],[[25809,8684],[-3,-3],[-182,-9],[-181,-10]],[[24880,8692],[47,47],[16,29],[30,86],[31,72],[29,98],[19,117],[3,127],[-8,64],[-22,99],[-28,92],[-30,117],[28,79],[-2,65],[-19,76],[28,48],[32,82],[15,111],[-5,59],[6,26],[31,14],[10,24],[8,61],[29,26],[26,-3],[20,26],[26,58],[27,17],[-4,-97],[-12,-67],[26,-11],[18,13],[26,97],[-1,82],[5,41],[42,40],[69,19],[-32,57],[0,50],[30,51],[7,35],[49,10],[90,-60],[31,-5],[37,-28],[17,-40],[39,-14],[32,-24],[106,-59],[22,-33],[0,-27],[21,-41],[-35,-58],[3,-32],[30,-46],[9,-94],[-7,-64],[-3,-92],[-36,-56],[-13,3],[-28,-107],[-21,-22],[-51,-26],[-12,-44],[-4,-70],[24,-37],[31,-19],[31,5],[13,31],[29,17],[31,75],[4,25],[35,32],[72,37],[43,-32],[34,-75],[18,-111],[22,-179],[15,-100],[12,-32]],[[24678,10382],[-13,7],[-19,61],[16,52],[-9,19],[-46,-8],[-2,23],[16,47],[5,45],[-4,46],[-31,37],[-52,28],[5,32],[-15,34],[-93,34],[-94,12],[-77,57],[-272,97],[-30,85],[-45,20],[-4,22]],[[23914,11132],[70,31],[43,32],[36,38],[45,23],[93,22],[55,54],[56,27],[32,52],[49,43],[5,-32],[32,-36],[16,-67],[-6,-64],[71,58],[120,-32],[38,-31],[35,-73],[31,-47],[10,-33],[51,-5],[37,15],[50,-37],[32,10],[26,-19],[29,36],[82,67],[78,23],[90,-6],[47,8],[44,23],[72,16],[5,-22],[-10,-118],[70,-19],[22,11],[46,-27],[28,27],[27,1],[26,-64],[7,-60],[-6,-32],[32,-3],[20,-20],[-1,-25],[40,-44],[-23,-20],[-99,19],[-49,-2],[-30,29],[-26,-48],[0,-49],[-21,7],[-46,58],[-35,29],[-43,21],[-60,10],[-37,-44],[-25,-17],[-35,-4],[-34,-20],[-44,15],[-41,-17],[-35,-68],[-35,-24],[-32,-65],[-12,39],[36,53],[-1,38],[-65,-26],[-10,-33],[-37,-21],[-5,60],[-13,1],[-13,-58],[-35,-62],[-50,-122],[-55,-105],[-1,-25]],[[24627,11545],[-62,-10],[1,-36],[-56,-55],[-4,-26],[-55,-75],[1,56],[-50,26],[15,50],[46,59],[77,42],[63,6],[23,-9],[1,-28]],[[24322,11810],[-16,-28],[-51,-20],[-11,16],[14,40],[106,76],[46,28],[-10,-59],[-78,-53]],[[25253,10711],[0,-49],[-32,-27],[10,65],[22,11]],[[25788,10844],[3,-33],[-58,-7],[-37,18],[71,59],[21,-37]],[[25624,11004],[-16,10],[-23,65],[36,21],[-6,-61],[9,-35]],[[24140,11843],[-15,-4],[-53,-47],[-100,-58],[-109,-49],[-68,-52],[-30,-31],[-117,-156],[-71,-83],[-114,-117],[-12,-29]],[[23451,11217],[-31,-46],[-19,-5],[0,-291],[-110,-89],[-49,-102],[-5,-56],[37,-26],[14,-25],[6,-57],[-22,-73],[2,-64],[-10,-19],[9,-70],[-10,-72],[8,-33],[34,-41],[76,-47],[14,-27],[71,-49],[38,-57],[20,-47],[34,-35],[55,-36],[36,-32],[17,-28],[11,-70],[9,-149]],[[22267,9571],[0,340],[0,227],[0,340],[-20,37],[-46,31],[-40,82],[8,32],[54,71],[16,73]],[[22239,10804],[-1,98],[-15,96],[-31,73],[-15,89],[-3,75],[7,62],[-13,27],[-11,310],[-35,117],[-6,54],[-19,59],[-13,66],[1,57],[-6,126],[5,73],[-27,158]],[[22767,6036],[0,253]],[[22569,7804],[-24,48],[-8,51],[-16,34],[-23,15],[-15,34],[-6,53],[-22,68]],[[24175,6036],[-10,2]],[[24146,6037],[1,-65],[-14,-11],[3,-46],[-22,-20],[11,-27],[-2,-32],[-17,-52]],[[24461,2942],[-79,-8],[-35,26],[-23,4],[-41,-24],[-46,-18],[-27,7],[-33,-62],[-21,-15]],[[23945,5279],[287,0],[291,-1]],[[20202,12344],[2,-385],[1,-289],[1,-289],[1,-289],[2,-288]],[[20209,10804],[1,-476],[-9,0]],[[20201,10328],[-238,0],[-239,0],[-239,0],[-179,0],[-179,0],[-239,0],[-299,0],[-298,0],[0,-253]],[[25571,5273],[11,107],[15,30],[46,13],[16,65],[36,51],[40,22],[58,13],[62,49],[44,47],[36,13],[28,76],[20,-2],[31,42],[24,9],[10,-38],[15,4],[41,60],[51,27],[22,-24],[14,9],[39,92],[27,24],[20,1],[8,90],[13,39]],[[26298,6092],[4,-8],[296,-4],[296,-4],[295,-5],[197,-3],[198,-2],[264,-4]],[[27848,6062],[-7,-39],[34,-104],[-4,-47],[-72,53],[-1,-68],[-21,-15],[-42,9],[-26,-53],[-22,-7],[-33,31],[-13,-60],[32,-7],[69,5],[41,19],[37,-10],[3,-47],[-7,-95],[23,16],[6,88],[34,32],[22,-30],[8,-68],[-8,-60],[-53,-70],[-37,-64],[-19,-13],[-59,23],[-27,-2],[-12,56],[-21,11],[-7,-38],[-29,-11],[-40,16],[-23,-3],[103,-65],[22,-30],[-25,-59],[-6,-40],[-41,-42],[9,-25],[79,24],[25,-26],[-42,-81],[-51,-12],[-8,-23],[-86,-5],[-23,5],[-32,-46],[-29,2],[-15,-14],[9,-33],[-37,-38],[-64,-84],[-48,-107],[-23,-82],[-107,3],[-43,-20]],[[27141,4712],[-113,176],[-184,288],[-308,16],[-3,68],[-38,76],[-30,-24],[-12,45],[-173,12],[-172,12],[-115,-59],[-92,-43]],[[27869,6062],[9,0]],[[27878,6062],[27,-163],[54,-176],[-6,-4],[-40,117],[-22,84],[-22,142]],[[22239,10804],[-127,0],[-253,0],[-127,0],[-254,0],[-127,0],[-317,0],[-191,0],[-254,0],[-254,0],[-126,0]],[[20206,8308],[-1,379],[-1,253],[0,378]],[[20204,9318],[188,0],[189,0],[188,0],[189,0],[189,0],[188,0],[189,0],[202,0],[81,-75],[49,-23],[22,22],[50,12],[77,4],[48,-11],[60,-50],[59,-30],[35,-32],[12,-35],[41,-29]],[[29274,9353],[-20,-97]],[[28802,9182],[-20,50],[-3,29],[21,76],[18,148],[6,100],[24,94],[20,38],[29,88],[15,76],[8,93],[45,28],[50,48],[20,32],[11,47],[-14,83],[30,69],[-2,50]],[[28091,7795],[45,48],[48,38],[-66,127],[-16,7],[-17,62],[-25,27],[4,91],[17,13],[10,42],[-13,58],[25,26],[32,52],[20,58],[39,44]],[[28194,8488],[214,-184]],[[28408,8304],[-4,-39],[-27,-80],[-55,-75],[-10,-40],[6,-37],[52,-13],[14,11],[12,-62],[-20,-129],[-24,-68],[-16,-107],[-22,-57],[-36,-64],[-11,-50],[-47,-48],[-12,-43],[-40,-104],[-44,-26],[16,99],[-42,33],[-24,-2],[-25,39],[-34,28],[-46,76],[0,56],[28,95],[18,20],[55,21],[21,57]],[[20484,6289],[0,-252]],[[20484,6037],[-11,-1],[-2,-425],[0,-284],[-1,-284],[-1,-284],[-1,-284],[-1,-284],[-1,-426],[-245,0],[-245,0],[-245,0],[-246,0],[27,-92],[33,-26]],[[19545,3647],[-241,3],[-240,3],[0,-227],[-228,-1]],[[15853,8813],[203,0],[203,0],[203,0],[204,0]],[[17476,8813],[0,-394],[0,-236],[0,-394],[0,-316],[0,-394],[0,-473],[0,-315]],[[27005,9461],[102,39],[63,13],[77,4],[91,-19],[59,-39],[26,-6],[70,13],[80,-8],[77,42],[31,47],[29,27],[50,15],[24,19],[2,35],[-11,78],[10,38],[27,47],[-27,59],[-22,-19],[-20,30],[31,39],[111,108],[9,47]],[[28560,10330],[-3,-73],[4,-43],[-6,-90],[17,-69],[-3,-97],[-24,-69],[-4,-47],[11,-82],[2,-57],[-10,-68],[29,-6],[16,-32],[-5,-348],[3,-56]],[[28485,8304],[-61,-81],[-26,4],[11,37],[-1,40]],[[28194,8488],[-16,34],[-49,26],[-24,39],[-12,81],[4,26],[-19,45],[-32,19],[-30,54],[-225,0],[-188,1],[-226,0],[-300,0],[-263,0],[0,139]],[[26814,8952],[68,58],[26,34],[57,49],[40,63],[50,60],[-11,76],[3,61],[-43,18]],[[28790,8301],[-2,-36],[62,55],[38,14],[13,-18],[-65,-61],[-59,-35],[-56,-24],[-118,-63],[-19,5],[-97,-32],[-39,-4],[-16,32],[-29,-29],[-23,-8],[-5,29],[42,77],[64,45],[19,19],[55,13],[26,-10],[45,16],[133,18],[68,67],[-37,-70]],[[28333,8067],[0,46],[24,23],[5,-37],[-29,-32]],[[26608,8808],[0,-307],[0,-371]],[[26608,8130],[-32,-22],[-6,-21],[15,-59],[1,-41],[-19,-90],[-36,-140],[-19,-87],[-1,-36],[-24,-47],[-48,-59],[-53,-47],[-29,15],[-30,-42],[-29,-18],[-22,-64],[-18,-30],[7,-57],[-35,-37],[-3,26],[-23,30],[-19,-24],[-36,-99],[7,-66],[-31,-43],[-7,-36],[-19,-24],[-30,-9],[-31,17]],[[25809,8684],[148,-104],[21,21],[36,-30],[-74,-23],[3,-18],[38,13],[49,1],[32,-31],[37,11],[73,33],[42,7],[66,2],[23,19],[52,62],[60,50],[98,62],[95,49]],[[22804,4596],[-61,26],[-44,35],[-25,38],[-51,47],[-25,10],[-10,-23],[-27,-19],[-40,7],[-27,22],[-39,-20],[-5,-14],[-43,10],[-54,-31],[-38,-33],[-10,-23],[-26,32],[-49,38],[-16,29],[-26,-27],[-28,15],[-12,26],[-23,2],[-24,-41],[-9,-54],[-20,3],[-12,71],[-41,-29],[-18,-1],[-38,39],[-24,38],[-50,-58],[-28,11],[-2,37],[-32,35],[-8,53],[-70,-2],[-17,-27],[-20,1],[-32,36],[-52,-3],[-46,27],[-59,14],[-4,42],[-19,44],[-25,25],[-12,-32],[-35,6],[-27,-8],[-64,78],[-36,22],[0,302],[0,242],[0,423],[-256,0],[-306,0],[-255,0]],[[14701,8813],[-34,62],[-18,130],[4,103],[-37,115],[24,101],[29,166],[16,35],[28,122],[10,20],[5,184],[8,140],[14,46],[-5,48],[6,65],[-4,65],[30,315],[-4,38],[10,51],[-9,134],[4,150],[9,21],[65,1],[42,20],[32,-34],[50,-1]],[[14976,10910],[27,13],[49,-33],[26,-69],[32,-153],[21,-24],[122,-18],[32,6],[55,44],[60,14],[71,-3],[43,-13],[16,-21],[38,-3],[90,26],[99,18],[53,18],[61,39],[137,48],[69,6],[43,28],[294,0],[284,0]],[[28091,7795],[-27,-45],[-40,-15],[-27,-25]],[[26892,7663],[-284,0],[0,467]],[[26608,8808],[50,25],[156,119]],[[29138,8665],[-27,28],[-26,-65],[-5,-91],[-21,-38],[-67,-24],[-20,2]],[[27141,4712],[-76,-77],[-21,-33],[-60,-128],[-15,-82],[-22,-55],[-39,-47],[-46,-21],[-8,-46],[-51,-61],[-35,11],[10,-41],[-12,-31],[-22,-24],[-57,-22],[-36,-38],[-57,17],[26,-56],[-3,-36],[-33,-30],[-11,52],[-32,-49],[19,-42],[-17,-37],[-25,-14],[-6,-43]],[[20204,9318],[-1,379],[-1,252],[-1,379]],[[25750,6089],[194,-1],[274,-2],[12,10],[68,-4]],[[22992,2743],[-13,1],[-29,-83],[16,-64],[-57,-9],[-178,-127],[61,65],[-21,10],[-35,-16],[-12,6],[14,54],[-4,48],[-25,1],[-16,-38],[-36,13],[8,-86],[29,-81],[-69,-103],[-3,-45],[-33,-58],[-104,-110],[-54,-53],[-46,-27],[-25,18],[-44,16],[-47,-30],[-21,27],[-14,34],[-17,-4],[31,-111],[28,-16],[-38,-45],[-31,-13],[-27,40],[-9,-101],[-23,-32],[-20,16],[-14,-13],[-38,-10],[4,-42],[29,17],[-10,-54],[-27,-54],[-22,-13],[-33,8],[-16,-17],[39,-84],[-25,-127],[-16,-46],[-23,-7],[-43,41],[-3,-54],[56,-25],[3,-61],[-22,-76],[25,-139],[8,-104],[9,-45],[51,-166],[18,-1],[1,-53],[-37,-10],[-21,-36],[-22,11],[-41,47],[-58,29],[-77,11],[-52,24],[-28,36],[-29,21],[-31,7],[-25,19],[-20,33],[-30,20],[-39,9],[-25,24],[-17,60],[-16,99],[-20,62],[-38,77],[1,68],[-19,86],[7,64],[-6,41],[-25,44],[-43,47],[-37,70],[-30,91],[-30,63],[-30,35],[-20,43],[-11,51],[1,37],[-18,57],[-42,104],[-23,77],[-6,48],[-26,58],[-72,110],[-6,30],[-72,88],[-21,54],[-16,18],[-85,3],[-65,6],[-47,15],[-29,23],[-20,-3],[-12,-29],[-24,-19],[-38,-9],[-33,-54],[-27,-101],[-16,-115],[-35,-43],[-19,-45],[-21,-22],[-24,1],[-45,35],[-66,69],[-51,42],[-38,17],[-33,31],[-28,48],[-51,47],[-28,53],[-32,89],[-16,69],[0,72],[-42,158],[-22,69],[-17,31],[-33,38],[-48,44],[-65,88],[-81,131],[-58,79],[-34,27],[-29,47],[-25,68],[-27,45]],[[22652,2323],[46,73],[2,-19],[-48,-54]],[[22115,1695],[9,57],[22,52],[16,-15],[-22,-36],[-25,-58]],[[22022,1391],[-6,14],[22,99],[12,9],[-28,-122]],[[22072,815],[-26,86],[-37,248],[-1,141],[15,-151],[40,-253],[9,-71]],[[18291,8813],[0,-505],[273,0],[273,0]],[[27558,7242],[-4,-57],[-12,-28],[-46,-61],[-15,-102],[23,-29],[33,15],[17,-7],[38,-80],[71,-32],[26,-20],[22,-42],[31,-24],[25,-35],[-11,-87],[-11,-23],[-40,3],[-82,129],[21,-64],[26,-28],[37,-64],[49,-29],[10,-39],[1,-69],[-37,14],[2,-43],[-28,-27],[26,-17],[20,-32],[12,-49],[-32,-31],[-55,77],[-8,39],[-19,-2],[-83,51],[1,-24],[21,-23],[69,-26],[10,-64],[36,-43],[4,-33],[24,-3],[43,32],[66,-21],[39,-182]],[[27869,6062],[-21,0]],[[26214,6561],[-2,-34],[16,-52],[23,-36],[46,-39],[26,-6],[59,63],[33,-36],[77,20],[20,19],[14,42],[27,-11],[58,34],[12,-18],[40,39],[20,56],[-11,28],[26,75],[59,112],[19,68],[24,43],[27,69],[13,71],[19,17],[25,-44],[30,-18],[42,11],[38,99],[17,59],[13,15],[34,-1],[30,46],[18,7],[30,39],[40,81],[31,147],[132,-154],[30,101]],[[28009,6806],[-60,-199],[3,-36],[-31,-22],[-31,-46],[-11,-65],[-22,-73],[-17,57],[6,68],[24,111],[26,69],[20,33],[16,67]],[[28050,6814],[-7,-2]],[[14976,10910],[-51,34],[-16,25],[-50,-1],[-11,16],[-56,-16],[-18,16],[-31,-11],[8,48],[-1,59],[28,-29],[16,115],[-50,42],[-11,60],[73,51],[-39,10],[-15,23],[-33,-7],[-10,98],[-30,99],[-18,128],[-23,63],[-44,61],[-11,35],[-11,90],[6,68],[-8,47],[21,-2],[55,-38],[69,-29],[55,-38],[186,-24],[46,15],[27,-35],[45,5],[22,24],[11,-63],[22,-67],[-33,-72],[-12,28],[-62,-122],[18,4],[62,71],[11,31],[35,47],[17,-43],[-29,-39],[21,-145],[-6,-57],[-64,18],[-29,-19],[-19,-60],[20,-20],[28,7],[30,-18],[27,29],[16,54],[34,19],[18,30],[-4,117],[-12,25],[8,32],[-3,53],[41,96],[-30,52],[-38,8],[-8,50],[31,33],[-22,40],[-49,46],[11,17],[36,3],[-5,83],[-13,55],[-25,-7],[-19,45],[-18,71]],[[15022,12344],[-10,0]],[[15095,12182],[-30,-40],[-25,17],[24,42],[31,-19]],[[15032,12095],[-34,4],[-8,39],[15,14],[27,-57]],[[15152,11921],[19,-38],[37,-48],[-12,-34],[-40,37],[-13,69],[-39,56],[22,56],[26,5],[8,-33],[-40,-27],[32,-43]],[[15201,11537],[-32,-19],[7,66],[25,-47]],[[23451,11217],[33,-32],[41,2],[92,32],[39,23],[38,38],[24,-10],[50,44],[29,4],[19,-27],[-26,-75],[-16,-79],[42,21],[9,22],[54,-45],[35,-3]],[[24678,10382],[-6,-38],[-56,-46],[-33,-90],[-21,-73],[17,-34],[49,51],[30,70],[20,28],[50,10],[23,-46],[-41,-135],[-11,-62],[-2,-106],[-34,-59],[-18,-76],[7,-93],[-20,-63],[-14,-79],[-15,-48],[-6,-52],[6,-71],[-4,-31],[29,-118],[-7,-157]],[[24829,10480],[-2,-101],[-33,-60],[-3,-31],[-25,-44],[-27,15],[-4,45],[22,39],[15,54],[37,44],[20,39]]],"transform":{"scale":[0.0036708928579972937,0.001980669819389381],"translate":[-178.19451843993753,24.544579564750705]},"objects":{"ne_50m_admin_1_states_provinces_lakes":{"type":"GeometryCollection","geometries":[{"arcs":[[0,1,2,3]],"type":"Polygon"},{"arcs":[[[-2,4,5,6,7,8,9,10,11,12]],[[13,14]],[[15,16]],[[17]],[[18]],[[19]],[[20]],[[21]],[[22]],[[23]],[[24]],[[25]],[[26]],[[27]],[[28]],[[29]],[[30]],[[31]],[[32]],[[33]],[[34]],[[35]]],"type":"MultiPolygon"},{"arcs":[[36,37,38,39,40,41]],"type":"Polygon"},{"arcs":[[[42,43,44,45,46]],[[47]],[[48]]],"type":"MultiPolygon"},{"arcs":[[[49,50]],[[51]],[[52]],[[53]],[[54]]],"type":"MultiPolygon"},{"arcs":[[[-43,55]],[[56]],[[57]]],"type":"MultiPolygon"},{"arcs":[[[-3,-13,58,59,60,61]],[[62]],[[63,64]]],"type":"MultiPolygon"},{"arcs":[[[-41,-61,65]],[[-64,66]],[[67]],[[68]],[[69]],[[70]],[[71]],[[72]],[[73]],[[74]],[[75]],[[76]],[[77]],[[78]],[[79]],[[80]],[[81]],[[82]],[[83]],[[84]],[[85]],[[86]],[[87]],[[88]],[[89]],[[90]],[[91]],[[92]],[[93]],[[94]],[[95]],[[96]],[[97]],[[98]],[[99]],[[100]],[[101]],[[102]],[[103]],[[104]],[[105]],[[106]],[[107]],[[108]],[[109]],[[110]],[[111]],[[112]],[[113]],[[114]],[[115]],[[116]],[[117]],[[118]],[[119]]],"type":"MultiPolygon"},{"arcs":[[[-37,120,121,122,123,124,125,126,127,128,129]],[[130,131]],[[132]],[[133]],[[134]],[[135]],[[136]],[[137]]],"type":"MultiPolygon"},{"arcs":[[138]],"type":"Polygon"},{"arcs":[[[-46,139,140,141,142,143]],[[-50,144,-122,145,-132,146]],[[147]],[[148]],[[149]],[[150]],[[151]]],"type":"MultiPolygon"},{"arcs":[[-4,-62,-40,152,153]],"type":"Polygon"},{"arcs":[[-12,154,155,-59]],"type":"Polygon"},{"arcs":[[[-11,-14,-17,156,-155]],[[157]],[[158]],[[159]],[[160]],[[161]],[[162]],[[163]],[[164]],[[165]],[[166]],[[167]],[[168]],[[169]],[[170]],[[171]],[[172]],[[173]],[[174]],[[175]],[[176]],[[177]],[[178]],[[179]],[[180]],[[181]],[[182]],[[183]],[[184]],[[185]],[[186]],[[187]],[[188]],[[189]],[[190]],[[191]],[[192]],[[193]],[[194]],[[195]],[[196]],[[197]],[[198]],[[199]],[[200]],[[201]],[[202]],[[203]],[[204]],[[205]],[[206]],[[207]],[[208]],[[209]]],"type":"MultiPolygon"},{"arcs":[[210,211,212,213,214]],"type":"Polygon"},{"arcs":[[215,216,217,218,219,220]],"type":"Polygon"},{"arcs":[[221,222,223,224,225]],"type":"Polygon"},{"arcs":[[[-223,226,227,228]],[[229]],[[230]],[[231]],[[232]]],"type":"MultiPolygon"},{"arcs":[[233,234,235,236,237,238]],"type":"Polygon"},{"arcs":[[239,240,241,242]],"type":"Polygon"},{"arcs":[[243,244]],"type":"Polygon"},{"arcs":[[245,246,247]],"type":"Polygon"},{"arcs":[[[-215,248,249]],[[250]],[[251]]],"type":"MultiPolygon"},{"arcs":[[-214,252,253,254,255,-249]],"type":"Polygon"},{"arcs":[[256,257,258,259,260,261]],"type":"Polygon"},{"arcs":[[-6,262,263,264,265,266,267]],"type":"Polygon"},{"arcs":[[-262,268,269,270,271,272]],"type":"Polygon"},{"arcs":[[-271,273,274,275,276]],"type":"Polygon"},{"arcs":[[-237,277,278,279]],"type":"Polygon"},{"arcs":[[-272,-277,280,281,282,283,284]],"type":"Polygon"},{"arcs":[[[-218,287,288,289]],[[290]]],"type":"MultiPolygon"},{"arcs":[[[-243,291,292,293,294,295]],[[296]],[[297]]],"type":"MultiPolygon"},{"arcs":[[-245,298,299,300,-246,301,302,303]],"type":"Polygon"},{"arcs":[[[-45,306,307,-140]],[[308]]],"type":"MultiPolygon"},{"arcs":[[[-128,309,310,-275,311]],[[312,313]],[[314]],[[315]],[[316]],[[317]],[[318]]],"type":"MultiPolygon"},{"arcs":[[-38,-130,319,320,-260,321,322]],"type":"Polygon"},{"arcs":[[-221,323,-279,324,-257,-273,-285,325,-287,326]],"type":"Polygon"},{"arcs":[[-212,327,-288,-217,328]],"type":"Polygon"},{"arcs":[[-1,-154,329,330,331,-263,-5]],"type":"Polygon"},{"arcs":[[[-254,332,333,334,335]],[[336,337]]],"type":"MultiPolygon"},{"arcs":[[-39,-323,338,-330,-153]],"type":"Polygon"},{"arcs":[[-236,339,340,-258,-325,-278]],"type":"Polygon"},{"arcs":[[-141,-308,341,-294,342]],"type":"Polygon"},{"arcs":[[343,344,345]],"type":"Polygon"},{"arcs":[[-226,-239,346,347,348]],"type":"Polygon"},{"arcs":[[-224,-229,349,-266,350]],"type":"Polygon"},{"arcs":[[[-126,351,-124,-143,352,-292,-242,353,-345,354,355]],[[356]],[[357]]],"type":"MultiPolygon"},{"arcs":[[358,359,-281,-276,-311,360]],"type":"Polygon"},{"arcs":[[-220,361,-347,-238,-280,-324]],"type":"Polygon"},{"arcs":[[-228,362,363,-267,-350]],"type":"Polygon"},{"arcs":[[-355,-344,364,-247,-301,365,-359,366]],"type":"Polygon"},{"arcs":[[-240,-296,367]],"type":"Polygon"},{"arcs":[[-255,-336,368]],"type":"Polygon"},{"arcs":[[-259,-341,369,-331,-339,-322]],"type":"Polygon"},{"arcs":[[-213,-329,-216,-327,-286,-326,-284,370,-333,-253]],"type":"Polygon"},{"arcs":[[[-219,-290,371,-348,-362]],[[372]],[[373]],[[374]],[[375]]],"type":"MultiPolygon"},{"arcs":[[-225,-351,-265,376,-234]],"type":"Polygon"},{"arcs":[[[-244,377,-337,378,-334,-371,-283,379,-299]],[[-303,380]]],"type":"MultiPolygon"},{"arcs":[[-142,-343,-293,-353]],"type":"Polygon"},{"arcs":[[[-7,-268,-364,382]],[[384]],[[385]],[[386]],[[387]]],"type":"MultiPolygon"},{"arcs":[[[-261,-321,388,-313,389,-269]],[[390]]],"type":"MultiPolygon"},{"arcs":[[-282,-360,-366,-300,-380]],"type":"Polygon"},{"arcs":[[-235,-377,-264,-332,-370,-340]],"type":"Polygon"}]}}}
},{}],20:[function(_dereq_,module,exports){
module.exports={"type":"Topology","arcs":[[[50499,34740],[-44,-3],[-55,-36],[-6,-28],[-64,2],[-122,-28],[-142,-12],[-127,-40],[-65,-48]],[[50907,35561],[-16,49],[-45,26],[14,50],[-57,78],[-6,90],[-51,36],[-43,93],[-48,61],[60,39],[25,52],[111,38],[87,45],[37,37],[65,13],[59,-40],[146,-12],[98,35],[34,58],[92,-21],[13,-88],[34,-69],[55,-44],[31,-81],[-1,-57],[-102,-46],[-14,-49],[21,-95],[-40,-78],[68,-42],[-19,-126],[-97,-36],[-39,41],[-58,-59],[-53,-7],[47,-93],[22,-19],[-20,-66],[-96,-8],[-15,-63],[-102,-19],[-87,-38],[-18,-55],[-80,11],[28,-133],[68,-53],[-26,-59],[-63,-3],[-64,-51],[-102,-8],[-82,31],[-42,-33],[-46,3],[-12,-58],[-80,2]],[[47052,20362],[35,16],[28,-89]],[[47388,28234],[-73,29],[-39,41],[25,88]],[[18304,27306],[18,-72]],[[18379,27488],[-30,-47]],[[17191,33219],[55,-21],[46,5],[90,33],[17,-8]],[[17124,33251],[42,-31]],[[17082,33308],[-25,40]],[[17050,33448],[5,-90]],[[17047,33491],[-9,31]],[[45716,17193],[3,38],[22,31]],[[70294,33019],[-4,-78],[-42,-39],[-24,-60],[9,-51]],[[71001,33448],[-59,-60]],[[63185,33652],[-65,-8],[-75,-80],[-79,-68],[-45,-20],[-87,-9],[-29,17]],[[62459,33889],[2,160],[-25,66],[-44,77],[-7,61],[31,45],[11,80],[-30,57],[6,39],[45,95],[-33,129],[1,39],[-69,143],[-81,22],[-159,-41],[-14,72]],[[62296,35475],[-16,-74],[-55,-89],[-2,-145]],[[46775,21372],[-7,-104]],[[12609,35720],[-14,-61]],[[12547,35752],[62,-32]],[[13186,37645],[82,41]],[[13290,37695],[25,2]],[[15483,38093],[-78,49],[-111,-47],[-89,-16],[-61,8],[-109,-13],[-170,6],[-106,19]],[[15806,38108],[63,58]],[[58827,33318],[-16,89],[-21,32],[-48,9]],[[18666,26615],[19,-67]],[[18643,26848],[-14,46],[14,90]],[[18934,27743],[-61,-65]],[[23382,33895],[2,84]],[[11827,27801],[16,-35],[43,-38],[-9,-22]],[[11778,28067],[-13,91]],[[11743,28284],[34,44],[44,0],[51,-39],[42,54],[57,-61],[29,-58]],[[12760,28902],[-24,-65],[-52,-71],[-52,-107],[-119,-51],[-11,-23]],[[13817,29651],[-9,-6]],[[10968,32432],[-24,-50],[2,-92]],[[10964,32805],[-2,-104],[45,-47]],[[11007,33094],[-3,-71],[21,-101],[-12,-18]],[[11097,33414],[-82,69]],[[10863,33560],[37,-5]],[[12300,37756],[69,-41],[19,-45]],[[11875,37999],[-45,-32],[135,-68],[99,-12],[226,-53],[16,-58]],[[42244,36200],[-26,-11]],[[42472,36316],[-162,-106],[-51,-7]],[[41806,36477],[92,-63],[6,-26]],[[41612,36613],[-118,82]],[[45509,31592],[-49,72]],[[46221,31831],[74,104],[32,86],[37,26],[68,17],[149,-5],[14,43],[-10,61]],[[46140,32573],[-71,43],[-66,25],[-167,83]],[[45733,32870],[6,-80]],[[45516,33269],[-6,-153],[8,-50]],[[48213,32084],[180,92],[38,99],[44,67],[68,74],[-12,50],[31,41],[67,25],[20,55],[-54,23],[-47,98]],[[46416,20338],[7,-47],[23,-40],[63,-72]],[[57806,32753],[55,-40],[52,-22],[88,-96],[-113,-58],[-43,-31],[-19,-73],[32,-101],[-17,-45],[196,-161]],[[47436,29202],[11,20],[-17,67]],[[65390,26076],[25,51],[25,103]],[[28101,9436],[-53,48]],[[27911,9558],[-55,26],[-68,85],[-80,27],[-57,1]],[[27477,9823],[77,-64],[36,-47]],[[26565,9883],[96,92]],[[26939,9923],[-227,40]],[[64343,27154],[-111,57]],[[61044,27898],[-51,31]],[[61154,27980],[-12,-41],[-38,-28]],[[55345,27828],[16,-120],[-33,-57]],[[46659,26747],[9,77],[23,67]],[[46699,27263],[-4,-59]],[[44790,38368],[30,-48],[-36,-50]],[[10073,30403],[26,-74]],[[56791,21797],[37,1],[12,22],[7,70],[17,21]],[[21950,33994],[-64,54]],[[21673,34050],[51,16]],[[21115,34078],[-61,-17],[-214,14],[-64,12],[-180,-14],[-87,3],[-88,40]],[[11062,36578],[104,-51],[92,-23],[193,-11],[115,36],[100,46],[117,14],[201,-31]],[[12589,31151],[-48,36]],[[57921,23531],[-34,51],[-31,20],[-36,-26],[-29,12]],[[16489,32020],[55,9]],[[16448,32036],[-31,-5]],[[15318,30900],[9,-32],[-24,-54],[27,-70],[120,-121],[55,-29]],[[12493,31677],[-6,35],[-36,64]],[[14990,31792],[35,-34],[16,-70],[-13,-65],[13,-51],[36,-72],[8,-48],[-15,-46],[31,-72],[-5,-86],[-43,-83],[-44,-8],[37,-98],[33,-14]],[[14845,32052],[-87,9],[-90,21],[-43,67],[-64,26],[1,52],[-60,33],[-71,-1],[-30,-25],[-75,18],[-33,-29]],[[13681,32221],[-45,-65],[-65,-38],[-118,-4],[-70,-28],[-29,-32],[-30,23],[-59,7]],[[38024,19305],[78,42],[47,-29],[46,-3]],[[72612,4634],[54,19]],[[72377,4665],[43,-14]],[[38196,19314],[-32,64],[-9,86]],[[55466,30174],[-86,4],[-37,-9],[-27,-31]],[[15961,33191],[1,108],[-25,99],[-157,204],[-61,99],[-44,167],[-19,156]],[[15656,34024],[10,64],[-34,57],[7,33],[59,31]],[[15658,34413],[26,-44]],[[15648,34493],[161,94]],[[15972,34858],[115,63],[51,12]],[[16233,34949],[151,18]],[[45640,36133],[54,34],[33,59]],[[39498,19709],[-22,123],[1,61],[28,60]],[[37497,21485],[-43,28],[-32,47]],[[46067,24340],[7,-59],[-8,-70],[10,-86],[-14,-40],[-63,-109],[-2,-23],[-37,-44],[-21,-6],[-44,26],[-41,-34],[-16,-34],[-81,-90],[-34,-52],[-27,-81],[-44,-78]],[[57727,34459],[-166,-116],[-104,-53],[-110,-91],[-13,-55]],[[44461,6458],[-70,-39],[-38,17]],[[44173,6644],[57,-93]],[[68356,11184],[-7,-99],[11,-86],[-5,-30]],[[20015,31830],[-63,100],[-29,115],[15,31]],[[19893,19459],[-29,17],[-25,37]],[[26683,10238],[-3,22],[31,60],[28,26],[67,24],[39,0]],[[27209,10460],[-169,5],[-36,-20]],[[27390,10516],[-79,-69],[-53,-18]],[[25724,8135],[32,69],[14,116],[23,111],[8,73],[-9,58],[12,21]],[[26346,9450],[-26,-29],[-4,-37],[-78,-146],[-38,-46],[-92,-81]],[[26564,10025],[1,-142]],[[26565,9883],[-40,-68],[-53,-44]],[[10003,34844],[-60,30],[-49,4],[-147,-40],[-137,3],[-39,24],[-87,100],[-88,121]],[[11956,35726],[26,-13]],[[10033,29848],[-64,-7]],[[20412,32069],[-20,37]],[[20443,32075],[67,-32],[24,3]],[[20117,32110],[-20,-9]],[[20396,32130],[-10,28]],[[16709,32416],[140,21]],[[17189,32486],[-28,24],[-44,-12]],[[40545,32106],[-58,-5]],[[40557,32106],[87,-52]],[[39849,31626],[23,42],[45,32],[102,-10]],[[15297,25347],[53,-97],[13,-65]],[[14929,26135],[-62,57]],[[13475,27266],[2,67]],[[13501,27352],[2,60]],[[20560,33334],[62,-16]],[[20936,33344],[72,-26],[18,15]],[[14524,23134],[3,-36]],[[14308,23340],[27,11]],[[13984,23630],[16,-14]],[[9969,29841],[-9,-13]],[[10590,28608],[23,20]],[[18712,20120],[-49,30],[-64,80],[-143,211]],[[15350,33915],[266,82],[40,27]],[[15095,34048],[38,-29],[43,-77],[29,-19],[72,6],[49,-23]],[[14375,34068],[-37,-8],[-72,-49]],[[14604,34176],[13,-5]],[[36785,18784],[-4,-43],[15,-56],[-7,-123],[4,-33]],[[19219,27629],[77,-88],[11,-28]],[[45270,16276],[34,154]],[[45716,17193],[-49,-88],[-50,-79],[-108,-145]],[[46620,11718],[5,-83]],[[61541,33312],[-27,21]],[[12074,30461],[-36,7]],[[12179,30513],[53,69]],[[12602,30643],[-33,43]],[[13577,33284],[38,-43],[-49,-55],[-20,-46],[-40,-33],[-118,-28],[-79,7]],[[45993,26371],[-46,80]],[[47503,36106],[104,-83],[69,-29]],[[57325,26571],[-22,20]],[[56272,26774],[-22,-43],[-55,47]],[[45727,36226],[42,57],[114,73],[63,8],[59,-19],[49,-37]],[[54781,29734],[-39,-35],[-84,-19]],[[54340,29954],[-95,19],[-13,30]],[[52783,31565],[-46,29]],[[27938,10227],[-20,82],[-6,69],[-31,69],[17,58]],[[28389,12820],[22,51],[64,65],[34,50],[31,69],[23,111],[24,23],[40,82],[50,24],[83,0],[70,17],[35,63],[23,15],[43,-11]],[[37028,29543],[79,-18],[103,5]],[[37787,29764],[36,65]],[[18217,27812],[-117,80],[-70,-15]],[[17966,28261],[4,89],[-36,171],[-39,94]],[[7707,36021],[-59,102],[-45,49],[-96,63],[-67,61]],[[14047,37156],[29,66]],[[16037,37560],[192,-47],[93,0],[144,-38],[47,11]],[[14749,37596],[64,27]],[[15844,37711],[-14,-16],[-142,9],[-141,-29],[-12,-22],[-100,-32],[-100,31],[-153,19],[-105,-22],[-91,10]],[[26869,15217],[11,-123],[22,-110],[36,-138],[26,-83],[38,-38],[41,11]],[[19085,27820],[74,-53],[6,-23]],[[52103,33591],[-23,-135],[-25,-30]],[[44489,7423],[-7,-28],[-46,31],[-56,-24]],[[45048,7624],[-87,40]],[[59700,34095],[-13,31],[23,82],[-66,131],[-77,79],[29,93],[67,27],[117,29],[77,38],[52,70]],[[46139,16615],[17,-16],[-56,-142],[-29,-93],[-35,-68],[-43,-43],[-62,-37],[-122,-48]],[[45989,16973],[33,-31],[35,-7],[14,-30]],[[45919,17021],[55,-38]],[[63903,37282],[53,-9],[76,18],[47,-19],[18,-65],[52,-32],[5,-39],[57,-24],[122,46],[71,53],[39,-41],[-31,-49],[51,-34],[113,53]],[[42796,32846],[30,2]],[[49040,32943],[39,44],[1,33],[-36,62],[15,31],[66,41],[-1,44],[-34,41],[-10,44],[33,68],[73,79],[43,22],[69,-9]],[[49830,34668],[52,-81],[-8,-40]],[[49874,34547],[-53,-26],[-29,-41],[-5,-57],[-22,-39],[-64,-42],[-11,-76],[87,-91],[30,-50],[10,-81],[47,-39],[162,-1]],[[48521,35059],[-21,19],[-60,139],[22,56],[-71,27],[-69,2]],[[47061,35085],[-63,9],[-84,-45],[-34,20]],[[47202,35261],[74,30],[52,62]],[[48035,35327],[-32,8],[-85,93],[-80,-15],[-113,-62]],[[47360,35499],[-28,44],[21,45],[69,-57]],[[38195,19315],[58,-98],[47,-54],[54,-43],[33,-40],[32,-68],[26,-81],[26,-43],[5,-38],[-9,-183],[-20,-44],[-1,-70]],[[41019,35949],[8,36]],[[41026,36041],[-9,79]],[[45509,36338],[-42,26],[-136,7],[-27,13]],[[45522,36356],[103,-12],[42,-35],[59,-82]],[[79299,3706],[60,80]],[[16423,32868],[24,-94],[-9,-44],[-59,-48],[-14,-55],[14,-84],[21,-46]],[[16094,32963],[9,40]],[[15961,33191],[49,-96],[22,-20]],[[64675,25994],[-24,-34],[-60,-29],[-45,7],[-83,41]],[[64695,26032],[18,79]],[[44705,10585],[59,104],[38,43],[46,82],[28,66],[41,54],[109,76],[84,44]],[[46028,11366],[-50,5],[-72,-27],[-61,-7],[-83,6],[-117,-18],[-62,20]],[[14290,35720],[28,-48]],[[12609,35720],[50,17],[171,118],[157,56],[66,10],[137,-2],[234,21],[127,-19],[220,-11]],[[13880,35871],[125,7]],[[14422,36701],[-53,-104]],[[14491,36910],[49,-45],[-79,-67],[-25,-46]],[[15042,36921],[111,-42],[55,-4],[51,28]],[[14900,36925],[100,-3]],[[14744,36976],[14,-68],[95,12]],[[15554,37007],[26,33],[126,35],[15,51]],[[14815,37023],[-31,83],[-87,2]],[[14552,37054],[6,-79],[-24,-35]],[[15761,37185],[40,117]],[[16129,37348],[-99,-14],[-151,8]],[[16217,37478],[12,35]],[[47067,20093],[-37,51],[-9,32],[16,96],[30,57],[-15,33]],[[47115,20289],[99,-113],[48,-34],[66,-23],[37,-25],[12,-87],[-8,-150],[-24,-59],[-62,-29],[-62,-44],[-42,-10],[-40,46],[-42,6],[-22,27],[-13,47],[-112,31],[-53,31],[-29,-73],[-27,-25],[-35,-5],[-33,-28],[-67,-13],[-63,46]],[[47388,28234],[83,17],[72,-11],[90,-32],[51,-33],[13,-33],[125,-136],[12,-49],[33,-59],[55,-57],[32,-77]],[[47264,28553],[7,-31],[41,-40],[8,-45],[-19,-45]],[[18322,27234],[41,-62],[6,-31],[-22,-39],[-48,-15],[-76,19],[-58,-41],[-3,-39],[-74,-47],[8,-39],[-27,-117],[-22,-48],[-38,-31],[-36,-66],[-1,-65],[-14,-73]],[[18349,27441],[-36,-83],[-9,-52]],[[18616,27708],[-121,-18],[-58,-63],[-49,-98],[-9,-41]],[[17166,33221],[25,-2]],[[17399,33228],[129,35],[62,26],[125,78],[34,-11],[76,29],[27,-19],[88,-6],[55,-20],[53,-39],[70,-14],[89,46],[56,10],[61,46],[20,55],[24,8],[57,-35],[90,-2],[50,-18],[115,-117],[62,-30],[65,13],[173,77],[70,46],[65,81],[33,25],[90,29],[59,41],[99,37]],[[17081,33308],[43,-57]],[[17055,33358],[2,-10]],[[17047,33491],[3,-43]],[[16961,33518],[77,4]],[[45867,17660],[-19,-29],[-50,-9],[-55,-125],[-24,-83],[10,-70],[21,-58],[-9,-24]],[[67154,35011],[-117,17],[30,122],[84,128],[-2,101],[62,22],[-22,85],[7,32],[56,71],[72,24],[34,87],[31,30],[106,45],[63,58],[7,92],[45,-12],[48,-43],[171,-83],[156,-48],[120,33],[81,-45],[52,33],[58,2],[21,-41],[57,-49],[119,0],[77,-15],[134,48],[82,-5],[3,67],[93,34],[123,3],[84,39],[27,35],[47,1],[100,64],[75,11],[33,55],[100,72],[155,156],[82,58],[106,-18],[29,19],[-48,78],[102,85],[33,79],[-6,68],[-23,25],[35,48],[-56,15],[-7,33],[62,27],[42,63],[-5,40],[32,33],[-78,137],[-88,-20],[-110,33],[-71,37],[-43,52],[-80,47],[-135,-13],[-18,28],[20,55],[-177,-11],[-105,3],[-157,22],[-78,-24],[-134,8],[-244,50]],[[20319,30192],[-33,-17],[-46,3],[-31,29],[-5,44],[-86,43],[-40,-27],[-19,-53],[-77,-36],[-48,-102],[-69,-35],[-2,-36],[32,-35],[32,-109],[-67,-100],[-37,-23]],[[26177,15966],[-49,-20],[-63,-50],[-49,-25],[-66,-17],[-47,-22],[-47,-60],[-5,-41],[-28,-30],[-66,-27],[-39,7]],[[25718,15681],[-13,77],[-41,32],[-65,-31],[-53,33],[-41,51],[-27,0],[-72,-26],[-35,-25],[-36,-43],[-45,-79],[-30,-34],[-34,-8],[-46,17],[-44,-7],[-69,48],[-39,-25],[-31,-59],[-34,-41],[-35,-20],[-23,-37],[-11,-52],[-20,-32],[-65,3],[-33,-44]],[[24776,15379],[-27,-2],[-41,25],[-93,28],[-82,35],[-38,-16]],[[24495,15449],[-14,-36],[-57,-29],[-48,8],[-44,-6],[-35,-26],[-25,-46],[-42,-23],[-56,-2],[-44,-20]],[[24130,15269],[-58,-68],[-41,-21],[-110,75],[-41,-44],[-74,-28],[-55,-46],[-43,51],[-68,23],[-65,3],[-102,70],[-38,35],[-77,99]],[[23358,15418],[-66,82],[-25,45],[-50,42],[-20,30],[-65,4],[-51,45],[-54,2],[-29,-25],[-31,32],[-24,-13],[-6,-36],[-25,-34],[-51,-12],[-77,36],[-30,-2],[-44,-51],[-38,-19],[-25,-30],[-20,-59],[-13,-71],[-65,6],[-77,-13],[-113,-34],[-63,-33],[-14,-70],[3,-37],[-17,-56],[10,-37],[-29,-45],[-67,-11],[-24,15],[-16,41],[-35,36],[-17,55],[-23,23],[-130,-48],[-37,-9],[-52,39],[-29,2],[-51,33],[-38,132],[-117,-47],[-106,15],[-33,-4],[-48,-31],[-19,-26],[-8,-46],[5,-104],[-33,-6],[-30,-71],[-6,-36]],[[54283,28615],[-34,-24],[-54,37],[-31,54],[-43,13],[-46,-19],[-58,41],[-56,12],[-83,-9],[-84,23],[-207,126],[-21,46],[-51,49],[-84,52],[-64,61],[-46,69],[-40,86],[-42,47],[-72,60],[-93,55],[-27,68],[-75,78],[-31,63],[-15,97],[-18,44],[-48,71],[-32,114],[-72,69],[-89,19],[-11,20],[-118,105],[-85,105],[-33,53],[-58,42],[-85,41],[-66,71],[-48,217]],[[52164,30670],[-1,1]],[[52163,30671],[-14,67],[51,44],[-10,57],[-30,54],[40,31],[67,27]],[[52163,30671],[-21,42],[3,102],[-104,71],[-14,83],[-45,74]],[[70233,32791],[4,-54],[-101,-74],[-27,-93],[-27,-32],[-96,-73],[-102,-6],[-54,-33],[-10,-54],[-49,-31],[-78,6],[-150,-43],[-87,-58],[-75,3],[-62,-42],[-40,-10],[-34,-58]],[[69245,32139],[-20,-18],[-76,1],[-138,-12],[-55,15],[-73,-12],[-20,50],[-43,58],[-5,35],[22,70],[-59,86],[15,56],[-15,29],[-83,10],[-37,26],[-57,74],[-68,63],[-88,-9],[-70,33],[-27,45],[-108,-13],[-55,3],[-68,33],[-48,67],[2,65],[18,45]],[[70942,33388],[-59,-11],[-53,-46],[-107,-18],[-198,-162],[-83,-47],[-121,-35],[-27,-50]],[[71001,33448],[16,66],[32,52],[24,79],[-19,60],[-38,21],[-17,41],[-70,61],[-11,36],[58,61],[79,-10],[63,-31],[52,-1]],[[78175,38477],[88,-41],[50,-75],[-139,-52],[-115,-67],[-150,-68],[-13,-126],[-149,-47],[-112,-72],[-15,-41],[55,-53],[-46,-69],[66,-42],[198,-2],[150,-25],[78,19],[12,44],[71,37],[0,71],[115,14],[140,-45],[163,80],[80,50],[103,-76],[-25,-86],[51,-14],[88,-75],[89,-31]],[[62806,33483],[-100,97],[-54,77],[-34,23],[-67,86],[-17,46],[-72,52],[-3,25]],[[62093,34933],[4,88],[-35,66],[-64,29],[16,40],[74,6],[107,-32],[28,37]],[[62296,35475],[28,58],[10,83],[-42,106],[-45,46]],[[62253,35765],[-6,3]],[[62247,35768],[-84,31]],[[62247,35768],[-84,31]],[[62163,35799],[-61,8],[-131,-5],[-69,-14],[-38,-46],[-88,-43],[6,-42],[-50,-50],[-61,21],[-54,41],[-66,19],[-53,35],[-65,0],[-63,-33],[-30,-40],[-154,-21],[-170,5],[-57,-10],[-37,30],[-155,17],[-75,-6],[-83,-32],[-103,-73],[-88,-20],[-114,-7],[-193,13],[-65,-20]],[[62163,35799],[-36,13]],[[26072,10620],[-22,146],[2,35],[60,160],[36,36],[40,91],[42,68],[4,45],[38,48],[95,62],[27,84],[15,75],[45,71],[66,23],[17,45],[32,123],[9,86],[16,34],[2,95],[14,64],[47,119],[8,77],[19,56],[-19,8],[-21,50],[6,124],[-9,58],[2,84],[-13,35],[11,32],[-9,100],[22,113],[3,102],[34,99],[27,116],[38,100]],[[26756,13284],[21,101],[24,66],[63,106],[48,63],[51,132],[26,143],[-38,71],[-1,46],[23,57],[10,56],[-3,52],[30,49],[64,46],[35,40],[30,73],[27,5],[-1,49],[22,31],[15,55],[-1,73],[24,34],[-47,64],[-85,23]],[[27093,14719],[-50,17]],[[66651,33958],[-97,-11],[-66,-37],[-65,-67],[-77,-59],[-6,-50],[106,-6],[41,-40],[-11,-32],[22,-63],[-15,-76],[-40,-41],[-103,-81],[-59,-109],[-16,-53],[-39,-52],[-30,-80],[-39,-47],[-20,-59],[35,-34],[-20,-69],[-65,-29],[-52,-5],[-71,-39],[-134,-108]],[[17245,27557],[-21,41],[-30,13],[-127,94],[-30,9],[-32,56],[-10,74],[-17,27],[-51,40],[-26,38],[-64,55],[-50,8],[-72,59],[-38,16],[-38,-5],[-63,31],[-58,-42],[-35,21],[-65,-43],[-102,62],[-40,49],[-11,48],[-78,11],[-22,36],[-33,14],[-51,60],[-41,7],[-13,27],[-53,25],[-46,68],[-42,-23],[-27,50],[41,41],[-16,77],[-33,38],[-63,195],[-82,77],[-138,91],[-80,47],[-51,-13],[-85,-83],[-59,-82],[-52,-35],[-68,16],[-183,62],[-82,12],[-41,-18],[-75,11],[-109,38],[-107,17],[-107,-4],[-134,-19],[-56,29],[-94,33],[-95,8],[-68,36],[-80,29],[-86,-20],[-42,8],[-60,62],[-10,62],[-57,138]],[[46696,21721],[-54,103],[-19,62],[-54,108],[-21,61],[-42,40],[-80,59],[-45,44],[-25,50],[-39,30]],[[46697,21722],[6,-94],[23,-81],[44,-99],[13,-55],[-8,-21]],[[11084,33582],[5,33],[-28,87],[-44,36],[-34,55],[-17,67],[41,80],[52,50],[116,86],[18,75],[161,106],[65,-2],[121,-24],[102,17],[74,42],[39,3],[68,47],[31,105],[50,122]],[[11905,34565],[66,20],[37,-19],[31,-50],[13,-81],[35,0],[43,45],[79,30],[-26,38],[33,88],[-8,86],[30,46],[-11,52],[24,72],[-30,67],[47,55],[118,12],[106,39],[22,34],[-15,74],[-29,48],[-7,75],[61,134],[4,118],[-31,42],[3,46],[95,23]],[[12547,35752],[12,52],[-42,19]],[[61097,24902],[-16,-93],[-37,-35],[1,-62],[-56,-94],[15,-132],[15,-71],[-45,-23],[-79,64],[-39,-8],[-44,-84],[-3,-37],[-45,-97],[-19,-96],[-6,-129],[-13,-99],[23,-142],[-3,-106],[-100,5],[-37,-47],[-24,-105],[-33,-37],[-64,-38],[-19,-27],[1,-39],[-19,-57],[5,-50],[22,-74],[23,-133],[43,-74],[11,-118],[-12,-131],[11,-103],[-27,-60],[78,-134],[12,-106]],[[60617,22330],[4,-37]],[[60621,22293],[19,-71],[21,-29],[-9,-61],[-1,-115]],[[60651,22017],[-28,-29]],[[60623,21988],[-35,-30],[-21,-66],[-67,-102],[-13,-9]],[[13118,37645],[68,0]],[[13268,37686],[21,9]],[[13315,37697],[30,11],[111,116],[-8,63],[39,48],[66,4],[72,-51],[174,-45],[145,-29],[63,-2],[207,42],[94,28],[100,53],[28,36],[228,24],[25,67],[70,37]],[[15483,38093],[119,0],[66,28],[77,-21],[61,8]],[[15869,38166],[121,-19],[135,97],[22,75]],[[35611,19947],[20,-20],[40,7],[24,80],[25,6],[10,85],[32,59],[32,21],[49,85],[52,7],[29,45],[36,84],[29,49],[22,13],[6,71],[18,47],[10,85],[-48,136],[-87,28],[16,51],[-36,106]],[[45523,19021],[19,87],[5,121],[7,38],[26,41],[13,76],[-5,37],[17,50],[58,94]],[[45867,17660],[-71,108],[-39,82],[-9,99],[30,105],[31,158],[5,95],[-11,70],[-66,182],[-30,73],[-42,82]],[[45665,18714],[-53,50],[-35,56],[-36,78],[-24,79],[5,44],[-25,89],[-23,34],[2,113],[-32,147],[28,70],[3,53],[30,59],[37,14],[89,-12],[28,-21]],[[43010,18946],[-24,7],[-8,35],[-48,37],[-13,35],[4,54],[-15,65],[-51,115],[0,38],[-41,40]],[[65124,16660],[8,-35],[-32,-34],[-33,-70],[-69,8],[-24,-41],[19,-28],[-31,-50],[-17,-59],[53,-95],[31,-23],[24,5],[25,34],[43,-15],[23,-53],[-18,-49],[7,-39],[-16,-45],[19,-51],[-21,-126],[18,-63],[-6,-50],[9,-78],[-5,-50],[-27,-98],[4,-64],[-31,-36],[-26,-88]],[[73409,7079],[-30,29],[-44,2],[-110,-18],[-38,11],[-79,-75],[-67,-79],[-21,-46],[-36,-111],[-49,-53],[-41,-63],[-64,-43],[-131,-37]],[[41195,19587],[-33,2],[-48,-29],[-76,-59],[-49,-55],[-23,-53],[-36,-51],[-89,-20],[-43,-31],[-56,-86],[-104,-68],[-45,-47],[-123,-62],[-86,4],[-141,62],[-95,10],[-62,-17],[-88,-57]],[[58742,33448],[-50,28],[-16,42],[-8,104],[-48,53],[-101,59],[-43,-31],[-53,25],[-146,-12],[-78,-48]],[[44762,30948],[117,73],[38,41],[-13,85]],[[26684,12262],[21,75],[67,100],[11,72],[33,19],[-4,72],[10,137],[-19,131],[4,185],[-2,80],[-20,75],[-18,13],[-11,63]],[[59288,24453],[-7,-121],[37,-80],[48,-48]],[[59366,24204],[45,-37]],[[59288,24453],[16,-52],[72,-66],[23,-64],[33,-22],[47,-56]],[[60604,25657],[11,-39],[-18,-26],[-59,-35],[-73,-69],[-34,-54],[-22,-63],[-27,-38],[-58,-42],[-40,-7],[-63,-44],[-81,2],[-56,-22],[-118,-30],[-25,-19],[-54,-8],[-81,-74],[-75,-35],[-48,-4],[-66,17],[-55,-17],[-78,17],[-52,-9],[-63,-36],[-45,-54],[-31,-72],[-14,-63],[-4,-105],[-15,-90],[26,-81],[2,-104]],[[23592,17408],[18,51],[-13,53],[-4,98],[28,79],[27,4],[34,-27],[48,0],[34,-79],[106,-65],[47,-17],[51,11],[27,25],[76,26],[51,-2],[58,33],[37,9],[100,0],[22,-45],[10,-56],[-12,-65],[-34,-74],[-34,-143],[-37,-57],[11,-48],[-40,-80],[-24,-82],[-18,-118],[-29,-97],[-41,-79],[-28,-107],[-16,-135],[-6,-84],[9,-52],[-13,-70],[2,-115],[-5,-81]],[[45330,31316],[-33,36],[-80,11],[-107,-56]],[[45110,31307],[28,-17]],[[45138,31290],[75,-33],[17,-24],[76,-30]],[[45138,31290],[185,-10]],[[15054,27264],[14,18],[76,25],[27,19],[2,81],[23,31],[92,48],[58,-8],[84,-102],[43,-36],[12,-48],[50,1],[8,-48],[32,2],[17,-43],[50,-24],[78,-86],[-11,-42],[52,-34],[10,-49],[36,-83],[68,-108],[15,-58],[58,-149],[33,-59],[28,-19],[25,-46],[45,-29],[-13,-52],[28,-164],[63,-21],[37,-31],[6,-68],[-11,-52],[42,-44],[12,-35]],[[43248,16382],[-33,-4],[-91,22],[-6,24],[-60,16],[-52,-15],[-70,12],[-46,-13],[-44,-34],[-17,13]],[[67880,40181],[9,75],[232,-36]],[[41567,19025],[35,25],[14,71],[52,56],[10,48],[-17,39],[9,72],[-10,35],[-29,28],[-55,107],[-35,6],[-50,-21],[-26,8],[-62,50],[-56,-33],[-28,10],[-27,45],[-34,20],[-63,-4]],[[21857,3762],[27,54],[-41,82],[-23,89],[-99,118],[-61,1],[-38,41],[-89,34],[-24,69],[-47,54],[-4,32],[-41,73]],[[20627,16910],[-1,-62],[27,-45],[37,-8],[57,13],[76,-57],[64,-24],[38,-31],[10,-39],[38,-26],[35,-78],[50,-46],[10,-55],[54,-8],[56,-56],[34,-23],[75,-24],[41,-36],[38,0],[28,-21],[61,12],[62,-15],[41,-20],[30,25],[60,-6],[35,-57],[56,-48],[48,-21],[16,-44],[58,17],[65,-59],[75,-56],[85,-10],[52,-26],[33,80],[21,20],[24,-13],[55,8]],[[62569,25875],[22,-13],[91,3],[37,13],[30,33]],[[62749,25911],[30,13],[55,-9],[24,-21],[60,26],[53,38],[48,-5],[19,39],[49,18],[41,64],[42,29],[32,74],[54,0],[46,57],[19,8],[36,-29],[25,6],[26,46],[41,5],[48,56],[109,176],[44,102],[60,15],[164,24],[43,0],[79,17],[114,-14],[82,-32],[93,-46],[40,-90],[40,-72],[67,44],[25,-32],[65,-15],[-8,-50],[35,-35],[11,-92],[69,20],[58,-36],[-8,-52],[48,-22]],[[42814,19372],[-61,26],[-35,41],[-56,91],[-103,132],[-39,63],[-24,71],[-48,72],[-33,20],[-39,51],[-92,-28],[-62,39],[-75,57],[-72,20],[-44,40],[-20,49],[-11,78],[-15,50],[0,91],[-9,40],[-49,62]],[[41927,20437],[-29,-14],[-60,34],[-40,50],[-75,28],[-5,29]],[[41927,20437],[-27,15],[-18,46],[-8,76],[-19,51],[-51,53],[-6,53]],[[18685,26548],[-54,-180],[1,-30],[33,-110]],[[18643,26848],[-6,-111],[6,-59],[23,-63]],[[18873,27678],[-36,-56],[-25,1],[-91,-124],[-45,-41],[-36,-56],[1,-46],[-26,-81],[7,-44],[33,-91],[20,-92],[-32,-63]],[[18950,27858],[25,-42],[-42,-73]],[[19095,27977],[-11,-48],[-44,-73]],[[55196,26700],[-8,-45],[-77,-142],[-22,-10],[-57,-66],[-50,-130],[-20,-110],[-20,-26],[-8,-80]],[[54934,26091],[-36,-34],[-40,-61],[-17,-4]],[[56427,27064],[-70,54],[-45,44],[-90,47],[-27,35],[-3,48],[-20,44],[-74,66],[-60,15],[6,-65],[-74,-9],[-80,33],[-73,-40],[-5,-46],[-48,-63],[-92,-82],[-80,-55],[-69,-32],[-49,-43],[-46,-79],[-38,-48],[-91,-44],[-58,-54],[-28,-81],[-17,-9]],[[17198,21399],[6,88],[48,7],[114,-27],[18,30],[-46,36],[15,43],[-41,46],[14,83],[66,39],[7,105]],[[21862,2183],[48,7],[25,54]],[[21937,2245],[20,98],[43,63],[11,45],[40,46],[38,3],[130,-16],[41,-51],[60,-97],[27,-62],[-4,-46],[41,-45],[31,-63],[101,-22],[71,6],[173,36],[68,-44],[133,18],[12,38],[151,106],[105,22],[69,-15]],[[59269,33036],[-59,30],[-106,3],[-43,16],[-58,58],[-68,32],[-84,72],[-24,71]],[[23397,33818],[-15,77]],[[23384,33979],[53,49],[135,-8],[114,-71],[80,-26],[158,-22],[40,-60],[34,-16],[90,-2],[155,96],[60,23],[90,3]],[[18230,20875],[100,71],[62,-57],[49,11],[18,28],[41,132],[47,41],[24,38],[45,39],[-2,34],[36,53],[12,47],[29,19],[75,-58],[33,19],[57,9],[20,22],[46,3],[105,70],[59,-5]],[[23315,3750],[87,-12],[86,-28],[180,-128],[85,-44],[91,-24],[42,-27],[61,-90]],[[22076,4903],[18,-99],[35,-73],[48,-16],[24,-34],[22,-97],[-2,-81],[-16,-53],[-12,-115],[41,-41],[92,-15],[35,-16],[48,-47],[53,-13],[44,-30],[114,-27],[41,-20],[8,-53],[-25,-67],[60,-63],[105,-18],[36,-54],[79,-42],[18,-35],[33,-18],[85,-5],[96,-22],[90,-7],[68,8]],[[11748,27205],[-16,-52],[-40,-47],[-15,-42],[-33,-39],[-7,-38],[61,-52]],[[11877,27706],[-38,-36],[-42,-78],[-3,-97],[-42,-62],[-1,-82],[55,-84],[-26,-58]],[[11780,27209],[-32,-4]],[[11778,28067],[-5,-146],[54,-120]],[[11743,28284],[12,-32],[10,-94]],[[12176,28371],[-76,-51],[-28,-42],[0,-53],[-23,-24],[-49,23]],[[12502,28585],[-26,-27],[-47,-82],[-20,-61],[12,-68],[-15,-41],[-60,7],[-41,39],[-40,10],[-9,49],[-81,-40]],[[13808,29645],[-62,-24],[-37,0],[-46,-24],[-26,-31],[-65,-7],[-60,-85],[-24,-13],[-116,-8],[-49,-22],[-69,-65],[-17,-44],[-56,-30],[-57,31],[-35,5],[-57,-54],[-39,-98],[-63,-62],[-35,-87],[-80,-104],[-55,-21]],[[13811,29723],[6,-71]],[[10020,31425],[-116,11],[-18,16]],[[9886,31452],[-27,99],[-23,45],[-65,13]],[[10946,32290],[-15,-48],[-31,-29],[-2,-48],[-67,25],[-45,1],[-46,50],[-40,6],[-66,-20],[-82,9],[-2,-45],[-23,-38],[-57,-39],[-14,-68],[48,-52],[14,-36],[16,-144],[24,-34],[97,1],[25,-36],[26,-99],[38,-20]],[[10744,31626],[25,-47],[-58,-39],[-58,-3],[-117,-32],[-98,-37],[-161,-28],[-83,24],[-111,-7],[-63,-32]],[[11007,32654],[43,-3],[9,-104]],[[11059,32547],[-79,-82],[-12,-33]],[[11013,32904],[-49,-99]],[[11481,32964],[-50,100],[-66,98],[-92,75],[-98,116],[-56,28],[-22,33]],[[11015,33483],[-115,72]],[[10863,33560],[-13,-17],[-21,-117],[48,-88],[41,-48],[19,-92],[70,-104]],[[44571,13762],[77,35],[4,31],[42,68],[-31,47],[6,63],[33,98],[10,78],[-21,108],[2,79],[-15,81],[-1,117],[7,68]],[[44071,16731],[-27,47],[-40,31],[-63,34],[-64,71],[-64,111],[-50,64],[-54,25],[-42,49],[-37,20],[-47,6],[-110,-7],[-97,-38],[-84,-15],[-122,-8],[-120,-35],[-177,-95],[-99,-115],[-65,-92],[-22,-83],[-37,-206]],[[42650,16495],[-17,-61],[-83,-140]],[[42550,16294],[-49,-74],[-61,-73],[-71,-46],[-29,-44]],[[42340,16057],[-72,-186],[-75,-111],[-6,-33],[6,-246],[-6,-54]],[[42187,15427],[-10,-88],[-37,-99],[-27,-56],[-63,-31],[-48,-71],[-30,-20],[-37,-51],[-76,-125],[-34,-13],[-33,11],[-60,-20],[-42,-26],[-44,-94],[-56,-8],[-19,-19],[12,-40],[-19,-86],[-49,-45],[-54,5]],[[5019,36961],[-27,56],[-99,17],[-124,-43],[-140,-92],[-11,-70],[50,-51],[166,-99],[-8,-53],[-56,-29],[-33,-70],[36,-66],[42,-39],[-36,-53]],[[12563,37652],[-56,-7],[-119,25]],[[12306,37776],[-6,-20]],[[11875,37998],[-60,74],[26,76],[-13,45],[-87,116],[-120,46],[-246,53],[58,110],[110,4],[-22,67],[55,82],[76,42]],[[42176,12325],[34,-23],[9,-32],[24,-202],[-4,-90],[-30,-159],[2,-97],[35,-49],[29,-19],[18,-39],[24,-114],[59,-89],[87,-38],[19,-52],[29,-46],[18,-55],[24,-130],[26,-72],[51,-104],[114,-202],[31,-44],[55,-33],[98,-7],[124,-18],[66,6],[46,-8],[82,-43],[86,21],[40,-3],[50,-32]],[[25725,11561],[58,29],[52,45],[39,-10],[72,65],[34,80],[32,124],[41,104],[27,101],[37,106],[11,70],[-18,64],[-40,82],[-48,65],[-14,38],[-33,42],[-16,56],[36,100]],[[61690,27297],[40,-54],[31,-20],[31,-49],[54,-47],[54,-90],[6,-66],[51,-57],[48,-12],[63,13],[43,-5],[37,-58],[6,-42],[-13,-86],[-28,-50],[-1,-59],[41,-88],[32,-158],[4,-64],[-8,-120],[6,-57],[18,-41],[27,-12],[50,23],[61,-36],[112,24],[69,-8],[-22,37],[37,45]],[[69078,11652],[-22,12],[-16,59],[-26,48],[-72,89],[-15,63],[-16,14],[-69,-25],[-35,63],[-62,87],[-19,6],[-16,44]],[[42243,36200],[16,3]],[[42472,36316],[10,44]],[[41904,36388],[31,-56],[55,-5],[62,-36],[92,-14],[25,-71],[49,-17]],[[41612,36613],[12,-42],[42,-6],[140,-88]],[[41234,36854],[46,-10],[28,-39],[186,-110]],[[45668,26362],[-19,56],[39,98],[-5,81],[45,85],[59,49],[38,81]],[[42423,32234],[64,-56],[65,-39],[90,-2],[141,8],[57,10],[29,-36],[-15,-75]],[[42854,32044],[-24,-38],[-17,-70],[12,-39]],[[42825,31897],[7,-114],[-21,-32],[10,-103],[-44,-85],[26,-58]],[[42803,31505],[26,-92]],[[42829,31413],[12,-51],[63,-40],[99,-15],[48,9],[87,-44]],[[43138,31272],[32,-86],[51,-12],[22,-44],[42,5],[104,39],[64,-38],[63,-15],[43,-40],[60,44],[35,9],[51,-33],[-51,-38],[7,-36],[57,-60],[59,-39],[-37,-42],[12,-38],[71,13],[71,-6],[162,-30],[47,10],[169,-35],[79,-5],[75,31],[26,34],[66,45],[64,25],[138,27],[42,-9]],[[44762,30948],[34,0],[104,41],[51,51],[-20,80],[-27,27]],[[44904,31147],[-11,17],[32,155],[26,48],[31,11],[56,-51],[72,-20]],[[72699,6596],[-36,1],[-49,31],[-92,33],[-45,-12],[-97,-7],[-31,-39],[-51,-39],[-52,-61],[-72,-23],[-35,-48],[-23,-12],[-42,-58],[-74,-41],[-9,-23],[-61,-46],[-21,-52],[-141,-80],[-26,-53],[-60,-54],[-15,-61],[-102,-80],[-31,-56],[-2,-147],[-16,-63],[9,-31],[47,-53],[-28,-76],[-104,-77],[-17,-82]],[[45881,35135],[56,-51],[-11,-46],[-43,-71],[-4,-66],[-45,-41],[-92,-9],[-24,-50],[9,-76],[-57,-18],[-27,-30],[-96,-24],[-100,-50],[-15,-29],[-91,-3],[-74,-25],[-28,15],[-93,98],[-45,23],[-148,31],[-54,65],[-92,-5],[-103,32],[-92,-7],[-48,28],[-17,51],[-64,97],[-83,64],[-124,-5],[-66,24],[-90,52],[-66,22],[-38,56]],[[22759,5851],[30,-90],[15,-207],[28,-94],[70,-165],[-15,-83],[13,-69],[34,-78],[21,-86],[-3,-95],[-14,-146],[-20,-94],[-27,-43],[-9,-61],[9,-133],[-7,-62],[10,-86],[75,-81],[42,-16],[26,-57],[81,-63],[0,-60],[-22,-61],[66,-67],[74,-57],[79,-47]],[[47640,29072],[87,-23],[51,-84],[4,-62],[41,-23],[125,-10],[43,-14],[113,-4],[38,-12],[19,-88],[49,-15],[51,-54],[15,-55]],[[71134,14882],[12,-57],[-32,-73],[-65,-93],[-22,-45],[3,-51],[17,-43],[-11,-63],[6,-55],[22,-97],[-15,-44],[-38,-42],[-15,-61],[-21,-29],[-35,16],[-33,-33],[-44,-16],[-80,21]],[[60362,26085],[39,43],[55,16],[32,30],[2,36],[38,20],[23,41],[29,-8],[23,-53],[-2,-92],[-74,-69],[-55,-87],[62,-140],[-18,-62],[1,-38],[62,-24],[25,-41]],[[46267,34796],[35,-53],[-28,-46],[-43,-24],[-3,-69],[-42,-121]],[[45509,31592],[3,-4]],[[43734,32630],[55,63],[71,31],[51,4],[124,-43],[33,-47],[73,-36],[60,-51],[31,-10],[61,-50],[63,-11],[35,-40],[52,-26],[154,-44],[65,35],[37,-7],[65,13],[237,-118],[8,-26],[49,-18],[8,-25],[84,-43],[18,-34],[1,-70],[32,-114],[48,-55],[41,-87],[37,-45],[101,-79],[32,-33]],[[45997,31760],[114,47],[110,24]],[[46585,32162],[-7,79],[28,38],[-20,64],[-163,81],[-80,55],[-163,85],[-40,9]],[[45837,32724],[-99,66]],[[45733,32870],[-86,8],[-85,65],[-32,64],[-12,59]],[[46186,34483],[-96,27],[-56,-28],[-52,-73],[-118,28],[-95,-24],[-204,-30],[-71,-39],[-26,-55],[24,-140],[-30,-32],[2,-92],[20,-74],[-20,-49],[-45,-15],[-13,-48],[25,-80],[59,-109],[93,-125],[-57,-111],[14,-108],[-24,-37]],[[48213,32084],[-22,-14],[-141,13],[-78,-17],[-34,13],[-78,-31],[-71,-68],[-40,-21],[-98,-11],[-92,-27]],[[48548,32708],[-39,-6],[-50,-52],[-77,47],[-15,33],[-246,14],[-103,-2],[-30,53],[-51,19],[-17,42],[-40,-4],[-40,58],[-42,-19],[-75,23],[9,112],[-65,85],[13,32],[-16,43],[-112,8],[-2,51],[-47,94],[-4,85],[22,106]],[[47521,33532],[8,21],[-41,54],[-20,63],[2,86],[42,78],[15,58],[-1,71],[-31,75],[-88,41],[-26,57],[-52,67]],[[40333,32230],[63,-39],[48,3],[66,37],[74,9],[49,41],[64,28],[41,33],[77,31],[98,64],[73,18],[153,20],[19,25],[111,48],[61,-31]],[[42423,32234],[-41,35],[-86,-2],[-40,17],[-65,55],[-81,0],[-83,12],[-56,-62],[-124,-10],[-76,31],[-67,12],[-56,44],[-184,75],[-40,47],[-94,32]],[[14783,27473],[56,-46],[40,-66],[34,-20],[64,-12],[77,-65]],[[41636,31916],[127,19],[33,20],[56,1],[92,62],[36,-3],[14,-89],[63,-126],[95,-24],[40,-31],[83,-27],[70,-48],[38,-50],[87,-77],[46,-13],[23,-28],[60,-19],[74,-2],[100,-62],[56,-6]],[[37762,30223],[92,-43],[20,-63],[-33,-27],[-101,-4],[-111,40],[-123,20],[-48,-20],[-52,3],[-72,-16],[-92,-34],[-54,-8],[-49,21],[-86,-18],[-67,25],[-53,-57],[-66,-53],[-20,-35],[-55,-9],[-54,36],[-95,-3],[-135,-29],[-45,-2],[-58,23]],[[37448,30579],[40,-50],[61,-14],[30,26],[41,-30],[102,-16],[44,-51],[42,-28],[57,-1],[93,-40],[137,-135],[69,-58],[84,-59],[104,-103],[80,-1],[29,34],[45,-26],[1,-24],[52,-12],[15,-33],[-27,-45],[3,-59],[24,-19],[56,-7]],[[45659,19567],[52,-22],[40,4],[53,92],[67,58],[47,78],[8,57],[-6,41],[-33,46],[95,130],[31,64],[6,102],[21,125],[-7,130],[12,98],[-4,109],[7,60],[-37,89],[-12,48],[0,70],[-50,82],[-39,124],[-7,66],[21,104],[25,39],[16,75],[12,165]],[[46643,19805],[-26,118],[-22,18],[8,132],[-5,54],[-39,36],[-50,16]],[[46416,20338],[-13,87],[4,86],[-33,55],[-10,82],[-36,57],[-16,117],[-56,40],[-23,63],[4,144],[11,32],[-21,77],[-54,122],[-46,146],[-21,27],[-47,32],[-68,92],[-14,4]],[[42046,33097],[71,-120],[-17,-72],[-32,-18],[-103,-5],[-34,38],[-107,23],[-56,71],[-81,42],[38,60],[-16,55],[-161,89],[-60,93]],[[41489,33353],[-59,49],[-45,72],[-138,11],[-36,11],[-70,53],[15,66],[52,55],[8,109],[-22,48],[-85,45],[-104,31],[-99,72],[-140,23]],[[58248,32213],[-62,-22],[-63,-48],[-86,-17]],[[57806,32753],[-55,13],[-100,87],[-55,24],[-72,65],[-88,3],[-68,32],[-91,19],[-154,-4],[-95,41],[-73,7],[-91,47],[-96,-3],[-128,67],[-41,45],[-29,76],[-49,65],[-74,45],[-16,69],[-27,38],[-57,36],[-64,119],[-9,63],[-27,60],[-46,54],[-115,44],[-81,70],[-38,68],[-94,84],[3,79],[-16,11]],[[59425,32099],[-17,-71],[-65,-53],[-132,-92],[-40,5],[-100,50],[-67,25],[-155,37],[-68,26],[-68,51],[-38,3],[-115,79],[-40,15],[-65,0],[-84,34],[-80,12],[-43,-7]],[[48993,26740],[18,-25],[135,-47],[84,-60],[94,18],[136,3]],[[47954,27764],[29,-13],[37,18],[73,1],[55,20],[38,-20],[30,-48],[39,-3],[59,-132],[48,-23],[5,-48],[61,-56],[144,-39],[42,-51],[46,-22],[52,-109],[-4,-67],[28,-153],[21,-57],[7,-72],[17,-38],[51,-55],[63,-13],[59,-39],[40,-5]],[[42821,31646],[-2,-104],[-16,-37]],[[8838,35105],[72,33],[2,54],[85,110],[8,69],[52,16],[104,-77],[49,-16],[56,-43],[104,-119],[26,-46]],[[47436,29202],[-73,-33],[-22,-43],[59,-43],[39,-5],[69,-31],[38,-78],[-61,-71],[-38,-79],[-48,-17],[-39,-33],[-54,-21],[-56,-1],[-26,-44],[2,-58],[30,-39],[8,-53]],[[48061,29650],[-29,-54],[-80,-23],[-137,-1],[-52,-16],[13,-44],[-24,-63],[-119,32],[-148,-23],[-114,-49],[-14,-33],[68,-66],[5,-21]],[[67988,10774],[-58,-31],[-117,7],[-38,-7],[-59,-75],[-45,-92],[-65,-95],[-86,-57],[-51,-13],[-43,12],[-36,33],[-29,56],[-44,36],[-42,-15],[-39,15],[-81,135],[-2,46]],[[10911,33697],[-66,99],[-79,52],[-112,9],[-57,15],[-134,95],[-101,87],[-80,55],[-131,50],[-73,69],[-91,33],[-51,-8],[-44,-88],[-10,-44],[8,-136],[-22,-44],[38,-37],[39,-72],[0,-96],[44,-125],[3,-88],[-32,-76],[33,-79],[8,-57],[87,-145],[-14,-67],[40,-68],[55,-145],[14,-84],[1,-85],[-10,-48],[-49,-46],[-111,-35],[-62,1],[-95,21],[-53,-25]],[[65025,24576],[28,10],[46,47],[4,56],[35,47],[8,77],[37,104],[-5,39],[-30,5],[-7,37],[5,95],[-45,70],[-5,39],[19,35],[46,42],[-1,76],[20,38],[-1,65],[13,38],[56,90],[17,59],[31,47],[44,35],[19,33],[-7,88],[37,99],[9,54],[-8,75]],[[59229,23779],[-14,72],[-29,41]],[[59186,23892],[-10,-57],[28,-127],[-5,-111]],[[58986,24012],[-22,-82],[0,-35],[-53,-72],[-7,-63]],[[58986,24012],[74,-23],[35,-84],[33,-126],[1,-60]],[[59366,24204],[13,-105],[23,-24],[31,-71],[17,-69]],[[59142,24342],[25,-61],[66,-77],[61,-181],[38,-59],[19,-51],[-5,-72]],[[59142,24342],[90,-49],[59,1]],[[58976,24472],[26,-42],[30,-82],[53,-51],[-6,-62],[-15,-22],[9,-44],[41,-43],[29,-96],[9,-76],[34,-62]],[[58976,24472],[76,-37],[11,-43],[36,3],[43,-53]],[[58865,24615],[5,-52],[54,-91],[15,-51],[-6,-51],[4,-81],[-10,-43],[41,-79],[-2,-41],[31,-66],[-11,-48]],[[58865,24615],[19,-54],[34,-57],[33,-2],[25,-30]],[[56987,26603],[5,-69],[-28,-44],[-63,-9],[-56,-71],[-68,-2],[-63,-30],[-34,-37],[-56,1],[-26,-66],[3,-40],[-33,-53],[-4,-38],[19,-102],[-1,-98],[43,-173],[34,-68],[99,-83],[72,-35],[21,-35],[44,-19],[34,-51],[33,-82],[44,-30],[40,-94],[33,-62],[71,-109],[9,-31],[41,-51],[51,-4],[34,-22],[42,-47],[36,-65],[69,-48],[67,-61],[42,-23],[92,-8],[39,-25],[41,15],[41,93],[23,15],[58,-24],[7,29],[35,24],[31,-29],[36,20],[50,47],[60,8],[47,-12],[45,12],[68,-21],[51,-45],[51,0],[53,-16],[73,-57],[32,-6],[77,47],[11,-43],[29,-12],[52,13],[35,-19],[77,44],[78,-66],[10,-43],[28,-35],[-6,-43]],[[38843,30452],[-4,35],[-41,41],[-7,95],[-30,88]],[[47438,18640],[-3,-50],[53,-107],[28,-23],[69,-34],[17,-51],[41,-29],[10,-44],[57,-48],[44,-20],[72,11],[43,29],[76,70],[40,-14],[62,-186],[28,-55],[17,-64],[23,-42],[65,-91]],[[49074,27132],[21,-23],[22,-98],[-6,-42],[28,-76],[15,-103],[58,-113],[-9,-49]],[[12555,27330],[-48,-21],[-45,4],[-58,24],[-98,71],[-45,10],[-59,-44],[15,-50],[-25,-18],[-64,12],[-172,-94],[-85,-28],[-55,26],[-36,-13]],[[13263,27485],[21,-54],[-31,-78],[-51,-32],[-36,-82],[-47,-55],[-57,17],[-59,82],[-63,-36],[-59,26],[-46,55],[-42,28],[-40,3],[-69,-45],[-28,-6],[-52,27],[-51,-3]],[[41092,37031],[-69,-63],[-123,-82],[10,-61],[44,-36],[-6,-82],[49,-91],[57,-62],[20,-67],[40,-59],[105,-107]],[[57428,22030],[27,-128],[-34,-50],[1,-21]],[[57428,22030],[35,-54],[34,-21],[71,-2]],[[55548,23037],[55,-8],[45,16],[79,-39],[76,-63],[104,-48],[64,-51],[48,-21],[39,5],[57,-28],[53,-11],[4,-38],[46,-14],[41,-37],[68,39],[44,-5],[29,15],[61,-55],[48,-31],[42,-1],[62,37],[80,-5],[15,14],[72,18],[23,-10],[35,-54],[46,-9],[76,-50],[46,39],[51,-35],[41,-5],[21,-73],[30,-67],[35,-31],[35,-55],[13,-92],[21,-31],[65,-26],[75,-43],[14,-43],[21,-111]],[[28087,9209],[34,36],[45,74],[12,45],[-26,53],[-51,19]],[[28048,9484],[-16,27],[-53,35],[-49,-13],[-19,25]],[[27651,9697],[-61,15]],[[27477,9823],[-90,93],[-46,-7],[-37,14],[-11,-25],[-44,-21],[-29,12],[-42,-19],[-105,-6],[-18,-69],[-35,60],[-37,-35],[-20,48],[1,53],[-25,2]],[[26712,9963],[-16,19],[-35,-7]],[[23597,10235],[33,-14],[16,31],[82,72],[58,109],[21,75],[19,130],[-3,66],[12,117],[-43,39],[-32,79],[-47,49],[-79,48],[-47,44],[-15,42],[-52,43],[-49,54],[-84,55]],[[23001,10711],[71,-94],[13,-29],[39,-19],[94,-119],[46,-83],[33,-2],[43,42],[53,-39],[47,1],[22,-22],[23,-55],[34,-43],[26,-59],[33,45],[21,0]],[[19713,28355],[-4,-45],[24,-76],[-19,-29],[44,-39],[46,-83],[3,-74],[34,-26],[0,-82],[24,-64],[3,-42],[47,-77],[29,-116],[32,-49],[50,-51],[-25,-73]],[[23301,12568],[-13,56],[-54,99],[-16,64],[15,72],[-2,43],[-25,101],[11,64]],[[41290,35616],[-38,-44],[-33,-102],[-69,-54]],[[65830,32711],[43,-64],[89,38],[71,-73],[113,-5],[51,25],[60,0],[95,-16],[93,10],[91,37],[167,87],[58,15],[95,56]],[[64343,27154],[29,-33],[41,-107],[32,-25],[43,-9],[20,-24],[9,-64],[38,-37],[12,-53],[-2,-61],[29,-29],[-8,-60],[26,-62],[-7,-44],[26,-45],[60,-15],[55,-33],[70,4],[37,50],[68,20],[31,-23],[44,-5]],[[63162,27318],[64,24],[37,-2],[45,-30],[38,6],[94,45],[49,-7],[91,-50],[21,-71],[29,-24],[23,-61],[33,-14],[25,30],[38,12],[32,33],[46,7],[53,42],[65,16],[87,-36],[78,6],[35,-9],[42,11],[45,-35]],[[68089,32939],[-45,29],[-21,51],[1,83],[-89,129],[-28,133],[-24,65],[-45,73],[-12,65],[-33,90],[-70,69],[8,43],[-41,43],[-34,6],[-18,50],[-134,52],[-60,-22],[-81,31],[-41,37],[-69,21],[-90,46],[-105,0],[-149,-29],[-41,9],[-217,-55]],[[52746,26707],[35,-63],[15,-69],[6,-97],[44,-64],[81,-31],[71,-8],[58,15],[60,34],[61,53],[48,22],[57,-11],[37,60],[23,116],[38,96],[14,95],[96,112],[30,60],[61,46],[32,54],[33,97],[34,60],[74,51],[41,68],[113,112],[56,20],[91,85],[8,55],[54,49],[55,29],[111,39],[38,-2],[26,21]],[[19491,28610],[-272,-131],[-18,-28],[-63,-40],[-105,-47],[-84,-60],[-51,-67]],[[63190,23128],[-54,16],[-32,77],[-29,30],[-44,126],[-26,24],[-93,49],[-48,52],[-33,57],[-66,89],[-149,180],[-133,156],[-74,72],[-127,76],[-115,53],[-53,61],[-68,114],[-32,94],[-26,56],[-26,11],[-99,152],[-40,6],[-67,42],[-39,97]],[[63289,24746],[17,-54],[32,-17],[28,-44],[20,-89],[57,-120],[-26,-49],[12,-33],[29,-8],[53,-79],[59,40],[44,50],[49,4],[47,-52],[18,-4],[20,-57],[63,19],[21,-11],[44,52],[14,-16]],[[61044,27898],[60,13]],[[61154,27980],[40,-1],[61,-69],[43,-25],[-8,-67],[28,-49],[45,-25],[34,35],[45,-15],[13,-42],[1,-72],[14,-54],[47,-55],[27,21],[65,-7],[79,-53],[31,-9],[32,51],[95,56],[47,-50],[42,17],[42,-33],[84,-16],[34,-22],[45,-73],[44,-6],[67,37],[-30,51],[-34,85],[-53,29],[-68,52],[-83,42],[-38,58],[-43,29],[-24,45],[-57,11],[-60,39],[-15,36],[-33,160],[20,88],[68,57],[46,60],[62,31],[77,-26],[67,28],[61,-21],[40,-66],[28,-14],[99,6],[79,-18],[53,2],[28,21],[-10,53],[31,33],[40,6],[101,-27],[18,50],[65,66],[84,65],[-5,33],[-85,71],[6,75],[79,61],[71,34],[152,19],[59,46],[31,83],[37,49],[33,72],[27,31],[43,132],[54,97],[-26,134],[2,90],[-12,84],[100,125],[21,64],[84,48],[34,33],[76,12],[94,-32],[68,-67],[126,-27],[57,13],[65,-13],[57,11],[52,-13],[86,-60],[75,-5],[30,-17],[65,-73],[-6,-73],[14,-45],[-15,-54],[-48,-29],[7,-52],[-44,-80],[-28,-144],[-67,-98],[-22,-48],[18,-76],[46,-88],[-27,-59],[-13,-71],[-42,-64],[-10,-65],[23,-129],[-6,-63],[7,-76],[20,-89],[1,-66],[-50,-98],[-20,-74],[-2,-122],[46,-13],[76,16],[71,49],[82,19],[91,73],[64,-17],[102,-49],[70,-21],[69,22],[56,32],[90,-23],[87,11],[74,-9],[48,6],[34,30]],[[60782,27994],[136,-25],[75,-40]],[[65124,27958],[16,68],[25,34],[76,56],[29,59],[94,65],[58,86],[107,97],[36,71],[68,39],[76,101],[79,57],[31,1],[70,29],[40,64],[150,92],[21,-5]],[[21199,30935],[-24,-78],[67,-66],[13,-50],[-8,-41],[20,-35],[38,-14],[-1,-76],[-25,-94],[-10,-105],[-40,-123],[-9,-146],[3,-106]],[[61211,32214],[-91,25],[-38,21],[-2,97],[62,78],[79,21],[115,-3],[88,-10],[59,-20],[95,29],[206,115],[33,7],[20,55]],[[27611,9811],[41,-10],[83,11],[16,47],[59,35],[34,124],[49,133],[38,40],[7,36]],[[72071,37056],[-63,80],[-221,112],[-243,82],[-124,116],[2,60],[55,45],[146,9],[17,65],[95,44],[16,78],[-59,86],[29,83],[-133,24],[10,76],[108,106],[18,96],[-36,87],[-73,78],[-42,173],[19,37],[231,60],[85,-10],[148,11],[-33,48],[120,28],[50,41],[28,60],[116,41],[26,64],[55,26],[170,-2],[48,75],[79,46],[-21,116],[57,48],[-50,75],[22,110],[-30,20],[-2,68],[37,50],[126,69],[146,31],[54,35],[144,29],[125,99]],[[72070,37057],[69,-85]],[[61802,16289],[26,-2],[41,-34],[39,6],[15,22],[44,-15],[52,1],[34,40],[33,8],[42,-17],[34,1],[64,59],[35,-7],[74,-43],[90,58],[41,-10]],[[54201,24312],[15,33]],[[54151,24318],[16,48]],[[54216,24345],[26,-51]],[[54216,24345],[-49,21]],[[54167,24366],[-12,13]],[[56962,27113],[-71,99]],[[56893,27210],[-22,47],[-70,82],[-94,25],[-91,62],[-22,47],[-145,185],[-131,60],[-71,69],[-63,46],[-39,49],[-3,44],[-50,62],[-30,61],[-69,20],[-54,66],[-117,44],[-41,43],[-5,-74],[-14,-33],[-54,-28],[-66,7],[-96,32],[-51,2],[-29,-22],[-5,-46],[-24,-42],[-44,-36],[-27,-50],[17,-43],[-30,-27],[-8,-34]],[[55328,27651],[-13,-17],[-82,-24],[-15,-49],[-45,-17],[-29,-83],[-36,-43],[-10,-86],[-16,-38],[-52,-52],[-17,-48],[-2,-68],[-16,-56],[-48,-103],[-48,-141],[-14,-88],[6,-125],[-13,-35],[12,-69],[-7,-128],[8,-68],[-8,-59],[-32,-90],[11,-64],[-21,-108]],[[54841,25992],[-50,-36],[-66,-78],[-93,-76],[-22,-43],[-53,-46],[-77,-37],[-30,-35],[-24,-91],[-28,27],[-86,-51],[-5,-60],[-36,-39],[3,-31],[-61,-97],[-3,-178],[39,-58],[73,-130],[4,-57],[-17,-58],[-3,-72],[10,-22],[-20,-49],[-47,-45],[-19,-107],[-31,-49],[-15,-51],[-37,-53],[-40,-12]],[[60443,21889],[9,130],[37,84],[17,62],[34,39],[43,96],[34,30]],[[60617,22330],[28,15],[49,-152],[27,-62],[12,-58],[32,-74],[23,-26]],[[60623,21988],[-36,-159],[-20,-64],[-24,-36],[-13,-61]],[[60651,22017],[12,-44],[51,-57],[-6,-38],[46,-45]],[[60621,22293],[-27,-26],[-15,-53],[12,-24],[-36,-115],[-3,-57],[11,-36],[-18,-61],[-30,-35],[-25,-71],[-58,-98]],[[55860,34177],[-101,136],[-101,61],[-92,35]],[[55564,34409],[-94,64],[-34,100],[-35,37],[94,70],[105,56],[92,16],[46,31],[56,74],[-7,67],[16,44],[-55,162],[-77,41],[-31,48],[-42,1],[-224,96],[-29,-5],[-103,40],[-90,5],[-80,18],[-48,27],[-68,-2],[-158,49],[-13,31],[-74,17],[-94,-16],[-52,19],[-141,-1],[-83,44],[-55,7]],[[54434,36518],[35,-35],[81,-19],[110,-111],[-13,-79],[22,-40],[-33,-29],[29,-37],[8,-63],[-44,-40],[-47,-5],[-50,-50],[-83,-41],[-27,-83],[16,-12],[-40,-122],[-32,-12],[-5,-75],[-34,-8],[-42,-56],[0,-53]],[[71688,30804],[-42,36],[-59,28],[-88,-51],[-72,-30],[-19,-72],[-35,-33],[-34,-65],[-48,36]],[[45685,26351],[33,52],[61,68],[53,21],[102,14]],[[22271,16071],[3,-71],[26,-52],[49,-49],[50,-21],[106,23],[62,25],[52,11],[78,-15],[63,-37],[51,-18],[78,10],[102,38],[89,-3],[76,-45],[49,-17],[41,0],[18,-52],[-8,-25],[37,-109],[20,-78],[16,-8],[9,-63],[23,-60],[-3,-37]],[[28329,10476],[-35,92],[19,57],[-31,17],[-13,58],[22,34],[48,0],[35,55],[39,41],[33,55],[30,28],[18,46],[34,34],[48,21],[39,0],[31,-21],[29,4],[54,46],[83,13],[33,29],[34,2],[76,79],[45,24],[82,-3],[77,45],[76,15],[44,-33],[37,29],[33,44],[32,-7]],[[61534,25695],[8,-91],[39,-136],[29,-59],[54,-83],[25,42],[39,91],[-3,51],[18,54],[31,11],[10,-88],[17,-30],[-19,-63],[5,-50],[-18,-82],[12,-185],[79,-18],[37,10],[56,48],[35,8],[56,70],[58,-13],[17,-31],[3,-65],[-19,-53],[16,-39],[49,25],[48,69],[58,33],[13,-38],[27,3],[56,44],[11,42],[-33,160],[1,99],[16,32],[41,8],[12,32],[49,65],[29,62],[-21,48],[34,60],[40,-2],[23,44],[-12,61],[9,34]],[[61335,26866],[60,-83],[-1,-58],[-37,7],[4,-51],[34,-53],[27,-65],[6,-95],[18,-123],[5,-83],[-6,-98],[19,-98],[11,-200],[18,-92],[41,-76]],[[46659,26746],[-5,-31]],[[46695,27204],[-7,-107],[5,-51],[-9,-83],[7,-72]],[[46761,27491],[-39,-54],[-23,-105],[0,-69]],[[48180,17892],[24,-39],[63,-41],[35,-45],[12,-48],[-10,-53],[8,-77],[-14,-39],[24,-69],[-42,-76],[-22,-91],[4,-57],[-37,-65],[32,-64],[2,-52],[28,-72],[0,-50],[31,-80],[11,-48],[4,-87],[27,-74],[-2,-147],[-34,-92],[12,-34]],[[24687,11709],[8,93],[-9,69],[8,19],[5,143],[14,84],[6,95],[4,165],[13,103],[30,86],[8,86],[18,82],[-14,57],[5,84]],[[41821,18398],[0,-122],[-29,-68],[-2,-58],[-38,-55],[-15,-48],[-3,-65],[13,-55],[28,-45],[47,-40],[48,20],[30,-37],[44,-9],[21,-24],[31,12],[37,-16],[29,-82],[53,-35],[35,-38],[18,-44],[11,-97],[22,-43],[-36,-122],[20,-108],[-26,-82],[15,-73],[-9,-47]],[[44737,12662],[6,-74],[43,-99],[47,-74],[20,-16],[38,3],[42,-18],[31,-37],[21,-55],[3,-62],[-12,-68],[-29,-62],[-70,-22],[-37,-54],[-9,-65],[-41,-79],[5,-50],[-21,-38],[-46,-34],[-108,-19],[-36,-25],[-48,-57],[-33,15],[6,-64],[-45,-22],[5,-30],[-9,-95],[5,-97],[14,-56],[55,-7],[39,28],[41,0],[68,31],[37,-3],[34,-30],[34,-7],[70,23],[51,-4],[64,-39],[73,-15],[74,-42],[40,1]],[[45303,15496],[55,-21],[45,-68],[-3,-36],[61,56],[43,66],[-3,73],[57,71],[33,71],[11,64],[-14,57],[-1,133],[-24,55],[-18,9],[-38,78],[12,22],[67,22],[73,-67],[66,-13],[39,27],[16,50],[29,23]],[[75075,34277],[27,-69],[31,12],[105,147],[81,56],[41,43],[56,91],[49,38],[11,49],[58,18],[26,36],[3,55],[-20,65],[65,88],[86,62],[90,-4],[255,-38],[61,-15],[97,4]],[[64805,16868],[-54,-77],[-65,-36],[-68,-7],[-95,-53],[-52,-12],[-52,-86],[-37,-21],[-34,-68],[-98,-3],[-57,-28],[-59,17],[-37,-41],[-56,3],[-40,-13],[-46,-57],[-50,-25],[-33,3],[-71,50]],[[43593,12785],[-14,118],[24,35],[5,127],[-24,117],[-26,58],[-54,88],[-10,42],[4,98],[18,156],[-2,115],[-22,76],[-5,80],[14,146],[-35,102],[-34,48],[-79,89],[-62,52],[-37,55],[-46,123],[13,37],[46,54],[47,123],[9,78],[-5,99],[-21,61],[-36,25],[-67,74],[-28,3],[-27,-29],[-89,75],[-115,66],[-155,109],[-54,9],[-51,47],[-37,15],[-78,7],[-102,62]],[[42458,15425],[-69,46],[-73,-9],[-53,6],[-55,-9],[-21,-32]],[[69378,12114],[16,-41],[6,-57],[-22,-46],[-38,10],[-31,-31],[-33,-122],[-91,-79],[-86,-94],[-21,-2]],[[65107,17052],[26,41],[90,41],[41,5],[26,-19],[33,30],[21,42],[41,-40],[42,26],[13,46],[30,59],[77,109],[23,-22],[53,26],[43,-18],[40,18],[28,37],[25,2]],[[45763,37418],[61,-34],[94,29],[28,43],[7,63],[-40,49],[29,57],[-44,81]],[[45873,37612],[36,-25]],[[45898,37706],[138,88],[102,29],[191,-42],[89,-10],[97,13]],[[45852,37744],[46,-38]],[[44784,38270],[-40,-39],[-182,15],[-162,33],[-72,-43],[-22,-56],[-118,-48],[-39,-61]],[[45115,38807],[34,-30],[-48,-107],[-101,-74],[-53,-76],[-130,-41],[3,-52],[-30,-59]],[[62870,40339],[-66,-9],[-226,-94],[-188,-37],[-111,-92],[-122,-38]],[[59865,39424],[235,-20],[180,-27],[101,6],[88,24],[216,91],[30,47],[-7,53],[-107,15],[-53,35],[-100,-6],[-135,31],[52,60],[115,6],[266,-33],[316,-11],[109,4],[71,52],[38,68],[185,-6],[96,85]],[[62157,40069],[-114,-45],[-117,-21],[-273,-31],[-73,-22],[-18,-51]],[[45484,17240],[16,107],[-18,54],[-19,17],[-63,15],[-77,52],[-25,4],[-18,51],[-58,43],[-25,51],[-24,20],[-51,-11],[-53,10]],[[9641,30010],[-49,80],[-24,6]],[[10099,30329],[-8,-37],[-48,-6],[-48,-42],[-88,-41],[-41,9],[-43,-22],[-52,13],[-28,-18],[-40,-92],[-5,-71],[-33,-29],[-24,17]],[[43672,36547],[-122,43],[-106,77]],[[72188,36758],[90,14],[114,-39],[58,-62],[73,-45],[107,14],[59,31],[10,92],[-36,57],[19,74],[88,9],[122,-73],[65,10],[54,-147],[63,14],[110,-9],[36,73],[95,59],[85,4],[66,36],[30,155],[132,-17],[77,-54],[58,33],[65,-10],[55,98],[-11,35],[29,125],[75,48],[58,70],[118,76],[92,19],[95,42],[-35,87],[-55,52],[-12,52],[-45,21],[-80,-11],[-205,38],[-150,74],[-97,73],[-157,171],[-30,65],[13,54],[65,24],[69,-2],[64,22],[62,72],[112,47],[-3,28],[148,124],[142,41],[-26,36],[29,39],[283,51],[114,9],[105,-12],[95,15],[28,39],[-47,32],[29,29],[95,-7],[178,50],[52,27],[27,61],[43,32],[2,52],[136,56],[27,-16],[114,14],[13,-37],[213,-42],[169,-1]],[[59206,29798],[-156,19],[-68,18],[-69,38],[-139,37],[-63,0],[-54,31],[-64,74],[-35,7],[-53,-24]],[[43290,18882],[-61,-43],[-18,-35],[-49,-25],[-40,11],[-54,47],[-58,109]],[[3354,38786],[-51,-141],[73,-83],[-15,-54],[-193,-90],[-37,-85],[-103,-34],[-79,5],[-313,-86],[-53,-37],[-74,-90],[-92,-42],[-237,72],[-75,-1],[-72,-93],[38,-63],[-70,-21],[-69,31],[-47,-24],[-110,0],[-8,-134],[34,-19],[-63,-57]],[[56864,21911],[80,36],[78,-4],[16,31],[47,-54],[61,-17],[45,-53],[26,-63],[14,-75],[-6,-62]],[[55581,22374],[39,-45],[36,-72],[1,-69],[26,-87],[52,-47],[18,-52],[29,-30],[-8,-32],[30,-24],[60,13],[26,-31],[25,27],[46,-28],[3,-29],[55,-9],[62,-34],[62,-49],[45,19],[68,47],[60,54],[24,-2],[47,-37],[82,-22],[53,-38],[82,-90],[39,38],[77,-9],[26,6],[45,55]],[[2447,37159],[-29,-20],[-122,-18],[-104,-64],[10,-60],[-58,-39],[-49,-65],[7,-86],[-29,-36],[-127,-59],[-86,3],[-149,55],[-69,-23],[-49,-70],[-74,-14],[-111,9],[-239,-17],[-75,-32],[-26,-45],[-78,-60],[-27,-44],[-88,-54],[-108,-46]],[[42793,12783],[-21,124],[-17,49],[4,40],[-36,65],[0,58],[36,53],[5,59],[-14,50],[-37,36],[-61,123],[-31,21],[-10,103],[-18,53],[-53,99],[-35,46],[-5,64],[-89,202],[-40,67],[-15,106],[-41,103],[-7,48],[1,76],[-26,66],[-21,213],[-1,76],[6,104],[-6,75],[6,98],[12,44],[66,38],[50,20],[23,27],[32,102],[-6,51],[16,28],[-2,55]],[[61209,33367],[-69,-12],[-94,-39],[-180,-4],[-99,12],[-74,-29]],[[22034,33987],[2,28],[-86,-21]],[[21886,34048],[-78,28],[-84,-10]],[[21673,34050],[-109,31],[-44,-4],[-37,28],[-111,49],[-147,-55],[-111,-21]],[[20421,34116],[-48,7],[-108,-18],[-70,-1],[-129,23]],[[60963,26956],[26,-89],[35,-75],[-4,-108],[34,-51],[10,-64],[28,-54],[108,-150],[46,-115],[34,-43],[25,-83],[48,-99],[-6,-38],[18,-78],[32,-96],[19,-121],[34,-112],[-7,-47],[19,-81],[22,-141],[-16,-137],[7,-49],[3,-152],[36,-89],[7,-108],[23,-39],[37,-29],[79,-114],[47,-23],[44,19],[42,-67],[8,-51],[-14,-92],[-23,-27],[-7,-57],[-40,-77],[-27,-112],[9,-41],[29,-42],[43,-127],[35,-44],[15,-76],[55,-114],[58,-36],[7,-76],[-15,-37],[-88,-71],[-21,-63],[-18,-22],[7,-55],[-37,0],[-37,-31],[-28,-117],[44,4],[27,-47],[28,-23],[10,-91],[13,-23],[47,17],[55,0],[91,-18],[50,22],[102,65],[22,-25],[-70,-90],[-25,-43],[14,-76],[-15,-63],[-3,-65],[8,-62],[-16,-74],[-66,-100],[-19,-47],[12,-71],[25,-27]],[[39662,33514],[-39,10],[-132,-32],[-73,7],[-46,25]],[[66475,36282],[-71,-3],[-63,-44],[-110,-3],[-57,-20],[-92,-74],[-65,-16],[-82,-63],[-204,86],[-74,48],[-57,81],[-131,41],[-122,69],[-68,19],[-83,0],[-94,-27],[-100,-106],[-66,-42],[-118,-108],[-65,-43],[-94,-35],[-45,-42],[-23,-60],[-35,-29],[-127,-28],[-122,-84],[-32,-53],[-45,-4],[-43,-33],[-122,-4],[-113,-35],[-24,-59],[-45,-42],[-97,11],[-84,-49],[-91,-95],[-54,-13],[-45,-49],[-8,-37],[-93,-76],[-35,-9],[-61,21],[-28,-19],[-17,-61],[-27,-14],[-114,-10],[-63,-58],[-52,-17],[115,-54],[-13,-68],[-118,-47],[30,-65],[-3,-65],[-46,-78],[-38,-100],[-87,-105],[8,-52],[29,-31],[-11,-98],[111,-142],[33,-5],[58,-49],[203,9],[57,22],[142,112],[75,-31],[-29,-64]],[[67880,40181],[111,-147],[19,-181],[60,-97],[-22,-84],[-145,-66],[-130,-5],[-100,-25],[51,-70],[-57,-55],[-133,-54],[-112,-165],[-146,-100],[-38,-152],[-36,-93],[-108,-84],[-37,-54],[-21,-79],[8,-76],[39,-72],[76,-215],[34,-39],[33,-86],[78,-107],[13,-121],[32,-83],[44,-50],[83,-55],[62,-17],[66,-70],[56,-29],[158,-57]],[[67818,37593],[100,-58],[102,-86],[183,-87],[78,-56],[180,-1],[44,-11],[48,-61],[4,-64],[53,-136],[52,-98],[-37,-78],[-24,-87],[-73,-86],[-169,-101],[-62,-18],[-137,-13],[-112,-28],[-253,-23],[-173,-71],[-258,-59],[-75,-25],[-139,-8],[-292,0],[-263,15],[-108,-49],[-11,-21]],[[67151,30733],[14,-28],[-19,-124],[18,-118],[51,-83],[-14,-26],[-98,-58],[-73,-21],[-33,-42],[-50,-124],[-9,-84],[-42,-36],[-4,-52],[-47,-67],[-23,-53]],[[7876,36544],[12,-43],[86,-6],[77,-63],[120,-54],[169,-127],[81,-41],[29,-32],[83,-43],[89,-31],[64,36],[57,5],[42,-28],[0,-54],[193,-40],[110,-52],[81,0],[85,-44],[49,-4],[83,34],[75,103],[33,-1],[77,-55],[35,50],[22,104],[97,112],[-31,78],[37,94],[60,67],[170,69],[206,57],[56,125]],[[44691,8691],[29,52],[53,50],[43,19],[22,38],[29,11],[50,60],[22,53],[43,59],[39,33],[74,19],[58,46],[15,38],[28,22],[123,22],[56,-13],[67,-35],[63,-12],[106,12],[65,-18],[46,-30],[55,-23],[57,-76],[32,-61],[51,-68],[59,-165],[21,-99],[22,-66],[36,-52],[10,-35],[55,-55],[30,-12],[58,-54],[19,7],[15,-65],[-23,-78]],[[41832,18671],[5,-37],[-21,-63],[-29,-30],[-12,-85],[-24,-38],[-43,-37],[-7,-38],[-35,-49],[-24,-16],[-10,-42],[-49,-24],[-30,27],[-20,-24],[-4,-83],[-30,-38],[-25,-2]],[[62997,35516],[44,-28],[74,14],[78,-39],[122,-34],[104,4],[60,27],[135,103],[38,73],[-21,94],[12,60],[-17,73],[-87,117],[-19,57],[34,52],[-11,128],[-55,45],[36,163],[-21,76],[40,101],[-14,28],[25,56],[-20,38],[11,58],[41,34],[16,98],[-86,36],[33,46],[122,32],[-81,24],[0,29],[73,35],[-101,14],[-154,73],[-42,30],[-104,9],[7,61],[-67,3],[-127,-53],[-175,-1],[-36,33],[-224,34],[-81,19],[-39,53],[-85,37],[-103,67],[-106,-19],[-88,16],[-88,35],[-138,-23],[-202,48],[-52,-39],[-113,-9],[-307,7],[8,-53],[-53,-19],[17,-56],[-150,21],[-205,-3],[-197,68],[-82,7],[-104,37],[-159,4],[-101,20],[-114,46],[-91,22],[-135,-24],[-40,15]],[[58875,38043],[203,64],[88,-22],[187,-81],[-10,-84],[16,-93],[54,-33],[31,-71],[58,-8],[81,-55],[78,-29],[172,-35]],[[44526,12624],[-20,23],[-31,-3],[-34,19],[-24,68],[-42,56],[-11,43],[40,77],[-7,89],[-36,46],[-21,-6],[-9,114],[27,69],[4,38],[-13,74],[17,105],[39,-3],[55,45],[16,47],[-1,67],[45,55],[33,109],[18,6]],[[44684,14635],[17,76],[-17,162],[-12,50],[-53,60],[-59,50],[-77,94],[-12,49],[6,54],[-6,48],[-20,43],[-5,51],[10,111],[-5,72],[-12,19]],[[44439,15572],[-20,97],[-1,112],[12,161],[-5,111],[-31,89],[-34,65],[-26,89],[4,52],[20,59],[6,56],[-15,89],[-18,36],[-64,53],[-80,22],[-75,58],[-41,10]],[[11062,36578],[-9,20],[-101,33],[-130,-55],[-120,-18],[-129,19],[-210,86],[-9,74],[-131,23]],[[10223,36760],[-72,36],[-150,48],[-211,28],[-58,35],[43,131],[-25,132],[-62,55],[-90,167],[-50,53],[-49,16],[-13,76],[-55,33],[-52,64],[-8,59],[-49,39],[-108,23]],[[9214,37755],[-95,19],[-101,73],[-140,39],[-225,102],[-160,24],[-81,88],[-6,53],[99,59],[-173,116],[-127,50],[-54,48],[-27,103],[-93,27],[-58,40],[-112,4],[-263,-64],[-105,-11],[-202,92],[-72,63],[-58,116],[18,72],[57,33],[12,41],[-40,62]],[[7208,39004],[-69,21],[-69,64]],[[7208,39004],[-89,102],[105,77]],[[7070,39089],[-186,-24]],[[7070,39089],[-46,32],[-1,39],[-50,47]],[[23217,13067],[19,45],[1,125],[-25,61],[12,43],[105,-9],[20,55],[46,60],[48,17],[31,41],[48,39],[54,104],[50,21],[46,65],[22,56],[49,3],[57,47],[5,29],[-27,121],[0,74],[37,114],[31,23],[17,36],[43,23],[46,156],[33,19],[43,-6],[83,131],[38,4],[63,69],[58,28],[31,27],[20,51],[60,48],[50,77],[32,4],[58,153],[20,21],[49,17],[30,53],[7,36],[45,56],[11,51],[22,40],[71,84]],[[12672,31141],[-43,-24],[-40,34]],[[12541,31187],[-57,8],[-38,122],[52,150],[20,149],[-25,61]],[[21642,12480],[17,39],[-10,75],[5,30],[43,27],[76,-9],[89,-91],[73,-96],[36,-24],[13,-26],[60,-49],[48,-13],[18,11],[90,-29],[44,-2],[55,35],[55,-7],[57,9],[83,44],[62,96],[18,38],[43,36],[22,63],[30,22],[38,51],[35,-22],[34,11],[38,51],[35,4],[44,41],[87,41],[46,48],[67,47],[37,68],[87,68]],[[20633,17112],[76,-15],[91,72],[23,26],[40,122],[33,67],[20,68],[6,69],[25,88],[43,109],[26,80],[9,52],[-9,74],[21,140],[3,112],[26,140],[1,70],[14,109],[46,100],[15,66],[62,59],[35,75],[-1,101],[14,86],[-1,139],[14,26]],[[21263,19147],[8,59],[-18,163],[-24,54],[-65,36],[-32,36],[-46,22],[-24,37],[-29,16],[-4,82],[-12,39],[3,83],[-16,22],[5,48],[25,72],[13,62],[-8,95],[-19,44]],[[57169,23229],[9,54],[23,11],[30,57],[62,37],[26,45],[37,120],[49,1],[25,57],[60,33],[25,-5],[61,-40],[81,4],[134,-15]],[[57921,23531],[21,-22],[-1,-38],[-27,-92],[24,-64],[88,-4],[31,-39],[62,-38],[49,-44],[99,-40],[74,28],[39,1],[50,-44],[30,-72],[16,-13]],[[60693,33295],[-52,59],[-76,36],[-178,49]],[[70673,15428],[-24,11],[-10,68],[-24,32],[-10,65],[-44,77],[-69,65],[-25,64],[10,39],[-2,85],[14,54]],[[23387,11274],[-13,20],[-40,156],[-25,54],[0,58],[-19,74],[24,101],[0,39],[-20,80],[-34,55],[-22,73],[0,53],[26,108],[38,87],[-6,80],[-27,99],[0,75],[13,58],[19,24]],[[49390,9257],[-102,90],[-15,63],[-40,54],[-72,34],[-67,-24],[-32,-31],[-25,5]],[[49037,9448],[-53,-24],[-30,-56],[-114,-17],[-117,-40],[-31,5],[-42,-20],[-26,18],[-21,65],[-57,106]],[[57582,26447],[99,-13],[26,-37],[52,-44],[177,-120],[32,-57],[77,-20],[25,-24],[4,-38],[27,-27],[51,-15],[38,4],[49,31],[85,-36],[64,-3],[99,10],[103,-6],[171,-26],[44,8]],[[20607,13224],[-34,144],[7,70],[-23,47],[-56,59],[-35,55],[-15,51],[-24,38],[-46,111],[-7,82],[-23,92],[-44,137],[-24,90],[-12,88],[-36,77],[-36,44],[-18,49],[-27,20],[-35,57],[3,50],[23,66],[35,68],[28,180],[30,32],[136,90],[34,-5],[15,-45],[34,-35],[100,-43],[62,-6],[13,-20],[79,-14],[8,-18],[118,-20],[43,64],[43,28],[60,10],[36,19],[12,30],[40,40],[72,-58],[32,50],[29,9],[76,-34],[55,44]],[[62825,20275],[18,-45],[4,-79],[14,-109],[48,-113],[121,-175]],[[62825,20275],[54,-47],[23,-56],[-19,-98],[12,-50],[56,-97],[48,-45],[54,-36]],[[63053,19846],[48,-20]],[[62038,22318],[28,23],[21,47],[40,16],[36,42],[58,-54],[85,-51],[78,45],[9,36]],[[62393,22422],[0,1]],[[62393,22422],[40,46],[-3,31],[55,15],[71,-33],[36,1],[112,-204],[72,-78],[18,-53],[-13,-73],[-2,-140],[16,-60],[76,-122],[42,-18],[16,-69],[52,-42],[-1,-70],[-24,-41],[35,-29],[48,-79],[13,-68],[-12,-72],[-4,-176],[29,-40],[5,-111],[-20,-67],[25,-69],[-11,-168],[9,-83],[-7,-46],[-82,-2],[-21,-26],[-13,-58],[-23,-29],[-51,-6],[-39,-32],[-12,-76]],[[22114,5903],[39,-38],[32,24],[10,81],[21,53],[25,30],[31,68],[9,96]],[[21246,18006],[-9,-73],[18,-62],[21,-22],[2,-36],[37,-19],[79,-21],[46,-1],[46,107],[108,24],[22,21],[74,8],[23,29],[116,73],[63,48],[46,87],[87,98],[67,39],[47,52],[29,64],[32,34],[61,11],[35,24],[20,-16],[42,27],[35,5],[47,-14],[61,0],[77,14],[66,30],[30,-2],[71,-27]],[[62539,26160],[11,-7],[37,-96],[76,-74],[30,-49],[56,-23]],[[17641,25978],[27,57],[56,15]],[[17668,26035],[-7,7]],[[17661,26042],[25,52]],[[17661,26042],[-33,38],[-55,40],[-71,90],[-9,40],[18,34],[-74,9],[-77,34],[-62,-2],[-61,55],[-29,64],[-9,70],[-26,42],[-33,4],[-4,50],[-26,20],[14,46],[-4,55],[20,22],[6,57],[31,84],[50,64],[15,44],[33,44],[1,30],[-37,53],[-9,56],[12,73],[-14,51],[10,58],[-8,73],[-17,47],[7,43],[25,27]],[[17245,27557],[2,58],[29,18],[8,52],[76,104],[10,73],[51,56],[6,35],[36,36],[7,65],[33,63],[6,57],[55,75],[15,64],[-8,39],[17,43],[48,77],[45,13],[17,60],[-6,60]],[[17692,28605],[-46,14],[-31,58],[6,73],[-15,69],[-23,33],[-157,140],[-18,45],[58,171],[-22,22],[-102,29],[-18,64],[-42,67],[-68,71],[-50,73],[-17,131],[24,102],[69,70],[34,87],[-30,87],[33,47],[82,23],[56,45],[40,83],[2,71],[-42,37],[-56,76],[-23,54],[-89,58],[-21,103],[12,58],[-30,53],[-10,117],[-24,64],[-106,66],[-50,68],[-73,50],[-64,30],[-124,88],[-54,95],[-152,78],[-72,135],[-5,112],[78,114],[85,15],[67,72],[22,71],[-66,72],[-48,0],[-68,63]],[[16327,31938],[35,96],[55,-3]],[[16447,32036],[42,-16]],[[15505,30594],[97,-68],[41,-15],[63,23],[65,2],[142,-59],[40,-44],[43,-20],[15,-29],[21,-93],[49,-103],[6,-90],[29,-40],[20,-113],[5,-132],[38,-103],[33,-32],[37,-84],[57,-45],[31,6],[21,-36],[-42,-60],[69,-132],[63,-23],[40,26],[46,-15],[167,27],[41,39],[31,0],[39,-31],[7,-77],[59,-18],[36,-28],[48,-76],[40,-24],[33,4],[84,36],[89,-13],[47,-27],[38,-2],[83,50],[46,42],[43,-13]],[[15079,31045],[11,-20],[98,-42],[9,-34],[47,-9],[74,-40]],[[12450,31776],[-43,40],[-25,89],[87,93],[63,24],[26,42],[69,64],[75,35],[104,66],[25,-53],[41,-40],[64,-2],[157,13],[34,-39],[139,-24]],[[14845,32053],[9,-55],[69,-27],[34,-56],[6,-60],[27,-62]],[[13681,32221],[9,12],[95,-8],[144,22],[50,-5],[87,22],[51,-5],[162,-47],[14,11]],[[71782,11069],[-78,34],[-90,91],[-43,34],[-46,17],[-22,24],[-53,132],[-39,55],[-69,36]],[[37198,20127],[19,-5],[45,46],[26,6],[38,41],[41,2],[61,39],[35,54],[47,99],[46,120],[27,12],[23,103],[61,-25],[9,-46],[72,-42],[19,-72],[-16,-25],[8,-49],[-7,-48],[19,-12],[2,-45],[-15,-119],[12,-25],[-24,-144],[-1,-44],[32,-102],[6,-122],[-4,-82],[14,-41],[7,-82],[-16,-58],[44,-123],[45,-72],[12,-44],[33,-60],[39,-10],[48,26],[-25,50],[7,57],[37,20]],[[72836,4352],[37,78],[-4,60],[-40,165],[-49,29],[-83,-5],[-31,-26]],[[72613,4634],[-35,3],[-99,35],[-59,-21]],[[72377,4665],[-50,12],[-53,37],[-36,10],[-81,-9],[-29,-20],[-10,-60],[-58,6],[-33,31],[-34,57],[-50,60],[-104,74],[-41,58],[-38,35],[-12,111],[-84,31],[-52,31],[-39,-58],[-34,60],[-12,85],[-50,52],[-54,19]],[[71423,5287],[-39,2],[-40,-26],[-115,41],[-34,22],[-43,-32],[3,-31],[-56,-80],[-67,85],[-70,5],[-45,35],[-22,-31],[-19,-192],[4,-48],[-45,-2],[-21,-93],[29,-42],[-30,-34]],[[3155,37560],[-124,8],[-44,-30],[-43,-87],[-203,-105],[-127,-126],[-167,-61]],[[37938,21069],[6,-109],[-12,-37],[10,-66],[-13,-74],[4,-32],[74,4],[37,-7],[50,-28],[46,-9],[17,-88],[14,-9],[-2,-95],[11,-59],[63,-131],[11,-65],[44,-34],[30,-85],[16,-182],[-63,-31],[-37,32],[-26,-38],[-41,-9],[9,-80],[-12,-31],[21,-125],[-10,-27],[-50,-3],[-35,-54],[38,-120],[17,-13]],[[62586,24961],[-20,-51],[11,-126],[-4,-40],[-36,-56],[-52,-36],[-38,25],[-37,-23],[-26,-110],[12,-29],[-6,-65],[27,-165],[97,32],[2,78],[29,16],[56,78],[38,10],[47,30],[43,58],[34,-27],[60,32],[62,77],[42,-18],[45,-48],[35,-18],[77,-15],[21,31],[3,47],[134,86],[47,12]],[[57390,23936],[-111,66],[-29,-7],[-28,31],[-80,-34],[-5,-103],[-18,-9],[-28,54],[-54,8],[-22,26],[-6,64],[-24,21],[-96,4],[-76,-31],[-98,3],[-28,-7],[-50,-39],[-82,-31],[-24,18],[-99,-72],[-133,-59],[-62,-69],[-37,20],[-35,-4],[-55,-24],[-107,-26],[-242,-41],[-24,5],[-40,-26],[-43,-6],[-64,-25],[-31,4],[-42,36],[-61,-21],[-34,-55]],[[55316,30138],[-71,-56],[-30,-40],[-39,-85],[-43,-42],[-63,-27],[-84,-11],[-86,-38],[-119,-105]],[[56498,30207],[65,-14],[-59,-69],[-41,-29],[-279,-15],[-69,-14],[-49,4],[-61,-20],[-129,15],[-88,-1],[-118,-14],[-54,39],[-14,58],[-54,29],[-82,-3]],[[39827,33484],[-68,36],[-97,-6]],[[22608,3694],[60,-25],[179,-26],[106,10],[104,-37],[74,-16],[31,-23],[46,-86],[15,-48],[31,-29],[78,-27],[40,-37],[83,-16],[46,-30],[62,-62],[77,-59],[44,-58],[49,-11],[73,-67],[15,-33]],[[22068,17149],[152,33],[65,34],[138,-15],[48,47],[50,-13],[26,17],[38,77],[80,75],[40,-21],[89,-125],[0,-37],[23,-74],[53,-227],[2,-138],[-60,-16],[-22,-39],[-21,-120],[10,-82],[18,-48],[29,-43],[70,-30],[55,-40],[38,-13],[79,18],[32,17],[104,-15],[44,-28],[59,-18],[48,2],[151,44],[82,-1],[38,-21],[76,-66],[66,-39],[34,-56],[34,-37],[70,-52],[128,-80]],[[24034,16019],[38,-4],[42,-42],[24,-93],[79,-133],[12,-46],[32,-64],[78,-121],[27,-34],[111,-43],[18,10]],[[15698,34209],[32,52],[1,80],[-47,28]],[[15658,34413],[-24,40],[14,40]],[[15809,34587],[70,22],[43,59],[50,157],[0,33]],[[16137,34933],[96,16]],[[16384,34967],[101,9],[82,42],[93,78],[77,40],[107,6]],[[21959,4496],[36,-180],[-12,-67],[17,-24],[100,-21],[37,-64],[25,-18],[28,-134],[29,-85],[54,-66],[54,-19],[60,55],[45,21],[98,-82],[41,-50],[37,-68]],[[45640,36133],[-57,-46],[-104,45],[-45,-14]],[[20039,30666],[-1,-56]],[[39879,18111],[-31,-48],[-12,-44],[-6,-95],[5,-42]],[[39879,18111],[1,-183]],[[39923,18233],[-35,-15],[-50,-63],[-82,22],[-60,-13]],[[39923,18233],[-25,-35],[-19,-53],[0,-34]],[[39947,18292],[11,-51],[-38,-106],[20,-31],[42,-119],[2,-87]],[[35927,19455],[35,87],[-18,69],[4,61],[-41,163],[35,27],[13,31],[58,33],[44,-18],[28,32],[58,31],[34,5],[28,33],[45,94],[34,10],[18,47],[7,59],[55,28],[47,70],[50,61],[16,42],[44,69],[29,97],[20,29],[30,3],[40,36],[90,117],[76,73],[93,14],[43,13],[51,54],[29,45],[55,12],[44,46],[37,25],[52,15],[51,48],[135,52],[49,78],[20,72],[39,63],[7,36],[-14,68],[38,50],[-25,23],[30,33],[20,70],[48,75]],[[39498,19709],[4,-36],[-13,-93],[16,-52],[40,-50],[53,-3],[55,-53],[72,-35],[46,-46],[51,4],[63,-14],[79,-145],[33,-86],[1,-70]],[[39998,19030],[-26,-269],[-1,-76],[19,-177],[-8,-74],[-35,-143],[-24,-58]],[[37422,21560],[-8,37],[35,80],[52,19],[64,62],[51,16]],[[37608,21736],[8,38]],[[37616,21774],[57,81]],[[37673,21855],[8,-51],[-48,-62],[-25,-6]],[[37673,21855],[33,65],[156,61],[72,11],[47,29],[105,40],[46,7],[46,-16],[45,0],[35,-19],[58,1],[26,-27],[33,-97],[27,-49],[21,-84],[33,-61],[59,-69],[35,-87],[44,-33],[-19,-119],[33,-39],[22,-58],[42,-83],[67,-71],[24,-44],[11,-54],[37,-77],[105,-89],[43,-58],[25,-157],[24,-63],[-1,-25],[47,-10],[100,-132],[49,-78],[22,-9],[56,-76],[77,-45],[50,-73],[64,-48],[40,-51],[4,-40],[-11,-75]],[[45977,21601],[10,133],[22,71],[36,53],[34,30],[54,19],[18,20],[65,36],[44,68],[7,99],[44,99],[6,49]],[[46317,22278],[4,31],[-16,159],[-43,58],[-10,68],[-19,26],[-4,99],[-8,43],[-30,44],[-7,38],[-84,28],[-42,-42],[-28,-47],[-59,-51],[-37,-63],[-36,-9],[-26,-39],[-8,-45],[-48,-51],[-36,-75],[-39,-48],[-77,-24],[-54,27],[-39,56],[-17,80],[-26,53],[-17,107],[-7,158],[-25,104],[5,23],[48,36],[-1,106],[-35,22],[-24,55],[18,92],[26,11],[42,78],[2,33],[63,72],[29,69]],[[46067,24340],[8,112],[-5,69],[7,76],[-12,88],[-56,65],[-20,37],[-11,71],[49,62],[13,41],[0,50],[-20,43],[-73,-32],[-31,3],[-17,31],[-32,10],[-55,112],[-17,7],[-86,162],[-86,98],[-23,47],[5,39],[-7,120],[-27,59],[36,206],[57,113],[17,47],[16,165],[-12,110]],[[61205,25795],[73,-96],[11,-96],[-30,-53],[-41,10],[11,-70],[49,-118],[15,-53],[-2,-89],[-45,-220],[-15,-54],[-54,-74],[-51,22],[-29,-2]],[[11187,33531],[45,18],[23,-25],[71,9],[22,19],[11,65],[30,33],[49,20],[111,0],[109,-28],[28,36],[-54,141],[2,38],[45,24],[41,63],[134,40],[73,-21],[60,19],[89,78],[66,81],[68,51],[125,-3],[106,-50],[37,-32],[67,-6],[135,31],[95,-59],[160,-39],[16,-40],[95,-39],[56,-74],[128,-87],[4,-27]],[[13233,33767],[53,-40],[191,-92],[53,11],[84,131],[127,133],[62,-3],[99,24],[91,2]],[[60882,26459],[60,-67],[62,-34],[40,-36],[20,-37],[14,-62],[26,-63],[63,-55],[30,-47],[50,-117],[42,-50],[18,-70],[12,-92],[49,-102],[6,-62],[25,-46],[12,-79],[6,-224],[-9,-158],[6,-264],[-6,-65],[17,-113],[23,-50],[-7,-66],[21,-52],[-38,-46],[-105,-49],[-12,-30],[41,-88],[5,-64],[24,-58],[-49,-24],[0,-76],[27,-23],[-18,-70],[18,-90],[-9,-31],[8,-52],[28,-37],[-32,-53],[-8,-37],[12,-46],[-7,-55],[20,-48],[10,-104],[-50,-22]],[[59702,26949],[6,-42],[34,-39],[24,-71],[39,-18],[25,-31],[36,-8],[48,59],[102,25],[37,-9],[74,3],[110,-14],[67,7],[98,-18],[1,-71],[51,-56],[66,35],[62,-5],[61,-61],[92,-2],[88,-47],[37,-101],[22,-26]],[[58199,33668],[-108,-4],[-137,-29],[-63,27],[-11,58],[1,161],[40,30],[-5,59],[-142,49],[-116,-26],[-79,40],[-35,45],[-94,35],[-96,6],[-20,25]],[[54781,36576],[145,-1],[51,-31],[116,21],[81,-4],[98,11],[256,-30],[101,-7],[119,12],[191,-29],[88,-37],[65,9],[73,-35],[74,-12],[137,-72],[90,-4],[81,-80],[80,-5],[12,-24],[148,-103],[0,-44],[43,-39],[27,-68],[89,-47],[70,5],[43,-46],[85,-15],[67,-44],[105,3],[74,-105],[47,-20],[30,-41],[62,-20],[53,-40],[177,-99],[58,-115],[85,-25],[51,-129],[66,-17],[29,-48],[7,-65]],[[58055,35136],[-8,-39],[-78,-109],[-3,-47],[-38,-91],[-16,-67],[-41,-67],[-74,-50],[-102,-24],[-12,-90],[51,-70],[-7,-23]],[[54433,36518],[167,52],[104,-10],[77,16]],[[54433,36518],[-63,-9],[-56,27],[-5,89],[-127,56],[-45,29],[-68,123],[-58,38],[-124,39],[-28,33],[-93,20],[-72,46],[-98,84],[-25,83],[75,70],[19,50],[49,28],[5,147],[49,57],[-18,54],[-53,55],[-1,40],[-44,32],[6,85],[-52,12],[-5,41],[64,77],[26,89],[49,82],[-13,67],[111,90],[36,52],[49,18],[63,-15],[88,20]],[[54074,38327],[206,-34],[106,70]],[[54074,38327],[58,39],[162,16],[92,-19]],[[54386,38363],[83,12]],[[42519,32785],[76,-34],[46,30],[38,52],[-13,109],[-51,73],[-35,92],[-121,67],[-39,46],[-75,53],[-100,16],[8,47],[-15,57],[-131,49],[-38,87],[-60,12],[-155,10],[-47,93],[15,55],[-25,39],[-90,76],[56,74],[58,138],[-8,37]],[[18615,29162],[-51,13],[-13,-66],[-43,-60],[-49,-43],[-24,-66],[-33,-7],[-42,45],[-34,3],[-28,-73],[-51,14],[-47,-7],[-39,-44],[-54,39],[-53,8],[-19,-20],[-52,5],[-32,-37]],[[17951,28866],[-20,-29],[5,-55],[-19,-22],[-71,-33],[5,-91]],[[17851,28636],[-93,32],[-41,-1],[-25,-62]],[[19823,29738],[-37,24],[-45,54],[-68,-27],[8,-79],[-15,-58],[-69,-201],[-70,-61],[-41,2],[-51,-40],[-33,-61],[5,-37],[-68,-3],[-31,-64],[6,-43],[-48,-67],[-55,6],[-31,30],[-29,55],[-28,1],[-83,-28],[-42,24],[-35,-17],[-60,41],[-65,30],[-31,60],[-81,26],[-30,-15],[-6,-49],[15,-37],[-90,-42]],[[44175,9501],[-12,27],[-57,48],[-18,31],[8,78],[22,87],[-66,63],[-54,-13],[-137,79],[36,53],[8,66]],[[43853,10062],[52,-42]],[[43893,10073],[12,-53]],[[43825,10133],[28,-71]],[[43825,10133],[67,-60]],[[43660,10263],[4,-43],[53,-97],[28,-34],[54,-23],[54,-4]],[[43660,10263],[13,-32],[103,-23],[49,-75]],[[43422,10552],[31,-12],[50,-79],[35,-27],[46,-93],[28,-39],[48,-39]],[[67880,40181],[-50,19]],[[67830,40200],[-106,-16],[-138,43],[-188,75],[-235,42],[-97,-6],[-101,43]],[[62827,38631],[89,-18],[75,20],[54,-30],[156,-38],[174,-32],[108,-1],[249,-47],[276,24],[66,64],[147,58],[-11,72],[99,44],[72,-32],[75,13],[87,-8],[-15,53],[-68,40],[52,29],[-5,93],[153,-3],[163,-18],[117,29],[151,-13],[20,-37],[169,12],[102,-55],[54,-156],[61,-3],[46,36],[83,12],[4,34],[105,31],[30,45],[-86,34],[-35,55],[51,30],[104,-41],[77,47],[-8,49],[92,-13],[88,42],[-198,36],[-5,30],[97,22],[5,86],[-44,81],[195,40],[65,61],[-1,47],[-91,76],[151,24],[68,48],[243,-12],[285,13],[91,23],[-50,61],[74,166],[94,67],[168,80],[-10,90],[-64,41],[-278,117],[-203,97],[-142,26],[-127,5],[-105,28]],[[47033,19475],[-1,-32],[36,-48],[-4,-81],[57,-97],[38,-44],[4,-115],[-15,-76],[-27,-73],[-8,-45],[34,-34],[-30,-95],[-34,-56],[-45,-32],[-41,27],[-29,1],[-40,-25],[-86,-72],[-29,1],[-31,-81],[9,-47],[-18,-209],[16,-15],[51,24],[-8,-69],[-24,-14],[2,-82],[-19,-68],[13,-55]],[[44353,6436],[-69,46],[-54,69]],[[44173,6644],[-57,52],[-45,58],[-44,19],[-42,47],[-69,120],[-30,-48],[-20,6],[-92,-74],[-38,-42],[-9,-36],[-41,58],[-67,48],[-53,140],[-50,30],[-57,127],[-41,13],[-58,-31],[-62,-72],[-85,-11],[-92,77],[-35,20],[-42,0],[-48,-18],[-44,-51],[-39,-24],[-7,-37],[-28,-30],[-75,23],[-56,5],[-67,-11],[-49,5],[-60,31],[-92,26],[-23,65],[9,48],[-17,41],[-25,10],[-35,55],[-28,-13],[-29,-64],[-13,-61],[-30,-12],[-42,-43]],[[45179,6938],[-10,-44],[-48,-47],[-20,-62],[11,-83],[-29,-71],[-44,-9],[-31,15],[-69,-16],[-24,-17],[-20,-49],[-42,-41],[-58,16],[-42,-57],[-34,-7],[-22,-29],[-44,-22],[-50,9],[-42,28],[-69,21],[-31,-15]],[[68086,10695],[82,53],[53,-20],[56,-59],[91,89],[24,57],[-13,69],[-32,38],[8,47]],[[68356,11184],[-8,33],[11,76],[-9,49],[-23,29],[-62,23],[-21,64]],[[23672,17253],[-47,-20],[-49,-51],[-33,-2],[-67,37],[-45,2],[-127,61],[-41,43],[-12,75],[-34,49],[-22,50],[-29,11],[-77,3],[-43,-45],[-66,22],[-33,48],[-21,60],[-36,61],[-34,33],[-16,62],[-23,17],[-85,11],[-59,34],[-15,97],[-13,41],[9,139],[-2,112],[8,34],[34,60],[3,50],[37,72],[0,62],[11,27]],[[22745,18508],[46,60],[25,57],[13,96],[-2,60],[26,38],[54,38],[38,43],[24,47],[38,29],[52,9],[37,19],[34,42],[30,17],[65,-4],[59,-19],[51,-51],[69,30],[18,21],[-5,51],[57,28],[41,-27],[28,-3],[33,39],[68,17],[72,42],[68,4],[38,16],[84,72],[65,29]],[[23971,19308],[63,-19],[61,5],[52,57],[69,55]],[[23991,19437],[6,-68],[-26,-61]],[[24068,19633],[-27,-113],[-17,-45],[-33,-38]],[[23937,19684],[-27,-74],[6,-46],[26,-58],[49,-69]],[[20541,31523],[18,-4],[37,-78],[59,-48],[68,3],[50,-27],[65,9],[177,63],[67,-3],[61,-33]],[[20432,31579],[50,-50],[59,-6]],[[20360,31620],[72,-40]],[[20359,31620],[1,0]],[[20015,31830],[35,-73],[54,-69],[75,-29],[181,-39]],[[41985,18593],[18,-39],[22,-6],[96,55],[39,-72],[47,10],[32,64],[39,46],[30,8],[52,-62],[80,-8],[76,22],[50,44],[42,122],[29,39],[35,73],[-30,41],[-38,-6],[-47,30],[-3,87],[-12,55],[-2,75],[40,110],[15,61],[11,87],[24,15],[26,47],[6,38]],[[55537,28753],[-29,9],[-63,-21],[-73,-46],[-33,-31],[-22,-47]],[[19893,19459],[32,-58]],[[19839,19513],[2,31]],[[55317,28617],[-70,-16],[-48,-26],[-82,-68],[-47,13],[-38,129],[16,156],[18,86],[-6,22],[-54,-10],[-9,35],[11,53],[-17,45],[-48,36],[-40,13],[-61,-40],[-46,-85],[-47,-50],[8,-85],[-48,-38],[-24,17],[-74,-2],[-31,-14],[-21,-50],[3,-63],[-25,-30],[-29,14],[-51,55],[-68,-22],[-66,-40],[-40,-37]],[[27948,8916],[-130,-99],[-96,-8],[18,58],[88,81],[59,30],[68,75],[36,28],[63,3],[54,34],[48,-24],[38,3],[131,99],[65,33]],[[28392,9229],[199,125],[68,37],[31,6],[91,-20],[58,-23],[45,34]],[[27653,10188],[-57,105],[10,19],[-4,75],[-38,69],[-39,23],[-47,2],[-39,54],[-19,10],[-30,-29]],[[27257,10429],[-48,31]],[[27004,10445],[-17,-6],[-36,-59],[-19,-6],[-32,31],[-55,-35]],[[24661,5389],[58,-26],[80,-50],[37,15]],[[24381,5596],[124,-118],[71,-49],[44,-8],[41,-32]],[[24381,5596],[64,-24],[67,-1],[37,-21],[78,-95],[34,-66]],[[24380,5596],[1,0]],[[24796,7518],[-64,-65],[13,-68],[-9,-83],[-40,-56],[-13,-48],[0,-74],[-27,-170],[-75,-71],[-29,-157],[9,-125],[-10,-74],[4,-110],[-99,-196],[-55,-94],[-34,-36],[-50,-34],[-10,-63]],[[24307,5994],[7,-57],[-19,-125],[14,-86],[36,-90],[36,-40]],[[25724,8135],[-5,-188],[-10,-100],[-35,-112],[-25,-17],[-47,-75],[-67,-25],[-33,-46],[-6,-46],[-28,-41],[-82,35],[-41,-19],[-30,-57],[-54,23],[-26,-8],[-37,17],[-55,-10],[-65,13],[-98,38],[-83,14],[-101,-13]],[[25916,8823],[-28,-24],[-21,-53],[-12,-83],[-18,-51],[-33,-29]],[[25804,8583],[47,20],[17,55],[31,136],[17,29]],[[26108,9111],[-36,-61],[-85,-69],[-32,-49],[-17,-65],[-22,-44]],[[26472,9771],[-32,-38],[-37,-126],[-41,-79],[-16,-78]],[[26683,10238],[-10,-22],[-43,-24],[-30,-50],[-5,-41],[-23,-25],[-8,-51]],[[25616,11779],[-83,-12],[-85,-63],[-87,-43],[-33,-33],[-42,-84],[12,-44],[40,-52],[3,-36],[40,-79],[-8,-55],[54,-72],[-17,-42],[-1,-55],[-71,-40],[-26,-26],[-55,-101],[7,-64],[-46,-105],[-68,-114],[-67,-46],[-15,-56],[11,-90],[-5,-58],[-33,-55],[-28,-103],[76,-41],[8,-81],[-23,-44],[-1,-33],[-45,-44],[-53,-82],[-74,-68],[28,-70],[11,-78],[25,-28],[-6,-57],[19,-40],[-11,-101],[-14,-48],[3,-107],[-16,-98],[3,-41],[37,-79],[5,-44],[-26,-53],[26,-75],[53,-76],[20,-45],[4,-49],[32,-64],[5,-91],[30,-53],[21,-89],[-10,-105],[-22,-43],[2,-35],[-67,-49],[-32,-44],[-2,-27]],[[25019,8169],[17,-67],[-43,-54],[-31,-102],[-45,-43],[-28,-79],[3,-82],[-100,-164],[4,-60]],[[10003,34844],[41,-15],[52,38],[107,49],[75,0],[152,-39],[152,20],[137,44],[26,-6],[38,-92],[98,-27],[88,5],[77,24],[108,87],[34,82],[-44,100],[43,153],[18,95],[-15,140],[38,24],[84,5],[32,48],[73,51],[87,7],[139,-26],[111,28],[141,72],[61,15]],[[11982,35713],[103,36],[76,59],[58,63],[177,-27],[-29,-35],[120,-12],[30,26]],[[52207,36824],[-118,-20],[-79,25],[-47,-2],[-76,-51],[-7,-23],[-108,-40],[-94,3],[-39,29],[-11,55],[-36,38],[20,80],[-23,82],[-86,40],[-13,55],[-62,90],[113,51],[28,42],[-17,32],[57,73],[54,25],[122,24],[31,59],[-10,35],[50,70],[24,108],[-87,81],[-101,52],[-5,90],[-63,86],[14,94],[-81,94],[-48,22],[-121,-7],[-116,-53],[-40,-76],[-59,-51],[-58,-22],[-64,-56],[-260,-60],[-122,8],[-120,36],[20,96],[62,66],[-13,129],[13,39],[-36,117],[7,93],[-58,67],[24,37],[-25,42],[116,46]],[[50659,38674],[112,54]],[[50659,38674],[92,-5],[114,74],[129,56],[39,61]],[[50771,38728],[-19,48],[55,85]],[[50958,38848],[-97,-80],[-90,-40]],[[11059,32547],[63,4],[20,-22],[-15,-68],[26,-74],[9,-64],[27,-30],[84,-20],[45,23]],[[22973,10123],[34,67],[74,-66],[35,8],[12,30],[45,43],[27,-8],[50,-60],[65,-15],[31,-57],[37,-114],[21,-86],[4,-59],[17,-34],[30,-10],[22,-38],[34,-102],[36,-34],[12,-38],[86,-41],[100,-74],[34,-54],[38,-137],[44,-79],[58,-49],[37,-57],[67,-146],[30,-40],[68,-69],[69,-54],[41,-65],[45,-34],[134,-51],[35,1],[51,-27],[66,-57],[55,-61],[151,-110],[47,-18],[36,-38],[54,-6],[73,-52],[41,-63]],[[10420,30151],[-20,-69],[-106,-33],[-33,-35],[-4,-87],[-109,-1],[-53,-12],[-62,-66]],[[40075,31134],[51,-13],[68,60],[13,63],[35,49],[35,10],[89,-7],[66,-48],[38,0],[97,28],[110,-17],[61,14],[21,-23],[102,-41],[89,45],[74,-2],[90,-53],[118,38],[102,-17]],[[5914,38002],[46,111],[127,22],[77,106],[155,127],[132,77],[-67,97],[0,48],[-205,69],[-143,-35],[-126,9],[-123,-34],[-174,-20],[-119,-39],[-29,-42],[-102,19],[-90,-55],[-69,-18],[-138,-3],[-136,-63],[-124,-13],[-52,-28],[-152,-28]],[[20210,29473],[34,15],[69,-24],[34,-85],[60,-43],[0,-36],[41,-19],[63,-54]],[[71787,14527],[34,-116],[0,-31],[30,-49],[56,-26],[82,43],[31,-28],[61,-76],[17,-37],[36,-25],[48,-12],[26,-45],[-11,-64],[-32,-17],[2,-42],[21,-48],[-5,-51]],[[21582,12918],[48,-18],[71,25],[26,35],[98,57],[57,62],[23,57],[5,45],[20,43],[33,41],[21,50],[23,20],[33,80],[39,54],[44,27],[58,-6],[20,18],[34,-4],[46,31],[22,30],[72,44],[71,-19],[28,-69],[27,-17],[68,-2],[154,94],[33,32],[-8,84],[41,116],[8,43],[64,88],[48,-1],[40,16],[40,38],[49,-42],[114,-17],[52,31],[19,33],[51,55],[84,24],[47,67],[-12,58],[10,34],[58,66],[19,105],[-9,73],[-36,52],[38,88],[46,-13],[39,15],[27,-32],[36,-5],[7,-28],[57,9],[44,35],[16,-20],[33,80],[35,26],[-3,79],[22,61],[31,9],[33,33],[19,41],[42,-1],[41,41],[3,48],[39,40],[40,25],[56,90],[8,46],[-14,51]],[[15498,23683],[28,-6],[59,45],[18,-15],[40,50],[43,23]],[[63890,24278],[23,-36],[21,-61],[42,-23],[44,-9]],[[20560,32056],[-26,-10]],[[20443,32075],[-31,-6]],[[20097,32101],[-41,-33],[-118,8]],[[20392,32106],[4,24]],[[20386,32158],[-37,-31],[-48,18],[-45,-33],[-50,25],[-33,54],[-50,-13],[-6,-68]],[[16711,32417],[-77,-30],[-81,43],[-132,28],[-21,39]],[[17117,32498],[-30,-34],[-131,3],[-107,-30]],[[17379,32650],[-32,-80],[-112,-37],[-46,-47]],[[40457,31780],[155,52]],[[40610,31832],[25,24],[-7,87],[29,72],[-13,39]],[[40545,32106],[12,0]],[[40487,32101],[-42,11],[-50,-37],[-67,11],[-62,-15],[-84,14],[-7,55],[11,95],[42,114],[17,74],[91,108]],[[40336,32531],[44,81],[18,61],[-27,70],[13,51],[-45,84],[-66,-11],[-36,17],[-46,59],[-11,58],[-50,29],[-81,120],[-7,47],[-64,100],[-1,73],[-21,37],[-56,45],[-73,32]],[[39827,33484],[-128,2],[-56,-21],[-62,3]],[[39496,30808],[-41,-20],[-8,-66]],[[39849,31626],[-54,-27],[-21,-36],[-5,-54],[-27,-57],[-76,56],[-37,-16],[-55,4],[-30,-32],[7,-49],[-16,-40],[18,-173],[-19,-50],[-25,-135],[0,-40],[40,-70],[-49,-61],[-1,-62],[37,-77]],[[40341,31716],[-69,-48],[-70,-2],[-131,-58],[-52,82]],[[15363,25185],[15,-38],[79,-35],[39,-34],[52,-18],[68,-39],[65,-7],[51,-19],[53,-37],[49,29]],[[14929,26135],[5,-19],[83,-110],[26,-81],[50,-104],[9,-57],[68,-92],[26,-60],[32,-44],[58,-60],[-2,-68],[17,-56],[-4,-37]],[[13643,26908],[31,-18],[42,-72],[103,-99],[69,-86],[56,-56],[69,-53],[33,-65],[36,-103],[0,-47],[41,-103],[92,-95],[103,-59],[95,-68],[39,14],[46,57],[13,74],[24,66],[27,36],[33,5],[31,31],[81,-22],[128,-6],[32,-47]],[[13475,27266],[12,-55],[61,-52],[31,-60],[16,-68],[48,-123]],[[13501,27352],[-23,-18]],[[13444,28862],[83,-11],[33,21],[90,-42],[80,-11],[73,-39],[34,-63],[-4,-71],[13,-153],[-15,-103],[-57,-67],[-22,-60],[-39,-62],[-51,-126],[-48,-94],[-5,-101],[-40,-208],[3,-72],[-13,-59],[-45,-83],[-11,-46]],[[45685,26351],[-17,11]],[[45668,26362],[-44,30],[-28,36],[-15,102],[-6,92],[-11,39],[-34,35],[-40,94]],[[42829,16403],[-34,32],[-145,60]],[[21026,33333],[53,22],[162,-42],[18,-92],[40,2]],[[20936,33344],[-9,27],[-56,6],[-249,-59]],[[20116,33361],[64,-33],[69,-11],[104,-3],[94,8],[87,19],[26,-7]],[[46700,12931],[-79,36],[-14,-60],[39,-173],[19,-60],[27,-4],[51,49],[29,-1],[40,-27],[25,-44],[77,-14],[36,10],[22,24],[49,15],[55,-40],[36,-8],[39,12],[43,31],[24,32],[21,62],[60,5],[32,-11],[41,-33],[71,60],[45,21],[42,0],[35,14],[56,47],[60,22],[40,30],[110,117]],[[14524,23134],[-2,19],[-65,25],[-28,72],[-24,24],[-25,87],[-28,13],[-17,-22]],[[14308,23340],[-31,17],[-21,42],[-34,29]],[[14222,23428],[-40,5],[-26,25],[-38,67],[-32,86],[-27,19],[-51,-30],[-8,16]],[[13984,23630],[-83,-49]],[[19921,29071],[43,72],[30,88],[69,73],[59,106],[88,63]],[[10152,29151],[16,-13],[0,-68],[-14,-71],[-40,-37]],[[9960,29828],[-6,-25],[44,-56],[22,-52],[-8,-33],[28,-104],[23,-49],[-6,-160],[32,-56],[13,-56],[50,-86]],[[22085,31619],[13,43],[45,59],[28,3],[38,64],[109,125],[49,9],[52,31],[120,49],[50,-26],[89,-95],[2,-131],[24,-30],[24,-151],[131,-39],[58,24],[37,-29],[66,-21],[41,-67],[-28,-97]],[[23107,8536],[11,-68],[-8,-37],[22,-86],[61,-83],[27,-63],[21,-7],[38,37],[36,-2],[18,-39],[40,-27],[60,31],[36,47],[26,4],[45,-118],[12,-70],[40,-95],[3,-68],[21,-159],[-8,-91],[5,-43],[34,-135],[-6,-61],[32,-96],[57,-83],[30,-68],[65,-68],[35,-119],[29,-55],[79,-66],[62,-2],[26,20],[168,-97],[75,-66],[8,-53],[0,-86],[11,-40],[1,-76],[11,-39],[-50,-45],[3,-78],[-18,-69],[23,-54],[5,-54],[-10,-89],[24,-26]],[[61327,23345],[-1,-41],[36,-136],[5,-46],[-30,-40],[-95,-53],[-46,-9],[-40,-36],[-56,-171],[-17,-157],[-24,-3],[-3,-68],[26,-90],[34,-19],[28,-87],[-5,-67],[13,-86],[-30,-152],[-4,-110],[14,-72]],[[10590,28608],[-19,-35],[-99,-30],[-43,20],[-26,52],[-82,85],[-81,121]],[[10814,28639],[-17,35],[-88,83],[-50,-95],[-46,-34]],[[10240,28821],[-33,105],[-42,27]],[[22281,6217],[79,-5],[32,-73],[27,-3],[87,19],[19,-34],[5,-102],[18,-62],[42,0],[81,-57],[68,-20],[20,-26]],[[18712,20120],[98,-45],[41,-66],[84,-12],[35,31],[4,25]],[[18401,20478],[55,-37]],[[41474,18092],[-52,20],[-63,-19],[-57,-67],[-24,-2],[-35,-49],[-52,-11],[-39,-39],[-33,-8],[-30,23],[-40,-19],[-50,-60],[-52,-40],[-52,-18],[-56,-2],[-37,-61],[-27,-15],[-42,-56],[-40,-6]],[[42165,17017],[13,-55],[48,-39],[45,-60],[11,-41],[15,-113],[15,-63],[71,-42],[34,-60],[-1,-107],[-7,-76],[19,-69],[0,-61],[-29,-58],[-44,-60],[-15,-56]],[[21685,5],[144,5],[79,-9],[40,12],[142,-13],[20,17],[75,1],[40,28],[121,0],[37,32]],[[15326,33906],[24,8]],[[13993,33933],[38,-26],[154,38],[49,22],[31,45]],[[14375,34068],[87,51],[-1,37],[35,45],[109,-25]],[[14617,34171],[79,-27],[72,6],[63,30],[45,-43],[52,-1],[21,-56],[136,-12],[10,-20]],[[36793,18529],[13,-26],[58,-70],[15,-58],[36,-84],[48,-72],[25,-95],[25,-25]],[[36740,19408],[3,-92],[28,-162],[9,-173],[17,-30],[-20,-88],[8,-79]],[[36756,19660],[-15,-67],[6,-104],[-7,-79]],[[19307,27513],[48,-46],[16,-44],[55,-86],[28,-19],[21,-43],[26,-108],[47,-75],[23,-88],[47,-26]],[[19165,27744],[54,-115]],[[39517,32076],[-29,53],[-4,50],[-37,62],[-45,36],[-83,107],[-69,3],[-121,-44],[-51,5],[-91,80],[-6,52],[-68,48],[-91,11],[-52,31],[-71,73],[-110,57],[-65,-3]],[[61837,32629],[84,27],[71,35],[76,-21],[146,-12],[60,10],[92,-2],[48,-13],[61,8],[132,40],[40,-3],[135,71],[101,37],[93,50],[129,97],[-16,40],[43,61],[7,121],[88,53],[90,109],[76,61],[38,80],[-22,89],[-26,6],[-139,-29],[-32,29]],[[63212,33573],[0,40],[32,17]],[[63212,33573],[-27,30],[0,49]],[[45360,16081],[-13,19],[-50,-20],[-31,1],[-15,34],[19,161]],[[45304,16430],[10,81],[-24,90],[39,93],[30,34],[50,30],[32,36],[31,69],[37,18]],[[47215,21024],[-43,32],[-64,131],[-23,27],[-29,-16],[-16,-38],[-89,6],[-58,-19],[-49,10],[-31,-13],[-45,28],[-6,33],[6,63]],[[49460,26629],[54,-70],[33,-79],[52,-18],[77,-88],[8,-50],[26,-24]],[[57405,26787],[-32,-16],[-49,44],[-23,44],[5,53],[32,42],[-38,120],[-37,31],[-59,29],[-56,15],[-106,-9],[-80,-27]],[[46625,11635],[3,-53],[-9,-45],[-34,-81],[-34,-22],[-37,-57],[-13,-46],[12,-100],[22,-46],[20,-86],[59,-81],[12,-55],[2,-294]],[[46620,11718],[-5,34]],[[61397,33076],[59,-24],[34,64],[27,19],[36,79],[-13,99]],[[61514,33333],[-71,20],[-173,-10],[-61,24]],[[8527,35192],[67,-57],[101,-62],[53,-47],[55,-83],[-80,-10],[-47,-20],[10,-85],[-11,-47],[42,-23],[30,-116],[-36,-79],[-103,-29],[-35,-35],[18,-52],[-40,-79],[-236,-107],[-33,-4]],[[12517,35823],[4,112],[-13,110],[-42,79],[-39,27],[-63,8],[-60,67],[-54,37],[26,21],[-23,43],[-44,16],[-56,69],[-63,36],[27,47],[-19,61],[-114,2]],[[12038,30468],[-35,-32],[-90,-18],[-64,11],[-123,49],[-22,55],[-62,37],[-92,-10],[-118,21],[-123,113],[-41,65],[-65,56],[25,76],[0,92],[-43,15],[-14,48],[59,109],[92,225],[17,61],[-32,52],[-34,24],[-39,60],[-3,58],[-30,71],[-100,85],[-124,-40],[-82,1],[-37,-21],[-65,-77],[-49,-28]],[[12179,30513],[-63,-54],[-42,2]],[[12569,30686],[-55,73],[-109,55],[-35,-7],[-7,-86],[-55,-72],[-76,-66]],[[12699,30849],[-45,-107],[16,-77],[-68,-22]],[[67411,31378],[64,24],[48,-44],[56,10],[36,23],[78,20],[51,50],[100,13],[90,84],[130,-22],[18,27],[116,-31],[57,14],[54,-13],[76,9],[40,-17],[49,47],[87,137],[49,40],[29,47],[52,36],[36,0]],[[69245,32139],[-34,-67],[-95,-94],[-70,-14],[-132,-10],[-76,-18],[-111,-103]],[[42854,32044],[8,-29],[-28,-38],[-9,-80]],[[13577,33284],[-38,41],[-11,48],[67,70],[5,75],[42,73],[65,80],[34,57],[22,85],[93,35],[65,56],[72,29]],[[11352,33417],[35,-52],[96,-77],[103,-38],[24,-22],[49,4],[75,31],[70,-1],[102,-64],[23,-49],[84,-17],[63,16],[71,-26],[45,7],[66,-20],[49,-83],[34,-27],[16,-51],[96,-38],[18,-72],[43,11],[55,39],[123,1],[15,114],[47,19],[11,41],[61,41],[17,56],[26,18],[165,21],[61,0],[81,-26],[87,-82],[46,-5]],[[10315,30378],[-129,27],[-79,31],[-34,-33]],[[19259,30581],[-17,-92],[-36,-59]],[[21051,31232],[-34,0],[-75,-35],[-143,-131]],[[8630,35261],[99,32],[25,33],[-30,71],[-58,5],[-106,77],[-178,3],[-59,-11],[-171,53],[-98,16],[-78,-33],[-111,-72],[-72,-88],[-24,-95],[25,-34],[-25,-146],[-87,-6],[-35,-23]],[[45993,26371],[1,-68]],[[45947,26451],[-13,55]],[[45934,26506],[6,32],[-2,186]],[[47676,35994],[97,-26],[31,-48],[132,39],[55,-10],[147,28],[60,50],[98,24],[90,92],[205,69],[176,93],[110,32],[106,11],[96,30],[118,17]],[[47503,36106],[-29,23],[-64,-3]],[[4315,36971],[-46,46],[-51,19],[-91,-5],[-92,40],[-144,13],[-59,-13],[-99,16],[-117,-19],[-82,-60],[-27,-61],[0,-96],[14,-56],[-49,-96],[-58,-62],[-1,-59]],[[57344,26560],[-19,12]],[[57303,26591],[-109,88],[-73,-18],[-72,56],[-45,63],[-70,24],[-70,43],[-91,21],[-40,46],[-53,-39],[-59,-58],[-77,4],[-42,-13],[-61,-52],[-52,-29],[-42,2],[-55,46],[-20,-1]],[[56195,26778],[-14,-22],[37,-69],[-2,-46],[-106,3],[-126,-6],[-52,13],[-75,47],[-59,-28],[-90,-65],[-77,-123],[-48,-51],[-98,-44],[-122,-88],[-70,-43],[-58,-12],[-43,-77],[-21,-15],[-66,-11],[-58,-35],[-67,14],[-46,-29]],[[49197,36395],[40,53],[23,77],[-4,46]],[[49256,36571],[-79,61],[-163,80],[-61,41],[-53,95],[-55,52],[-149,29],[-293,146],[-102,73],[-95,31],[-45,136],[-101,58],[42,49],[53,91],[0,28],[-211,57],[-98,45],[-28,65]],[[46053,36307],[41,10],[69,55],[68,24],[78,70],[39,-13],[177,33]],[[54658,29680],[-33,10],[-63,-22],[-45,8],[-28,66],[-4,45],[-88,99],[-21,40],[-36,28]],[[54232,30003],[17,26],[0,163],[51,154],[-18,82],[4,79],[-32,65],[-58,42],[-11,73],[-23,-4],[-41,42],[-3,25],[-76,118],[-44,36],[-50,-15],[-49,17],[-42,84],[-45,52],[-39,10],[-53,88],[-69,30],[-61,58],[-76,-1],[-102,41],[-63,92],[-73,38],[-1,40],[-41,28],[-85,-27],[-49,-1],[-52,28],[-39,-9],[-39,40],[-155,-35],[-46,42],[28,38],[-14,23]],[[52737,31595],[-70,-20],[-39,3],[-50,26]],[[27898,10505],[17,32],[14,73],[28,55],[0,66],[13,55],[30,71],[-11,63],[-51,149],[12,55],[-2,107],[29,10],[42,54],[42,27],[67,84],[28,48],[9,41],[33,85],[40,182],[34,101],[15,92],[25,74],[9,180],[48,102],[6,45],[-11,51],[5,59],[-26,96],[-5,79],[27,78],[24,101]],[[28931,13379],[34,-11],[61,25],[21,34],[11,49],[45,18],[40,36],[5,41],[34,21],[13,35],[56,41],[191,-111],[18,22],[26,-57],[23,-2],[26,-111],[66,-62],[176,-116],[38,-21],[32,-65],[64,-46],[15,-38],[30,-26]],[[35890,20992],[-28,26],[-51,110],[-16,49],[-41,49],[-61,35],[-46,16],[-29,28],[-56,28],[-56,73],[-45,30],[-29,37],[-16,69],[-28,14],[-48,67],[-21,89],[-49,56],[-72,2],[-25,53],[-76,88],[-54,25],[-59,-3],[-47,10],[-31,-24],[-92,-16],[-24,-16],[-78,-5],[-49,19],[-31,-32],[-42,-116],[-13,-85]],[[37028,29543],[-77,-29],[-132,-3],[-52,-24]],[[37787,29764],[-53,-83],[-62,-61],[-79,-16],[-115,-52],[-60,-13],[-102,46],[-73,-59],[-34,4]],[[37979,29788],[-55,65],[-23,4],[-78,-28]],[[5399,36814],[-14,47],[62,45],[98,97],[-53,99],[-88,65],[-126,40],[-33,37],[-82,30],[-174,-3],[-26,74],[-181,47],[-55,80],[-158,41],[-101,11],[-111,28],[-39,-7],[-60,42],[-35,65],[-56,42],[-70,23],[-86,1],[-133,-54],[-133,-8],[-11,101],[-73,2],[-79,-31],[-57,6],[-137,47],[-173,-24],[-87,73],[-60,13]],[[24783,12875],[16,102],[18,48],[25,25],[21,55],[15,83],[5,72],[-7,64],[-22,44],[-60,56],[-7,39],[10,72],[49,69],[-9,76],[17,87],[5,156],[25,68],[18,102],[-55,103],[-17,51],[19,73],[54,51],[33,12],[58,47],[23,78],[37,58],[68,73],[37,53],[5,35],[46,88],[17,65],[28,51],[54,4],[26,53],[4,32],[26,33],[86,68],[51,77],[40,87],[38,114],[20,98],[16,118],[19,65],[46,12],[37,-11]],[[57230,29744],[83,22],[70,57],[34,4],[100,55],[85,33],[61,5],[78,-25],[41,15],[94,12],[54,-4],[81,61],[54,4],[118,-22],[98,-47],[49,3],[98,25],[34,-10],[70,6],[-27,60]],[[58505,29998],[-48,36],[-67,25],[-21,47],[10,45],[89,35],[84,-10],[42,14]],[[71043,15258],[-82,26],[-35,26],[-50,7],[-20,19],[-65,9],[-58,62],[-60,21]],[[58005,37014],[84,82],[154,62],[87,74],[-101,77],[-159,56],[-136,-39],[-90,52],[-47,45],[-149,-40],[-141,67],[-67,83],[-24,79],[-58,27],[27,51],[-29,92],[23,44],[84,26],[49,38],[-19,26],[37,56],[52,31],[-28,95],[-36,32],[20,56],[-32,27],[46,37],[1,82],[-87,80],[-87,36],[-186,36],[-56,-15],[-198,45],[-117,49],[-126,72]],[[36767,29487],[-78,-1],[-46,-15],[-77,-47],[-91,-7],[-67,-72],[-30,-56]],[[18030,27877],[-26,-3],[-88,81],[-28,57],[29,59],[36,13],[11,35],[-6,84],[8,58]],[[18898,28237],[-87,-44],[-46,41],[-27,-20],[-25,-56],[-58,-78],[-39,-85],[-24,-23],[-81,-4],[-48,-81],[-64,-90],[-42,-41],[-86,61],[-54,-5]],[[17895,28615],[-44,21]],[[7914,35850],[-71,-1],[-65,18],[-32,35],[-39,119]],[[7440,36296],[-242,210],[-83,87],[-22,54],[-40,17]],[[37952,33423],[46,-11],[91,17],[35,18],[50,-90],[82,21],[49,-44],[106,28],[110,-11]],[[13789,36922],[65,18],[207,16],[9,42],[-20,72],[-3,86]],[[14076,37222],[82,164],[105,67],[49,82],[169,33],[67,-21],[21,-32],[113,23],[67,58]],[[16514,37486],[132,-34]],[[14813,37623],[173,36]],[[15844,37711],[132,-49],[12,-46],[49,-56]],[[48276,28628],[42,-76],[78,-46],[-29,-50],[93,-67],[21,-42],[20,-85],[-19,-156],[71,-153],[3,-57],[47,-124],[16,-85],[25,-31],[63,-16],[40,-25],[3,-77],[-25,-68],[12,-61],[43,-72],[76,-66],[82,-48],[17,-38],[52,-35],[67,-18]],[[49074,27132],[76,43],[48,-30],[81,-32],[10,-40],[-5,-56],[26,-24],[27,-52],[46,-31],[19,-64],[33,-56],[17,-56],[1,-95]],[[43114,31644],[-26,-40],[17,-118],[-16,-48],[35,-66],[20,-78],[-6,-22]],[[43386,32223],[-39,-7],[-28,-56],[-44,-48],[-76,-58],[-34,-42],[9,-25],[-52,-42],[1,-46],[-32,-75],[26,-51],[-3,-129]],[[44062,32320],[18,-29],[-38,-91],[-39,-7],[-68,16],[-93,50],[-177,-19],[-52,26],[-24,84],[-109,-75],[-52,-3],[-44,-49]],[[52585,33417],[24,-57],[81,-17],[53,23],[-3,165],[-16,47],[53,93],[131,28],[67,27],[57,59],[117,69],[58,48],[65,91],[3,47],[-23,70],[-3,67],[26,55],[43,49],[123,33],[51,25],[40,69],[-4,93],[53,78],[16,65],[43,30],[31,50],[41,123],[40,24],[16,49],[44,19],[35,57],[-2,49],[58,108],[83,64],[52,112],[-2,78],[38,31],[172,26],[38,85]],[[26987,11252],[-18,52],[-22,14],[-46,81],[-10,42],[18,61],[37,50],[24,55],[57,71],[-2,97],[10,29],[41,39],[46,-10],[50,52],[8,46],[31,12],[27,75],[-20,40],[5,120],[-4,68],[-24,90],[7,82],[-38,46],[-24,66],[-29,129],[0,33],[31,49],[6,42],[-12,48],[3,43],[26,81],[-9,99],[26,183],[-11,61],[2,93],[46,160],[13,63],[5,84],[-13,107],[35,23],[16,64],[51,55],[53,115],[-52,44],[42,64],[15,76],[16,123],[2,114],[-11,106],[-18,80],[-57,32],[-44,50],[-50,5],[-49,23],[-33,-16],[-47,-44]],[[26869,15217],[9,155],[-3,89],[-13,116],[19,28]],[[59296,34066],[-88,-69],[-51,1],[-62,70],[-98,12],[-59,-9],[-152,48],[-75,-9],[-19,26],[74,80],[10,96],[-9,57],[-95,48],[6,93],[-93,46],[-65,-11],[-40,74],[-69,33],[-63,67],[-116,13],[-53,29],[5,49],[-19,105],[1,92],[-111,129]],[[60045,27656],[23,31],[143,15],[70,86],[111,71],[41,8],[47,-16],[59,-53],[41,-15],[52,-43],[91,-139],[43,-48],[61,-15],[80,-73],[13,-47],[48,-13],[56,-109],[31,-31],[18,-52],[49,-44],[136,-93],[55,-128],[2,-51],[20,-31]],[[67830,40200],[-7,51],[146,153],[195,51],[81,45]],[[44141,15433],[-3,90],[-33,43],[-16,96],[-23,55],[-43,71],[-27,87],[-12,102],[-17,37],[-41,43],[-42,20],[-75,22],[-48,27],[-44,52],[-50,17],[-66,66],[-114,54],[-48,38],[-86,23],[-77,34],[-28,-28]],[[19040,27856],[44,-36]],[[67830,40200],[-38,43],[9,135],[81,60],[-34,84]],[[59534,27722],[83,7],[115,-13],[91,-38],[63,24],[43,-6],[69,-43],[47,3]],[[42746,17681],[6,-123],[-30,-125],[-65,-167],[-33,-131],[-3,-157],[-36,-137],[-4,-85],[9,-105],[-8,-100],[-27,-94],[-11,-73],[6,-90]],[[21758,11445],[15,57],[2,75],[32,80],[-14,140],[-51,73],[-14,49],[-27,27],[-39,5],[-164,115],[-38,11],[-63,-4],[-60,81],[-37,136],[-65,121],[-26,73],[-6,53],[15,32],[-9,63],[-53,151],[13,34],[51,-11],[61,-28],[-19,193],[-47,57],[-6,35],[6,108],[-16,37],[16,42],[-34,21],[-14,88],[-39,44],[-3,57],[-34,113],[7,38],[24,37],[32,22],[-54,35],[-2,41],[-31,62],[2,33],[-43,48],[-31,2],[-16,168],[-51,38],[-6,40],[26,25],[14,74],[-3,78],[-17,44],[0,98],[37,30],[69,25],[39,75],[57,7],[27,44],[10,75],[14,23],[9,89],[35,-8],[43,56],[3,38],[32,24],[18,83]],[[45069,17653],[-36,18],[-45,-4],[-22,30],[-71,-19],[-65,-9],[-29,17],[-48,51],[-32,-101],[-49,-9],[-10,23],[-11,91],[-31,49],[-98,-109],[-74,8],[-35,-53],[-34,11],[-71,-15],[-66,-42],[-19,20],[-70,5],[-43,18],[-43,39],[-45,19],[-48,0],[-77,44],[-57,-3],[-24,15],[-81,-10],[-42,20],[-57,73],[-120,39],[-51,9],[-35,-11],[-43,21],[-28,-3],[-64,37],[-92,16],[-39,73],[-39,47],[-52,38],[-45,47],[-71,13],[-42,-19],[-59,-59],[-55,-120],[-62,-87],[16,-100],[-9,-90]],[[52055,33426],[28,-88],[-9,-46],[-124,-70],[-88,8],[-128,58],[-80,26],[-92,15],[-67,24],[-119,73],[-137,24],[-110,-14],[-140,-37],[-234,-49],[-73,12],[-50,-36],[-66,10],[-50,-10],[-143,-79],[-50,-106],[-16,-169],[50,-62],[-2,-54],[45,-87],[101,-160],[-25,-91],[8,-36],[-23,-76],[-35,-23],[-13,-114],[17,-49],[-4,-99],[12,-45],[47,-59],[13,-64],[-49,-85]],[[52240,34352],[-54,-64],[-23,-76],[34,-158],[-39,-26],[-11,-52],[25,-35],[-53,-108],[-1,-41],[37,-128],[-53,-74]],[[26858,7125],[-77,22],[-39,-8],[-80,27],[-37,49],[-89,141],[-86,80],[-36,21],[-53,-3],[-26,43],[-36,31],[-64,15],[-52,-7],[-91,40],[-46,-12],[-18,25],[-35,-2],[-37,-35],[-24,17],[-50,-1],[-63,-42],[-29,-44],[-67,-18],[-50,-24],[-29,-64],[-34,-40],[-57,-29],[-30,-43],[-58,-38],[12,-46],[-43,-3],[-56,-98],[-52,-66],[-35,-78],[-108,-165],[-67,-61],[-42,-83],[-56,-64],[-15,-51],[-37,-68],[14,-88],[-20,-38],[7,-42],[-43,-128],[11,-62],[-38,-62],[-2,-81],[12,-66],[-18,-49]],[[17399,21849],[-51,39],[-6,40],[-36,51],[-69,62],[-26,45],[-45,50],[-20,63],[3,34],[26,41],[-27,37],[-34,-21],[-34,36],[-71,36],[-28,82]],[[16981,22444],[-58,37],[-44,5],[-14,92]],[[16981,22444],[-12,51],[-37,56],[-14,44]],[[44380,7402],[-60,-30],[-87,-102],[-50,-26],[-19,-61],[1,-36],[-27,-45],[-37,9],[-29,40],[-43,-51],[-28,-96],[-25,-35],[-60,-28]],[[44961,7664],[-28,50],[-103,-18],[-28,-33],[-74,-21],[-41,15],[-57,-33],[-58,-61],[-4,-73],[-47,-20],[-32,-47]],[[45429,7823],[13,-58],[-33,-49],[-48,-33],[-86,-19],[-47,-21],[-26,-30],[-70,-2],[-38,18],[-46,-5]],[[15091,23727],[81,-56],[49,-80],[38,-9],[37,15],[32,-16],[31,9],[30,33],[59,29],[18,33],[32,-2]],[[60387,33439],[-92,-43],[-98,-12],[-204,4],[-98,40],[-73,67],[8,35],[42,34],[45,62],[-31,86],[-126,16],[-76,19],[-15,40],[29,71],[-8,40],[75,31],[38,55],[-74,57],[-29,54]],[[60046,35526],[52,-66],[20,-85],[-10,-72],[45,-100],[-17,-40],[83,-130],[-2,-38],[-76,-83],[-60,-46],[-197,-54],[0,-36],[39,-32],[-14,-69]],[[46071,16905],[27,-165],[17,-72],[24,-53]],[[45974,16983],[15,-10]],[[45716,17193],[89,43],[46,4],[80,-44],[13,-64],[-25,-111]],[[67818,37593],[-1,-38],[-75,-75],[-149,-52],[-79,-2],[-120,42],[-229,23],[-91,-49],[-260,-17],[-100,-27],[-92,-4],[-171,-47],[-17,-42],[-185,14],[-182,-28],[-50,-55],[-152,-16],[-89,-36],[-38,-83],[-13,-108],[60,-113],[-3,-40],[-67,-38],[-78,-17],[-115,20],[-93,47],[1,34],[-74,37],[-16,65],[-83,21],[-42,75],[-7,54],[-102,1],[-24,55],[-73,-13],[-83,9],[-55,24],[-91,-68],[17,-55],[-117,60],[-51,-31],[-53,20]],[[62492,38115],[59,0],[116,-86],[126,-28],[159,39],[82,-59],[27,-55],[96,2],[103,43],[138,21],[93,-17],[15,-60],[43,-29],[66,6],[99,-67],[77,5],[89,-61],[-40,-72],[54,-31],[-42,-46],[9,-41],[-68,-23],[3,-41],[69,1],[-4,-47],[32,-68],[-63,-63],[73,-56]],[[42812,32753],[-21,57],[5,35]],[[42826,32848],[61,36],[90,-14],[61,2],[156,43],[50,41],[128,41],[89,65],[26,34],[10,52],[-5,128],[31,32],[-1,35],[-47,49],[-41,18],[-73,89],[4,25],[-66,101],[-56,38],[-119,-9],[-60,10],[-75,42],[-102,31],[-57,38],[-44,55],[-53,36],[-62,17],[-17,28],[61,69],[64,44],[31,85],[-7,87],[32,100]],[[49546,31800],[64,-18],[31,-65],[101,-130]],[[49040,32943],[-66,-76],[-5,-52],[-59,-104],[-30,-105],[-90,-153],[2,-60],[103,-36],[69,-5],[26,13],[81,-25],[65,-54],[16,-74],[61,-72],[48,-33],[60,-75],[31,-14],[19,-62],[88,-60],[78,-84],[47,-88],[-1,-44],[-62,-103],[-40,-87],[21,-44]],[[50026,34004],[32,-1],[25,-49],[-45,-39],[-55,5],[-68,-16],[-72,25],[-138,-41],[9,-61],[-90,-135],[-64,-125],[-60,-27],[-32,4],[-91,-45],[-79,-56]],[[48521,35059],[52,-58],[119,-74],[21,-30],[72,-31],[166,17],[37,-11],[85,32],[31,-14],[128,58],[85,-43],[72,-18],[154,-2],[33,-50],[56,-39],[139,-29],[45,-30],[14,-69]],[[47060,35084],[77,60],[18,57],[47,60]],[[46002,35248],[30,-100],[129,-15],[91,-44],[71,-84],[-28,-29],[131,-58],[57,-12],[55,42],[23,79],[51,12],[26,29],[79,44],[55,10],[101,-36],[7,-17]],[[48036,35327],[86,-15],[106,12],[94,-22]],[[47725,35351],[-51,67],[-63,46],[-189,67]],[[47328,35353],[42,99],[-10,47]],[[38446,18553],[-1,-49],[14,-25],[90,-30],[22,-41],[8,-45]],[[41019,35949],[-52,-72]],[[41026,36041],[1,-56]],[[41219,36321],[-3,-83],[-87,-14],[-59,-41],[-53,-63]],[[45509,36338],[13,18]],[[45116,36544],[8,-24],[108,-92],[72,-44]],[[50906,36690],[-35,52],[-48,-8],[-69,27],[-105,21],[-118,-5],[-69,-30],[-110,-2],[-100,-32],[-27,81],[-104,11],[-13,73],[-91,-16],[-59,-28],[-103,14],[-62,-11],[-117,-72],[-38,-55],[-101,-71],[-135,-48],[-146,-20]],[[51303,36924],[-137,48],[-56,-29],[-62,-7],[-21,-48],[43,-101],[116,-49],[-3,-54],[-175,-32],[-69,2],[-34,36]],[[39581,33468],[-94,-34],[-140,42]],[[18715,29709],[-9,68],[-84,61],[-78,39],[-43,-3],[-102,-33],[-62,-4],[-55,-30],[-46,-54],[-94,-52],[-52,-68],[15,-105],[-17,-112],[-37,-58],[-5,-60],[24,-65],[0,-54],[-25,-66],[-73,-120],[-33,-99],[12,-28]],[[79286,3662],[13,44]],[[79359,3786],[57,42],[-70,53],[-52,4],[-23,39],[-6,102],[-47,23],[-49,55],[-24,97],[-32,47],[-54,-18]],[[78069,2046],[-54,-105]],[[78015,1941],[-5,-13]],[[78010,1928],[6,-76],[64,-67],[103,-38]],[[73733,7050],[-44,-27],[-34,-70],[-31,-34],[-38,80],[-42,45],[-45,33],[-90,2]],[[10033,30713],[-64,31],[-20,75],[-117,92],[-37,41],[-24,119],[34,55],[-15,55],[27,49],[10,88],[61,13],[27,40],[-29,81]],[[16423,32868],[-13,46],[-77,2],[-48,65],[-44,10],[-77,-65],[-60,-7],[-10,44]],[[16103,33003],[-71,72]],[[65675,30441],[4,132],[29,22],[92,41],[99,27],[70,2],[174,-15],[358,63]],[[64705,23949],[14,-21],[39,-120],[36,-61]],[[64705,23949],[62,-31],[33,-33],[25,-4]],[[64692,24030],[108,5]],[[64652,24065],[20,-76],[33,-40]],[[64652,24065],[24,-12],[42,-67],[58,-6]],[[64302,24174],[30,-10],[26,-66],[38,-28],[100,-24],[31,17],[53,-21],[34,29],[38,-6]],[[66501,30713],[70,29],[67,2],[89,-14],[69,16],[52,44],[87,46],[61,4],[89,-53],[31,-6],[35,-48]],[[25995,12722],[33,48],[7,40],[-6,49],[31,57],[15,48],[7,74],[23,75],[-9,32],[16,49],[82,123],[67,3],[57,101],[-6,114],[-33,46],[-27,59],[-54,40],[-45,174],[6,87],[29,118],[24,77],[20,36],[85,52],[12,32],[-12,56],[-56,37],[-31,63],[-18,66],[-5,61],[-31,111],[-50,87],[-6,47],[12,73],[31,99],[10,62],[-8,46],[19,61],[5,84]],[[26190,15210],[61,86],[55,105],[37,-55],[10,-33],[28,-6],[22,38],[8,97],[-42,42],[-26,44],[-21,130],[-26,73],[-16,74],[-7,97],[13,35]],[[64020,24150],[26,45],[31,4],[67,-35],[28,5],[29,-22],[30,0],[71,27]],[[19454,28304],[53,10],[105,39],[52,11],[49,-9]],[[56638,26612],[-62,-58],[-20,-56],[-78,-37],[-18,-21],[-11,-62],[-53,-57],[-26,-51],[-27,-110],[17,-64],[3,-90],[21,-138],[19,-42],[37,-47],[7,-53],[-8,-44],[15,-72],[31,-36],[-3,-49],[22,-75],[59,-56],[97,-58],[12,-21],[51,-16],[74,-62],[59,-95],[45,-5],[28,-57],[28,-29],[78,-27],[38,-50],[73,-16],[24,-49],[62,-8],[57,-63],[8,-32],[52,-39],[94,20],[23,21]],[[64727,26136],[93,120],[39,39],[69,33],[-20,48],[48,43],[40,80]],[[64996,26499],[26,31],[61,-27],[29,5],[23,-58],[26,-11],[17,-39],[42,-46],[29,-53],[30,-29],[112,-43],[53,12]],[[65440,26230],[4,11]],[[65444,26241],[56,28],[35,45],[40,30],[28,57],[1,42],[72,50],[8,40],[84,34],[27,35],[28,79],[29,1],[28,41],[35,12],[36,119],[69,138],[63,42],[158,11],[30,-10],[71,-80]],[[56448,28122],[26,20],[17,44],[-29,46],[27,112],[-80,79],[-79,1],[-54,-18],[-61,15],[-42,51],[-66,55],[-26,84],[11,36],[-5,60],[40,2],[20,70],[-10,69],[7,46],[61,39],[55,-32],[49,33],[71,75],[32,57],[5,84],[28,47],[26,101],[73,69],[97,70],[59,22],[64,8],[89,30],[193,130],[111,102],[73,15]],[[58805,26034],[38,64],[44,11],[39,-11],[284,6],[107,-10],[51,8],[72,-28],[82,22],[68,-8],[242,-12],[90,-16],[53,-43],[47,-10],[43,13],[40,-32],[41,42],[33,13],[84,-3],[68,23],[31,22]],[[57747,39464],[158,-105],[91,-39],[71,-7],[142,46],[74,-22],[125,-76],[23,-122],[-14,-91],[68,-57],[14,-69],[41,-50],[-20,-62],[11,-61],[-46,-88],[16,-56],[50,-23],[38,-188],[-104,-85],[106,5],[68,-15],[112,-78],[-47,-51],[105,-19],[39,-33],[-39,-39],[62,-63],[-49,-267],[-60,-93],[-3,-150],[12,-53],[-59,-98],[-19,-97],[176,-95],[96,-78],[37,-72],[98,-81],[15,-68],[57,-94],[88,-47],[33,-54],[65,-30],[-74,-29],[-40,-53],[27,-101],[60,-91],[38,-115],[75,-63],[47,-62],[43,-108],[128,-142],[126,-94],[27,-105],[79,-76],[25,-41],[107,-33]],[[63419,25111],[12,-57],[35,18],[20,87],[85,108],[44,7],[51,-65],[35,20],[35,-40],[48,11],[28,23],[56,70],[1,62],[53,31],[64,-9],[25,13],[39,77],[31,89],[5,43],[-28,51],[-32,31],[24,69],[31,35],[29,68],[37,30],[120,48],[29,-25],[30,9],[27,64],[37,9],[73,-9]],[[64675,25994],[20,38]],[[64713,26111],[14,25]],[[7052,36664],[10,19],[-21,93],[-106,37],[-118,7],[-89,25],[-11,50],[-53,50],[-180,116],[-87,13],[-155,-6],[-237,69],[-27,19],[7,93],[-53,52],[-3,116],[61,63],[-24,55],[-165,46],[-43,53],[-122,53],[-29,65],[-57,26],[-7,39],[-84,32],[-55,38],[-152,31],[-133,-12],[-84,23],[-101,61],[-6,99],[-133,98],[-89,78],[-104,44]],[[4602,38309],[-87,8]],[[462,36811],[4,-44],[-35,-18],[-104,-11],[-80,7],[-160,-50],[21,-37],[-95,-52],[-13,-41]],[[462,36811],[-150,38],[10,85],[-16,45],[-68,59],[-86,0]],[[4515,38317],[-195,-34],[-214,-73],[-181,-20],[-141,-48],[-73,-56],[-100,1],[-26,-51],[-116,-77],[-109,-24],[-141,-61],[-151,-31]],[[3068,37843],[-70,-14],[-118,26],[-110,-30],[-95,-9],[-89,-35],[-235,-36],[-103,-36],[-269,-30],[-70,43],[-171,18]],[[1738,37740],[-33,2],[-87,-70],[-91,-99],[1,-67],[-36,-47],[-75,-36],[-72,-102],[2,-45],[-44,-44],[-43,-87],[-45,-39],[-35,-109],[50,-100],[49,-22],[-9,-33],[-72,-46],[-227,-21],[-43,-18],[-34,-59],[-121,-32],[-50,94],[-91,47],[-38,-27],[-82,-5],[-50,36]],[[4515,38317],[-2,0]],[[60298,27347],[85,2],[14,14],[63,-8],[58,-32],[18,-39],[80,-31],[101,-92],[70,-9],[46,-42],[47,-92],[43,-22],[40,-40]],[[46802,10462],[24,-66],[-4,-100]],[[46802,10462],[107,-79]],[[45110,11053],[16,36],[11,127],[22,28]],[[45159,11244],[73,57],[47,26],[56,17],[248,1]],[[46027,11366],[51,-2],[104,-124],[45,-60],[97,-148],[189,-107],[23,-44],[47,-152],[45,-60]],[[46628,10669],[44,-71],[64,-52],[66,-84]],[[44065,12743],[14,22],[-5,56],[-18,13],[-48,-8],[-63,-67],[-64,-45],[-71,-72],[-60,-44],[-36,-43],[-48,-129],[7,-113],[40,-145],[63,-73],[47,-141],[-7,-113],[5,-63],[-43,-152],[-15,-96],[3,-65],[21,-36],[15,-58],[-8,-34],[20,-56],[24,-140],[47,-93],[82,-105],[53,-89],[38,-111],[162,-21],[37,-30],[38,-52],[45,-17],[44,7],[52,-42],[31,-6],[33,19],[45,-6],[57,-30],[47,-6],[56,26]],[[77951,1952],[59,-24]],[[77998,1990],[17,-49]],[[14345,35643],[-27,29]],[[14290,35720],[-50,39],[9,38],[-42,17],[-151,19],[-51,45]],[[13880,35871],[-32,33],[-77,6]],[[6567,36551],[-21,-103],[-60,-29],[-73,-9],[-37,-35],[-3,-56],[-37,-94],[31,-54],[-30,-61],[4,-47],[34,-74],[-83,-36],[-9,-75],[-68,-2]],[[9713,37838],[-62,-29],[-109,-11],[-280,-15],[-48,-28]],[[77666,1611],[-33,-80],[15,-71],[52,-79],[1,-59],[-37,-117]],[[21451,31617],[-41,-44]],[[14369,36597],[-5,-27],[-68,-76],[10,-33]],[[14436,36752],[-14,-51]],[[15259,36904],[243,-16],[85,0],[1,42],[-34,77]],[[14853,36920],[11,-2]],[[15000,36922],[42,-1]],[[14534,36940],[-43,-30]],[[14744,36976],[71,47]],[[14697,37108],[-80,-44],[-65,-10]],[[15721,37126],[40,59]],[[15801,37302],[78,40]],[[16129,37348],[47,47],[41,83]],[[63177,19851],[-18,11],[-58,-8]],[[63101,19854],[-48,-8]],[[63101,19854],[52,-31],[25,-53],[44,-33]],[[14900,36925],[-36,-7]],[[38761,30711],[-22,90],[-65,112],[-52,31],[-81,29],[-44,26],[-15,42],[-59,46],[-81,24],[-44,53],[-3,67]],[[39139,31881],[-39,99],[1,45],[-30,55],[-100,70],[-50,22],[-70,7],[-159,-137],[-40,-16],[-130,-28],[-52,-31],[-57,1],[-41,38],[-47,12],[-148,-9],[-134,-57],[-26,5]],[[39139,31881],[39,-11],[18,-34],[38,-5],[54,-63],[6,-27],[55,-31],[27,-131],[-16,-28],[42,-76],[10,-88],[-40,-78],[-40,-40],[4,-81],[63,-15]]],"transform":{"scale":[0.004300731962112824,0.0030495790212123936],"translate":[-165.24393917525268,-50.24013722075609]},"objects":{"ne_50m_rivers_lake_centerlines":{"type":"GeometryCollection","geometries":[{"type":"LineString","arcs":[0]},{"type":"LineString","arcs":[1]},{"type":"LineString","arcs":[2]},{"type":"LineString","arcs":[3]},{"type":"MultiLineString","arcs":[[4],[5]]},{"type":"MultiLineString","arcs":[[6],[7],[8],[9],[10]]},{"type":"LineString","arcs":[11]},{"type":"MultiLineString","arcs":[[12],[13]]},{"type":"MultiLineString","arcs":[[14],[15],[16]]},{"type":"LineString","arcs":[17]},{"type":"MultiLineString","arcs":[[18],[19]]},{"type":"MultiLineString","arcs":[[20],[21],[22],[23]]},{"type":"LineString","arcs":[24]},{"type":"MultiLineString","arcs":[[25],[26],[27]]},{"type":"LineString","arcs":[28]},{"type":"MultiLineString","arcs":[[29],[30],[31],[32],[33]]},{"type":"MultiLineString","arcs":[[34],[35],[36],[37],[38]]},{"type":"MultiLineString","arcs":[[39],[40]]},{"type":"MultiLineString","arcs":[[41],[42],[43],[44]]},{"type":"LineString","arcs":[45]},{"type":"MultiLineString","arcs":[[46],[47],[48],[49]]},{"type":"LineString","arcs":[50]},{"type":"LineString","arcs":[51]},{"type":"LineString","arcs":[52]},{"type":"LineString","arcs":[53]},{"type":"LineString","arcs":[54]},{"type":"MultiLineString","arcs":[[55],[56],[57],[58],[59]]},{"type":"LineString","arcs":[60]},{"type":"MultiLineString","arcs":[[61],[62]]},{"type":"LineString","arcs":[63]},{"type":"MultiLineString","arcs":[[64],[65]]},{"type":"LineString","arcs":[66]},{"type":"LineString","arcs":[67]},{"type":"LineString","arcs":[68]},{"type":"MultiLineString","arcs":[[69],[70],[71]]},{"type":"LineString","arcs":[72]},{"type":"LineString","arcs":[73]},{"type":"LineString","arcs":[74]},{"type":"MultiLineString","arcs":[[75],[76]]},{"type":"MultiLineString","arcs":[[77],[78],[79],[80],[81]]},{"type":"LineString","arcs":[82]},{"type":"MultiLineString","arcs":[[83],[84]]},{"type":"LineString","arcs":[85]},{"type":"LineString","arcs":[86]},{"type":"MultiLineString","arcs":[[87,88],[89],[90],[91],[92]]},{"type":"LineString","arcs":[93]},{"type":"MultiLineString","arcs":[[94],[95]]},{"type":"LineString","arcs":[96]},{"type":"LineString","arcs":[97]},{"type":"MultiLineString","arcs":[[98],[99]]},{"type":"LineString","arcs":[100]},{"type":"LineString","arcs":[101]},{"type":"LineString","arcs":[102]},{"type":"MultiLineString","arcs":[[103],[104],[105]]},{"type":"MultiLineString","arcs":[[106],[107],[108,109]]},{"type":"MultiLineString","arcs":[[110],[111]]},{"type":"LineString","arcs":[112]},{"type":"MultiLineString","arcs":[[113],[114],[115],[116]]},{"type":"MultiLineString","arcs":[[117],[118]]},{"type":"MultiLineString","arcs":[[119],[120]]},{"type":"LineString","arcs":[121]},{"type":"MultiLineString","arcs":[[122],[123],[124],[125]]},{"type":"MultiLineString","arcs":[[126],[127]]},{"type":"MultiLineString","arcs":[[128],[129],[130]]},{"type":"LineString","arcs":[131]},{"type":"LineString","arcs":[132]},{"type":"LineString","arcs":[133]},{"type":"MultiLineString","arcs":[[134],[135],[136],[137]]},{"type":"LineString","arcs":[138]},{"type":"LineString","arcs":[139]},{"type":"MultiLineString","arcs":[[140],[141]]},{"type":"LineString","arcs":[142]},{"type":"LineString","arcs":[143]},{"type":"MultiLineString","arcs":[[144],[145],[146]]},{"type":"LineString","arcs":[147]},{"type":"LineString","arcs":[148]},{"type":"LineString","arcs":[149]},{"type":"MultiLineString","arcs":[[150],[151]]},{"type":"LineString","arcs":[152]},{"type":"MultiLineString","arcs":[[153],[154],[155]]},{"type":"MultiLineString","arcs":[[156],[157]]},{"type":"MultiLineString","arcs":[[158],[159]]},{"type":"MultiLineString","arcs":[[160],[161]]},{"type":"LineString","arcs":[162]},{"type":"MultiLineString","arcs":[[163],[164],[165],[166]]},{"type":"LineString","arcs":[167]},{"type":"LineString","arcs":[168]},{"type":"LineString","arcs":[169]},{"type":"MultiLineString","arcs":[[170],[171]]},{"type":"LineString","arcs":[172]},{"type":"MultiLineString","arcs":[[173],[174],[175]]},{"type":"LineString","arcs":[176]},{"type":"LineString","arcs":[177]},{"type":"MultiLineString","arcs":[[178],[179,180],[181],[182],[183],[184],[185]]},{"type":"LineString","arcs":[186]},{"type":"MultiLineString","arcs":[[187],[188]]},{"type":"MultiLineString","arcs":[[189],[190]]},{"type":"LineString","arcs":[191]},{"type":"MultiLineString","arcs":[[192],[193],[194]]},{"type":"MultiLineString","arcs":[[195],[196]]},{"type":"MultiLineString","arcs":[[197],[198]]},{"type":"MultiLineString","arcs":[[199],[200],[201]]},{"type":"MultiLineString","arcs":[[202],[203],[204],[205],[206],[207],[208],[209],[210],[211],[212]]},{"type":"LineString","arcs":[213]},{"type":"LineString","arcs":[214]},{"type":"MultiLineString","arcs":[[215],[216]]},{"type":"MultiLineString","arcs":[[217],[218],[219]]},{"type":"MultiLineString","arcs":[[220],[221],[222],[223],[224],[225]]},{"type":"LineString","arcs":[226]},{"type":"LineString","arcs":[227]},{"type":"LineString","arcs":[228]},{"type":"LineString","arcs":[229,230,231,232,233,234]},{"type":"MultiLineString","arcs":[[235],[236],[237],[238]]},{"type":"MultiLineString","arcs":[[239,240],[241],[242]]},{"type":"LineString","arcs":[243]},{"type":"MultiLineString","arcs":[[244],[245],[246],[247],[248],[249],[250],[251]]},{"type":"LineString","arcs":[252,253,254]},{"type":"LineString","arcs":[255]},{"type":"LineString","arcs":[256]},{"type":"MultiLineString","arcs":[[257],[258]]},{"type":"MultiLineString","arcs":[[259],[260],[261]]},{"type":"LineString","arcs":[262,263,264,265,266]},{"type":"MultiLineString","arcs":[[267],[268],[269],[270],[271]]},{"type":"LineString","arcs":[272]},{"type":"LineString","arcs":[273]},{"type":"LineString","arcs":[274]},{"type":"LineString","arcs":[275]},{"type":"LineString","arcs":[276]},{"type":"LineString","arcs":[277]},{"type":"LineString","arcs":[278]},{"type":"LineString","arcs":[279]},{"type":"LineString","arcs":[280]},{"type":"LineString","arcs":[281]},{"type":"LineString","arcs":[282]},{"type":"MultiLineString","arcs":[[283,284],[285],[286]]},{"type":"LineString","arcs":[287]},{"type":"LineString","arcs":[288,289]},{"type":"LineString","arcs":[290]},{"type":"LineString","arcs":[291]},{"type":"LineString","arcs":[292]},{"type":"LineString","arcs":[293]},{"type":"LineString","arcs":[294]},{"type":"LineString","arcs":[295]},{"type":"LineString","arcs":[296]},{"type":"LineString","arcs":[297]},{"type":"LineString","arcs":[298,299]},{"type":"MultiLineString","arcs":[[300],[301],[302]]},{"type":"MultiLineString","arcs":[[303],[304],[305],[306]]},{"type":"LineString","arcs":[307]},{"type":"MultiLineString","arcs":[[308,309],[310]]},{"type":"LineString","arcs":[311]},{"type":"MultiLineString","arcs":[[312],[313]]},{"type":"LineString","arcs":[314]},{"type":"MultiLineString","arcs":[[315],[316]]},{"type":"LineString","arcs":[317]},{"type":"MultiLineString","arcs":[[318],[319],[320],[321,322],[323],[324],[325],[326],[327],[328]]},{"type":"MultiLineString","arcs":[[329,330],[331,332],[333,334],[335],[336],[337],[338]]},{"type":"MultiLineString","arcs":[[339],[340,341,342,343,344]]},{"type":"LineString","arcs":[345]},{"type":"MultiLineString","arcs":[[346],[347],[348]]},{"type":"LineString","arcs":[349]},{"type":"LineString","arcs":[350]},{"type":"LineString","arcs":[351]},{"type":"LineString","arcs":[352]},{"type":"MultiLineString","arcs":[[353],[354],[355],[356],[357]]},{"type":"LineString","arcs":[358]},{"type":"LineString","arcs":[359,360,361,362,363,364,365,366]},{"type":"LineString","arcs":[367]},{"type":"LineString","arcs":[368]},{"type":"LineString","arcs":[369]},{"type":"LineString","arcs":[370]},{"type":"LineString","arcs":[371]},{"type":"LineString","arcs":[372]},{"type":"LineString","arcs":[373]},{"type":"MultiLineString","arcs":[[374],[375]]},{"type":"MultiLineString","arcs":[[376],[377],[378],[379],[380]]},{"type":"MultiLineString","arcs":[[381],[382],[383]]},{"type":"MultiLineString","arcs":[[384],[385]]},{"type":"LineString","arcs":[386]},{"type":"LineString","arcs":[387]},{"type":"LineString","arcs":[388]},{"type":"LineString","arcs":[389]},{"type":"LineString","arcs":[390]},{"type":"MultiLineString","arcs":[[391],[392]]},{"type":"MultiLineString","arcs":[[393],[394]]},{"type":"MultiLineString","arcs":[[395],[396]]},{"type":"LineString","arcs":[397]},{"type":"MultiLineString","arcs":[[398],[399]]},{"type":"LineString","arcs":[400]},{"type":"LineString","arcs":[401]},{"type":"MultiLineString","arcs":[[402],[403]]},{"type":"LineString","arcs":[404]},{"type":"LineString","arcs":[405]},{"type":"LineString","arcs":[406]},{"type":"MultiLineString","arcs":[[407],[408],[409],[410],[411],[412],[413],[414],[415],[416],[417],[418]]},{"type":"LineString","arcs":[419]},{"type":"LineString","arcs":[420]},{"type":"LineString","arcs":[421]},{"type":"MultiLineString","arcs":[[422],[423]]},{"type":"LineString","arcs":[424]},{"type":"MultiLineString","arcs":[[425],[426],[427]]},{"type":"MultiLineString","arcs":[[428],[429],[430],[431],[432],[433],[434]]},{"type":"LineString","arcs":[435]},{"type":"LineString","arcs":[436]},{"type":"LineString","arcs":[437]},{"type":"LineString","arcs":[438]},{"type":"MultiLineString","arcs":[[439],[440]]},{"type":"LineString","arcs":[441]},{"type":"LineString","arcs":[442]},{"type":"LineString","arcs":[443]},{"type":"LineString","arcs":[444]},{"type":"LineString","arcs":[445]},{"type":"MultiLineString","arcs":[[446],[447],[448]]},{"type":"LineString","arcs":[449]},{"type":"LineString","arcs":[450]},{"type":"LineString","arcs":[451]},{"type":"LineString","arcs":[452]},{"type":"MultiLineString","arcs":[[453],[454]]},{"type":"LineString","arcs":[455]},{"type":"MultiLineString","arcs":[[456],[457],[458],[459],[460],[461],[462],[463,464]]},{"type":"MultiLineString","arcs":[[465,466],[467],[468],[469]]},{"type":"MultiLineString","arcs":[[470],[471],[472]]},{"type":"LineString","arcs":[473]},{"type":"LineString","arcs":[474]},{"type":"LineString","arcs":[475]},{"type":"LineString","arcs":[476]},{"type":"MultiLineString","arcs":[[477],[478]]},{"type":"MultiLineString","arcs":[[479],[480],[481]]},{"type":"LineString","arcs":[482]},{"type":"LineString","arcs":[483]},{"type":"LineString","arcs":[484]},{"type":"LineString","arcs":[485]},{"type":"LineString","arcs":[486]},{"type":"LineString","arcs":[487]},{"type":"LineString","arcs":[488]},{"type":"LineString","arcs":[489,490]},{"type":"LineString","arcs":[491]},{"type":"LineString","arcs":[492]},{"type":"MultiLineString","arcs":[[493],[494],[495],[496]]},{"type":"MultiLineString","arcs":[[497],[498]]},{"type":"LineString","arcs":[499]},{"type":"MultiLineString","arcs":[[500],[501]]},{"type":"LineString","arcs":[502]},{"type":"MultiLineString","arcs":[[503],[504]]},{"type":"LineString","arcs":[505]},{"type":"LineString","arcs":[506]},{"type":"LineString","arcs":[507]},{"type":"LineString","arcs":[508]},{"type":"LineString","arcs":[509]},{"type":"MultiLineString","arcs":[[510],[511]]},{"type":"LineString","arcs":[512]},{"type":"LineString","arcs":[513]},{"type":"LineString","arcs":[514]},{"type":"MultiLineString","arcs":[[515],[516],[517],[518]]},{"type":"LineString","arcs":[519]},{"type":"LineString","arcs":[520]},{"type":"MultiLineString","arcs":[[521],[522,523]]},{"type":"LineString","arcs":[524]},{"type":"LineString","arcs":[525]},{"type":"LineString","arcs":[526]},{"type":"LineString","arcs":[527]},{"type":"MultiLineString","arcs":[[528],[529]]},{"type":"LineString","arcs":[530]},{"type":"MultiLineString","arcs":[[531],[532]]},{"type":"MultiLineString","arcs":[[533,534,535],[536],[537],[538],[539]]},{"type":"LineString","arcs":[540]},{"type":"MultiLineString","arcs":[[541],[542]]},{"type":"LineString","arcs":[543]},{"type":"MultiLineString","arcs":[[544],[545]]},{"type":"MultiLineString","arcs":[[546],[547]]},{"type":"LineString","arcs":[548]},{"type":"LineString","arcs":[549]},{"type":"LineString","arcs":[550]},{"type":"LineString","arcs":[551]},{"type":"LineString","arcs":[552]},{"type":"LineString","arcs":[553]},{"type":"LineString","arcs":[554]},{"type":"MultiLineString","arcs":[[555],[556,557],[558],[559],[560]]},{"type":"LineString","arcs":[561]},{"type":"LineString","arcs":[562]},{"type":"LineString","arcs":[563]},{"type":"MultiLineString","arcs":[[564],[565],[566],[567,568,569],[570],[571]]},{"type":"MultiLineString","arcs":[[572],[573],[574],[575],[576]]},{"type":"LineString","arcs":[577]},{"type":"LineString","arcs":[578]},{"type":"MultiLineString","arcs":[[579],[580],[581,582]]},{"type":"LineString","arcs":[583]},{"type":"LineString","arcs":[584]},{"type":"LineString","arcs":[585]},{"type":"LineString","arcs":[586]},{"type":"MultiLineString","arcs":[[587],[588]]},{"type":"LineString","arcs":[589]},{"type":"MultiLineString","arcs":[[590],[591,592]]},{"type":"MultiLineString","arcs":[[593],[594],[595],[596],[597]]},{"type":"LineString","arcs":[598]},{"type":"LineString","arcs":[599]},{"type":"LineString","arcs":[600]},{"type":"MultiLineString","arcs":[[601],[602],[603],[604],[605],[606],[607,608],[609],[610],[611],[612],[613]]},{"type":"MultiLineString","arcs":[[614,615],[616]]},{"type":"LineString","arcs":[617]},{"type":"MultiLineString","arcs":[[618],[619]]},{"type":"LineString","arcs":[620]},{"type":"LineString","arcs":[621]},{"type":"MultiLineString","arcs":[[622],[623,624]]},{"type":"LineString","arcs":[625]},{"type":"MultiLineString","arcs":[[626],[627],[628],[629]]},{"type":"LineString","arcs":[630]},{"type":"MultiLineString","arcs":[[631,632,633],[634]]},{"type":"MultiLineString","arcs":[[635],[636],[637],[638],[639],[640],[641],[642]]},{"type":"LineString","arcs":[643,644]},{"type":"LineString","arcs":[645]},{"type":"LineString","arcs":[646]},{"type":"MultiLineString","arcs":[[647],[648],[649]]},{"type":"MultiLineString","arcs":[[650],[651]]},{"type":"MultiLineString","arcs":[[652,653],[654],[655],[656],[657]]},{"type":"MultiLineString","arcs":[[658],[659],[660],[661],[662]]},{"type":"LineString","arcs":[663]},{"type":"LineString","arcs":[664]},{"type":"MultiLineString","arcs":[[665],[666]]},{"type":"LineString","arcs":[667]},{"type":"MultiLineString","arcs":[[668],[669]]},{"type":"MultiLineString","arcs":[[670],[671],[672]]},{"type":"MultiLineString","arcs":[[673],[674],[675],[676],[677,678],[679],[680,681],[682],[683],[684],[685,686]]},{"type":"MultiLineString","arcs":[[687],[688]]},{"type":"MultiLineString","arcs":[[689],[690],[691],[692],[693]]},{"type":"LineString","arcs":[694]},{"type":"LineString","arcs":[695]},{"type":"LineString","arcs":[696]},{"type":"LineString","arcs":[697]},{"type":"LineString","arcs":[698]},{"type":"LineString","arcs":[699]},{"type":"LineString","arcs":[700]},{"type":"LineString","arcs":[701]},{"type":"LineString","arcs":[702]},{"type":"LineString","arcs":[703]},{"type":"MultiLineString","arcs":[[704],[705],[706],[707],[708]]},{"type":"MultiLineString","arcs":[[709],[710],[711]]},{"type":"MultiLineString","arcs":[[712],[713],[714]]},{"type":"LineString","arcs":[715]},{"type":"LineString","arcs":[716,717]},{"type":"MultiLineString","arcs":[[718],[719],[720]]},{"type":"MultiLineString","arcs":[[721],[722],[723],[724],[725],[726]]},{"type":"LineString","arcs":[727,728]},{"type":"LineString","arcs":[729]},{"type":"MultiLineString","arcs":[[730],[731],[732]]},{"type":"LineString","arcs":[733]},{"type":"MultiLineString","arcs":[[734],[735],[736],[737]]},{"type":"LineString","arcs":[738]},{"type":"MultiLineString","arcs":[[739],[740]]},{"type":"LineString","arcs":[741]},{"type":"LineString","arcs":[742]},{"type":"LineString","arcs":[743]},{"type":"MultiLineString","arcs":[[744],[745],[746]]},{"type":"MultiLineString","arcs":[[747],[748],[749]]},{"type":"LineString","arcs":[750]},{"type":"LineString","arcs":[751]},{"type":"LineString","arcs":[752]},{"type":"MultiLineString","arcs":[[753],[754],[755],[756]]},{"type":"MultiLineString","arcs":[[757],[758],[759]]},{"type":"MultiLineString","arcs":[[760],[761]]},{"type":"LineString","arcs":[762]},{"type":"MultiLineString","arcs":[[763],[764],[765]]},{"type":"MultiLineString","arcs":[[766],[767]]},{"type":"LineString","arcs":[768]},{"type":"LineString","arcs":[769]},{"type":"LineString","arcs":[770]},{"type":"MultiLineString","arcs":[[771],[772]]},{"type":"MultiLineString","arcs":[[773],[774]]},{"type":"LineString","arcs":[775]},{"type":"LineString","arcs":[776]},{"type":"MultiLineString","arcs":[[777],[778],[779],[780]]},{"type":"MultiLineString","arcs":[[781],[782]]},{"type":"LineString","arcs":[783]},{"type":"MultiLineString","arcs":[[784],[785]]},{"type":"LineString","arcs":[786]},{"type":"LineString","arcs":[787]},{"type":"LineString","arcs":[788]},{"type":"LineString","arcs":[789]},{"type":"MultiLineString","arcs":[[790],[791,792]]},{"type":"MultiLineString","arcs":[[793],[794]]},{"type":"LineString","arcs":[795]},{"type":"MultiLineString","arcs":[[796],[797],[798]]},{"type":"LineString","arcs":[799,800]},{"type":"LineString","arcs":[801]},{"type":"MultiLineString","arcs":[[802],[803],[804]]},{"type":"MultiLineString","arcs":[[805],[806]]},{"type":"LineString","arcs":[807]},{"type":"MultiLineString","arcs":[[808],[809],[810]]},{"type":"LineString","arcs":[811]},{"type":"LineString","arcs":[812]},{"type":"LineString","arcs":[813,814]},{"type":"LineString","arcs":[815]},{"type":"LineString","arcs":[816]},{"type":"LineString","arcs":[817]},{"type":"MultiLineString","arcs":[[818],[819],[820]]},{"type":"MultiLineString","arcs":[[821],[822]]},{"type":"LineString","arcs":[823]},{"type":"MultiLineString","arcs":[[824],[825],[826],[827],[828]]},{"type":"LineString","arcs":[829,830]},{"type":"LineString","arcs":[831]},{"type":"MultiLineString","arcs":[[832],[833]]},{"type":"LineString","arcs":[834]},{"type":"MultiLineString","arcs":[[835],[836]]},{"type":"LineString","arcs":[837]},{"type":"LineString","arcs":[838]},{"type":"LineString","arcs":[839]},{"type":"LineString","arcs":[840]},{"type":"LineString","arcs":[841]},{"type":"LineString","arcs":[842]},{"type":"LineString","arcs":[843]},{"type":"LineString","arcs":[844]},{"type":"LineString","arcs":[845]},{"type":"LineString","arcs":[846]},{"type":"MultiLineString","arcs":[[847],[848]]},{"type":"LineString","arcs":[849]},{"type":"MultiLineString","arcs":[[850],[851],[852]]},{"type":"MultiLineString","arcs":[[853],[854],[855]]},{"type":"LineString","arcs":[856]},{"type":"MultiLineString","arcs":[[857],[858]]},{"type":"MultiLineString","arcs":[[859],[860],[861]]},{"type":"MultiLineString","arcs":[[862],[863]]},{"type":"MultiLineString","arcs":[[864],[865]]},{"type":"MultiLineString","arcs":[[866],[867],[868],[869],[870],[871],[872],[873],[874]]},{"type":"LineString","arcs":[875]},{"type":"MultiLineString","arcs":[[876],[877],[878]]},{"type":"MultiLineString","arcs":[[879],[880]]},{"type":"MultiLineString","arcs":[[881],[882]]},{"type":"LineString","arcs":[883]},{"type":"LineString","arcs":[884]},{"type":"MultiLineString","arcs":[[885],[886]]},{"type":"LineString","arcs":[887,888,889]},{"type":"LineString","arcs":[890]},{"type":"LineString","arcs":[891]},{"type":"MultiLineString","arcs":[[892],[893]]},{"type":"LineString","arcs":[894]},{"type":"MultiLineString","arcs":[[895],[896],[897],[898],[899],[900]]},{"type":"LineString","arcs":[901]},{"type":"MultiLineString","arcs":[[902],[903]]},{"type":"LineString","arcs":[904]},{"type":"LineString","arcs":[905]},{"type":"LineString","arcs":[906]},{"type":"MultiLineString","arcs":[[907,908],[909],[910]]},{"type":"LineString","arcs":[911]},{"type":"LineString","arcs":[912]},{"type":"LineString","arcs":[913]},{"type":"MultiLineString","arcs":[[914],[915],[916]]},{"type":"MultiLineString","arcs":[[917,918],[919],[920],[921,922,923],[924]]},{"type":"LineString","arcs":[925]},{"type":"MultiLineString","arcs":[[926],[927],[928,929],[930,931],[932]]},{"type":"LineString","arcs":[933]},{"type":"LineString","arcs":[934]},{"type":"MultiLineString","arcs":[[935],[936],[937]]},{"type":"LineString","arcs":[938]},{"type":"LineString","arcs":[939]},{"type":"LineString","arcs":[940]},{"type":"LineString","arcs":[941]},{"type":"MultiLineString","arcs":[[942],[943],[944],[945],[946],[947],[948],[949],[950],[951],[952]]},{"type":"MultiLineString","arcs":[[953],[954],[955]]},{"type":"LineString","arcs":[956]},{"type":"LineString","arcs":[957]},{"type":null},{"type":"MultiLineString","arcs":[[958],[959]]}]}}}
},{}],21:[function(_dereq_,module,exports){
module.exports={"type":"Topology","arcs":[[[40717,16964],[-28,-64],[20,-69],[-58,-14],[-24,22],[-22,-43],[-34,-18],[-14,82],[25,-3],[32,45],[31,13],[21,30],[51,22],[0,-3]],[[36552,17932],[-7,-80],[33,-78],[21,-184],[-5,-69],[24,-139],[23,-26],[28,37],[30,16],[20,-60],[72,-104],[-2,-122],[21,-134],[-3,-80],[53,-152],[26,-130],[-11,-147],[34,-64],[-4,-66],[21,-55],[74,-72],[32,3],[19,-56],[37,-57],[51,-50],[27,-13],[18,-33],[-9,-67],[12,-36],[47,-89],[29,-147],[15,-113],[12,-42],[43,63],[17,-48],[35,-43],[25,-4],[2,-138],[7,-64],[75,-128],[21,-6],[24,-37],[26,-91],[37,-47],[22,-92],[29,-53],[1,-59],[28,-63],[-9,-78],[9,-156],[-6,-49],[52,-240],[4,-141],[-29,-100],[-9,-139],[-27,-154],[2,-79],[-12,-121],[-18,-81],[-26,-59],[-10,-91],[-34,-73],[-40,-33],[-37,-102],[-29,-157],[-30,-62],[-31,-202],[-14,-8],[-42,-143],[-15,-164],[-12,-68],[3,-95],[-6,-62],[-71,-63],[-118,-7],[-43,-23],[-55,-66],[-61,-102],[-72,-27],[-32,-28],[-16,54],[-45,31],[17,32],[-28,36],[-12,-34],[-32,10],[24,58],[-26,44],[-40,-41],[15,-31],[-64,-58],[-33,-54],[-31,-28],[-35,18],[-70,67],[-48,12],[-52,29],[-34,-21],[-47,66],[-43,11],[-27,30],[-20,59],[-49,91],[-5,43],[14,91],[-35,130],[-58,94],[33,60],[-15,18],[-76,-68],[-38,7],[24,66],[12,69],[-2,60],[-45,137],[-19,-66],[-6,-63],[-20,-95],[-26,3],[-66,-25],[15,74],[43,1],[11,72],[0,101],[18,69],[32,64],[-9,87],[16,25],[-9,74],[-53,-97],[-24,-100],[-68,-62],[-23,-30],[-88,-209],[-23,82],[-30,155],[-33,64],[-11,69],[-21,32],[-35,5],[-22,96],[15,46],[-77,84],[-39,0],[-51,52],[-62,-11],[-55,71],[-66,46],[-41,-25],[-73,6],[-134,-29],[-100,-83],[-84,-46],[-61,-10],[-73,12],[-25,-9],[-51,-60],[-80,-75],[-44,-17],[-27,-39],[-30,-108],[-24,-55],[-50,-35],[-49,23],[-71,-23],[-10,26],[-74,11],[-67,-10],[-45,-21],[-64,-3],[-24,-29],[-21,-61],[-63,-26],[-86,-116],[-63,-26],[-120,26],[-60,44],[-30,62],[-50,51],[-31,11],[-1,170],[21,-29],[38,25],[19,78],[-8,121],[10,23],[4,156],[-5,44],[-59,204],[-20,136],[-4,180],[-12,67],[-26,62],[-10,75],[-42,107],[-8,123],[-8,45],[-36,115],[-54,139],[18,32],[22,-104],[17,-9],[14,60],[-30,54],[-22,89],[16,20],[39,-93],[7,-53],[31,-7],[-1,101],[-25,70],[-36,131],[-29,124],[1,70],[15,92],[23,72],[4,91],[-13,88],[32,161],[23,-88],[24,-5],[25,92],[29,47],[68,57],[63,106],[78,86],[49,3],[30,-18],[32,18],[58,61],[63,26],[40,61],[54,-9],[21,18],[48,11],[77,56],[34,43],[36,86],[38,144],[20,19],[40,82],[-14,16],[-10,97],[1,54],[31,78],[33,42],[28,81],[19,-98],[44,-143],[15,106],[13,36],[-32,87],[27,68],[14,-43],[81,61],[-22,82],[52,137],[46,50],[-7,53],[56,75],[39,-25],[10,88],[33,22],[19,-33],[38,96],[44,-44],[43,-60],[24,-67],[35,-61],[-3,-67],[35,59],[111,-33],[33,33],[-16,52],[-26,39],[9,41],[28,54],[15,93],[33,29],[14,33],[-14,36],[4,43],[26,62],[25,10],[6,55],[40,15],[22,38],[25,-23],[50,10],[55,-1],[24,29],[10,106],[-22,36],[-39,-2],[-36,53],[17,21],[34,-6],[47,-69],[20,26],[19,-14],[18,-58],[39,-25],[42,-4],[38,-40],[34,-12],[22,17],[34,-48],[21,-6],[65,72],[24,-65],[9,-54],[20,-2],[1,69],[31,40],[20,-62],[26,-28],[-47,-100],[7,-50],[-15,-51],[-19,20],[-42,-38],[7,-116],[-12,-79],[-54,-139],[14,-56],[77,-93],[9,-38],[33,-31],[23,-43],[25,3],[33,-43],[46,-38],[24,-56],[37,-56],[43,-13],[44,-28],[27,-99],[22,-12],[67,-74],[54,18],[36,48],[15,92],[29,85],[22,132],[4,107],[20,126],[-12,135],[8,73],[-14,81],[20,124],[-4,73],[27,84],[-19,20],[14,94],[17,42],[21,143],[3,76],[33,54],[30,-69],[13,-68],[4,-118],[34,-31]],[[40522,16687],[36,-64],[8,-99],[-38,-12],[-19,-27],[-36,-1],[-59,41],[-6,49],[15,54],[25,39],[64,34],[10,-14]],[[38930,16027],[27,-8],[61,-92],[24,-20],[53,-115],[73,-85],[44,-72],[28,-32],[3,-53],[-22,-13],[-67,66],[-10,31],[-96,99],[-36,54],[-53,106],[-32,80],[3,54]],[[34805,20874],[17,-4],[10,48],[22,27],[-1,36],[31,46],[30,13],[1,-107],[-46,-53],[-4,-33],[40,-42],[9,-44],[-80,26],[-11,-40],[0,-57],[10,-68],[31,-111],[-23,5],[-18,65],[-22,42],[2,122],[-18,46],[-4,100],[-11,75],[13,50],[11,86],[30,68],[-1,-70],[14,-30],[0,-85],[-41,-73],[9,-38]],[[28184,9395],[39,31],[1,-65],[35,8],[53,38],[27,-33],[-30,-51],[-43,17],[-7,-45],[-27,-24],[-84,47],[-11,87],[47,-10]],[[27526,39221],[-29,61],[256,21],[-94,-69],[-133,-13]],[[26914,39088],[-64,4],[-36,81],[159,-6],[-59,-79]],[[17482,36974],[3,-62],[-92,-45],[-44,22],[-146,-22],[22,101],[82,-8],[131,42],[44,-28]],[[38427,18548],[25,-37],[43,3],[31,-39],[22,-62],[-20,-16],[-37,25],[-53,7],[-26,60],[15,59]],[[33627,22618],[4,63],[61,121],[26,24],[102,232],[30,50],[-3,63],[23,103],[7,-79],[18,-92],[-10,-33],[-35,-36],[-11,-44],[-46,-33],[-1,-34],[-38,-118],[-34,-35],[-17,-52],[-53,-78],[-23,-22]],[[34213,23347],[49,-19],[26,1],[-16,-96],[-24,-29],[-4,-38],[-28,-31],[-37,-16],[-23,-37],[-3,95],[10,51],[6,126],[21,29],[23,-36]],[[34285,22761],[-15,-1],[-15,60],[-34,37],[-17,49],[6,61],[44,29],[-5,96],[19,88],[31,25],[29,-17],[6,-32],[-46,-211],[-1,-59],[18,-57],[-20,-68]],[[34448,23285],[40,7],[11,-37],[-1,-98],[20,-46],[7,-74],[-26,-53],[-27,31],[-5,62],[5,79],[-14,41],[-24,-8],[-13,140],[27,-44]],[[34523,23557],[10,-47],[24,-30],[-9,-55],[6,-94],[13,-95],[-44,3],[-23,46],[-17,99],[-23,57],[-33,51],[-10,75],[30,-10],[66,10],[10,-10]],[[34011,23776],[56,-11],[36,-69],[-5,-68],[7,-46],[-16,-77],[-18,-19],[-36,67],[-18,106],[-29,66],[-12,59],[35,-8]],[[37998,19142],[-27,-40],[-43,32],[-15,45],[-1,50],[-21,23],[-26,55],[-5,98],[29,1],[53,-139],[40,-54],[16,-71]],[[37659,19586],[-32,59],[2,50],[-11,66],[-46,112],[93,-138],[12,-33],[-18,-116]],[[36027,18778],[-27,-31],[-48,6],[-21,27],[37,143],[32,47],[54,11],[25,-71],[-23,-83],[-29,-49]],[[35034,20020],[26,-25],[45,-3],[21,-33],[12,-60],[21,-41],[-6,-66],[-89,88],[-20,34],[-37,0],[-5,-29],[-57,29],[-12,21],[-25,-45],[-25,4],[-16,40],[-22,13],[13,73],[97,6],[7,-22],[35,33],[37,-17]],[[34707,19969],[41,-69],[1,-56],[-62,-44],[-53,50],[-18,43],[-4,57],[60,24],[35,-5]],[[34246,18700],[-36,-26],[-42,-5],[-50,-35],[-27,19],[-10,-23],[-33,-4],[-55,30],[-48,6],[-24,-19],[-12,37],[13,58],[35,36],[43,11],[75,-54],[19,-23],[59,28],[34,-39],[19,6],[19,45],[28,23],[-7,-71]],[[33932,18525],[32,-63],[24,-6],[37,-84],[-15,-38],[-29,-21],[-34,22],[-24,54],[-37,44],[-59,15],[-14,43],[26,31],[48,8],[45,-5]],[[33732,18768],[28,10],[56,-15],[2,-95],[-59,-24],[-12,42],[-21,-39],[-128,-59],[-30,21],[5,109],[37,38],[46,-13],[27,-66],[18,-4],[30,32],[-47,57],[-7,44],[41,6],[14,-44]],[[33551,18700],[-30,-73],[-34,42],[-2,71],[39,54],[37,-42],[-10,-52]],[[33416,18805],[29,-58],[-42,-48],[-11,-42],[-45,85],[-25,12],[-16,59],[59,-13],[18,25],[33,-20]],[[32353,20295],[13,-45],[5,-74],[18,-64],[51,-25],[-16,-30],[-15,-80],[-62,53],[-6,76],[-18,71],[-26,24],[-24,-11],[-24,19],[28,53],[-1,36],[25,30],[36,5],[16,-38]],[[22917,34148],[76,-32],[-61,-51],[-68,-33],[-31,33],[23,50],[61,33]],[[22802,38737],[50,-4],[29,-80],[130,-18],[-38,-37],[127,-21],[-45,-92],[-118,-51],[-84,52],[-114,-14],[63,109],[-64,33],[-92,96],[156,27]],[[27390,39252],[-11,-50],[-115,-46],[-198,24],[33,68],[226,17],[65,-13]],[[25724,39256],[-15,-57],[-132,-36],[-115,47],[262,46]],[[26045,39273],[161,-55],[-84,-34],[-233,-39],[13,-37],[-156,-22],[-112,43],[115,46],[172,13],[7,61],[117,24]],[[26044,36575],[-73,-75],[-108,-29],[-44,71],[5,52],[35,39],[67,17],[118,-75]],[[27196,36747],[-2,-48],[-97,3],[-71,36],[-54,86],[60,45],[164,-122]],[[28352,37473],[-85,-2],[8,63],[107,36],[78,-78],[-108,-19]],[[30841,38988],[-176,50],[119,32],[184,-32],[-127,-50]],[[33183,37773],[-68,-70],[-129,41],[50,63],[147,-34]],[[36307,37681],[-57,-18],[-46,61],[85,20],[18,-63]],[[23056,28844],[36,13],[16,-53],[47,14],[85,-27],[29,10],[1,-48],[63,29],[-6,-51],[-164,-26],[-10,36],[-28,16],[-102,31],[1,55],[32,1]],[[23007,29630],[6,-25],[46,-37],[25,-4],[20,-104],[33,-16],[-7,-41],[-33,34],[-19,62],[-32,3],[-58,92],[19,36]],[[24256,28857],[-59,-69],[13,-70],[-40,-4],[-21,-38],[-45,-25],[-20,-28],[-55,35],[-15,52],[38,52],[25,0],[7,48],[59,-12],[113,59]],[[22516,33967],[-30,-29],[-12,-107],[-47,-37],[-21,44],[4,65],[45,63],[61,1]],[[21780,33496],[-40,-73],[22,-41],[-36,-23],[-5,-86],[-35,23],[-9,63],[-42,4],[-35,119],[39,7],[17,44],[45,5],[39,35],[41,-13],[-1,-64]],[[21562,33456],[20,-66],[-4,-44],[-38,-19],[-52,26],[-14,81],[88,22]],[[20714,29821],[36,-21],[-25,-72],[-19,-19],[-34,25],[-16,33],[-29,13],[46,56],[41,-15]],[[6723,3732],[32,-54],[-81,-38],[-8,-41],[-138,-16],[-65,14],[26,87],[-55,65],[289,-17]],[[6070,3840],[67,-27],[84,-98],[98,-1],[33,-90],[-117,4],[-181,91],[-3,26],[-103,42],[16,49],[106,4]],[[12874,2280],[-83,38],[-8,79],[109,-13],[58,-76],[-76,-28]],[[13816,5961],[67,-10],[-8,-32],[-111,-18],[18,50],[34,10]],[[13213,5866],[17,-33],[-70,-44],[-59,51],[94,72],[18,-46]],[[12022,4449],[-56,-66],[-223,-48],[-37,44],[146,39],[23,37],[110,13],[37,-19]],[[11877,4658],[62,-51],[-46,-53],[-99,19],[5,49],[78,36]],[[9264,4156],[70,-21],[68,35],[84,-10],[59,-40],[4,-79],[-54,-38],[-136,13],[-131,-6],[-202,63],[-234,31],[16,35],[231,34],[28,-48],[66,17],[72,48],[59,-34]],[[11949,3883],[-36,-40],[11,-78],[-150,64],[-13,57],[65,21],[7,40],[102,-17],[14,-47]],[[12669,5175],[-67,-59],[-59,30],[2,46],[85,138],[73,-79],[-55,-46],[21,-30]],[[127,37125],[59,-11],[93,-62],[-33,-44],[-190,-34],[-56,16],[0,24],[0,85],[0,17],[127,9]],[[13702,8899],[39,-9],[4,-46],[56,28],[16,-74],[-96,-52],[4,-38],[-47,-17],[-66,9],[9,48],[58,55],[-4,45],[27,51]],[[16162,8256],[45,-10],[43,-32],[60,-118],[-33,-24],[-47,82],[-26,10],[-108,90],[34,21],[32,-19]],[[12015,10719],[-39,-2],[-27,21],[16,89],[6,152],[14,96],[57,-23],[7,-131],[-36,-27],[40,-81],[-36,-51],[-2,-43]],[[13540,8854],[115,8],[-74,-124],[-37,-4],[-12,-35],[-37,-11],[-31,30],[81,65],[-5,71]],[[11935,9386],[-13,-197],[-80,36],[-1,82],[25,51],[-11,43],[27,28],[17,59],[26,-12],[10,-90]],[[12104,10390],[-42,-12],[-10,44],[27,70],[50,-49],[-25,-53]],[[13369,31959],[-48,-4],[-93,34],[-59,40],[-36,59],[-63,39],[35,27],[95,-29],[49,-26],[118,-97],[2,-43]],[[13447,31235],[34,-47],[38,48],[-4,59],[73,-10],[-17,-65],[-42,-51],[-33,-15],[-69,-4],[-24,84],[71,197],[29,41],[32,-48],[-28,-137],[-60,-52]],[[966,35301],[48,-13],[68,26],[47,-47],[67,-33],[80,-9],[-15,-32],[-105,-18],[-85,66],[-60,15],[-38,-24],[-45,19],[38,50]],[[1568,34553],[46,-16],[16,-92],[-68,-32],[-107,54],[-34,45],[68,2],[41,39],[38,0]],[[1869,33311],[13,-56],[-25,-25],[-56,-1],[-44,-45],[-40,-2],[-8,44],[41,62],[75,37],[44,-14]],[[11620,24920],[40,-16],[23,-30],[40,-25],[16,-54],[-36,-11],[-37,24],[-24,-16],[-16,-43],[-17,27],[-46,10],[-35,72],[-25,7],[9,52],[39,17],[69,-14]],[[12879,24917],[57,-15],[1,-32],[-40,-61],[-49,-6],[-90,10],[5,116],[109,-5],[7,-7]],[[3065,33965],[44,6],[33,-63],[-52,-53],[-67,-28],[-59,-54],[-31,33],[-31,-49],[-41,95],[19,52],[29,18],[52,-10],[70,52],[34,1]],[[3120,34088],[35,-46],[-99,-38],[-45,21],[80,84],[29,-21]],[[5544,33428],[-11,-51],[-83,-10],[-4,44],[23,39],[2,57],[41,28],[34,-65],[-2,-42]],[[5252,33623],[41,-4],[12,-53],[38,-22],[37,-70],[2,-47],[50,-63],[-1,-93],[-91,56],[-75,172],[24,28],[-33,52],[-4,44]],[[5274,33776],[47,-28],[-1,-47],[-28,-49],[-44,0],[-9,77],[-21,62],[56,-15]],[[5093,33856],[40,-145],[-2,-96],[-79,119],[15,48],[-66,66],[49,50],[43,-42]],[[5126,34042],[50,-4],[-3,-59],[37,-70],[1,-56],[-58,-65],[-17,8],[12,90],[-20,36],[-11,82],[9,38]],[[5007,34061],[16,-22],[66,-22],[-9,-124],[-47,34],[-13,-41],[-33,-8],[-19,52],[-56,69],[7,28],[88,34]],[[5819,32842],[5,-96],[-65,69],[-11,45],[16,61],[46,-39],[9,-40]],[[5355,33115],[108,-32],[-35,-166],[-58,-16],[-26,27],[-34,107],[-1,84],[46,-4]],[[5457,32901],[13,-63],[-10,-38],[-50,6],[-20,82],[67,13]],[[13141,31356],[15,21],[62,-32],[110,15],[15,-15],[-59,-59],[2,-43],[-39,5],[-27,58],[-18,-23],[-41,18],[-20,55]],[[12157,30096],[20,-21],[-97,-55],[-65,-15],[-12,46],[34,30],[120,15]],[[37541,19691],[22,20],[33,-30],[-1,-89],[-16,-52],[-27,-11],[-3,-33],[17,-48],[-31,-48],[-40,3],[-32,-85],[-48,-34],[-43,-47],[-87,-4],[-31,49],[-29,-11],[-36,48],[-46,35],[-6,51],[74,14],[28,-21],[66,11],[31,22],[22,-24],[37,3],[37,22],[21,72],[21,35],[39,18],[-1,56],[-13,67],[5,33],[37,-22]],[[34756,18743],[-43,-67],[-28,-15],[-32,-40],[-25,-9],[-86,-65],[-46,-66],[-4,-23],[-60,-111],[-29,-8],[-36,-37],[-24,8],[8,53],[-14,26],[13,81],[31,55],[25,27],[42,23],[51,68],[17,49],[71,35],[93,8],[39,33],[37,-25]],[[36763,11307],[27,4],[61,-58],[56,-32],[29,5],[47,34],[53,-2],[19,36],[39,15],[37,-39],[1,-251],[-32,-50],[-11,-72],[7,-134],[-23,-15],[-23,86],[-19,-11],[-19,-74],[-19,-14],[-23,-77],[-36,24],[-58,-9],[-7,62],[-20,8],[-35,73],[-25,87],[-4,121],[-52,138],[-15,101],[20,59],[25,-15]],[[31276,39307],[148,-64],[-92,-60],[17,-68],[-295,-26],[-78,-27],[-189,39],[-77,41],[119,40],[78,60],[-23,45],[310,67],[82,-47]],[[31406,39096],[183,-32],[87,-55],[-43,-105],[28,-83],[-56,-29],[-213,-2],[-136,44],[-192,25],[-49,64],[-130,21],[217,138],[304,14]],[[31995,38889],[89,-24],[88,-62],[79,-14],[19,-74],[-55,-34],[-126,-21],[-284,-15],[-127,-50],[-90,14],[83,69],[30,78],[148,155],[146,-22]],[[36198,38102],[102,-30],[10,66],[51,34],[110,-53],[139,-9],[177,-64],[-42,-73],[-97,-52],[-196,-44],[-34,32],[-195,-33],[-75,23],[-57,-67],[-114,33],[-123,100],[32,26],[-2,92],[26,37],[149,67],[139,-85]],[[36440,37657],[131,-75],[-17,-80],[-181,21],[-105,32],[48,96],[124,6]],[[36961,37996],[186,10],[4,-32],[185,-13],[61,-63],[-119,-39],[-170,12],[-220,86],[73,39]],[[39956,12652],[65,-21],[54,-51],[25,-77],[-17,-44],[46,-123],[2,-70],[-11,-52],[65,-35],[28,-48],[-9,167],[47,-111],[12,-110],[14,-48],[75,-56],[64,-23],[76,100],[30,-3],[29,-29],[-16,-61],[-14,-136],[-33,-39],[-1,-97],[-63,14],[-32,-24],[-19,-42],[17,-70],[-30,-111],[-52,-117],[-45,-126],[-76,-91],[-17,44],[-29,-3],[-31,33],[60,153],[10,77],[-27,77],[-50,32],[-25,39],[-47,31],[-19,44],[9,41],[62,40],[23,41],[13,129],[24,96],[-22,82],[6,116],[-35,1],[-10,49],[7,64],[-22,46],[-35,17],[-94,223],[9,43],[-55,129],[29,-18],[35,-92]],[[39938,11194],[13,-1],[58,68],[64,-7],[-8,-67],[-16,-45],[21,-54],[-83,-179],[-37,-105],[-57,-65],[10,-83],[31,-13],[-18,-49],[-49,38],[-11,-24],[-82,-66],[-25,-4],[-24,-89],[-13,-119],[-20,-41],[-27,-107],[9,-43],[-41,-16],[-82,-140],[-67,-18],[-38,14],[-43,-9],[-22,56],[-39,-1],[-18,40],[-36,-11],[-72,10],[12,97],[-11,53],[52,134],[67,84],[65,118],[50,21],[41,49],[56,43],[72,111],[50,41],[58,107],[26,154],[39,32],[20,48],[15,114],[15,44],[42,55],[14,-59],[25,-17],[14,-99]],[[34610,22820],[21,-10],[15,-99],[-21,-50],[26,-33],[10,-65],[-2,-85],[15,-35],[1,-99],[-16,-54],[-40,-64],[-29,138],[-48,-124],[-3,-24],[24,-52],[9,-112],[-24,-72],[-19,-7],[-13,95],[-16,-40],[-45,29],[-48,54],[-15,39],[-10,136],[25,92],[-27,62],[-35,35],[-19,-2],[-12,-92],[-24,28],[-38,0],[-5,44],[-34,-11],[-36,-154],[-24,-8],[-4,53],[13,38],[10,102],[23,50],[65,30],[12,55],[37,40],[10,30],[28,-19],[19,-43],[3,-56],[40,19],[25,77],[32,-10],[16,94],[30,-24],[8,37],[33,-3],[4,29],[-14,121],[12,21],[40,-57],[15,-44]],[[34055,24956],[85,-74],[51,25],[-14,-78],[-3,-70],[10,-74],[32,-72],[-11,-69],[-33,-147],[-39,-24],[-22,-33],[2,-61],[-25,-79],[35,-135],[-6,-59],[17,-84],[47,-43],[-1,50],[33,40],[41,-16],[42,-116],[23,51],[33,-19],[-13,-82],[24,-59],[-4,-37],[40,-17],[-9,-108],[-21,28],[5,58],[-33,-7],[-36,31],[-16,91],[-31,35],[-34,72],[-10,-57],[19,-61],[-18,-31],[-12,54],[-38,68],[-33,34],[-31,-22],[-7,-30],[-27,-16],[-41,56],[-23,-19],[-2,89],[34,70],[-4,51],[-38,11],[4,-65],[-17,-6],[-21,78],[-19,13],[-17,133],[-4,94],[-14,39],[7,73],[37,-64],[26,40],[-10,71],[12,97],[2,111],[-8,46],[28,200],[24,22],[32,3]],[[35455,20512],[56,-4],[33,-142],[-18,-83],[6,-109],[34,-146],[3,57],[16,11],[9,-91],[38,-92],[50,-3],[43,80],[17,60],[26,35],[17,71],[51,18],[37,39],[-5,43],[77,82],[95,-71],[129,-128],[41,0],[53,-22],[14,-37],[50,-5],[116,-105],[78,-54],[69,-26],[57,-81],[52,-9],[69,-121],[28,-8],[49,-100],[10,-149],[62,-34],[71,-69],[39,-8],[22,-24],[26,-59],[5,-80],[-15,-14],[-67,0],[-18,-47],[26,-104],[61,-114],[45,-52],[14,-104],[23,-32],[15,-82],[23,-9],[35,17],[17,-13],[-5,-77],[31,-41],[56,-17],[-24,-33],[13,-48],[89,-56],[-20,-39],[-1,-48],[-52,13],[-30,52],[-92,22],[-26,21],[-37,-5],[-69,27],[-25,36],[-16,55],[-44,66],[-12,68],[-25,15],[-58,179],[-17,39],[-34,28],[-30,6],[-56,29],[-22,37],[-30,17],[-32,-45],[-37,20],[5,-62],[-37,-58],[-59,-26],[2,-37],[32,-76],[-3,-36],[-82,-85],[-47,37],[-69,-10],[-24,14],[-31,-16],[-53,86],[-20,52],[-44,73],[-18,46],[-63,-20],[-21,24],[-30,-41],[-8,27],[26,128],[-39,78],[-16,72],[29,18],[-48,118],[-21,151],[-18,-5],[-3,56],[-35,48],[-75,76],[-54,24],[-72,67],[-60,20],[-29,-2],[-49,59],[-9,26],[-60,65],[-20,-5],[-28,56],[-8,49],[-46,-161],[-31,-7],[-25,90],[13,35],[-13,59],[-46,73],[-33,14],[-9,29],[29,25],[56,-25],[35,69],[18,12],[47,-26],[34,36],[1,67],[-49,-29],[-59,-10],[-52,12],[-22,-6],[-32,58],[-11,99],[-72,39],[-20,-15],[16,138],[62,35],[37,57],[30,23],[26,-1],[51,-35],[45,-49]],[[34484,20907],[-22,-39],[-30,-81],[-24,-21],[-53,-17],[-55,5],[-21,36],[-90,-1],[-50,-10],[-47,13],[-47,-12],[-35,17],[-40,-15],[-25,-65],[-13,-83],[10,-105],[19,-57],[28,-32],[17,-76],[42,-9],[55,127],[51,-18],[35,41],[69,0],[32,43],[23,-18],[0,-82],[-37,30],[-28,-21],[-34,-88],[-39,-57],[-34,-22],[-16,-39],[-42,-19],[57,-89],[26,-96],[24,-36],[12,-68],[-9,-17],[-7,-79],[40,-70],[9,-37],[21,-5],[3,-56],[-79,-33],[-32,-81],[-37,20],[-12,41],[15,112],[-35,40],[-45,84],[-1,35],[17,54],[0,95],[-7,19],[-38,1],[-36,-47],[-9,-42],[15,-68],[5,-82],[-8,-87],[6,-122],[-16,-122],[13,-56],[-9,-35],[-41,-8],[-26,-26],[-39,61],[-2,26],[18,100],[11,104],[1,90],[-23,128],[-48,-14],[-14,32],[-6,55],[5,53],[-9,37],[35,65],[9,79],[17,48],[-1,120],[22,115],[23,52],[7,46],[-6,90],[17,30],[-7,46],[28,104],[24,64],[28,-35],[27,50],[18,58],[19,9],[41,-26],[17,-34],[31,5],[64,-16],[52,-39],[44,18],[64,-21],[48,43],[30,47],[7,38],[21,17],[23,54],[28,-46],[-39,-117]],[[32503,19298],[33,-48],[53,-18],[18,7],[25,-47],[18,-74],[70,-17],[23,13],[37,-21],[68,-11],[18,33],[17,77],[27,8],[20,-54],[44,5],[62,-56],[51,-8],[12,-68],[17,-19],[0,-57],[51,-38],[57,4],[36,16],[36,-32],[6,-28],[-6,-118],[22,-83],[-101,67],[-49,44],[-65,-28],[-64,19],[-68,4],[-103,36],[-64,60],[-85,42],[-62,9],[-32,-30],[-61,17],[-42,42],[-29,16],[-74,13],[-24,40],[11,43],[-65,45],[-52,18],[19,71],[15,3],[9,78],[24,47],[85,-42],[25,44],[37,-24]],[[31272,21880],[40,11],[39,-16],[40,0],[41,-75],[13,-56],[26,-51],[6,-74],[40,-38],[14,-39],[42,-40],[66,-92],[44,-123],[37,-90],[26,-31],[16,28],[25,3],[29,-56],[20,-73],[35,-17],[42,-83],[9,-62],[26,-48],[43,-15],[25,-51],[46,-3],[22,-44],[13,-56],[-40,-54],[0,-79],[11,-51],[22,-30],[54,-39],[18,4],[24,-203],[31,-39],[-20,-67],[7,-43],[33,49],[42,-5],[21,-26],[52,-141],[-22,-117],[9,-61],[-11,-65],[7,-90],[0,-104],[-8,-154],[-23,-28],[-30,58],[-31,-45],[-49,51],[-5,-88],[-60,119],[-27,71],[-103,139],[-43,74],[-47,127],[-62,99],[-31,101],[-21,32],[-30,102],[1,48],[-42,146],[-20,109],[-51,118],[-30,95],[-49,57],[-41,263],[-27,94],[-63,77],[-34,28],[-13,112],[-22,29],[-48,138],[-19,31],[-40,25],[-108,217],[-32,120],[2,64],[19,15],[50,-26],[33,-51],[41,-15]],[[33570,22216],[12,-7],[24,71],[14,-12],[5,-55],[40,-47],[-4,-122],[23,0],[21,26],[13,-44],[38,-23],[16,-39],[41,-38],[30,-4],[3,-49],[-13,-22],[-52,-32],[-33,13],[-22,-44],[42,-75],[-7,-32],[-55,-26],[-31,20],[-31,-60],[30,-64],[-11,-37],[-53,-10],[25,-52],[-4,-40],[30,-29],[3,-55],[16,-18],[32,-99],[-23,-80],[25,-61],[63,-88],[39,-78],[-10,-22],[-41,-16],[-38,14],[-21,37],[-10,-47],[-20,-24],[-25,-113],[-7,-129],[11,-102],[-36,-36],[-27,-60],[-26,-18],[-25,-48],[-15,-71],[0,-62],[16,-57],[-4,-48],[-17,-16],[-5,-74],[-43,-159],[-142,-132],[-8,13],[-2,94],[-8,48],[-51,49],[-33,-40],[-18,8],[-12,54],[-21,-14],[-35,72],[-7,-59],[-42,-49],[-36,19],[-37,-48],[-15,-1],[0,109],[-52,29],[-36,-28],[-39,8],[-14,29],[-39,-7],[-1,54],[-17,170],[-13,20],[7,106],[-26,89],[-38,33],[-9,53],[-24,33],[-5,53],[16,68],[-36,74],[-4,101],[19,161],[35,98],[28,25],[10,-39],[30,-33],[41,1],[37,-40],[30,-7],[24,51],[13,92],[-3,68],[23,-12],[0,58],[32,50],[44,14],[70,36],[35,31],[46,123],[54,116],[14,80],[28,3],[79,96],[3,-29],[40,8],[16,27],[-11,83],[15,36],[27,-7],[30,79],[9,57],[40,90],[29,105],[7,-66]],[[28023,38195],[-168,-58],[-296,-72],[-261,-82],[-128,-119],[-111,-23],[-80,-45],[-11,-85],[-76,-26],[12,-38],[-103,-111],[-60,-16],[-161,36],[-80,-24],[-61,96],[100,44],[106,154],[104,75],[66,98],[79,29],[124,80],[246,60],[17,38],[200,-11],[168,33],[159,54],[22,30],[167,60],[107,-17],[47,-82],[-128,-78]],[[26615,37522],[126,-24],[-39,-95],[-41,-5],[-42,-69],[19,-56],[-25,-66],[84,-136],[96,-96],[83,-45],[-54,-32],[-165,23],[-123,-2],[-99,30],[-16,55],[34,33],[-54,33],[-5,46],[-175,-11],[-43,66],[17,67],[75,14],[42,39],[34,84],[-44,16],[81,54],[-6,40],[63,34],[119,22],[58,-19]],[[36505,33177],[25,-59],[-7,-79],[34,-115],[12,-76],[0,-81],[-17,-60],[-2,-61],[18,-115],[16,-26],[40,-273],[27,-89],[25,-134],[-17,-24],[-44,24],[-71,-26],[-60,-259],[-2,-77],[13,-46],[38,-72],[31,-153],[-51,11],[-29,22],[-19,-34],[-22,-107],[-28,-17],[-15,101],[23,158],[-8,103],[25,98],[-6,63],[-30,106],[17,75],[14,113],[1,149],[-9,68],[15,163],[-22,68],[-32,50],[-7,123],[16,65],[8,111],[-6,69],[36,36],[44,-11],[16,120],[-38,72],[48,26]],[[36625,30815],[32,-3],[43,-35],[36,-2],[62,66],[-28,-107],[28,-135],[-36,-33],[-45,-20],[-38,9],[-37,-24],[-61,-106],[-22,-96],[-82,59],[-75,74],[-50,-8],[-48,-47],[-31,49],[-26,1],[-18,-51],[24,-47],[23,-3],[47,-72],[-17,-16],[-38,18],[-32,-68],[-26,-22],[-18,36],[13,77],[-32,109],[8,60],[37,34],[30,58],[0,66],[33,-28],[59,-4],[13,42],[-2,60],[28,87],[15,160],[-22,101],[10,56],[30,25],[54,-88],[33,-70],[66,-98],[60,-64]],[[36332,30185],[26,7],[-7,-71],[7,-111],[38,-74],[21,-102],[-1,-96],[-8,-73],[-28,-31],[-21,-131],[-41,-16],[-20,-89],[12,-111],[-8,-107],[-20,-35],[-18,-79],[-1,-102],[29,-77],[-47,-49],[-5,-56],[-40,-53],[-30,-18],[3,78],[30,56],[-29,27],[-21,-58],[-45,-30],[-15,-41],[-15,-92],[-25,0],[4,64],[-25,25],[-44,-112],[-73,15],[-30,25],[-44,8],[-22,41],[-18,-70],[39,-57],[-3,-25],[-59,-34],[-47,-141],[-25,-17],[-27,15],[-32,79],[-8,90],[32,49],[-3,35],[-36,-5],[-34,31],[-60,-16],[-27,-39],[-59,-22],[-35,-30],[-55,-13],[-26,25],[-21,-29],[-17,-86],[-39,46],[-58,-24],[-35,6],[-4,66],[13,30],[40,5],[57,72],[120,180],[39,10],[13,-22],[94,18],[28,21],[81,27],[17,-51],[40,-5],[47,61],[-10,51],[78,173],[3,105],[14,42],[19,-125],[37,-16],[20,42],[91,61],[44,79],[20,65],[58,69],[20,104],[25,65],[28,134],[1,66],[-17,62],[13,68],[-11,65],[41,57],[12,88],[27,-8],[6,-69],[49,-3],[17,51],[-7,32],[-46,-24],[16,84],[33,-30]],[[35555,28550],[31,-7],[12,-94],[-41,-48],[-22,-83],[-25,46],[-37,14],[-39,-34],[-35,-119],[-55,17],[2,86],[-17,51],[32,40],[16,70],[23,22],[16,-37],[51,21],[14,50],[42,28],[32,-23]],[[35195,28399],[46,12],[14,-35],[-19,-52],[40,-5],[9,-94],[-27,-58],[-31,-163],[0,-48],[-14,-62],[-27,-34],[-35,3],[-23,-21],[-43,26],[13,72],[-15,38],[1,74],[30,49],[20,72],[-29,109],[-29,3],[6,-58],[-40,-29],[8,68],[-36,49],[9,30],[80,62],[14,46],[26,21],[27,-13],[25,-62]],[[34045,25876],[-13,-55],[-6,-105],[-17,25],[-12,74],[-29,43],[-29,140],[10,128],[53,177],[47,127],[62,56],[38,-69],[-12,-35],[1,-66],[-22,-93],[-18,-162],[-21,-105],[-32,-80]],[[32900,25272],[14,-77],[-21,-23],[-28,-88],[-14,-98],[-106,-121],[-55,34],[-37,38],[-8,86],[7,99],[55,77],[4,39],[49,32],[29,-5],[30,21],[47,-17],[11,36],[23,-33]],[[29405,22933],[30,-3],[52,-99],[75,-215],[7,-72],[26,-77],[24,-114],[-2,-89],[-25,-109],[-29,-42],[-74,-61],[-52,8],[-19,32],[-27,156],[-8,173],[-9,111],[11,-4],[15,148],[-1,47],[21,99],[-15,111]],[[25961,17822],[46,-147],[26,-222],[7,-159],[24,-96],[4,-54],[-31,-132],[-21,36],[-15,79],[-27,-25],[7,-117],[14,-41],[-8,-129],[-25,-50],[-12,-72],[5,-127],[-13,-100],[-35,-180],[-29,-192],[-22,-116],[-29,-205],[-50,-256],[-19,-177],[-21,-147],[-21,-79],[-22,-131],[-16,-44],[-34,-39],[-65,-18],[-73,-77],[-45,5],[-34,48],[-53,25],[-35,53],[-21,105],[-18,41],[-5,142],[9,48],[-18,104],[-19,44],[-15,116],[0,76],[19,93],[8,67],[33,40],[13,72],[37,112],[19,105],[5,114],[-24,81],[-1,77],[-21,104],[-7,206],[50,158],[6,111],[49,10],[29,44],[49,-2],[11,41],[72,23],[16,46],[46,65],[13,-3],[32,69],[60,92],[-4,40],[25,94],[-11,53],[16,31],[24,-29],[17,42],[45,62],[12,76],[-1,110],[35,87],[37,-81]],[[2369,2342],[-141,-9],[-177,33],[-267,87],[11,75],[99,63],[182,-41],[222,-115],[71,-93]],[[12404,2371],[-135,13],[-29,41],[59,88],[81,41],[346,120],[85,-5],[-285,-190],[-54,-94],[-68,-14]],[[15244,2571],[169,-2],[86,-151],[-3,-75],[-61,-90],[-640,-103],[-66,-33],[-496,-20],[-21,70],[98,87],[134,11],[221,152],[-45,45],[59,156],[128,128],[157,48],[98,12],[147,-22],[135,-48],[60,-45],[-16,-71],[-127,-15],[-17,-34]],[[13602,2218],[-10,-89],[-86,-49],[-266,41],[-150,4],[-194,84],[391,3],[137,15],[-43,40],[39,56],[87,34],[99,-31],[-4,-108]],[[12436,4781],[71,-90],[37,-117],[86,-176],[10,-222],[-45,-89],[-64,-74],[-357,-31],[-82,60],[135,66],[89,17],[32,32],[-112,19],[-68,37],[-58,-59],[-110,-53],[-100,20],[-64,42],[7,61],[125,53],[177,-1],[-50,60],[80,12],[72,-16],[59,37],[105,4],[-3,30],[-85,17],[-2,52],[83,40],[0,46],[-89,-8],[-83,52],[12,126],[-44,71],[30,39],[164,35],[35,-39],[7,-53]],[[12536,8577],[42,21],[58,-93],[-24,-64],[39,-13],[16,-57],[80,-111],[120,-111],[55,-28],[65,-5],[-19,-46],[-56,-8],[-76,-28],[-48,25],[-158,25],[-85,-22],[-46,12],[-26,34],[-35,-24],[-53,1],[-107,44],[62,65],[22,-34],[46,0],[-64,86],[-1,52],[38,59],[17,-83],[-28,-4],[52,-86],[49,17],[-28,45],[-19,80],[91,65],[-33,33],[-51,-19],[-42,49],[8,41],[71,68],[30,55],[38,-41]],[[11102,26001],[64,-8],[11,-16],[60,10],[32,-37],[33,0],[45,-46],[14,-38],[31,-39],[67,-9],[61,-59],[32,-48],[41,-7],[96,-135],[44,-13],[21,-22],[40,-5],[14,-26],[-5,-60],[86,-20],[42,-61],[39,-21],[-10,-49],[-43,-5],[-54,-36],[-50,-2],[-68,23],[-95,-23],[-81,-8],[18,52],[38,50],[12,44],[-13,35],[-87,16],[-56,78],[-17,106],[-21,24],[-51,-13],[-72,41],[-36,30],[-30,44],[-53,-2],[-54,30],[-42,4],[-30,43],[37,18],[-17,49],[-94,2],[-73,-107],[-59,-12],[-15,-53],[-47,-34],[-6,31],[20,34],[-4,70],[36,66],[89,69],[67,17],[64,34],[29,-6]],[[22722,39117],[167,-36],[9,71],[499,-55],[38,-59],[-135,-66],[-41,-49],[-192,-48],[-142,31],[-207,16],[-242,83],[-43,69],[228,78],[61,-35]],[[22257,39039],[90,-6],[113,-72],[58,-95],[154,-6],[5,-34],[101,-55],[-184,-28],[-150,-137],[-24,-115],[-68,-29],[-105,-188],[-264,151],[-41,48],[104,56],[-128,38],[-12,40],[148,28],[165,81],[-205,-20],[-60,-33],[-108,-8],[-156,110],[31,61],[-95,34],[-53,90],[14,63],[214,-5],[71,-67],[189,44],[86,2],[110,52]],[[21448,30072],[20,-88],[-19,-53],[7,-41],[-16,-212],[-57,16],[-20,-74],[-26,3],[-26,64],[3,118],[-5,45],[8,87],[-32,82],[1,50],[30,-8],[40,27],[52,57],[38,-42],[2,-31]],[[21431,30514],[8,-155],[-18,-92],[-24,-79],[-42,46],[-12,117],[-16,28],[17,76],[87,59]],[[22120,29460],[-38,-100],[-16,-75],[22,-102],[-20,-75],[-69,26],[-41,70],[-27,-1],[-111,108],[-25,0],[-30,57],[13,54],[21,30],[19,-34],[29,36],[22,-3],[37,-41],[93,10],[39,29],[30,-4],[52,15]],[[14394,36750],[78,-32],[-8,-73],[-165,-54],[-24,58],[-115,28],[16,77],[-19,45],[52,36],[113,-22],[72,-63]],[[5972,32314],[56,-29],[138,-45],[95,-204],[73,-47],[20,-33],[37,-115],[-11,-57],[-124,61],[-49,37],[22,58],[-60,-16],[-42,34],[-1,43],[-46,38],[-34,-6],[-16,106],[-48,0],[-40,66],[-45,-9],[-5,73],[-41,38],[3,50],[37,8],[81,-51]],[[12220,25210],[25,3],[20,37],[37,-11],[31,15],[36,-32],[23,-1],[35,-31],[19,9],[25,-86],[47,7],[0,-30],[72,-68],[34,-53],[3,-31],[-34,-73],[-31,43],[-95,8],[-46,-44],[-34,-8],[-19,28],[-43,-15],[-9,-46],[-33,-98],[-26,29],[-5,42],[-44,63],[-50,-2],[-42,-16],[-58,23],[-41,-14],[-15,-34],[-23,40],[-38,30],[4,64],[18,9],[72,-33],[91,-20],[50,56],[-53,91],[13,85],[-40,39],[-39,11],[0,34],[32,22],[54,-1],[47,-36],[30,-5]],[[14745,20626],[58,15],[60,-16],[24,-27],[-18,-103],[-30,-121],[-14,-36],[-35,-27],[-31,20],[-9,-46],[-25,-18],[-23,19],[-56,-18],[-11,21],[-18,105],[-1,136],[14,86],[45,36],[70,-26]],[[20007,34124],[-12,-45],[-88,-83],[14,-83],[53,25],[36,-8],[114,7],[33,-48],[-54,-145],[-38,-69],[-7,-56],[-59,-55],[15,-29],[51,17],[51,-29],[55,-76],[41,-183],[71,-62],[66,-89],[-14,-22],[64,-198],[-35,-58],[26,-23],[32,37],[56,-2],[68,-47],[10,-66],[-21,-88],[-32,-55],[-36,-9],[-8,-54],[-23,-42],[83,-6],[-2,-41],[-50,-59],[-85,-38],[-47,12],[-65,-11],[-57,21],[-26,-25],[-58,-5],[-1,-28],[-58,3],[-51,23],[-46,-20],[-20,-71],[-24,-21],[-42,37],[-63,-23],[-56,-62],[-14,52],[86,121],[7,47],[38,48],[39,10],[80,-6],[-18,42],[-30,6],[-94,75],[-58,-26],[-25,18],[4,67],[80,46],[45,79],[-13,86],[-61,-6],[60,76],[118,47],[16,71],[-28,90],[-45,78],[15,114],[-59,-42],[-63,-5],[-70,23],[47,130],[-20,63],[6,81],[-46,-28],[-35,-124],[-24,-6],[18,159],[25,106],[-62,21],[42,144],[-29,47],[24,102],[31,82],[37,5],[-7,60],[40,-3],[188,29],[-6,-27]],[[19547,33329],[26,29],[92,8],[75,-165],[-16,-52],[-37,-43],[-46,-23],[23,-94],[13,-149],[-22,-89],[-39,-81],[-25,5],[-84,-26],[-15,-27],[-130,-91],[-55,-20],[-61,-3],[29,47],[-90,41],[47,62],[6,61],[36,36],[33,124],[-56,68],[-10,45],[3,95],[-16,62],[77,6],[82,-16],[-4,27],[45,37],[-61,40],[44,47],[12,59],[53,13],[63,35],[32,-27],[-24,-41]],[[23,36519],[104,-43],[21,-48],[50,-1],[95,-68],[233,-130],[65,-206],[61,-35],[26,34],[-50,75],[87,9],[117,-50],[96,1],[49,-58],[132,-110],[47,-13],[-86,-64],[-14,-56],[-175,-45],[0,-86],[-81,-79],[50,-37],[-59,-57],[-82,15],[-66,63],[-82,41],[-40,-3],[-52,48],[-27,106],[-109,32],[-49,-25],[-93,-4],[-22,59],[-42,45],[-4,56],[-80,-12],[-12,-75],[42,-66],[-40,-76],[-33,-27],[0,45],[0,11],[0,56],[0,57],[0,56],[0,56],[0,56],[0,57],[0,56],[0,56],[0,56],[0,57],[0,56],[0,56],[0,56],[0,56],[0,57],[23,-10]],[[18601,35896],[90,-41],[-2,-68],[129,-62],[-21,-30],[27,-66],[-111,-143],[-53,-34],[-118,-43],[-66,-56],[-133,-35],[-14,-41],[-80,-30],[-175,35],[-119,87],[-22,-18],[-131,-5],[11,36],[76,49],[-13,71],[-38,16],[-21,48],[-176,16],[152,38],[77,107],[-102,20],[-113,-40],[-36,18],[78,156],[103,-38],[-40,75],[62,24],[91,-78],[36,-49],[-18,-56],[26,-49],[56,55],[43,14],[0,72],[37,4],[49,-63],[24,69],[69,21],[98,-6],[37,-39],[48,53],[46,-17],[36,34],[-12,39],[63,15],[50,-65]],[[39256,2913],[155,-15],[102,-32],[-117,-36],[-105,8],[-75,-47],[-56,70],[96,52]],[[31780,5588],[-49,-5],[-16,54],[54,20],[38,-22],[-27,-47]],[[5534,3542],[-87,9],[9,51],[88,-21],[-10,-39]],[[6837,3651],[-29,68],[96,1],[-67,-69]],[[2209,1932],[-373,33],[148,31],[225,-64]],[[3480,2959],[-139,10],[67,46],[72,-56]],[[20072,4419],[50,-13],[-31,-51],[-49,14],[30,50]],[[22193,4593],[-69,16],[10,38],[62,16],[-3,-70]],[[23396,4507],[-96,2],[-3,40],[80,20],[19,-62]],[[16520,2454],[-288,6],[107,46],[181,-52]],[[18028,3702],[21,-99],[-48,-28],[-50,103],[77,24]],[[14078,6164],[-96,-61],[-1,47],[97,14]],[[13346,4659],[-56,100],[77,-21],[-21,-79]],[[10048,3864],[-58,128],[98,3],[4,-88],[-44,-43]],[[26814,39141],[0,-59],[-99,-4],[-20,60],[119,3]],[[11818,9474],[-16,41],[15,59],[-5,60],[20,11],[26,-93],[0,-45],[-40,-33]],[[12475,8060],[91,-23],[70,-54],[-2,-61],[-65,16],[-2,43],[-46,22],[-13,-64],[-64,69],[31,52]],[[11864,9458],[-18,-54],[-25,-3],[-15,45],[58,12]],[[26439,23582],[48,6],[37,-25],[-43,-44],[-60,-4],[-32,43],[13,30],[37,-6]],[[26714,26864],[-27,-50],[-22,25],[36,40],[13,-15]],[[39393,16925],[-31,5],[16,54],[15,-59]],[[39410,16648],[9,-59],[-31,4],[3,53],[19,2]],[[39218,17272],[7,-76],[37,5],[7,-80],[-42,-35],[-22,53],[2,45],[-14,83],[25,5]],[[39293,16980],[48,-81],[-44,-24],[-15,101],[11,4]],[[38441,18718],[-26,37],[-79,76],[-40,65],[16,35],[42,-68],[37,-29],[46,-69],[4,-47]],[[38171,18994],[-5,-21],[-39,23],[-22,33],[-24,62],[-10,62],[48,-58],[13,-51],[39,-50]],[[11925,9514],[-41,-8],[-32,111],[27,51],[27,-47],[19,-107]],[[12111,8390],[50,-24],[31,-50],[-19,-40],[-56,-34],[12,60],[-61,-18],[-44,119],[87,-13]],[[12772,8006],[-78,-24],[-34,9],[-5,67],[98,-11],[19,-41]],[[11870,9122],[-45,-11],[14,77],[51,-23],[-20,-43]],[[11896,8816],[5,-45],[28,-38],[-36,-64],[-29,110],[32,37]],[[12284,8264],[42,-18],[-14,-61],[-37,33],[-57,6],[-27,36],[24,38],[69,-34]],[[12019,10478],[-28,-23],[29,-143],[-14,-51],[-28,3],[-29,96],[-30,60],[13,40],[46,20],[22,58],[19,-60]],[[11954,10180],[-43,7],[43,113],[0,-120]],[[33235,19046],[-43,-26],[-48,2],[-32,16],[12,55],[125,6],[-14,-53]],[[3470,16614],[-29,-10],[-4,47],[26,2],[7,-39]],[[867,17585],[18,-51],[-6,-27],[-35,3],[-23,49],[46,26]],[[39511,16327],[-40,16],[18,55],[22,-71]],[[39292,15817],[-30,23],[18,39],[12,-62]],[[38520,17968],[-65,50],[9,20],[56,-70]],[[38540,18768],[28,-68],[-6,-43],[24,-37],[11,-82],[-24,2],[-19,35],[-32,180],[18,13]],[[38649,18292],[44,-15],[21,-59],[-19,-26],[-38,24],[-28,35],[20,41]],[[38202,18785],[15,-61],[-35,14],[-3,41],[-24,-2],[10,67],[37,-59]],[[39376,9911],[11,-47],[-38,-26],[-31,9],[37,103],[21,-39]],[[35821,17507],[4,-62],[17,-50],[-60,14],[7,84],[32,14]],[[35114,17995],[-30,5],[-5,42],[21,32],[14,-79]],[[35132,18065],[46,9],[27,34],[31,-57],[-28,-63],[-39,-50],[-34,43],[-25,68],[22,16]],[[36637,11461],[-7,87],[26,23],[5,-72],[-24,-38]],[[37098,11544],[33,-53],[3,-36],[-24,-27],[-25,57],[13,59]],[[35921,12467],[51,-4],[-4,-35],[-39,3],[-25,-41],[-78,9],[-24,33],[11,33],[78,36],[30,-34]],[[37672,14762],[-11,46],[7,82],[35,30],[-31,-158]],[[34312,22850],[-6,9],[8,110],[37,116],[24,113],[14,-9],[0,-78],[-11,-62],[-35,-68],[-17,-99],[-14,-32]],[[34450,22927],[-26,-36],[-48,-1],[-13,44],[40,73],[18,6],[28,-31],[1,-55]],[[34168,22156],[-16,-3],[-17,45],[29,29],[26,-23],[-22,-48]],[[34147,23790],[23,-18],[-12,-59],[-20,28],[9,49]],[[34351,23502],[22,-27],[15,-47],[-7,-34],[-63,86],[-18,47],[23,24],[28,-49]],[[34423,23811],[-20,-23],[-15,30],[9,73],[12,22],[21,-47],[-7,-55]],[[36992,20228],[16,-51],[-52,8],[-10,40],[46,3]],[[37373,20067],[-30,0],[-14,34],[21,30],[22,-20],[1,-44]],[[37384,18531],[39,-38],[2,-36],[-52,10],[11,64]],[[35671,20529],[57,-17],[44,-78],[-19,-35],[-31,22],[-51,108]],[[35681,20313],[82,-15],[-1,-47],[-79,44],[-2,18]],[[35102,20290],[7,-64],[-19,-18],[-41,14],[12,53],[41,15]],[[35175,20376],[-25,14],[-12,68],[25,16],[20,-18],[-8,-80]],[[35154,20678],[52,-34],[-20,-41],[-39,-26],[-7,34],[-51,19],[6,25],[59,23]],[[34887,21150],[-18,-4],[-9,61],[13,39],[31,30],[10,-29],[-8,-57],[-19,-40]],[[34787,20605],[13,-34],[-24,-40],[-20,33],[31,41]],[[34853,20297],[-67,-15],[-19,19],[7,44],[33,21],[46,-69]],[[35599,19367],[1,-106],[-14,-45],[-22,7],[-32,63],[16,12],[5,68],[31,63],[15,-62]],[[35575,19198],[-20,-85],[-18,-22],[-16,32],[12,137],[42,-62]],[[34290,19633],[-14,-45],[-9,-90],[26,-31],[-24,-28],[-20,-63],[-26,29],[21,77],[10,136],[20,51],[16,-36]],[[34230,19468],[-9,-27],[-32,15],[13,58],[-3,69],[37,34],[7,-72],[-17,-47],[4,-30]],[[34612,20267],[-7,-29],[-54,-5],[4,32],[57,2]],[[34493,20287],[40,-24],[-36,-31],[-19,11],[-48,-25],[-4,73],[32,12],[35,-16]],[[34294,20410],[38,-39],[-20,-39],[-25,47],[-29,-65],[-9,36],[11,57],[34,3]],[[35212,18841],[-27,31],[19,90],[34,77],[15,-70],[-7,-43],[-34,-85]],[[34700,18917],[-37,-65],[-71,31],[15,35],[27,-10],[44,31],[22,-22]],[[34264,18172],[13,49],[35,51],[5,-41],[-53,-59]],[[34448,18808],[62,-15],[-3,-34],[-74,-21],[-1,60],[16,10]],[[34416,18765],[-25,-51],[-11,45],[36,6]],[[34375,18778],[-58,-72],[-10,34],[23,40],[28,17],[17,-19]],[[33513,19790],[-24,-43],[-9,82],[12,82],[16,-12],[5,-109]],[[32767,20407],[-32,-14],[6,59],[30,-12],[-4,-33]],[[32609,21526],[-25,38],[-10,44],[22,36],[22,-49],[-9,-69]],[[31862,21156],[7,-31],[-19,-52],[-15,10],[-7,60],[34,13]],[[32011,20850],[-55,20],[1,68],[51,-62],[3,-26]],[[32187,20958],[8,-59],[-9,-29],[-36,42],[12,38],[25,8]],[[32175,20602],[1,-64],[-26,34],[25,30]],[[32597,19990],[-17,-53],[-13,23],[-37,-19],[-2,94],[8,54],[46,-7],[24,-54],[-9,-38]],[[31384,21015],[51,-113],[-6,-79],[-22,-7],[-9,54],[-23,26],[-18,110],[27,9]],[[31574,20270],[-33,26],[-30,107],[8,53],[22,12],[45,-163],[-12,-35]],[[31269,21221],[-50,54],[-24,14],[-1,60],[70,-92],[5,-36]],[[22951,34195],[-43,-31],[-15,39],[20,51],[37,-21],[1,-38]],[[26897,39415],[23,-67],[-260,-15],[156,81],[81,1]],[[26547,39316],[164,-20],[173,-54],[-224,-36],[-186,41],[73,69]],[[26369,37048],[-35,-36],[-61,62],[55,10],[41,-36]],[[29139,37289],[-83,17],[96,61],[70,-34],[-83,-44]],[[29350,37388],[-95,39],[57,46],[43,-40],[-5,-45]],[[39496,36666],[-97,19],[-63,40],[64,41],[115,-31],[-19,-69]],[[38866,34144],[14,94],[99,23],[-2,-52],[-38,-11],[-73,-54]],[[37994,32236],[-15,-23],[-62,-24],[-3,46],[52,24],[24,65],[26,-3],[-22,-85]],[[37289,31166],[-3,45],[77,86],[-17,-64],[-57,-67]],[[37166,31091],[-39,-23],[-39,-52],[-16,46],[94,29]],[[35960,33337],[6,-62],[-38,-39],[-25,83],[16,42],[41,-24]],[[36895,30903],[-11,-58],[-39,-70],[-11,13],[50,115],[11,0]],[[34646,28312],[2,55],[40,20],[20,-48],[-62,-27]],[[30846,23329],[-22,78],[19,78],[5,130],[7,23],[6,110],[24,-24],[-21,-72],[12,-92],[-14,-24],[1,-48],[-14,-54],[5,-27],[-8,-78]],[[2762,25047],[-34,13],[-1,72],[-18,85],[26,60],[-1,60],[71,-64],[15,-57],[30,-51],[-28,-47],[-29,-14],[-31,-57]],[[2659,25488],[24,5],[33,-45],[-14,-26],[-37,-6],[-6,72]],[[2511,25609],[-35,-32],[-19,61],[36,27],[18,-56]],[[2333,25718],[-27,-5],[-20,30],[24,42],[25,-1],[-2,-66]],[[39904,32859],[71,-37],[-57,-23],[-14,60]],[[1514,33063],[19,-57],[-37,-40],[-41,11],[18,101],[41,-15]],[[23508,28934],[-15,56],[23,40],[26,-32],[-34,-64]],[[23346,29715],[19,-50],[-13,-32],[-63,38],[0,35],[29,19],[28,-10]],[[23236,29865],[-10,-40],[-33,10],[0,34],[43,-4]],[[18414,28231],[56,-26],[-16,-25],[-38,17],[-2,34]],[[18511,27199],[-9,-52],[-28,-33],[-28,76],[44,17],[21,-8]],[[18753,27151],[22,123],[20,-27],[-12,-77],[-30,-19]],[[18617,27146],[1,-63],[-36,-20],[-12,48],[15,37],[32,-2]],[[18341,27225],[-18,61],[29,7],[-11,-68]],[[17707,24127],[-30,-11],[1,72],[29,-61]],[[25275,17944],[-27,34],[17,32],[10,-66]],[[24826,19260],[-29,-23],[-6,79],[20,-5],[15,-51]],[[26879,15972],[-30,-4],[0,63],[21,53],[25,-49],[-16,-63]],[[26669,15776],[-49,15],[-15,49],[9,36],[40,-1],[20,-53],[-5,-46]],[[21347,21542],[23,-30],[-27,-93],[-26,10],[30,113]],[[19705,33625],[-45,-13],[8,82],[39,-36],[-2,-33]],[[19657,34089],[-25,-79],[-50,-20],[-25,57],[96,74],[4,-32]],[[19664,33891],[1,-43],[-71,29],[58,48],[12,-34]],[[21909,33174],[-39,-31],[-21,68],[28,36],[32,-73]],[[21644,33291],[42,-19],[-32,-41],[-47,33],[37,27]],[[22629,29739],[-22,7],[-26,60],[32,5],[16,-72]],[[14136,37404],[-57,-52],[-54,50],[111,2]],[[11575,26250],[16,-87],[-29,-27],[-9,66],[-17,42],[39,6]],[[12099,25548],[-16,-49],[-56,-10],[16,59],[56,0]],[[11565,26356],[0,-56],[-34,-41],[-28,107],[30,90],[32,-100]],[[11021,25635],[-33,-29],[-26,80],[10,34],[31,-12],[18,-73]],[[13479,24009],[-27,-6],[-18,77],[34,-16],[11,-55]],[[13458,23007],[-66,-16],[15,109],[-15,48],[77,21],[-14,-39],[8,-79],[-5,-44]],[[13427,24182],[-22,64],[18,14],[4,-78]],[[13393,24357],[-20,12],[-4,55],[28,-7],[-4,-60]],[[3649,34422],[-9,20],[79,106],[15,-25],[-58,-84],[-27,-17]],[[5204,33740],[15,-61],[-38,-38],[-21,97],[44,2]],[[5416,33571],[-2,-38],[-60,31],[32,96],[35,-58],[-5,-31]],[[14496,20349],[4,53],[13,27],[1,53],[14,47],[34,25],[10,-29],[-13,-81],[-41,-74],[-22,-21]],[[14733,20740],[-11,-60],[-19,-8],[-38,17],[8,43],[44,17],[16,-9]],[[14630,20648],[-31,-45],[-8,51],[39,-6]],[[14766,20653],[-30,-7],[0,36],[24,16],[6,-45]],[[10141,20501],[-23,22],[8,37],[23,7],[8,-40],[-16,-26]],[[10035,20684],[11,-57],[23,-44],[2,-41],[18,-36],[-12,-43],[-53,-18],[-10,50],[36,41],[-26,77],[-6,60],[17,11]],[[7668,27343],[-9,-54],[-26,18],[5,73],[21,17],[9,-54]],[[10565,36995],[107,32],[97,-14],[14,114],[-67,20],[-74,72],[120,96],[-56,4],[-34,81],[43,41],[-22,35],[50,53],[266,90],[142,-22],[29,-73],[62,-45],[37,-86],[-45,-43],[-30,-80],[93,1],[90,41],[36,-28],[81,20],[-22,35],[54,41],[126,7],[193,-65],[32,-77],[72,-15],[-6,-44],[68,-29],[27,-55],[65,46],[163,-49],[76,-83],[26,-84],[69,24],[88,-18],[204,-168],[12,-72],[-96,8],[-46,-41],[52,-23],[79,-1],[67,-31],[-4,-42],[-77,3],[-30,-30],[-33,-112],[144,-29],[-30,-22],[91,-85],[105,-4],[31,-35],[56,6],[42,-44],[17,-63],[47,7],[62,-33],[-24,-48],[99,-28],[47,26],[75,-85],[-63,-77],[-76,-20],[61,-44],[-75,-91],[-80,-22],[-4,-97],[-85,-12],[-58,23],[-84,137],[11,56],[-62,18],[-70,44],[-35,64],[-60,-58],[15,-62],[-42,-26],[-80,4],[108,-105],[-2,-29],[52,-61],[72,-54],[52,0],[56,-51],[30,-143],[37,4],[19,-175],[-63,-1],[45,-78],[-56,2],[-11,-50],[-81,66],[-74,18],[-62,52],[-60,76],[-14,-43],[-74,64],[-48,-6],[88,-121],[52,-19],[79,-90],[38,-23],[41,-68],[25,-90],[-22,-10],[-134,65],[-106,19],[-80,39],[-52,78],[-78,5],[-115,64],[-82,89],[59,19],[-31,44],[-56,1],[-138,164],[-121,57],[-63,-49],[-47,19],[-32,-34],[-108,-35],[-121,29],[-43,58],[9,73],[73,52],[14,67],[95,-20],[74,-32],[96,35],[161,37],[-21,52],[-78,85],[130,123],[28,13],[28,69],[64,51],[-95,196],[-30,37],[-149,101],[-31,64],[-80,-22],[-88,-45],[-17,73],[68,5],[19,66],[-152,76],[19,40],[-76,17],[-24,80],[-50,-4],[-88,86],[-38,-38],[67,-60],[-9,-46],[-82,-19],[-168,45],[-11,-25],[-107,-32],[-107,37],[-87,-9],[-178,34],[-101,8],[-34,58],[-56,3],[-88,-37],[-108,62],[-40,55],[-18,70],[172,-28],[-3,61],[-137,18],[-86,47],[-21,105],[41,48],[-22,58],[57,91],[11,60],[63,77],[111,74],[108,25],[187,-6],[16,-26],[-126,-100],[-62,-89],[35,-94],[-2,-78],[35,-81],[114,-96],[-43,-29],[-134,-49]],[[12499,39753],[533,-32],[148,-24],[43,-48],[192,-28],[13,-37],[-102,-55],[-271,-71],[-103,-103],[-356,-139],[-87,-59],[-84,-6],[-141,-146],[-208,-27],[-27,-34],[-203,-17],[19,-50],[-131,-45],[206,-65],[-162,-164],[-125,-19],[-118,4],[-7,-98],[-90,-43],[-193,-4],[34,-38],[106,2],[13,-51],[106,10],[34,-52],[-34,-42],[-138,-60],[-134,-31],[-54,78],[-307,-14],[-142,-34],[-131,43],[-43,-27],[-194,7],[-133,20],[8,77],[124,64],[173,22],[-103,67],[-35,52],[122,40],[195,-77],[-29,85],[-76,35],[-151,22],[7,55],[77,82],[179,30],[221,-31],[24,39],[-142,8],[14,58],[-101,89],[-151,53],[13,109],[292,-20],[211,-117],[96,-7],[-221,135],[48,22],[282,35],[188,59],[-28,60],[-269,-90],[-270,14],[-77,-36],[-139,-7],[-409,64],[-116,84],[110,85],[-302,33],[131,44],[442,71],[190,61],[267,-47],[47,82],[253,77],[272,-13],[197,38],[273,-11],[58,19],[322,8],[54,-23]],[[9692,37168],[119,-70],[35,-46],[-4,-94],[67,-46],[93,-109],[-134,-106],[95,-42],[153,-7],[-30,-49],[40,-98],[-12,-89],[37,-48],[42,59],[24,113],[69,58],[115,-105],[26,-93],[-61,-27],[17,-120],[105,-134],[83,77],[19,70],[47,57],[25,129],[52,25],[46,77],[-59,36],[-14,141],[144,-3],[65,-31],[118,-2],[-3,-53],[149,-78],[-70,-48],[64,-14],[11,-44],[-82,-45],[-61,-4],[65,-134],[80,-92],[-23,-90],[-46,-19],[-86,-92],[-134,-68],[-74,-26],[-70,34],[-39,49],[-135,-2],[6,-48],[66,-31],[-5,-38],[-141,-153],[-77,-1],[-201,135],[-42,-13],[129,-123],[220,-34],[-28,-85],[-93,-148],[-62,-40],[-85,7],[-84,-13],[15,-41],[-149,-69],[66,-35],[6,-56],[-20,-38],[-66,-32],[-103,3],[18,-53],[-39,-9],[-53,-105],[-40,-37],[-6,-51],[-72,-90],[-1,-41],[-65,-166],[-16,-106],[0,-158],[8,-83],[48,-43],[114,9],[95,-340],[-19,-55],[41,-7],[129,54],[58,-4],[91,-56],[95,-30],[99,-88],[58,-94],[213,-105],[70,-73],[40,-6],[90,13],[151,-38],[40,-77],[-22,-105],[31,-124],[-2,-127],[-11,-70],[74,-121],[-7,-32],[79,-74],[35,-49],[33,-98],[60,-36],[38,91],[87,-17],[-28,54],[55,119],[-26,87],[0,52],[-20,43],[-24,155],[13,66],[-28,22],[-48,137],[86,43],[115,81],[64,73],[76,127],[15,138],[-6,109],[-36,133],[-30,59],[-97,88],[-57,65],[10,50],[76,110],[3,68],[45,28],[2,56],[-41,89],[-22,82],[-33,15],[38,54],[10,81],[25,27],[-49,47],[-21,80],[8,57],[78,50],[56,-11],[120,-48],[54,-1],[74,-30],[120,63],[175,-167],[50,-26],[-23,-35],[47,-70],[80,-24],[50,3],[72,-87],[-13,-69],[14,-130],[-6,-111],[32,-1],[-1,-57],[27,-44],[54,2],[42,-71],[85,-88],[109,77],[28,56],[36,7],[49,66],[20,153],[29,30],[28,80],[48,3],[65,-143],[44,-72],[40,-109],[67,-88],[-9,-36],[42,-79],[46,-26],[-8,-57],[13,-53],[58,-85],[-4,-76],[-54,-22],[38,-45],[48,-120],[29,21],[73,-110],[-30,-64],[59,13],[106,-25],[24,-71],[44,-15],[22,21],[90,-68],[-33,-42],[-52,-8],[-54,-70],[-24,-1],[-74,-48],[-37,0],[-36,-55],[36,-21],[57,32],[46,53],[67,41],[15,40],[88,-15],[30,-72],[-39,-49],[18,-38],[56,62],[42,6],[40,-41],[36,-86],[-4,-207],[15,-40],[-36,-48],[-109,-108],[-89,-7],[-29,-24],[-55,-6],[-74,-119],[-104,-121],[-137,-12],[-49,-22],[-22,29],[-90,16],[-59,-13],[-56,14],[-128,-7],[-46,10],[-96,-27],[-42,3],[-51,-50],[-48,-149],[-103,-35],[-74,-84],[-50,-102],[-45,-66],[-15,-61],[-61,-94],[-30,-65],[-63,-79],[-69,-25],[-91,-92],[-34,-18],[-56,-108],[-27,-6],[-36,-49],[7,-36],[47,19],[45,134],[48,40],[27,39],[67,64],[73,29],[84,63],[106,182],[12,33],[55,71],[85,80],[76,53],[157,82],[74,12],[77,-17],[65,-62],[1,-85],[-58,-74],[-56,-48],[-75,39],[-32,-46],[41,-18],[20,-49],[50,26],[59,-20],[-23,-82],[-46,-61],[55,-9],[-9,-40],[21,-50],[21,-99],[70,-17],[-15,-33],[86,-62],[67,-3],[24,-27],[60,57],[9,-36],[44,-6],[23,-57],[0,-47],[-139,-90],[-90,-46],[-27,3],[-24,-39],[-27,37],[-32,-22],[-11,-57],[-55,-100],[-52,-46],[-18,-34],[-29,9],[-16,54],[-28,5],[-7,75],[10,52],[51,90],[85,83],[51,31],[24,-33],[17,69],[-58,-1],[-34,50],[-111,-93],[-25,8],[-38,-37],[-49,-6],[-56,-86],[-61,-55],[-42,12],[-38,-45],[-39,28],[-12,-77],[-18,-26],[-64,-41],[-44,-9],[-3,-33],[-36,-63],[-34,-121],[10,-35],[-35,-78],[35,-24],[38,-138],[-29,-21],[-87,23],[-11,-59],[-84,-20],[-66,-6],[-83,-56],[-41,-42],[-32,-53],[-1,-35],[32,-29],[-25,-139],[-29,-77],[-53,-56],[-51,45],[-5,-57],[41,-136],[-36,-100],[-28,-100],[-33,-3],[26,77],[-62,94],[-6,55],[11,58],[-13,108],[-14,-5],[-13,-93],[2,-84],[27,-147],[0,-123],[-16,-33],[14,-37],[32,-32],[12,-59],[-29,-87],[-37,-57],[79,-30],[0,-58],[-45,-67],[-38,-19],[-27,-76],[2,-54],[-46,0],[-69,-97],[-30,-86],[-64,-8],[-39,-50],[-40,-120],[-84,-118],[-24,-10],[-70,-107],[-8,-38],[-31,-56],[2,-40],[-39,-150],[31,-249],[39,-171],[43,-129],[-14,-70],[48,-224],[11,-29],[10,-118],[-10,-169],[-20,-49],[-7,-66],[-42,-40],[-42,-4],[1,41],[-30,118],[-40,35],[-17,104],[-19,26],[-16,67],[-30,50],[-22,107],[-3,44],[-20,30],[20,147],[1,92],[-13,38],[-59,92],[-46,109],[-39,40],[-30,-8],[-8,-37],[-48,-31],[-58,-21],[-4,45],[-46,67],[-47,38],[-9,37],[-48,-21],[-35,7],[-21,24],[-12,-45],[-58,-11],[-27,89],[-12,-72],[-63,-3],[-24,14],[-47,-16],[-30,-41],[-72,49],[-21,-55],[26,-26],[21,7],[46,-33],[6,-33],[-21,-35],[23,-46],[57,-50],[-16,-42],[-23,12],[-59,89],[-20,-6],[-5,-57],[-32,26],[-42,-38],[-61,36],[-4,52],[-21,9],[-18,44],[-25,21],[-41,-64],[-47,9],[-57,42],[-59,-2],[-45,-23],[-90,-69],[-43,-93],[-66,-75],[-48,4],[-19,-14],[-13,-53],[-28,-33],[-39,-18],[-2,-61],[-24,-105],[-16,-25],[-3,-79],[17,-131],[30,-85],[-10,-102],[-32,-131],[-18,-143],[-11,-249],[2,-84],[-12,-73],[11,-137],[13,-98],[16,-47],[43,-180],[46,-98],[29,-73],[19,-121],[34,-66],[19,-69],[72,-13],[44,-42],[28,-78],[40,4],[72,54],[76,9],[20,33],[83,24],[19,-56],[31,-4],[29,39],[-8,63],[44,58],[25,46],[5,87],[21,42],[2,148],[15,104],[61,61],[106,32],[47,36],[38,10],[99,-39],[25,34],[23,-39],[6,-63],[-10,-61],[-40,-86],[-28,-94],[4,-46],[-29,-60],[30,-13],[-39,-261],[-10,-41],[-23,59],[-2,71],[-31,-109],[25,-2],[3,-53],[-19,-117],[7,-22],[-10,-74],[-2,-129],[-28,-78],[-15,-10],[-23,-82],[33,-29],[8,19],[46,-34],[29,41],[29,7],[14,-27],[23,10],[43,-17],[62,5],[65,51],[34,-24],[58,21],[80,-39],[57,-138],[24,17],[20,-13],[24,-56],[-16,-57],[13,-93],[-26,-79],[-17,-156],[6,-86],[0,-123],[-16,-28],[-13,-82],[13,-66],[-24,-79],[4,-39],[21,-49],[22,-104],[37,-97],[39,-86],[46,-55],[14,-91],[19,-22],[34,5],[48,-41],[58,25],[33,44],[48,30],[62,89],[63,-20],[10,-19],[48,-5],[48,-39],[29,-39],[15,-41],[40,-58],[55,-125],[9,50],[-13,71],[23,17],[46,68],[16,64],[26,38],[30,4],[-4,64],[10,61],[-9,46],[25,96],[68,114],[45,-27],[11,-41],[17,104],[16,17],[30,-15],[55,5],[66,96],[51,40],[16,69],[47,53],[25,3],[27,-22],[14,-66],[-40,-76],[-43,-21],[-9,-49],[40,-174],[-8,-49],[-33,-77],[-18,-67],[16,-60],[24,-50],[9,-61],[44,12],[23,52],[4,82],[-23,100],[-15,28],[-17,118],[8,43],[23,8],[51,48],[66,38],[21,36],[18,-20],[8,53],[-42,-11],[-11,60],[9,48],[23,19],[19,-42],[14,-99],[9,-20],[45,9],[46,-20],[48,-62],[19,-136],[41,-22],[100,32],[84,5],[44,-86],[82,-44],[32,7],[75,83],[40,-3],[-5,50],[43,-5],[90,25],[52,-12],[-16,-35],[-35,4],[-17,-34],[20,-50],[17,-5],[22,-91],[23,47],[40,-39],[28,7],[23,-49],[34,-18],[25,-45],[-30,-61],[-28,-157],[57,42],[36,-10],[16,18],[41,-22],[16,-36],[19,-3],[53,-66],[65,-122],[17,-50],[-1,-66],[21,-37],[36,-21],[42,-77],[30,-41],[32,-77],[10,13],[57,-13],[60,-33],[12,38],[113,7],[111,-48],[45,-50],[62,-32],[51,-93],[18,-18],[42,-96],[15,11],[14,-79],[23,20],[27,-51],[17,-97],[2,-89],[45,-265],[22,-69],[31,-10],[26,-29],[9,-55],[-2,-59],[-45,-75],[-32,-96],[-27,-57],[-53,-59],[-13,-71],[-34,-84],[-2,-59],[-23,-38],[-11,-50],[-24,9],[4,-64],[28,12],[74,84],[43,24],[25,-148],[33,-58],[45,42],[32,-22],[46,45],[-27,-181],[11,3],[28,135],[24,20],[32,79],[30,-8],[0,86],[38,94],[34,17],[18,-12],[29,21],[30,-28],[36,-7],[22,-44],[45,-14],[65,-73],[21,-2],[20,-78],[23,53],[33,-60],[15,-4],[31,-143],[-16,-9],[-23,-169],[32,45],[15,94],[23,10],[20,-21],[54,19],[8,29],[50,-21],[39,-45],[39,-30],[42,11],[27,-31],[36,-13],[96,32],[57,-15],[108,-121],[37,-61],[24,-14],[48,-115],[48,-86],[37,-28],[14,-46],[35,-12],[31,-31],[69,10],[49,-17],[18,-28],[17,-72],[16,-143],[12,-47],[21,-206],[-8,-105],[5,-52],[-37,-220],[-21,-69],[-48,-109],[-71,-179],[-42,-43],[-19,-34],[-54,-156],[-31,-138],[-38,-113],[-25,-58],[-29,-26],[-3,44],[-30,-4],[-5,-84],[-22,-50],[-7,-51],[7,-93],[10,-9],[-13,-143],[7,-138],[13,-140],[-28,-207],[-10,-127],[7,-89],[-38,-66],[-18,-60],[-10,-89],[4,-147],[-16,-85],[-18,-21],[-33,-127],[-11,-63],[-45,-78],[-29,-138],[7,-95],[-16,-38],[-66,-51],[-31,-63],[5,-47],[-12,-37],[-103,-4],[-65,-19],[-38,32],[-32,-23],[-55,-11],[-3,-64],[-33,-11],[-53,-69],[-62,-26],[-101,-101],[-31,-59],[-83,-116],[-6,-39],[-29,-33],[-25,9],[0,-75],[-16,-50],[-15,-204],[14,-113],[-10,-83],[2,-119],[-20,-115],[-53,-68],[-54,-113],[-32,-101],[-31,-143],[-36,-109],[-60,-133],[-1,83],[53,95],[-3,64],[-50,14],[-22,-71],[-12,-87],[-58,-76],[-25,-116],[8,-64],[-24,-63],[-36,-160],[-30,-61],[-56,-87],[-42,-138],[-43,-66],[-83,-61],[-87,36],[-51,-30],[-83,53],[-36,52],[-74,-6],[-65,130],[1,-65],[-15,-23],[27,-89],[59,-49],[52,-67],[17,-73],[-22,-50],[10,-97],[64,-65],[3,-97],[-47,-137],[-35,-69],[-17,-78],[-72,-81],[-93,-54],[-93,-38],[-145,-36],[-56,-1],[-52,18],[-31,-53],[32,-51],[-9,-104],[-17,-16],[-13,-69],[1,-60],[16,-50],[-17,-50],[-63,-50],[-92,-9],[-40,32],[-83,45],[-31,-15],[-3,-52],[16,-106],[-5,-93],[18,-44],[43,-31],[11,-31],[37,9],[34,64],[23,-67],[-11,-92],[-50,-12],[-22,66],[-37,9],[-42,-52],[45,-34],[22,-35],[-63,-54],[-33,-77],[5,-96],[-14,-99],[-32,-42],[1,-80],[-63,10],[-44,-50],[-41,-17],[-69,-164],[-1,-86],[20,-48],[69,-104],[88,-20],[29,-58],[-22,-110],[14,-27],[-66,-92],[-72,-65],[-49,-75],[-25,-68],[-11,-140],[-41,-54],[-90,-66],[-33,-131],[23,-114],[0,-44],[72,-154],[-96,24],[-36,-50],[-140,-80],[-17,-105],[-5,-127],[-34,-24],[-100,58],[-23,37],[8,50],[52,-8],[51,41],[28,61],[-27,18],[-84,-69],[-37,-43],[-10,-48],[-57,50],[19,86],[49,13],[72,63],[-94,2],[-42,-83],[-36,-37],[-77,95],[-32,120],[53,-27],[35,44],[-26,57],[-25,2],[9,104],[-104,62],[-32,88],[46,4],[-10,44],[17,64],[33,50],[5,85],[-9,47],[3,185],[-26,74],[-1,63],[82,-10],[37,-36],[-21,125],[-25,-63],[-33,-9],[-48,61],[19,63],[37,57],[1,48],[-97,53],[-16,33],[-57,-3],[82,104],[-16,65],[69,8],[34,17],[14,80],[55,-15],[21,116],[32,14],[39,42],[8,69],[-68,61],[5,62],[25,62],[-11,40],[25,94],[10,171],[33,74],[-37,16],[17,58],[-39,26],[-25,-54],[-32,-4],[-33,64],[-17,92],[33,238],[30,67],[20,130],[-28,138],[-6,59],[7,75],[-22,79],[7,117],[44,5],[11,110],[27,68],[11,95],[22,50],[-4,40],[45,113],[19,110],[6,104],[19,79],[19,38],[-9,129],[17,28],[16,72],[4,62],[-11,41],[-1,93],[-16,146],[-5,124],[4,69],[31,43],[9,113],[-4,69],[-15,35],[-4,62],[25,59],[13,67],[15,150],[15,32],[3,92],[28,197],[1,77],[-8,48],[25,94],[5,47],[-15,121],[10,198],[11,50],[-22,45],[3,72],[26,48],[27,311],[0,55],[-12,122],[6,211],[-15,124],[-10,199],[-65,108],[-45,57],[-7,60],[-15,29],[-65,67],[-41,68],[-37,21],[-68,72],[-48,34],[-37,56],[-108,115],[-49,117],[-45,61],[-40,115],[-1,43],[13,99],[-73,268],[-37,67],[-7,91],[-47,86],[-12,104],[-51,172],[-29,165],[-15,49],[-21,124],[-28,93],[-41,86],[-44,177],[-39,96],[-80,84],[-37,52],[-4,26],[33,42],[-6,66],[-26,70],[7,32],[-26,83],[6,80],[55,135],[33,54],[54,62],[12,32],[22,117],[-45,12],[-18,-41],[-19,18],[-58,101],[22,26],[-7,102],[4,57],[-12,70],[40,53],[11,60],[-3,50],[39,83],[12,96],[-7,86],[40,45],[95,52],[1,73],[-15,23],[26,51],[24,-17],[-4,123],[36,46],[20,-5],[51,90],[48,162],[9,64],[-30,46],[16,152],[-6,26],[-7,133],[-15,28],[32,56],[-25,91],[13,75],[-19,44],[-31,102],[-42,93],[-28,119],[35,73],[-41,15],[-5,42],[-64,85],[-41,2],[-27,-36],[-7,-58],[-43,-57],[-27,-14],[-11,-48],[44,-95],[7,-38],[-31,-17],[-17,-35],[-46,-12],[-22,112],[-26,-19],[-27,22],[-25,95],[-108,43],[-21,-18],[-29,25],[-4,54],[-30,7],[-14,-40],[-21,39],[14,44],[-3,53],[-29,55],[-87,85],[-6,58],[-19,-18],[-26,-55],[-23,53],[-35,21],[-19,53],[-4,61],[18,80],[-27,35],[18,43],[-24,56],[-57,94],[-33,96],[-42,64],[-61,107],[29,52],[-9,52],[-25,1],[-25,-41],[-65,1],[-41,23],[-46,44],[-60,19],[-18,29],[-72,56],[-61,-1],[-27,15],[-49,55],[-51,78],[-61,131],[-126,210],[-51,53],[-33,-19],[-27,50],[-26,-55],[-38,-47],[-85,-65],[-33,-9],[-34,17],[-42,42],[-65,13],[-43,55],[-43,23],[-28,53],[-105,42],[-38,46],[-93,65],[-85,103],[-28,63],[-41,8],[-55,24],[-84,60],[-53,116],[-56,60],[-60,50],[-19,58],[-43,95],[-21,95],[46,44],[-26,45],[29,79],[3,86],[-25,29],[-25,85],[1,78],[-17,69],[-69,130],[-60,158],[-67,110],[2,29],[-50,29],[-48,133],[-43,53],[-31,13],[-40,57],[-5,70],[26,63],[-9,52],[-23,41],[-31,-2],[-21,87],[-38,19],[-23,36],[-15,77],[9,48],[-44,6],[-23,18],[-39,95],[-24,20],[-29,81],[-25,45],[-7,58],[-59,164],[-9,73],[-33,113],[8,89],[-66,38],[-1,28],[-35,37],[-23,-28],[-45,52],[-33,15],[-5,-148],[14,-45],[14,-104],[-2,-62],[32,-93],[22,-22],[47,-85],[29,-101],[34,-29],[12,-66],[25,-20],[16,-138],[46,-69],[15,-78],[21,-50],[50,-60],[27,-134],[4,-77],[29,-60],[16,-88],[25,-82],[-7,-45],[21,-87],[60,-9],[35,-85],[4,-33],[29,-41],[-5,-61],[-52,-76],[-19,28],[-10,77],[-21,60],[-30,31],[-46,84],[-43,52],[-74,117],[-7,47],[9,102],[-14,97],[-23,70],[-31,23],[-41,62],[-21,62],[-45,-31],[-27,57],[-69,57],[-10,49],[-52,70],[48,11],[31,21],[13,32],[15,95],[-11,41],[-96,179],[-19,10],[-58,76],[-16,125],[-20,26],[-7,89],[-27,39],[-4,53],[-37,83],[4,66],[-25,34],[-45,153],[-9,100],[-16,45],[-70,98],[-37,5],[-11,63],[-37,1],[-35,20],[-14,34],[-38,37],[-99,12],[-19,25],[2,114],[-26,31],[0,35],[-55,99],[-60,123],[-5,56],[15,36],[-11,48],[-32,12],[-26,50],[-12,77],[4,68],[-57,59],[-2,41],[-80,145],[-14,106],[4,58],[-11,55],[-50,90],[-5,55],[27,110],[7,95],[-20,93],[4,45],[-16,32],[-6,106],[-15,53],[22,121],[22,81],[10,249],[15,184],[-2,129],[-14,32],[20,88],[-20,19],[-35,210],[-27,56],[-10,57],[8,52],[74,-48],[127,-18],[23,-48],[11,-109],[-29,-18],[13,-50],[39,60],[-4,103],[17,44],[-62,226],[-30,0],[-32,107],[-59,5],[-32,71],[-39,4],[-34,49],[-31,91],[-152,26],[1,63],[-27,25],[-48,-11],[-74,65],[8,76],[-52,71],[-29,85],[36,37],[-27,27],[19,97],[-39,4],[-22,84],[-72,35],[-32,-33],[-87,109],[33,94],[-44,66],[24,71],[-71,18],[-22,80],[21,45],[-7,88],[-41,63],[-57,-31],[-9,-53],[-40,17],[43,101],[-39,60],[-71,143],[-74,30],[-10,120],[-54,87],[-3,32],[-81,55],[-82,172],[51,-205],[-29,-12],[-63,74],[-21,48],[-59,0],[53,-65],[-3,-35],[-51,-25],[-109,78],[-55,82],[-37,34],[-163,105],[24,70],[-71,-28],[-71,5],[-94,51],[-144,29],[-97,-19],[-124,74],[-113,32],[-28,42],[-82,80],[-82,-6],[-85,-29],[17,-158],[-26,-41],[-47,-9],[-67,12],[-86,-100],[-44,0],[-43,-74],[-88,-13],[-12,50],[75,53],[-65,15],[3,66],[32,49],[15,124],[99,65],[50,-21],[6,83],[-75,2],[-127,-90],[-2,-36],[-59,-52],[-43,-67],[-12,-73],[-39,-10],[-111,-118],[-5,-56],[54,-18],[36,-42],[-78,-84],[-25,-76],[-86,-33],[-167,-158],[1,-52],[-125,-108],[-95,-45],[16,-55],[-58,-48],[-83,-41],[-50,-4],[-56,-54],[-72,-35],[-73,-3],[-38,-58],[-89,-43],[-11,53],[91,123],[107,69],[41,-57],[35,54],[27,71],[58,57],[57,29],[71,80],[38,60],[81,70],[-1,103],[19,28],[-9,57],[62,72],[-19,33],[-99,-55],[-32,1],[-23,47],[-40,-30],[10,-45],[-35,-12],[-67,101],[-28,-21],[-50,54],[-113,-88],[-44,-13],[-5,114],[-26,40],[22,70],[-47,136],[-37,-44],[-73,-33],[-78,-9],[-85,114],[-78,55],[62,81],[61,-39],[95,8],[-37,48],[-28,-24],[-133,23],[-44,32],[-27,89],[-34,36],[3,36],[52,15],[-10,53],[58,86],[45,34],[-6,41],[49,95],[34,11],[85,-50],[84,51],[38,61],[35,-18],[99,24],[36,60],[-24,99],[-45,44],[52,32],[5,49],[-34,30],[-65,-25],[-118,-101],[-103,49],[-67,-1],[-67,-28],[-141,28],[-38,33],[7,46],[-46,39],[18,55],[-83,19],[-78,55],[115,50],[37,40],[110,62],[124,54],[70,13],[32,-19],[-29,-72],[38,-31],[213,-6],[29,45],[150,39],[-4,35],[-68,22],[-60,-27],[-52,34],[27,64],[-39,16],[-117,-5],[-82,39],[-46,94],[-142,101],[-65,26],[-48,61],[20,106],[79,-4],[136,16],[50,23],[79,81],[24,85],[121,131],[100,-6],[191,132],[36,-18],[112,10],[91,57],[47,52],[183,-50],[43,-47],[67,-22],[108,31],[114,-29],[-27,-36],[59,-46],[82,-8],[63,21],[58,-15],[100,13],[158,-45],[19,-20],[212,-13],[71,-35],[224,25],[205,-106],[162,-8],[90,-24],[119,-84],[159,-60],[96,6],[20,91],[63,36],[142,9],[81,-17],[37,50],[102,24],[163,89],[96,-33],[-148,-84],[-103,-19],[-93,-67],[-78,-86],[76,-16],[-1,61],[90,48],[15,33],[88,-1],[74,52],[123,46],[33,-33],[137,117],[-46,45],[47,23],[70,-64],[59,-117],[62,-60],[77,-27],[3,63],[67,84],[71,-77],[-32,-63],[92,-1],[45,38],[21,62],[108,1],[125,-36],[76,-55],[161,-37],[88,-50],[109,-31],[120,-12],[49,27],[138,-72],[51,-60],[-15,-30],[-76,1],[-75,-80],[34,-24],[140,-26],[165,-5],[163,25],[107,43],[47,-57],[117,-32],[69,-77],[-53,-45],[147,-41],[-92,195],[23,69],[105,22],[119,92],[-76,-3],[-17,-36],[-115,-6],[-13,-36],[-59,-6],[-48,29],[42,77],[99,18],[144,53],[71,-46],[12,-60],[140,-98],[81,19],[89,-70],[128,-27],[125,34],[148,-27],[82,24],[109,-48],[28,54],[-139,49],[25,58],[107,31],[33,-46],[107,-22],[-38,-139],[66,-82],[59,15],[-44,109],[22,66],[81,11],[123,109],[-3,79],[-90,-34],[-10,39],[54,49],[-21,76],[-209,96],[-48,104],[36,50],[-36,61],[16,107],[136,145],[73,17],[28,-49]],[[11228,32798],[-32,3],[-94,46],[-1,52],[81,5],[42,-65],[4,-41]],[[12032,36253],[-108,10],[22,61],[107,-28],[-21,-43]],[[11550,35262],[39,-42],[-46,-40],[-60,56],[67,26]],[[11084,35143],[-13,-69],[-102,-102],[-77,-11],[-23,73],[59,98],[156,11]],[[11362,35019],[29,-38],[-4,-51],[-44,-95],[-64,57],[5,67],[41,58],[37,2]],[[8533,38403],[-75,9],[-96,125],[54,3],[118,-92],[-1,-45]],[[9559,38553],[247,-12],[-25,-61],[-279,2],[-18,53],[75,18]],[[10198,38432],[-149,28],[0,50],[111,-3],[38,-75]],[[9779,37918],[-1,-78],[-130,-11],[-178,64],[-33,41],[101,99],[141,14],[76,-56],[24,-73]],[[7938,38750],[48,-72],[-116,-10],[-85,21],[-227,-23],[3,28],[310,80],[67,-24]],[[6975,38044],[-120,9],[87,66],[113,14],[-80,-89]],[[8450,37433],[-101,88],[-84,41],[35,49],[125,17],[97,-35],[8,-60],[-80,-100]],[[9048,37669],[95,-49],[179,36],[65,-51],[-30,-65],[-59,-23],[8,-57],[44,-23],[32,-69],[69,-34],[-40,-40],[20,-110],[-109,-47],[-68,7],[-2,-51],[-53,-32],[-57,15],[-64,90],[-97,91],[-69,38],[-59,-1],[-112,108],[24,51],[66,12],[64,-72],[98,6],[23,80],[-132,72],[54,38],[7,48],[103,32]],[[10754,35673],[35,43],[59,-53],[77,-26],[84,-74],[69,-31],[30,-50],[8,-93],[101,15],[59,-75],[-84,-69],[-151,56],[-10,51],[-116,39],[-26,-65],[-66,-52],[-37,-64],[-97,-39],[-25,119],[-167,3],[34,58],[72,49],[-14,99],[47,261],[46,50],[35,-29],[1,-65],[36,-58]],[[11799,36377],[59,-20],[3,-161],[-73,-58],[-139,-4],[-34,101],[64,121],[120,21]],[[11363,37602],[141,3],[122,-39],[104,-96],[-13,-60],[-162,18],[-188,-33],[-50,25],[-27,83],[-73,36],[-3,78],[50,10],[99,-25]],[[9308,38248],[11,-141],[22,-69],[-53,-31],[9,-69],[-87,-25],[-188,0],[-54,91],[90,64],[-297,-39],[33,61],[92,49],[-80,67],[198,76],[155,0],[149,-34]],[[8661,38903],[87,-51],[88,-2],[256,-116],[-16,-64],[70,-43],[-3,-59],[-161,8],[-55,69],[-188,41],[-180,-24],[91,115],[-40,31],[-188,30],[17,67],[222,-2]],[[9966,39320],[141,-124],[226,-47],[109,-4],[79,-109],[-19,-54],[186,-30],[24,-74],[-195,-71],[-201,-162],[-206,-10],[-146,20],[-88,35],[-168,139],[183,66],[-289,2],[-71,29],[49,54],[-132,45],[-35,66],[131,60],[-58,68],[97,71],[230,29],[153,1]],[[9693,38351],[55,0],[124,-72],[160,18],[51,-56],[177,-31],[-109,-57],[232,-124],[106,24],[135,-26],[165,34],[77,38],[78,-15],[123,18],[207,-46],[75,-42],[17,-130],[-80,-30],[-7,-39],[-173,-24],[-195,25],[-101,-18],[-191,7],[-167,-16],[-94,6],[-41,51],[-120,-39],[-106,34],[-82,-10],[-52,33],[-13,75],[-39,52],[40,71],[-13,44],[-102,116],[-65,-18],[-182,-2],[-58,62],[-97,38],[-13,62],[127,21],[151,-34]],[[9478,38722],[148,-42],[-50,-28],[28,-53],[-215,-43],[-143,143],[4,86],[184,-27],[44,-36]],[[7865,38624],[78,-33],[-48,-100],[-199,-41],[-137,42],[-5,86],[311,46]],[[7289,38454],[-83,-60],[49,-30],[-41,-83],[-204,-38],[-152,-51],[-62,-84],[-56,-7],[-69,45],[-101,3],[-98,44],[47,40],[85,9],[285,199],[123,17],[66,-14],[120,51],[91,-41]],[[8110,38154],[131,-41],[84,42],[89,-30],[14,-43],[-43,-127],[-146,-61],[-136,4],[-54,28],[-76,-42],[-76,-10],[-86,-46],[-179,-51],[-112,3],[-86,40],[16,37],[146,48],[-182,26],[-145,-27],[-84,47],[-127,23],[50,46],[21,62],[66,15],[-25,61],[116,83],[157,4],[42,-57],[127,-1],[134,-90],[52,-60],[223,-10],[16,45],[-97,37],[35,65],[-88,63],[105,78],[80,-39],[57,-81],[-19,-41]],[[9821,37718],[107,-43],[128,8],[83,-32],[-143,-193],[-56,-63],[-139,11],[-48,-31],[24,-56],[-54,-91],[-131,0],[-7,109],[-33,64],[-10,210],[67,76],[63,20],[149,11]],[[9338,36680],[129,-68],[69,-140],[-81,-61],[-121,17],[-47,29],[-154,44],[-35,39],[125,69],[-11,49],[32,48],[54,24],[40,-50]],[[7406,37358],[115,23],[33,66],[52,1],[168,-62],[55,-42],[65,30],[-52,79],[74,-5],[101,-59],[47,-52],[56,-172],[32,-26],[70,59],[-44,58],[-62,220],[59,51],[69,-31],[71,1],[121,-94],[20,-86],[100,-225],[-25,-76],[70,-78],[96,-58],[40,2],[86,-54],[95,-32],[24,-99],[-80,-4],[-61,34],[-49,-68],[73,-33],[13,-88],[-117,-46],[-65,-3],[-100,27],[-85,-2],[9,37],[-126,19],[-59,64],[-88,-101],[-105,-15],[-66,-41],[-114,-29],[-153,-21],[-206,-11],[-55,79],[-9,82],[-71,17],[-94,-1],[-153,37],[-67,87],[-4,68],[288,49],[221,-20],[96,11],[-37,40],[-186,56],[-252,-24],[-181,10],[-76,59],[50,62],[180,47],[28,25],[-250,-8],[51,52],[-128,6],[-8,68],[82,64],[-18,62],[93,69],[111,51],[219,73],[66,-69],[-53,-110]],[[6816,37707],[200,35],[87,-21],[191,-127],[7,-42],[-397,-173],[-51,-74],[-87,-34],[-30,-134],[-20,-29],[-293,-88],[-90,128],[-146,68],[-53,37],[88,104],[17,112],[118,159],[-53,43],[-49,90],[382,42],[86,-41],[113,-27],[-20,-28]],[[11431,33607],[-27,-88],[-27,12],[-40,-30],[-27,29],[46,72],[23,7],[14,75],[35,-42],[3,-35]],[[13026,34778],[-64,46],[30,28],[52,-21],[-18,-53]],[[11375,36714],[-14,-37],[-69,-21],[-56,49],[139,9]],[[11400,36497],[69,121],[42,-16],[-71,-89],[-40,-16]],[[10564,36242],[-30,17],[-11,67],[61,19],[4,-59],[-24,-44]],[[8796,38144],[-40,-54],[-83,-3],[-117,96],[142,29],[58,-6],[40,-62]],[[8593,38275],[110,-31],[-138,-28],[28,59]],[[8582,37919],[-87,26],[62,65],[79,-61],[-54,-30]],[[14086,32520],[-40,-79],[-7,-66],[-93,-194],[-17,-68],[10,-26],[58,101],[54,-30],[-14,-87],[24,-34],[34,13],[40,-43],[50,43],[71,-11],[38,-28],[5,-41],[-31,-76],[6,-66],[88,11],[-13,-39],[-45,-39],[-19,-72],[11,-59],[36,77],[37,7],[-25,-86],[59,-24],[-29,-103],[-18,-97],[-62,0],[4,64],[-67,-18],[35,117],[-14,84],[-23,24],[-27,-91],[-48,-18],[-52,-110],[-53,-9],[-6,47],[39,21],[34,66],[-20,49],[-23,-45],[-33,14],[2,60],[-42,-27],[-83,-22],[-110,23],[-46,0],[-89,-24],[-28,68],[11,29],[63,76],[25,43],[-25,20],[35,112],[34,-2],[-13,71],[23,36],[81,266],[19,20],[43,112],[82,67],[38,-22],[26,15]],[[22338,36670],[57,-21],[-6,-71],[-111,-14],[3,52],[57,54]],[[22536,36778],[38,-23],[-91,-89],[-49,43],[102,69]],[[22141,36432],[17,-40],[-111,-19],[55,71],[39,-12]],[[22078,36519],[2,-75],[-96,22],[94,53]],[[15284,39539],[-213,60],[-59,71],[99,8],[186,-66],[-13,-73]],[[18323,38005],[46,-96],[-122,2],[-21,73],[97,21]],[[16971,39879],[453,-63],[167,-95],[220,-20],[113,-45],[-107,-48],[-157,-23],[-647,-29],[4,-50],[193,25],[304,-11],[171,-60],[60,62],[200,13],[22,-86],[-66,-79],[110,-9],[84,58],[315,-30],[189,92],[209,-10],[171,-32],[76,-52],[-218,-89],[-117,-14],[-2,-49],[-255,-43],[93,-41],[-105,-46],[-136,-6],[-154,20],[-82,-56],[2,-48],[96,-28],[-14,-77],[38,-45],[-213,-137],[-89,-176],[100,28],[155,-44],[-135,-23],[51,-58],[106,-33],[108,-2],[-19,-100],[-176,34],[-99,-9],[-106,-72],[35,-65],[166,-17],[67,-106],[9,-121],[-119,19],[-57,-55],[113,-23],[86,-113],[-16,-45],[-199,-40],[-135,44],[-2,-55],[223,-50],[-16,-81],[-118,-14],[-72,-37],[-194,29],[119,-72],[92,-38],[-4,-119],[-25,-65],[-201,87],[-67,-14],[203,-107],[97,-65],[29,-48],[9,-143],[17,-76],[-47,-19],[-157,2],[-52,24],[-67,136],[-144,88],[-10,-80],[-110,-53],[-110,10],[-97,-119],[60,-17],[143,17],[67,-53],[71,24],[267,-54],[9,-54],[-70,-58],[-208,-134],[-17,-36],[-73,-43],[-329,-93],[-70,3],[-125,-57],[-50,15],[-81,59],[-26,-39],[10,-67],[-100,-68],[-106,-207],[-61,-67],[-112,-68],[-85,-71],[-63,-11],[-93,-43],[-28,26],[-71,-19],[-150,-16],[43,-49],[-122,-98],[54,-100],[-68,-59],[19,-21],[7,-93],[-56,-49],[-7,-55],[-59,-69],[-59,-93],[-19,-96],[24,-68],[-35,-73],[-34,-177],[-50,-100],[4,-62],[-22,-31],[-84,-1],[-85,29],[-64,35],[-24,73],[-52,22],[-10,37],[14,74],[-97,-65],[-104,3],[-25,-24],[-42,64],[-12,51],[-49,12],[-41,72],[-11,69],[-27,25],[8,54],[-87,55],[-2,82],[-46,51],[-82,135],[2,61],[-34,86],[-38,31],[-19,170],[-57,100],[-49,1],[10,88],[-57,41],[-3,59],[68,95],[-74,37],[-22,34],[10,65],[43,36],[-18,60],[48,85],[176,-35],[44,63],[-112,-25],[-66,21],[-68,-3],[39,72],[83,21],[59,-35],[60,44],[20,108],[65,167],[-199,28],[-87,56],[-112,27],[-55,54],[38,38],[199,-31],[100,-57],[45,123],[-28,33],[-78,1],[-119,44],[41,43],[-82,32],[-55,-41],[-82,-22],[-102,43],[31,128],[-39,26],[7,53],[92,43],[5,53],[-105,31],[42,59],[-91,68],[20,77],[-42,62],[-46,13],[-11,104],[-204,164],[10,69],[-302,108],[-255,42],[-107,-8],[-121,-40],[-56,25],[-99,-62],[-149,22],[-139,61],[14,82],[-174,42],[-8,70],[213,5],[95,42],[190,3],[-26,66],[-146,-37],[-73,19],[-82,-31],[-260,99],[-132,68],[47,71],[385,81],[171,60],[171,4],[62,47],[71,152],[-232,17],[-17,78],[153,62],[128,73],[121,3],[51,44],[176,-17],[38,61],[-12,80],[87,40],[131,-8],[63,48],[452,59],[220,-7],[213,-98],[136,18],[-115,61],[-32,51],[241,-12],[404,-133],[86,5],[33,119],[-38,40],[-215,94],[371,70],[214,-36],[107,54],[182,0],[104,35],[548,26],[343,-8]],[[13909,6118],[-140,-90],[-65,-97],[-7,-47],[-76,-32],[-44,23],[-79,-69],[-46,-58],[-76,-12],[-46,-68],[17,-59],[-35,-58],[52,-64],[65,32],[73,11],[-32,-60],[-110,-24],[-91,14],[-3,-102],[-62,96],[-56,-4],[-37,-86],[4,-63],[-73,28],[-21,-50],[6,-65],[-70,-5],[3,-76],[-26,-108],[73,-64],[11,-61],[176,-32],[-13,-53],[114,-124],[5,-56],[51,-68],[-48,-56],[43,-30],[5,-93],[73,8],[33,-69],[-22,-91],[20,-59],[-102,-19],[71,-48],[67,9],[-1,-81],[39,-134],[48,-5],[-164,-101],[60,-33],[-18,-117],[-60,-17],[63,-65],[-137,-6],[37,-48],[-92,-5],[-138,-62],[79,-34],[-15,-68],[-340,-130],[-130,-22],[-202,-55],[-89,-70],[-193,-22],[-236,12],[-157,27],[-217,-11],[3,-47],[104,-101],[97,-46],[287,-27],[-45,-70],[-177,-65],[-331,55],[-324,43],[-95,-31],[406,-107],[48,-19],[-50,-71],[-274,-18],[-334,107],[-58,-24],[57,-72],[180,-75],[87,-92],[84,54],[413,-14],[53,-69],[-59,-63],[96,-111],[53,-96],[250,-26],[227,68],[88,-54],[621,-153],[256,-8],[-6,-28],[-244,6],[63,-50],[-201,-8],[1,-64],[112,-43],[257,24],[-10,-60],[93,-104],[37,-83],[234,-27],[245,142],[175,83],[206,64],[377,52],[259,19],[262,-62],[-54,-68],[167,4],[155,39],[253,223],[301,96],[177,-42],[141,40],[72,55],[347,44],[331,66],[576,25],[94,47],[-209,22],[-531,38],[-41,100],[-124,4],[-179,-19],[-285,68],[-82,36],[6,70],[124,145],[61,25],[53,70],[166,65],[62,-4],[208,115],[146,66],[204,49],[121,46],[203,44],[204,24],[116,-5],[107,44],[113,-11],[110,35],[-26,34],[164,208],[98,3],[86,-20],[124,108],[-133,-1],[-21,43],[-66,26],[10,50],[95,76],[66,11],[78,-18],[36,29],[-34,61],[124,-12],[52,36],[141,50],[61,117],[-134,66],[-23,52],[77,13],[95,-51],[35,59],[71,62],[74,-34],[67,-114],[105,29],[11,97],[-29,42],[55,35],[62,-13],[102,30],[39,-33],[-55,-68],[15,-43],[169,3],[84,-11],[101,13],[91,-26],[126,19],[40,-83],[79,72],[77,45],[201,69],[101,13],[182,43],[259,37],[31,32],[95,-27],[70,67],[122,-75],[112,-50],[56,-9],[86,91],[53,36],[101,-30],[76,10],[148,-5],[88,27],[14,-44],[160,-33],[35,57],[87,0],[-20,-88],[27,-52],[98,-4],[107,17],[30,80],[41,55],[75,-50],[-15,-39],[84,-38],[49,20],[45,73],[49,-10],[59,-102],[137,-33],[140,29],[55,32],[137,43],[124,65],[156,17],[140,52],[42,86],[-48,127],[8,47],[55,41],[121,-3],[-43,-91],[88,1],[72,-125],[168,-3],[45,-36],[51,21],[66,-25],[40,-53],[41,12],[95,133],[18,105],[144,90],[144,52],[47,52],[212,51],[33,31],[159,34],[-9,50],[64,29],[59,-31],[-32,-37],[42,-36],[53,16],[78,-31],[7,133],[-28,39],[88,24],[86,-54],[68,6],[-9,76],[-25,15],[3,71],[176,98],[202,37],[144,-14],[86,-37],[74,-68],[56,-11],[37,-44],[-82,-38],[-34,-111],[83,47],[83,10],[126,-49],[58,-59],[139,23],[94,-36],[156,-24],[114,32],[99,-26],[217,-34],[77,0],[131,-28],[123,34],[64,-166],[-51,-63],[11,-113],[-82,-33],[32,-65],[-12,-48],[-102,6],[-103,-100],[76,-34],[148,-2],[-43,-139],[-94,-82],[-69,-140],[-61,-217],[-36,-60],[85,-21],[73,49],[0,77],[51,54],[128,32],[120,123],[75,50],[21,108],[50,103],[97,71],[7,69],[57,57],[70,24],[71,-31],[96,-3],[53,73],[147,85],[71,31],[76,102],[27,80],[63,36],[309,99],[145,20],[22,37],[112,75],[139,-11],[46,25],[71,4],[139,57],[195,-8],[96,20],[53,31],[119,21],[75,-26],[71,9],[68,-24],[99,43],[79,-30],[100,7],[77,23],[70,-24],[61,32],[68,-47],[58,5],[171,67],[67,93],[79,0],[56,20],[145,-28],[275,-97],[89,-14],[111,-39],[25,-29],[123,-33],[131,89],[-4,48],[36,55],[62,26],[186,37],[69,-30],[71,-99],[110,-46],[37,-48],[-41,-59],[-145,-42],[4,-53],[304,91],[47,-16],[114,7],[129,-26],[-16,-40],[266,65],[79,53],[202,57],[71,-23],[56,16],[34,50],[53,13],[91,-30],[23,-60],[76,-67],[100,-16],[92,18],[57,131],[42,41],[95,23],[99,-10],[118,13],[64,22],[84,-45],[12,-45],[77,33],[43,47],[95,-38],[41,-34],[38,23],[98,-15],[44,-32],[124,-7],[75,-32],[113,-9],[43,-18],[78,6],[40,-33],[64,-12],[63,29],[102,-33],[31,-28],[-77,-150],[52,0],[82,39],[96,0],[96,-78],[-3,-71],[187,-45],[296,19],[8,-87],[48,16],[111,-8],[60,-29],[71,29],[9,87],[47,-15],[89,-98],[60,-42],[110,-36],[63,1],[49,-30],[77,23],[184,-72],[38,-49],[10,-55],[50,-24],[43,-55],[44,-117],[98,-33],[-28,79],[20,67],[51,7],[77,-73],[73,-3],[46,29],[164,-31],[87,-5],[107,-33],[68,-86],[71,-21],[154,-82],[79,-54],[-72,-19],[-30,-104],[-128,-34],[114,-41],[-32,-74],[-271,-22],[52,-44],[-70,-47],[-94,2],[-6,-46],[-151,-68],[60,-129],[-128,16],[-114,-18],[-64,-51],[-7,-79],[-97,-12],[128,-144],[-42,-71],[32,-16],[4,-120],[-35,-48],[72,-18],[59,-81],[1,-45],[92,-106],[201,-94],[-10,-27],[-168,-8],[-303,-98],[-86,12],[-80,-51],[-31,-81],[39,-105],[-202,-65],[-22,-30],[222,0],[42,-222],[216,-109],[97,-74],[-284,-51],[238,-13],[167,25],[281,-137],[155,-26],[-56,-68],[352,-26],[79,-39],[846,-129],[68,-33],[0,-1298],[-40717,0],[0,96],[0,66],[0,162],[0,163],[0,162],[0,162],[0,162],[0,162],[0,163],[182,3],[756,-47],[542,-64],[450,-19],[657,-65],[16,86],[-745,65],[-126,86],[68,48],[-98,32],[-294,-1],[-183,83],[-554,126],[254,12],[385,-87],[249,6],[351,-30],[104,-43],[311,-6],[228,46],[-28,49],[251,35],[277,122],[-193,116],[86,56],[-108,43],[-240,44],[57,36],[508,31],[442,29],[-277,126],[40,47],[215,18],[16,68],[-261,53],[-181,70],[-18,32],[-144,-4],[-191,35],[-143,76],[199,54],[-71,48],[-218,-2],[-115,57],[-24,40],[32,143],[175,-13],[52,24],[158,-4],[102,-23],[51,-44],[197,-2],[226,-86],[156,57],[-46,48],[206,19],[60,-49],[90,4],[-42,89],[25,73],[-180,72],[-143,-13],[-119,30],[151,59],[275,-73],[28,23],[-84,52],[45,48],[100,2],[173,73],[141,16],[93,-45],[205,109],[235,32],[120,-15],[37,70],[98,33],[213,-37],[161,20],[256,-29],[218,40],[284,1],[165,-14],[266,6],[211,22],[151,62],[144,-21],[222,4],[35,108],[51,14],[94,-37],[-44,-128],[-4,-78],[253,44],[-12,120],[71,19],[98,-41],[1,-79],[-131,-97],[40,-14],[208,2],[165,-32],[135,-5],[194,53],[361,-3],[76,-68],[280,54],[-193,84],[17,71],[-128,6],[-53,107],[-119,32],[41,61],[161,-32],[136,8],[66,35],[-93,38],[-132,20],[-168,-1],[-46,31],[-33,84],[100,20],[9,-58],[179,1],[318,-12],[210,-63],[104,20],[112,-21],[61,22],[166,9],[149,-33],[56,18],[15,61],[108,20],[24,45],[68,-2],[37,-36],[-25,-85],[115,9],[41,-38],[120,34],[113,-68],[316,-81],[99,27],[5,97],[93,86],[104,-34],[126,-107],[176,22],[-32,-82],[96,3],[113,38],[146,-14],[121,58],[295,40],[170,39],[130,62],[41,51],[54,120],[-71,129],[-5,108],[-44,146],[-58,93],[-8,86],[-27,49],[65,19],[87,-15],[44,58],[-47,69],[68,142],[-40,50],[50,114],[-100,21],[11,74],[61,30],[39,-52],[27,138],[68,-3],[17,112],[109,62],[66,72],[4,86],[31,33],[72,21],[62,43],[107,35],[62,66],[25,56],[66,34],[87,24],[72,62],[74,23],[94,51],[41,-33]],[[40588,36952],[-26,51],[30,42],[125,71],[0,-126],[-129,-38]],[[22392,34310],[60,19],[-25,-69],[-94,-58],[-91,-97],[13,-51],[-24,-92],[11,-72],[-20,-54],[6,-45],[-20,-82],[-40,-112],[-76,-17],[-62,3],[-64,-76],[15,-70],[-19,-30],[-42,7],[-55,-19],[-49,15],[10,65],[-43,102],[23,29],[10,81],[-53,67],[-77,186],[-4,59],[-29,33],[-20,82],[-8,139],[-41,24],[-22,47],[-52,-87],[-64,-15],[-124,-148],[-40,-36],[-79,-28],[-52,1],[-39,54],[-31,2],[-77,58],[-21,47],[4,57],[34,-3],[13,52],[-14,38],[-55,-43],[-26,14],[19,95],[53,5],[-4,46],[-69,-7],[14,169],[-22,32],[-14,231],[68,138],[67,31],[26,38],[105,-2],[-31,81],[64,33],[60,-2],[33,102],[51,25],[54,-17],[42,12],[44,-44],[77,7],[-4,38],[-91,-24],[-24,79],[97,127],[74,60],[14,52],[92,63],[13,89],[79,104],[15,119],[58,59],[128,162],[5,71],[59,60],[84,30],[59,61],[8,37],[83,39],[36,69],[45,13],[22,74],[78,19],[24,43],[98,31],[124,2],[27,65],[145,46],[41,-62],[155,138],[27,68],[126,-34],[-56,-69],[22,-48],[118,131],[82,8],[41,33],[170,-53],[126,-73],[43,-4],[39,-58],[-152,-68],[22,-57],[78,15],[110,-24],[60,61],[108,-47],[-63,-62],[121,-12],[26,-27],[246,-27],[212,-115],[79,-78],[45,-7],[84,-58],[53,-13],[39,-42],[66,-27],[19,-103],[25,-13],[-19,-88],[-75,-87],[-140,-72],[-101,-16],[-124,37],[-202,39],[-78,50],[-152,35],[-71,59],[-44,-24],[125,-137],[94,-44],[44,-61],[-19,-81],[-24,-27],[45,-94],[7,-98],[19,-27],[69,-14],[57,-44],[24,-43],[122,-43],[71,49],[-13,67],[-87,21],[-74,100],[40,75],[73,-15],[54,-52],[198,-70],[59,97],[-54,81],[2,57],[58,35],[41,49],[89,37],[83,91],[44,-23],[71,-1],[98,-42],[44,101],[-22,98],[-58,35],[51,171],[-21,127],[117,7],[92,-23],[90,-117],[0,-28],[-131,-21],[-67,-64],[0,-29],[71,-38],[48,-77],[79,-6],[91,24],[25,26],[18,124],[114,34],[222,142],[193,56],[15,32],[127,61],[19,-35],[-24,-49],[23,-93],[68,10],[96,68],[71,19],[123,-22],[118,77],[84,26],[51,-61],[-30,-67],[90,-5],[0,65],[67,5],[50,59],[-30,29],[-56,110],[96,60],[85,-20],[180,-20],[94,-33],[95,-60],[156,-74],[67,-18],[82,-82],[88,-38],[36,50],[35,88],[-67,4],[-48,62],[-13,56],[-106,49],[-16,74],[25,47],[16,119],[-70,37],[28,88],[71,30],[82,62],[38,76],[49,151],[69,52],[209,4],[149,-51],[-21,-141],[-39,-88],[-42,-34],[19,-55],[57,-36],[17,-76],[-23,-149],[10,-158],[-10,-93],[32,-61],[83,-55],[-46,-91],[1,-72],[-63,-64],[-84,-133],[-55,-11],[21,-46],[-69,-49],[-63,13],[-57,55],[-98,-13],[21,-48],[129,-54],[91,5],[88,-27],[45,20],[11,52],[105,57],[51,43],[32,96],[78,81],[-3,71],[-40,79],[22,76],[61,26],[151,27],[88,-117],[-7,-159],[93,110],[-39,148],[-187,77],[-65,0],[-69,-34],[-117,25],[13,51],[-38,66],[31,108],[54,109],[-91,141],[-47,41],[66,92],[149,69],[12,63],[-35,90],[66,-3],[42,-115],[-56,-111],[29,-37],[-19,-71],[88,-29],[167,-11],[-12,33],[-132,65],[-32,72],[95,28],[77,-43],[77,25],[-86,55],[114,47],[106,-4],[152,-67],[85,-79],[166,1],[-7,-75],[-68,-37],[2,-104],[60,34],[18,-86],[-39,-84],[92,30],[27,46],[-66,128],[48,113],[-42,64],[-63,12],[-52,71],[-154,58],[-19,51],[21,55],[-38,31],[9,111],[139,21],[196,2],[344,52],[70,-4],[-72,84],[-65,17],[80,47],[-6,37],[87,60],[28,43],[186,78],[99,28],[146,13],[234,47],[-42,24],[158,45],[152,-9],[109,-42],[201,72],[17,54],[171,0],[88,52],[-12,65],[242,147],[100,23],[212,-55],[19,-23],[-130,-49],[91,-40],[139,8],[55,-25],[-89,-81],[111,-15],[46,48],[349,1],[170,-78],[20,-58],[71,-4],[51,-59],[-29,-130],[-78,-63],[-120,-64],[-229,-96],[-4,-35],[-83,-31],[-99,-78],[-105,-17],[-90,-110],[144,7],[238,69],[-11,58],[57,67],[90,-18],[124,-53],[73,9],[143,-37],[145,27],[222,-24],[130,-2],[-3,-79],[150,-61],[141,-10],[100,6],[91,-17],[53,15],[53,55],[-36,78],[122,51],[139,-54],[160,6],[80,-15],[130,-67],[41,-111],[-15,-51],[33,-39],[-14,-51],[-105,-1],[159,-222],[158,-87],[46,36],[75,171],[38,58],[97,-93],[77,-26],[201,58],[162,-62],[38,-35],[53,82],[100,-2],[23,-26],[88,11],[-33,48],[18,103],[-80,42],[52,38],[96,-1],[71,22],[-31,70],[142,-40],[165,-5],[222,-36],[126,-55],[-95,-83],[-71,-3],[7,-53],[69,12],[64,40],[86,85],[144,4],[125,-33],[52,-40],[-10,-34],[-108,-30],[190,-44],[111,-54],[105,-104],[145,10],[238,50],[242,-13],[149,-57],[42,-32],[32,-79],[-32,-101],[12,-19],[122,-41],[14,-117],[60,-44],[-3,109],[95,62],[201,16],[40,-26],[142,-5],[119,-20],[92,56],[59,-38],[17,-70],[114,-44],[34,-68],[105,9],[52,51],[-47,124],[-11,120],[234,-32],[82,-33],[76,11],[244,-3],[92,-49],[218,-60],[130,-92],[0,-844],[0,-56],[-62,-56],[-105,-51],[-88,26],[-59,34],[-110,-5],[54,-50],[84,26],[-4,-67],[59,-51],[48,9],[30,-68],[7,-102],[70,-76],[-8,-42],[35,-74],[-50,-84],[-125,52],[-76,9],[-73,-19],[-42,-36],[-120,-56],[-85,-64],[-101,-25],[-56,-71],[-31,14],[-52,-93],[-164,-124],[-38,-20],[-29,-108],[-42,23],[-41,85],[-44,37],[-124,-6],[-103,-38],[-29,-23],[-72,-100],[-22,25],[24,116],[-143,-89],[-7,-54],[-73,45],[-68,-4],[-41,-47],[-16,-123],[-30,-65],[-63,-71],[-35,-60],[-20,-85],[9,-37],[48,-35],[28,42],[49,-25],[7,-35],[-49,-76],[3,-125],[51,-28],[9,-105],[-32,-43],[-59,49],[-50,-39],[-35,-100],[-5,-68],[42,-139],[-43,-51],[-108,2],[-79,-81],[-26,-93],[6,-92],[-107,-77],[-38,-40],[-15,-57],[-2,-73],[-41,-114],[-64,-76],[-37,-61],[-52,-56],[-28,116],[-1,101],[-14,137],[-16,26],[-14,88],[-22,212],[-32,215],[-8,111],[18,166],[31,143],[84,105],[28,72],[-21,65],[75,12],[24,51],[62,1],[93,91],[63,88],[29,74],[79,97],[35,18],[67,92],[34,29],[29,61],[109,84],[84,30],[-11,45],[44,53],[-23,27],[25,58],[21,133],[44,41],[-39,52],[-104,-33],[-37,-175],[9,-51],[-78,22],[-154,-161],[-30,-48],[-62,18],[18,44],[-49,24],[-7,37],[49,116],[-15,25],[-72,-42],[-53,46],[-114,-38],[-68,10],[-90,-73],[-5,-48],[-65,-64],[-23,-50],[-100,-92],[-77,-125],[-16,-70],[49,3],[66,-41],[0,-39],[-52,-12],[-28,18],[-50,-33],[-36,34],[-38,-2],[-54,-63],[-55,20],[-35,-27],[-86,-8],[-23,48],[98,18],[-57,83],[-116,16],[-93,41],[-65,-32],[16,-33],[-150,-23],[-41,-27],[-121,37],[-20,-48],[-38,-6],[-43,50],[-121,-9],[-70,8],[-76,-9],[-69,-30],[-63,-55],[-31,-59],[-86,-75],[-35,-47],[-27,-92],[-50,-29],[-63,-82],[-30,-16],[-58,-68],[-110,-189],[-102,-103],[-37,-27],[-33,-51],[-71,-55],[-32,-39],[0,-49],[67,-34],[107,9],[-13,-159],[53,-25],[25,106],[31,-96],[-36,-75],[109,12],[41,35],[4,133],[71,-30],[44,20],[53,-52],[19,-55],[75,-73],[44,-71],[-35,-66],[19,-13],[-10,-109],[36,-43],[-14,-60],[-49,-69],[-27,-89],[-24,-158],[6,-95],[-8,-50],[6,-73],[-23,-124],[-16,-122],[-19,-46],[-71,-100],[-42,-116],[-47,-75],[-28,-118],[-74,-167],[-100,-149],[-74,-156],[-31,-27],[-44,-124],[-40,-71],[-49,-54],[-112,-105],[-62,-31],[-51,41],[-46,2],[4,81],[-35,-27],[-34,19],[-64,-132],[-37,-24],[-12,-45],[-30,-6],[-44,-59],[-36,-76],[-8,-50],[9,-44],[-6,-103],[-42,-30],[-26,-54],[-46,-40],[-45,-64],[-38,-10],[-45,-49],[-2,-74],[-18,-58],[44,-28],[67,-106],[28,-103],[49,-114],[32,-93],[15,-122],[-8,-159],[18,-24],[-16,-103],[-23,-73],[-26,-18],[-54,0],[-7,-53],[-83,19],[-36,-49],[-8,-64],[-33,33],[-59,-59],[-13,134],[-12,51],[23,80],[0,33],[26,28],[-3,56],[-17,35],[0,72],[-37,90],[41,35],[33,-11],[-25,109],[-9,80],[-23,19],[-48,-1],[-20,25],[-36,-58],[-29,70],[-57,21],[43,99],[27,25],[-17,47],[29,96],[-6,55],[-124,91],[-18,-19],[-85,-18],[-36,-21],[-79,-70],[-24,-50],[-54,-65],[-37,-14],[-24,26],[65,44],[12,67],[-58,-4],[-1,37],[28,21],[0,48],[33,24],[44,94],[9,42],[-47,69],[-77,13],[-46,-71],[-33,-83],[-100,-75],[-23,-35],[-19,-79],[-35,-55],[-33,2],[-37,-25],[-29,37],[-29,-22],[-25,-117],[23,-72],[28,-29],[105,-32],[17,-79],[-16,-85],[18,-30],[39,-17],[35,7],[14,45],[55,80],[43,31],[53,-54],[47,-32],[79,-13],[30,5],[-9,-110],[-20,-27],[-47,30],[-106,-83],[-35,-59],[0,-40],[-33,-26],[-34,15],[-4,-60],[-71,-128],[-24,-67],[15,-60],[27,-39],[70,-58],[13,-36],[21,-122],[41,-143],[-2,-82],[62,-66],[1,-38],[48,-71],[-17,-43],[-37,34],[-34,-32],[69,-92],[24,-93],[-39,-17],[-60,-65],[-20,-47],[-42,8],[20,-59],[31,7],[29,32],[31,-17],[28,-57],[30,-20],[2,-92],[-5,-82],[-23,28],[-10,-112],[-17,-30],[14,-63],[-28,-31],[-29,14],[-16,-58],[-42,-105],[-5,-61],[-30,-50],[-21,-98],[-35,-29],[11,-49],[-47,-64],[21,-31],[-12,-75],[-32,-29],[-3,-60],[-23,4],[-14,-68],[-39,-80],[-54,11],[-9,-33],[6,-54],[-48,-94],[-19,0],[-68,-89],[-37,-62],[-7,-54],[-70,-33],[-24,12],[-16,-31],[-35,22],[-38,-46],[-23,32],[-18,-52],[-26,3],[1,-57],[-42,31],[-32,99],[-32,12],[25,-73],[-1,-85],[-32,-31],[-29,7],[-48,-79],[-31,-8],[-41,25],[-39,-67],[-57,-17],[-37,-22],[-65,-81],[-2,-41],[38,-78],[-7,-37],[-37,-21],[-52,150],[2,49],[28,80],[-44,13],[-44,-25],[-14,46],[-80,15],[-9,-31],[-49,-14],[-45,-49],[-6,-53],[-77,-21],[9,-52],[-23,-48],[-4,-55],[-40,-68],[-20,-12],[-20,-81],[-2,-67],[-19,-76],[30,-106],[39,-65],[30,-63],[-2,-52],[29,-81],[43,-72],[7,-36],[61,-94],[13,-38],[22,2],[73,-195],[16,-24],[30,-152],[1,-38],[24,-160],[-4,-132],[18,-74],[-24,-72],[-2,-211],[-18,-31],[-6,-59],[-19,-4],[-62,-88],[-20,-9],[-28,-45],[-66,-69],[-20,23],[-31,-15],[-3,-73],[-21,-59],[-3,-50],[-46,-56],[-75,-70],[-20,-67],[-24,-40],[-25,-10],[-8,50],[3,185],[6,48],[22,20],[-8,38],[-25,31],[-16,-8],[-45,86],[-37,10],[-31,-18],[-15,22],[22,66],[-22,59],[-20,-52],[-23,-2],[-7,69],[4,57],[-41,127],[-25,22],[-57,97],[-35,36],[-32,-16],[-66,22],[12,165],[-35,20],[-48,-8],[-24,-30],[8,-71],[-14,-82],[3,-119],[-41,-163],[-16,-131],[-37,-131],[0,-135],[26,-119],[38,23],[20,-46],[6,-102],[36,-93],[20,-189],[-3,-59],[18,-3],[49,-72],[54,1],[34,-90],[34,-53],[27,-16],[22,-71],[51,-78],[24,-60],[25,-95],[6,-105],[-12,-143],[10,-57],[-1,-135],[11,-37],[31,-44],[46,-197],[8,-56],[-13,-27],[-21,21],[-33,-1],[-25,-28],[-14,50],[-71,71],[-20,43],[-46,47],[-41,75],[-29,25],[-26,47],[1,84],[-66,164],[9,13],[-21,81],[0,64],[-16,90],[-13,125],[-1,91],[-26,105],[-48,100],[3,52],[-21,52],[-20,10],[-6,43],[-29,75],[-29,40],[-24,65],[-18,-38],[-20,56],[0,79],[15,120],[21,126],[-11,193],[24,71],[8,58],[-1,81],[-11,34],[6,112],[-13,215],[-38,131],[-15,-5],[-1,104],[-33,160],[-11,233],[-14,33],[4,119],[-28,-3],[-22,124],[-24,59],[-18,-124],[-21,-50],[-69,-71],[-28,-19],[-10,-44],[-33,-59],[-24,24],[-26,-2],[-6,37],[-34,71],[-33,-64],[-10,27],[7,90],[9,28],[27,214],[-18,145],[-20,93],[-21,66],[-41,30],[33,96],[-44,77],[-28,64],[-41,4],[-24,32],[12,29],[-18,51],[-14,-21],[-34,74],[-29,88],[-5,117],[-36,189],[-24,67],[-37,-42],[-23,-4],[-36,114],[-19,-9],[-3,-70],[18,-72],[-1,-41],[-48,-102],[-8,47],[-32,-8],[-27,-50],[-13,12],[-46,-38],[-90,-7],[-34,46],[-31,-39],[-65,-35],[-27,-61],[-2,-30],[15,-93],[-26,-89],[-28,-33],[-25,-58],[-80,-52],[-13,47],[-24,-32],[-2,-57],[-52,-88],[-35,-100],[-40,-91],[-51,-51],[-52,-106],[-68,-77],[-26,-41],[0,-72],[-12,-51],[-56,-53],[-41,8],[-18,-23],[-12,-70],[-18,-47],[-37,32],[-40,-43],[-22,-88],[-5,-58],[12,-114],[-6,-84],[15,-101],[-20,-38],[31,-56],[-12,-155],[-10,-54],[-32,-107],[-12,-95],[11,-87],[-2,-201],[-50,-3],[-15,-63],[-30,-80],[-9,-53],[7,-43],[-63,-37],[-26,-50],[-15,-116],[-61,-70],[-24,15],[-38,60],[-47,114],[-26,126],[15,20],[-13,88],[-47,198],[-23,133],[-22,79],[-38,81],[-28,117],[-20,117],[-11,136],[-20,87],[-13,103],[-49,133],[-2,74],[-18,41],[-33,110],[-16,94],[-20,267],[-19,109],[-14,134],[-6,142],[-17,122],[25,168],[-8,128],[-14,14],[-17,119],[-2,64],[8,69],[-29,-1],[-10,-56],[-24,-47],[25,-37],[0,-30],[-27,-86],[-70,-66],[-42,-30],[-35,0],[-26,23],[-41,59],[-43,94],[-63,112],[-20,47],[61,48],[71,38],[17,56],[-9,36],[-66,-49],[-49,20],[-67,79],[-25,88],[-21,5],[-22,58],[-42,-8],[-41,83],[-15,134],[-53,24],[-1,84],[-30,80],[-85,-51],[-31,5],[-71,-16],[-13,-28],[-58,34],[-48,12],[-26,-40],[-104,10],[-29,-27],[-103,-8],[-84,42],[-73,23],[-46,4],[-18,18],[-46,-14],[-28,31],[-114,23],[-52,32],[-26,133],[-14,123],[-19,42],[-51,25],[-80,-51],[-26,-47],[-15,3],[-45,-53],[-28,-10],[-45,43],[-61,6],[-29,50],[-53,46],[-33,42],[-24,67],[-51,48],[-41,4],[-44,66],[-21,88],[-3,49],[-22,33],[1,44],[-23,20],[-3,62],[-54,116],[-21,66],[-49,-42],[-72,85],[0,-64],[-52,-41],[-54,16],[-1,-90],[-28,-46],[37,-14],[15,-87],[18,-50],[32,-144],[17,-40],[15,-76],[30,-44],[26,-59],[66,-81],[5,-69],[21,-51],[-21,-45],[28,-125],[20,-33],[12,-77],[32,-46],[-9,128],[27,123],[29,40],[32,-58],[-6,-87],[13,-86],[-2,-42],[-31,-144],[52,-19],[15,-62],[39,2],[45,37],[141,-18],[57,46],[40,122],[39,53],[48,105],[48,68],[25,95],[30,28],[-6,-120],[-1,-181],[31,-126],[30,-73],[35,-53],[45,-27],[57,-20],[23,-21],[29,4],[22,-29],[29,-89],[57,-127],[33,-16],[-3,-66],[-48,-166],[-54,-88],[-48,-163],[-30,4],[-12,34],[-27,-75],[-17,-146],[11,-136],[-72,-26],[-39,-34],[-20,-39],[-11,-96],[-32,-50],[-90,-25],[-25,-59],[5,-47],[-27,-78],[-33,-17],[-23,15],[-56,-6],[-52,-56],[-35,-9],[-81,-58],[-29,-41],[-17,-77],[4,-69],[-69,-73],[-67,-45],[-55,-24],[-40,-43],[-30,-5],[-63,-44],[-34,-42],[-43,-93],[-44,-12],[-33,11],[-40,-43],[-26,-45],[-70,-45],[-66,-10],[-62,-20],[-29,-62],[-26,-16],[-15,-42],[-47,1],[-29,-34],[-49,-12],[-51,51],[-28,98],[6,86],[-6,50],[-16,35],[-8,125],[-34,253],[14,86],[-6,97],[-10,65],[-36,89],[-10,72],[-61,103],[-36,130],[-23,52],[-13,93],[-41,155],[-76,117],[-22,6],[-31,52],[-38,105],[-21,77],[6,48],[-18,83],[12,118],[-17,112],[-54,191],[-20,46],[-42,63],[-42,24],[-41,121],[9,33],[-13,76],[-26,80],[-18,25],[-9,67],[-18,15],[-31,112],[-75,193],[-18,69],[-27,69],[-64,26],[18,83],[23,224],[-28,-49],[-13,-118],[-20,-92],[-5,-78],[-20,-58],[-52,65],[-58,120],[-20,116],[-40,103],[-17,104],[-24,-79],[24,-56],[7,-91],[30,-98],[35,-82],[39,-71],[0,-67],[39,-128],[7,-92],[42,-144],[27,-76],[71,-280],[49,-94],[-14,-66],[2,-77],[20,-114],[24,-48],[36,-25],[21,-54],[47,-68],[11,-117],[32,-93],[-7,-18],[9,-126],[-4,-100],[6,-124],[13,-113],[13,-62],[51,-61],[31,-70],[42,-41],[39,-148],[26,-161],[17,-185],[24,-90],[40,-14],[20,-31],[19,-74],[48,-30],[19,-32],[43,-28],[54,-147],[67,-90],[17,-87],[14,3],[31,-83],[32,-9],[37,-145],[-15,-51],[-54,-53],[-18,-41],[59,1],[32,-50],[46,-129],[61,-82],[63,2],[44,49],[54,42],[85,-20],[46,41],[49,57],[35,-14],[35,6],[74,42],[26,-15],[54,20],[29,25],[53,18],[40,45],[20,50],[17,9],[53,-35],[-20,-114],[7,-156],[-9,-51],[-19,-42],[-8,-189],[-45,-134],[-37,-148],[-28,-55],[-10,-69],[-22,-84],[-25,-70],[-29,-134],[-5,-54],[-45,-156],[-47,-124],[-29,-105],[-124,-278],[-94,-187],[-25,-38],[-103,-114],[-66,-97],[-98,-176],[-103,-217],[-65,-150],[-28,-105],[-40,-100],[-44,-20],[-11,-73],[-29,-62],[-27,-4],[-20,-30],[-12,-129],[-29,-75],[-5,-48],[-37,-159],[-30,-43],[-28,-196],[-18,-83],[6,-104],[68,-126],[5,-56],[-26,-91],[17,-113],[-15,-99],[17,-115],[21,-57],[-2,-50],[18,-116],[23,-57],[46,-44],[25,-71],[-14,-24],[6,-69],[-16,-61],[20,-300],[1,-169],[-5,-21],[11,-213],[20,-17],[2,-75],[-25,-74],[-6,-83],[-51,-116],[-30,-105],[-68,-82],[-18,-41],[-51,-24],[-55,-38],[-102,-114],[-35,-59],[-4,-31],[-41,-89],[-15,-58],[-32,-17],[-53,-74],[-33,-75],[-47,-74],[-22,-2],[-6,-134],[32,-92],[15,-89],[2,-46],[15,-59],[7,-89],[-1,-82],[25,4],[-5,-64],[9,-70],[-23,-192],[19,-6],[-11,-80],[-32,-85],[-62,-64],[-88,-57],[-55,-44],[-63,-89],[-22,-82],[38,-57],[-9,-190],[-22,-122],[-14,-136],[-28,-97],[-29,-50],[-28,-22],[-50,-101],[-36,-121],[-83,-245],[-36,-81],[-26,-34],[-69,-123],[-31,-66],[-112,-175],[-89,-107],[-73,-55],[-50,11],[-38,-32],[-2,-37],[-72,9],[-19,-45],[-27,-1],[-46,25],[-68,18],[-36,-22],[-81,16],[-34,-13],[-52,-70],[-61,-8],[-22,10],[-60,-23],[-57,-74],[-44,8],[-35,41],[-25,51],[-31,-3],[-2,59],[-31,5],[-21,34],[9,51],[-62,173],[-3,31],[45,40],[8,34],[-2,88],[-11,87],[-60,167],[-55,211],[-27,160],[-24,90],[-29,85],[-87,155],[-42,133],[-24,138],[4,43],[-23,65],[-13,136],[-1,159],[-38,191],[-3,211],[-8,72],[13,61],[-23,118],[-39,97],[-15,68],[-44,128],[-32,168],[-14,36],[-66,253],[-41,89],[-37,123],[-6,123],[11,176],[-2,165],[-5,35],[24,46],[35,228],[12,138],[14,65],[5,73],[39,94],[10,58],[49,59],[30,91],[12,71],[7,175],[-15,96],[-20,49],[-37,165],[-7,72],[-17,79],[40,83],[3,73],[-33,135],[-26,126],[-4,64],[-34,83],[-25,115],[12,24],[-29,81],[6,38],[-19,100],[-29,82],[0,26],[-85,169],[-9,38],[-35,61],[-33,88],[-67,114],[-56,159],[8,67],[-30,54],[-25,89],[-3,47],[36,31],[18,51],[13,127],[37,-27],[4,25],[-34,39],[-16,57],[33,39],[-13,77],[-2,73],[16,30],[20,84],[1,140],[15,124],[-8,53],[-24,52],[-27,89],[-49,39],[-11,79],[-16,33],[-21,-12],[-16,50],[-16,-43],[-128,-10],[-16,-27],[-29,-12],[-77,-12],[-56,81],[-23,111],[1,63],[-18,19],[-42,124],[-48,74],[-35,15],[-153,-8],[-108,-27],[-57,-26],[-23,-22],[-18,-55],[-23,-12],[-55,0],[-69,-59],[-51,-63],[-30,-11],[-65,-46],[-41,-50],[-45,38],[-96,43],[-89,26],[-173,-32],[-39,-18],[-106,-78],[-65,-68],[-27,-2],[-68,51],[-98,107],[-130,235],[-57,54],[-8,35],[-50,51],[-29,60],[-43,54],[-63,47],[-2,64],[-42,44],[-9,66],[-22,16],[-14,50],[11,33],[-10,137],[-34,85],[-14,-1],[-1,72],[-82,91],[-40,157],[-30,2],[-20,50],[-20,14],[-9,44],[-3,72],[8,51],[-44,-41],[-15,45],[-21,-10],[-34,66],[-31,34],[-6,183],[-6,44],[26,80],[-17,35],[-6,58],[-42,126],[-28,34],[30,30],[34,86],[35,125],[0,103],[22,147],[30,142],[5,79],[-6,145],[-14,111],[-34,82],[26,98],[8,101],[-25,98],[-22,-4],[-28,104],[-14,-11],[8,198],[15,59],[32,40],[17,60],[17,116],[46,149],[-11,22],[36,53],[45,93],[17,16],[20,74],[7,123],[31,119],[12,76],[52,54],[43,56],[45,212],[26,59],[109,50],[63,58],[39,76],[67,81],[72,171],[21,69],[2,77],[-25,63],[7,162],[15,66],[37,86],[12,112],[46,80],[27,61],[34,43],[83,61],[75,76],[16,37],[46,148],[48,232],[60,33],[33,-106],[30,-43],[58,-27],[72,27],[55,-10],[26,39],[15,-64],[93,-10],[39,23],[38,41],[48,70],[55,44],[27,-9],[27,18],[12,38],[41,45],[84,60],[151,18],[43,42],[62,3],[29,23],[111,0],[49,-50],[26,-1],[72,44],[48,51],[50,-39],[77,19],[34,-33],[103,32],[36,46],[62,33],[58,-31],[-1,-39],[25,-70],[18,34],[55,45],[8,-46],[-37,-87],[-31,-39],[-6,-35],[13,-66],[47,-58],[13,-90],[-29,-82],[-37,-78],[-38,-46],[-17,-66],[13,-48],[63,-77],[41,11],[8,-44],[28,-33],[47,-31],[52,-54],[54,-13],[60,26],[41,-26],[67,-27],[31,-39],[75,-28],[21,-53],[15,-116],[24,-53],[47,-37],[75,-11],[64,-31],[95,-70],[54,-83],[30,-29],[40,0],[48,46],[34,72],[16,63],[-22,110],[-4,60],[22,92],[57,83],[50,45],[41,6],[24,32],[62,-5],[102,-69],[2,-66],[20,-27],[58,-12],[38,-35],[62,2],[39,-31],[23,-80],[17,-5],[58,25],[153,-56],[33,-38],[49,-26],[62,-11],[63,-51],[58,42],[39,50],[33,7],[19,46],[35,-13],[34,14],[22,28],[37,-29],[41,19],[-13,-57],[37,-46],[25,37],[32,-43],[60,23],[60,-9],[58,40],[35,66],[22,70],[27,165],[38,179],[40,147],[4,50],[38,68],[-9,70],[5,86],[-20,80],[21,98],[-16,71],[43,81],[-16,57],[-58,-72],[-41,9],[-41,38],[-24,-4],[-34,-41],[-31,-60],[-38,-37],[-101,-34],[-47,34],[-41,81],[-75,61],[-81,15],[-18,-127],[-28,-1],[-62,-35],[-61,56],[-12,68],[-36,1],[-27,24],[-77,-16],[37,68],[-88,-2],[19,52],[-34,32],[-16,69],[18,67],[-63,50],[-44,18],[16,38],[49,-16],[-11,78],[28,41],[-38,93],[17,62],[-81,-22],[8,120],[15,8],[48,86],[61,13],[22,-31],[92,19],[56,30],[66,63],[-36,46],[5,39],[26,11],[67,-17],[49,10],[52,-25],[51,5],[22,49],[72,61],[24,33],[122,66],[155,-14],[29,24],[33,-77],[29,-21],[56,11],[14,-59],[38,-38],[30,24],[33,-42],[79,-24],[70,-35],[118,42],[43,-29],[52,-5],[78,58],[52,48],[32,65],[7,61],[-31,158],[-41,39],[-33,54],[-35,13],[-74,81],[-61,98],[-69,90],[-61,30],[-37,64],[-50,8],[-24,55],[-66,49],[37,24],[79,20],[-3,43],[36,100],[28,22],[-24,120],[150,112],[-92,16],[-34,-22],[-72,-1],[-23,-36],[-120,-60],[-51,-8],[-49,-55],[-19,13],[-43,-57],[1,-47],[19,-66],[39,-80],[21,-9],[69,32],[46,-13],[-20,-76],[-59,-14],[-45,22],[-44,-68],[-42,1],[-49,-62],[-42,-35],[-52,39],[18,81],[-6,44],[-42,22],[-30,35],[-47,13],[72,79],[59,46],[-52,52],[-82,-21],[-50,41],[-3,39],[40,10],[-44,54],[-133,-36],[-33,-103],[-33,-55],[-69,-42],[14,-78],[-19,-115],[-58,-20],[-46,-106],[2,-71],[-11,-111],[-11,-26],[-38,1],[-23,-48],[-4,-100],[-46,-65],[25,-27],[31,-69],[1,-45],[41,-89],[80,-54],[-11,-51],[-20,-8],[-69,25],[-28,-21],[-48,-4],[-28,-66],[-55,-43],[-50,-86],[-8,44],[61,71],[-78,-3],[-28,53],[-85,35],[-35,-32],[-47,13],[-57,-85],[-24,-73],[-40,-2],[-77,64],[-4,-105],[37,-94],[8,-68],[-17,-17],[14,-52],[-22,-53],[54,-32],[49,-71],[32,-18],[10,-130],[-53,69],[-44,-12],[0,-77],[33,-42],[-37,-24],[-45,13],[35,-142],[-52,1],[-14,-48],[-46,106],[-21,-67],[-35,79],[12,52],[-14,53],[-27,30],[-21,57],[29,64],[29,-5],[19,35],[83,-49],[27,-31],[37,22],[-80,84],[-22,-20],[-29,14],[-56,-21],[-41,14],[-13,62],[-25,35],[1,47],[-54,73],[-34,88],[-17,76],[-52,56],[-6,77],[14,72],[-3,113],[16,83],[-44,37],[-34,69],[-121,126],[-104,137],[-57,34],[-46,-5],[-5,31],[-50,58],[-43,80],[17,29],[-43,75],[-4,86],[-53,60],[-48,-117],[-38,63],[-13,86],[30,33],[-65,33],[-20,-30],[-60,-41],[-25,-3],[-6,-47],[34,-63],[-31,-57],[17,-114],[33,-53],[99,-97],[27,-90],[23,-113],[20,-42],[77,-104],[34,-27],[90,1],[21,-42],[-27,-31],[12,-43],[123,-86],[42,-51],[54,-42],[42,-65],[18,-62],[-16,-65],[-30,27],[-24,79],[-44,8],[-29,39],[-33,-6],[-45,-138],[8,-51],[58,-59],[7,-88],[-47,-25],[-23,-40],[-1,-70],[-30,-37],[-25,-71],[-38,0],[-9,54],[20,29],[17,95],[25,10],[-19,137],[-38,146],[-45,18],[-39,39],[4,32],[-24,67],[-35,-8],[-17,41],[-23,3],[-42,94],[-55,12],[-18,-10],[-52,52],[-63,108],[-30,33],[-19,47],[-53,57],[-52,92],[-44,132],[-9,78],[-22,39],[-36,19],[-50,50],[-59,23],[-24,-17],[-53,-98],[-102,-60],[-78,-113],[-43,-22],[-34,6],[-46,30],[-39,50],[-49,-13],[-65,47],[-90,-92],[-23,-64],[22,-223],[-28,-41],[-78,-69],[-26,-41],[-119,-52],[-36,-55],[-39,-115],[-50,-102],[-28,-82],[13,-105],[46,-70],[-66,-74],[-30,-65],[-19,-132],[-58,-4],[-53,-75],[-35,-105],[-55,7],[-22,-22],[-53,10],[-125,-9],[-35,-48],[-56,-19],[-22,-67],[-30,-25],[-47,38],[-26,93],[-13,10],[-12,73],[-44,55],[-69,-6],[-39,-38],[-86,27],[-45,-21],[21,92],[2,69],[-10,52],[9,52],[-13,75],[-34,-14],[-4,48],[-25,17],[11,139],[26,47],[35,132],[-6,15],[27,195],[-27,247],[23,39],[-34,71],[-6,53],[-22,37],[6,46],[35,36],[38,1],[33,24],[-1,32],[63,42],[49,-39],[181,-3],[129,-38],[104,24],[63,-34],[19,19],[61,-29],[80,25],[17,29],[15,105],[12,124],[6,138],[14,148],[-9,116],[-73,47],[-30,68],[4,52],[-85,110],[-44,41],[-84,33],[-46,-3],[-36,60],[27,65],[-37,-1],[0,40],[74,39],[39,1],[55,30],[26,-12],[35,-58],[28,25],[114,-1],[-15,38],[-2,91],[-30,110],[67,-1],[14,-67],[110,-21],[50,39],[-10,55],[48,36],[71,32],[35,67],[3,103],[38,58],[60,17],[129,89],[102,197],[22,115],[77,65],[19,40],[60,32],[85,8],[53,55],[82,2],[60,-31],[9,73],[32,-9],[1,98],[-29,31],[35,33],[-32,58],[-7,144],[-55,42],[8,88],[-9,36],[5,107],[13,56],[38,60],[77,8],[29,20],[46,80],[34,8],[31,-39],[-1,-47],[-26,-56],[-1,-87],[63,-23],[1,-52],[-60,-19],[-16,-79],[-67,-86],[6,-108],[65,-141],[16,7],[74,-59],[44,-55],[45,46],[35,6],[53,68],[51,-13],[48,-62],[31,3],[16,-69],[70,-41],[26,79],[166,63],[42,60],[80,41],[120,25],[39,-94],[35,-19],[49,9],[39,36],[25,87],[61,17],[76,-14],[6,119],[-20,83],[-6,102],[7,130],[31,44],[13,70],[30,58],[93,35],[66,-92],[17,-54],[41,-27],[72,46],[14,35],[-12,105],[26,131],[-47,-20],[-39,22],[-30,68],[-1,124],[66,22],[34,41],[120,12],[8,27],[126,-20],[39,-24],[104,-8],[64,92],[56,3],[22,40],[110,-29],[-45,73],[-74,0],[-48,42],[1,54],[-134,-34],[-25,17],[-80,-29],[-100,-18],[-56,-32],[-163,-52],[-104,25],[-25,65],[-129,55],[-9,85],[28,143],[-39,92],[11,54],[-29,91],[11,39],[51,57],[-12,38],[65,9],[14,45],[114,105],[116,142],[32,66],[92,48],[-3,108],[-86,63],[-101,17],[-67,-22],[-83,25],[-13,-56],[-77,-44],[1,-65],[-50,-73],[37,-98],[-50,-46],[-29,-72],[-118,-93],[-41,4],[-37,-55],[-48,-14],[-16,-65],[-64,-15],[-7,-75],[-53,-18],[22,-49],[-17,-108],[-32,-38],[14,-201],[79,-26],[68,-77],[43,-72],[4,-42],[-112,-99]],[[26091,31464],[-30,-14],[-47,-56],[-74,-27],[-30,-76],[-45,-23],[-60,-81],[-71,-31],[-4,-33],[-43,-126],[-40,-58],[-1,-54],[54,-52],[32,-157],[-5,-109],[86,-197],[18,-52],[32,-38],[50,-112],[13,-63],[62,-102],[25,-1],[35,-38],[-85,-50],[-26,-135],[5,-59],[-19,-64],[-21,14],[-18,-73],[2,-97],[10,-115],[14,-52],[44,-39],[74,-20],[46,-91],[66,-62],[73,-30],[48,2],[134,57],[73,19],[-19,165],[-3,59],[5,235],[-28,61],[-56,50],[13,40],[26,14],[-7,67],[-47,6],[-28,54],[-1,41],[24,147],[22,-49],[54,-2],[28,-39],[51,43],[44,13],[-13,70],[-56,74],[-28,132],[-25,10],[-52,-10],[-22,-27],[-16,-145],[-25,43],[-13,57],[-4,74],[18,75],[-3,76],[-85,35],[-26,57],[-36,4],[1,70],[-18,39],[-35,125],[-57,30],[9,69],[51,1],[36,-29],[-16,110],[43,86],[36,10],[90,0],[62,-21],[-35,61],[41,142],[-6,82],[10,28],[-29,65],[-49,9],[-33,-35],[-61,41],[-53,21],[-56,-39]]],"transform":{"scale":[0.00884151582876931,0.004352258509194726],"translate":[-180,-89.99892578125002]},"objects":{"ne_50m_land":{"type":"GeometryCollection","geometries":[{"arcs":[[0]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[1]],"type":"Polygon"},{"arcs":[[2]],"type":"Polygon"},{"arcs":[[3]],"type":"Polygon"},{"arcs":[[4]],"type":"Polygon"},{"arcs":[[5]],"type":"Polygon"},{"arcs":[[6]],"type":"Polygon"},{"arcs":[[7]],"type":"Polygon"},{"arcs":[[8]],"type":"Polygon"},{"arcs":[[9]],"type":"Polygon"},{"arcs":[[10]],"type":"Polygon"},{"arcs":[[11]],"type":"Polygon"},{"arcs":[[12]],"type":"Polygon"},{"arcs":[[13]],"type":"Polygon"},{"arcs":[[14]],"type":"Polygon"},{"arcs":[[15]],"type":"Polygon"},{"arcs":[[16]],"type":"Polygon"},{"arcs":[[17]],"type":"Polygon"},{"arcs":[[18]],"type":"Polygon"},{"arcs":[[19]],"type":"Polygon"},{"arcs":[[20]],"type":"Polygon"},{"arcs":[[21]],"type":"Polygon"},{"arcs":[[22]],"type":"Polygon"},{"arcs":[[23]],"type":"Polygon"},{"arcs":[[24]],"type":"Polygon"},{"arcs":[[25]],"type":"Polygon"},{"arcs":[[26]],"type":"Polygon"},{"arcs":[[27]],"type":"Polygon"},{"arcs":[[28]],"type":"Polygon"},{"arcs":[[29]],"type":"Polygon"},{"arcs":[[30]],"type":"Polygon"},{"arcs":[[31]],"type":"Polygon"},{"arcs":[[32]],"type":"Polygon"},{"arcs":[[33]],"type":"Polygon"},{"arcs":[[34]],"type":"Polygon"},{"arcs":[[35]],"type":"Polygon"},{"arcs":[[36]],"type":"Polygon"},{"arcs":[[37]],"type":"Polygon"},{"arcs":[[38]],"type":"Polygon"},{"arcs":[[39]],"type":"Polygon"},{"arcs":[[40]],"type":"Polygon"},{"arcs":[[41]],"type":"Polygon"},{"arcs":[[42]],"type":"Polygon"},{"arcs":[[43]],"type":"Polygon"},{"arcs":[[44]],"type":"Polygon"},{"arcs":[[45]],"type":"Polygon"},{"arcs":[[46]],"type":"Polygon"},{"arcs":[[47]],"type":"Polygon"},{"arcs":[[48]],"type":"Polygon"},{"arcs":[[49]],"type":"Polygon"},{"arcs":[[50]],"type":"Polygon"},{"arcs":[[51]],"type":"Polygon"},{"arcs":[[52]],"type":"Polygon"},{"arcs":[[53]],"type":"Polygon"},{"arcs":[[54]],"type":"Polygon"},{"arcs":[[55]],"type":"Polygon"},{"arcs":[[56]],"type":"Polygon"},{"arcs":[[57]],"type":"Polygon"},{"arcs":[[58]],"type":"Polygon"},{"arcs":[[59]],"type":"Polygon"},{"arcs":[[60]],"type":"Polygon"},{"arcs":[[61]],"type":"Polygon"},{"arcs":[[62]],"type":"Polygon"},{"arcs":[[63]],"type":"Polygon"},{"arcs":[[64]],"type":"Polygon"},{"arcs":[[65]],"type":"Polygon"},{"arcs":[[66]],"type":"Polygon"},{"arcs":[[67]],"type":"Polygon"},{"arcs":[[68]],"type":"Polygon"},{"arcs":[[69]],"type":"Polygon"},{"arcs":[[70]],"type":"Polygon"},{"arcs":[[71]],"type":"Polygon"},{"arcs":[[72]],"type":"Polygon"},{"arcs":[[73]],"type":"Polygon"},{"arcs":[[74]],"type":"Polygon"},{"arcs":[[75]],"type":"Polygon"},{"arcs":[[76]],"type":"Polygon"},{"arcs":[[77]],"type":"Polygon"},{"arcs":[[78]],"type":"Polygon"},{"arcs":[[79]],"type":"Polygon"},{"arcs":[[80]],"type":"Polygon"},{"arcs":[[81]],"type":"Polygon"},{"arcs":[[82]],"type":"Polygon"},{"arcs":[[83]],"type":"Polygon"},{"arcs":[[84]],"type":"Polygon"},{"arcs":[[85]],"type":"Polygon"},{"arcs":[[86]],"type":"Polygon"},{"arcs":[[87]],"type":"Polygon"},{"arcs":[[88]],"type":"Polygon"},{"arcs":[[89]],"type":"Polygon"},{"arcs":[[90]],"type":"Polygon"},{"arcs":[[91]],"type":"Polygon"},{"arcs":[[92]],"type":"Polygon"},{"arcs":[[93]],"type":"Polygon"},{"arcs":[[94]],"type":"Polygon"},{"arcs":[[95]],"type":"Polygon"},{"arcs":[[96]],"type":"Polygon"},{"arcs":[[97]],"type":"Polygon"},{"arcs":[[98]],"type":"Polygon"},{"arcs":[[99]],"type":"Polygon"},{"arcs":[[100]],"type":"Polygon"},{"arcs":[[101]],"type":"Polygon"},{"arcs":[[102]],"type":"Polygon"},{"arcs":[[103]],"type":"Polygon"},{"arcs":[[104]],"type":"Polygon"},{"arcs":[[105]],"type":"Polygon"},{"arcs":[[106]],"type":"Polygon"},{"arcs":[[107]],"type":"Polygon"},{"arcs":[[108]],"type":"Polygon"},{"arcs":[[109]],"type":"Polygon"},{"arcs":[[110]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"arcs":[[111]],"type":"Polygon"},{"arcs":[[112]],"type":"Polygon"},{"arcs":[[113]],"type":"Polygon"},{"arcs":[[114]],"type":"Polygon"},{"arcs":[[115]],"type":"Polygon"},{"arcs":[[116]],"type":"Polygon"},{"arcs":[[117]],"type":"Polygon"},{"arcs":[[118]],"type":"Polygon"},{"arcs":[[119]],"type":"Polygon"},{"arcs":[[120]],"type":"Polygon"},{"arcs":[[121]],"type":"Polygon"},{"arcs":[[122]],"type":"Polygon"},{"arcs":[[123]],"type":"Polygon"},{"arcs":[[124]],"type":"Polygon"},{"arcs":[[125]],"type":"Polygon"},{"arcs":[[126]],"type":"Polygon"},{"arcs":[[127]],"type":"Polygon"},{"arcs":[[128]],"type":"Polygon"},{"arcs":[[129]],"type":"Polygon"},{"arcs":[[130]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[131]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[132]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[133]],"type":"Polygon"},{"arcs":[[134]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[135]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[136]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[137]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[138]],"type":"Polygon"},{"arcs":[[139]],"type":"Polygon"},{"type":null},{"arcs":[[140]],"type":"Polygon"},{"type":null},{"arcs":[[141]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[142]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[143]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[144]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[145]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"arcs":[[146]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[147]],"type":"Polygon"},{"arcs":[[148]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[149]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"arcs":[[150]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[151]],"type":"Polygon"},{"arcs":[[152]],"type":"Polygon"},{"type":null},{"arcs":[[153]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[154]],"type":"Polygon"},{"type":null},{"arcs":[[155]],"type":"Polygon"},{"arcs":[[156]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[157]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[158]],"type":"Polygon"},{"arcs":[[159]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"arcs":[[160]],"type":"Polygon"},{"arcs":[[161]],"type":"Polygon"},{"type":null},{"arcs":[[162]],"type":"Polygon"},{"type":null},{"arcs":[[163]],"type":"Polygon"},{"arcs":[[164]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[165]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[166]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[167]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[168]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[169]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[170]],"type":"Polygon"},{"arcs":[[171]],"type":"Polygon"},{"arcs":[[172]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[173]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[174]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[175]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"arcs":[[176]],"type":"Polygon"},{"arcs":[[177]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[178]],"type":"Polygon"},{"type":null},{"arcs":[[179]],"type":"Polygon"},{"type":null},{"arcs":[[180]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[181]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[182]],"type":"Polygon"},{"type":null},{"arcs":[[183]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"arcs":[[184]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"arcs":[[185]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[186]],"type":"Polygon"},{"arcs":[[187]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[188]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[189]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[190]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[191]],"type":"Polygon"},{"arcs":[[192]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[193]],"type":"Polygon"},{"arcs":[[194]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[195]],"type":"Polygon"},{"type":null},{"arcs":[[196]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[197]],"type":"Polygon"},{"arcs":[[198]],"type":"Polygon"},{"arcs":[[199]],"type":"Polygon"},{"arcs":[[200]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[201]],"type":"Polygon"},{"arcs":[[202]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"arcs":[[203]],"type":"Polygon"},{"arcs":[[204]],"type":"Polygon"},{"type":null},{"arcs":[[205]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[206]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[207]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[208]],"type":"Polygon"},{"type":null},{"arcs":[[209]],"type":"Polygon"},{"arcs":[[210]],"type":"Polygon"},{"arcs":[[211]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[212]],"type":"Polygon"},{"type":null},{"arcs":[[213]],"type":"Polygon"},{"type":null},{"arcs":[[214]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[215]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[216]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[217]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"arcs":[[218]],"type":"Polygon"},{"type":null},{"arcs":[[219]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"arcs":[[220]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[221]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[222]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[223]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[224]],"type":"Polygon"},{"arcs":[[225]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[226]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"arcs":[[227]],"type":"Polygon"},{"arcs":[[228]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[229]],"type":"Polygon"},{"arcs":[[230]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"arcs":[[231]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[232]],"type":"Polygon"},{"arcs":[[233]],"type":"Polygon"},{"type":null},{"arcs":[[234]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[235]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[236]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[237]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[238]],"type":"Polygon"},{"type":null},{"arcs":[[239]],"type":"Polygon"},{"type":null},{"arcs":[[240]],"type":"Polygon"},{"arcs":[[241]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[242]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[243]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[244]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[245]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[246]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[247]],"type":"Polygon"},{"type":null},{"arcs":[[248]],"type":"Polygon"},{"type":null},{"arcs":[[249]],"type":"Polygon"},{"type":null},{"arcs":[[250]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[251]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"arcs":[[252]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"arcs":[[253]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[254]],"type":"Polygon"},{"arcs":[[255]],"type":"Polygon"},{"arcs":[[256]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[257]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[258]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[259]],"type":"Polygon"},{"type":null},{"arcs":[[260]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[261]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[262]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[263]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[264]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[265]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[266]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[267]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[268]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[269]],"type":"Polygon"},{"arcs":[[270]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[271]],"type":"Polygon"},{"type":null},{"arcs":[[272]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[273]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[274]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"arcs":[[275]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[276]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[277]],"type":"Polygon"},{"type":null},{"arcs":[[278]],"type":"Polygon"},{"arcs":[[279]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[280]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[281]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[282]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[283]],"type":"Polygon"},{"arcs":[[284]],"type":"Polygon"},{"arcs":[[285]],"type":"Polygon"},{"arcs":[[286]],"type":"Polygon"},{"arcs":[[287]],"type":"Polygon"},{"arcs":[[288]],"type":"Polygon"},{"arcs":[[289]],"type":"Polygon"},{"arcs":[[290]],"type":"Polygon"},{"arcs":[[291]],"type":"Polygon"},{"arcs":[[292]],"type":"Polygon"},{"type":null},{"arcs":[[293]],"type":"Polygon"},{"arcs":[[294]],"type":"Polygon"},{"arcs":[[295]],"type":"Polygon"},{"arcs":[[296]],"type":"Polygon"},{"arcs":[[297]],"type":"Polygon"},{"arcs":[[298]],"type":"Polygon"},{"arcs":[[299]],"type":"Polygon"},{"arcs":[[300]],"type":"Polygon"},{"arcs":[[301]],"type":"Polygon"},{"arcs":[[302]],"type":"Polygon"},{"arcs":[[303]],"type":"Polygon"},{"arcs":[[304]],"type":"Polygon"},{"arcs":[[305]],"type":"Polygon"},{"arcs":[[306]],"type":"Polygon"},{"arcs":[[307]],"type":"Polygon"},{"arcs":[[308]],"type":"Polygon"},{"arcs":[[309]],"type":"Polygon"},{"arcs":[[310]],"type":"Polygon"},{"arcs":[[311]],"type":"Polygon"},{"arcs":[[312]],"type":"Polygon"},{"arcs":[[313]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[314]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[315]],"type":"Polygon"},{"type":null},{"arcs":[[316]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"arcs":[[317]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[318]],"type":"Polygon"},{"type":null},{"type":null},{"arcs":[[319]],"type":"Polygon"},{"type":null},{"arcs":[[320]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[321]],"type":"Polygon"},{"arcs":[[322]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[323]],"type":"Polygon"},{"arcs":[[324]],"type":"Polygon"},{"arcs":[[325]],"type":"Polygon"},{"arcs":[[326]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[327]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[328]],"type":"Polygon"},{"type":null},{"arcs":[[329]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"arcs":[[330]],"type":"Polygon"},{"arcs":[[331]],"type":"Polygon"},{"arcs":[[332],[333]],"type":"Polygon"},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null},{"type":null}]}}}
},{}],22:[function(_dereq_,module,exports){
"use strict";

var World = _dereq_('./World'); //


var somehow = function somehow(obj) {
  return new World(obj);
};

module.exports = somehow;

},{"./World":9}],23:[function(_dereq_,module,exports){
"use strict";

function _typeof(obj) { if (typeof Symbol === "function" && typeof Symbol.iterator === "symbol") { _typeof = function _typeof(obj) { return typeof obj; }; } else { _typeof = function _typeof(obj) { return obj && typeof Symbol === "function" && obj.constructor === Symbol && obj !== Symbol.prototype ? "symbol" : typeof obj; }; } return _typeof(obj); }

function _templateObject() {
  var data = _taggedTemplateLiteral(["<circle ...", " ><title>", "</title></circle>"]);

  _templateObject = function _templateObject() {
    return data;
  };

  return data;
}

function _taggedTemplateLiteral(strings, raw) { if (!raw) { raw = strings.slice(0); } return Object.freeze(Object.defineProperties(strings, { raw: { value: Object.freeze(raw) } })); }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _possibleConstructorReturn(self, call) { if (call && (_typeof(call) === "object" || typeof call === "function")) { return call; } return _assertThisInitialized(self); }

function _assertThisInitialized(self) { if (self === void 0) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return self; }

function _getPrototypeOf(o) { _getPrototypeOf = Object.setPrototypeOf ? Object.getPrototypeOf : function _getPrototypeOf(o) { return o.__proto__ || Object.getPrototypeOf(o); }; return _getPrototypeOf(o); }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function"); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, writable: true, configurable: true } }); if (superClass) _setPrototypeOf(subClass, superClass); }

function _setPrototypeOf(o, p) { _setPrototypeOf = Object.setPrototypeOf || function _setPrototypeOf(o, p) { o.__proto__ = p; return o; }; return _setPrototypeOf(o, p); }

var Shape = _dereq_('./Shape');

var colors = _dereq_('spencer-color').colors;

var defaults = {
  fill: colors.blue,
  stroke: 'none'
};

var Dot =
/*#__PURE__*/
function (_Shape) {
  _inherits(Dot, _Shape);

  function Dot() {
    var _this;

    var obj = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : {};
    var world = arguments.length > 1 ? arguments[1] : undefined;

    _classCallCheck(this, Dot);

    obj = Object.assign({}, defaults, obj);
    _this = _possibleConstructorReturn(this, _getPrototypeOf(Dot).call(this, obj, world));
    _this._radius = obj.radius || 5;
    return _this;
  }

  _createClass(Dot, [{
    key: "build",
    value: function build() {
      var h = this.world.html;
      var point = this.world.projection([this.point[1], this.point[0]]);
      var attrs = Object.assign({}, this.attrs, {
        id: this._id,
        cx: point[0],
        cy: point[1],
        r: this._radius
      });
      return h(_templateObject(), attrs, this._title);
    }
  }]);

  return Dot;
}(Shape);

module.exports = Dot;

},{"./Shape":25,"spencer-color":6}],24:[function(_dereq_,module,exports){
"use strict";

function _typeof(obj) { if (typeof Symbol === "function" && typeof Symbol.iterator === "symbol") { _typeof = function _typeof(obj) { return typeof obj; }; } else { _typeof = function _typeof(obj) { return obj && typeof Symbol === "function" && obj.constructor === Symbol && obj !== Symbol.prototype ? "symbol" : typeof obj; }; } return _typeof(obj); }

function _templateObject2() {
  var data = _taggedTemplateLiteral(["<g>\n      ", "\n      <path ...", "></path>\n      ", "\n    </g>"]);

  _templateObject2 = function _templateObject2() {
    return data;
  };

  return data;
}

function _templateObject() {
  var data = _taggedTemplateLiteral(["<circle cx=\"", "\" cy=\"", "\" r=\"5\" fill=\"", "\"></circle>"]);

  _templateObject = function _templateObject() {
    return data;
  };

  return data;
}

function _taggedTemplateLiteral(strings, raw) { if (!raw) { raw = strings.slice(0); } return Object.freeze(Object.defineProperties(strings, { raw: { value: Object.freeze(raw) } })); }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _possibleConstructorReturn(self, call) { if (call && (_typeof(call) === "object" || typeof call === "function")) { return call; } return _assertThisInitialized(self); }

function _assertThisInitialized(self) { if (self === void 0) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return self; }

function _getPrototypeOf(o) { _getPrototypeOf = Object.setPrototypeOf ? Object.getPrototypeOf : function _getPrototypeOf(o) { return o.__proto__ || Object.getPrototypeOf(o); }; return _getPrototypeOf(o); }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function"); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, writable: true, configurable: true } }); if (superClass) _setPrototypeOf(subClass, superClass); }

function _setPrototypeOf(o, p) { _setPrototypeOf = Object.setPrototypeOf || function _setPrototypeOf(o, p) { o.__proto__ = p; return o; }; return _setPrototypeOf(o, p); }

var Shape = _dereq_('./Shape');

var colors = _dereq_('spencer-color').colors;

var d3Geo = _dereq_('d3-geo');

var data = _dereq_('../data');

var defaults = {
  fill: 'none',
  stroke: colors.blue,
  'stroke-width': 4
};

var Line =
/*#__PURE__*/
function (_Shape) {
  _inherits(Line, _Shape);

  function Line() {
    var _this;

    var obj = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : {};
    var world = arguments.length > 1 ? arguments[1] : undefined;

    _classCallCheck(this, Line);

    obj = Object.assign({}, defaults, obj);
    _this = _possibleConstructorReturn(this, _getPrototypeOf(Line).call(this, obj, world));
    _this._radius = obj.radius || 5;
    _this._data = [];
    return _this;
  }

  _createClass(Line, [{
    key: "color",
    value: function color(c) {
      this.attrs.stroke = colors[c] || c;
      return this;
    }
  }, {
    key: "toData",
    value: function toData() {
      return {
        type: 'Feature',
        geometry: {
          type: 'LineString',
          coordinates: [this._data[0].reverse(), this._data[1].reverse()]
        }
      };
    }
  }, {
    key: "from",
    value: function from(str) {
      this._data[0] = data.points[str] || str;
      return this;
    }
  }, {
    key: "to",
    value: function to(str) {
      this._data[1] = data.points[str] || str;
      return this;
    }
  }, {
    key: "makePoint",
    value: function makePoint(arr) {
      var h = this.world.html;
      var point = this.world.projection([arr[0], arr[1]]);
      console.log(point);
      return h(_templateObject(), point[0], point[1], this.attrs.stroke);
    }
  }, {
    key: "build",
    value: function build() {
      var h = this.world.html;
      var projection = this.world.projection;
      var toPath = d3Geo.geoPath().projection(projection);
      var geoJSON = this.toData();
      var d = toPath(geoJSON);
      var attrs = Object.assign({}, this.attrs, {
        id: this._id,
        d: d
      });
      return h(_templateObject2(), this.makePoint(this._data[0]), attrs, this.makePoint(this._data[1]));
    }
  }]);

  return Line;
}(Shape);

module.exports = Line;

},{"../data":10,"./Shape":25,"d3-geo":3,"spencer-color":6}],25:[function(_dereq_,module,exports){
"use strict";

function _templateObject() {
  var data = _taggedTemplateLiteral(["<path ...", "></path>"]);

  _templateObject = function _templateObject() {
    return data;
  };

  return data;
}

function _taggedTemplateLiteral(strings, raw) { if (!raw) { raw = strings.slice(0); } return Object.freeze(Object.defineProperties(strings, { raw: { value: Object.freeze(raw) } })); }

function _typeof(obj) { if (typeof Symbol === "function" && typeof Symbol.iterator === "symbol") { _typeof = function _typeof(obj) { return typeof obj; }; } else { _typeof = function _typeof(obj) { return obj && typeof Symbol === "function" && obj.constructor === Symbol && obj !== Symbol.prototype ? "symbol" : typeof obj; }; } return _typeof(obj); }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

var data = _dereq_('../data');

var d3Geo = _dereq_('d3-geo');

var topojson = _dereq_('topojson-client');

var colors = _dereq_('spencer-color').colors;

var defaults = {
  fill: 'none',
  stroke: 'grey'
};

var Shape =
/*#__PURE__*/
function () {
  function Shape() {
    var obj = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : {};
    var world = arguments.length > 1 ? arguments[1] : undefined;

    _classCallCheck(this, Shape);

    this.world = world;
    this.data = obj.data || [];
    this._id = obj.id; //|| fns.uid('input')

    this.attrs = Object.assign({}, defaults, obj);
    this.shape = data.shapes[obj.shape] || data.points[obj.point] || obj.shape;
    this.point = [];
    this.style = {};
  }

  _createClass(Shape, [{
    key: "at",
    value: function at(str) {
      if (typeof str === 'string') {
        str = str.toLowerCase().trim();
        this.point = data.points[str];
      }

      return this;
    }
  }, {
    key: "color",
    value: function color(c) {
      this.attrs.fill = colors[c] || c;
      return this;
    }
  }, {
    key: "stroke",
    value: function stroke(c) {
      this.attrs.stroke = colors[c] || c;
      return this;
    }
  }, {
    key: "fill",
    value: function fill(c) {
      this.attrs.fill = colors[c] || c;
      return this;
    }
  }, {
    key: "opacity",
    value: function opacity(n) {
      this.attrs.opacity = n;
      return this;
    }
  }, {
    key: "drawSyle",
    value: function drawSyle() {
      var _this = this;

      return Object.keys(this.style).map(function (k) {
        return "".concat(k, ":").concat(_this.style[k], ";");
      }).join(' ');
    }
  }, {
    key: "geoJSON",
    value: function geoJSON() {
      if (_typeof(this.shape) === 'object') {
        var key = Object.keys(this.shape.objects)[0];
        return topojson.feature(this.shape, this.shape.objects[key]);
      }

      return [];
    }
  }, {
    key: "build",
    value: function build() {
      var h = this.world.html;
      var projection = this.world.projection;
      var toPath = d3Geo.geoPath().projection(projection);
      var geojson = this.geoJSON();
      var d = toPath(geojson);
      var attrs = Object.assign({}, this.attrs, {
        id: this._id,
        d: d
      });
      return h(_templateObject(), attrs);
    }
  }]);

  return Shape;
}();

module.exports = Shape;

},{"../data":10,"d3-geo":3,"spencer-color":6,"topojson-client":7}],26:[function(_dereq_,module,exports){
"use strict";

function _typeof(obj) { if (typeof Symbol === "function" && typeof Symbol.iterator === "symbol") { _typeof = function _typeof(obj) { return typeof obj; }; } else { _typeof = function _typeof(obj) { return obj && typeof Symbol === "function" && obj.constructor === Symbol && obj !== Symbol.prototype ? "symbol" : typeof obj; }; } return _typeof(obj); }

function _templateObject2() {
  var data = _taggedTemplateLiteral(["<g transform=\"", "\" style=\"", "\">\n    <text ...", " >\n      ", "\n    </text>\n  </g>"]);

  _templateObject2 = function _templateObject2() {
    return data;
  };

  return data;
}

function _templateObject() {
  var data = _taggedTemplateLiteral(["<tspan x=\"0\" dy=\"1.2em\">", "</tspan>"]);

  _templateObject = function _templateObject() {
    return data;
  };

  return data;
}

function _taggedTemplateLiteral(strings, raw) { if (!raw) { raw = strings.slice(0); } return Object.freeze(Object.defineProperties(strings, { raw: { value: Object.freeze(raw) } })); }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

function _possibleConstructorReturn(self, call) { if (call && (_typeof(call) === "object" || typeof call === "function")) { return call; } return _assertThisInitialized(self); }

function _assertThisInitialized(self) { if (self === void 0) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return self; }

function _getPrototypeOf(o) { _getPrototypeOf = Object.setPrototypeOf ? Object.getPrototypeOf : function _getPrototypeOf(o) { return o.__proto__ || Object.getPrototypeOf(o); }; return _getPrototypeOf(o); }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function"); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, writable: true, configurable: true } }); if (superClass) _setPrototypeOf(subClass, superClass); }

function _setPrototypeOf(o, p) { _setPrototypeOf = Object.setPrototypeOf || function _setPrototypeOf(o, p) { o.__proto__ = p; return o; }; return _setPrototypeOf(o, p); }

var Shape = _dereq_('./Shape');

var defaults = {
  fill: 'grey',
  stroke: 'none',
  'stroke-width': 1,
  'stroke-linecap': 'round',
  'text-anchor': 'start'
};

var Text =
/*#__PURE__*/
function (_Shape) {
  _inherits(Text, _Shape);

  function Text(obj, world) {
    var _this;

    _classCallCheck(this, Text);

    var text = null;
    var textFn = null;

    if (typeof obj === 'string') {
      text = [obj];
      obj = {};
    } else if (typeof obj === 'function') {
      textFn = obj;
      obj = {};
    } else if (Array.isArray(obj)) {
      text = obj;
      obj = [];
    }

    obj = Object.assign({}, defaults, obj);
    _this = _possibleConstructorReturn(this, _getPrototypeOf(Text).call(this, obj, world));
    _this.textFn = textFn;
    _this.textLines = text || obj.text || [];
    return _this;
  }

  _createClass(Text, [{
    key: "build",
    value: function build() {
      var h = this.world.html;
      var textArr = this.textLines;

      if (this.textFn !== null) {
        textArr = this.textFn(this.world);
        textArr = typeof textArr === 'string' ? [textArr] : textArr;
      }

      var point = this.world.projection([this.point[1], this.point[0]]);
      var inside = textArr.map(function (str) {
        return h(_templateObject(), String(str));
      });
      var transform = "translate(".concat(point[0], " ").concat(point[1], ")");
      return h(_templateObject2(), transform, this.drawSyle(), this.attrs, inside);
    }
  }]);

  return Text;
}(Shape);

module.exports = Text;

},{"./Shape":25}]},{},[1]);
