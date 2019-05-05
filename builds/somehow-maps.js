/* somehow-maps v0.0.1
   github.com/spencermountain/somehow-maps
   MIT
*/

(function(f){if(typeof exports==="object"&&typeof module!=="undefined"){module.exports=f()}else if(typeof define==="function"&&define.amd){define([],f)}else{var g;if(typeof window!=="undefined"){g=window}else if(typeof global!=="undefined"){g=global}else if(typeof self!=="undefined"){g=self}else{g=this}g.somehowMaps = f()}})(function(){var define,module,exports;return (function(){function r(e,n,t){function o(i,f){if(!n[i]){if(!e[i]){var c="function"==typeof require&&require;if(!f&&c)return c(i,!0);if(u)return u(i,!0);var a=new Error("Cannot find module '"+i+"'");throw a.code="MODULE_NOT_FOUND",a}var p=n[i]={exports:{}};e[i][0].call(p.exports,function(r){var n=e[i][1][r];return o(n||r)},p,p.exports,r,e,n,t)}return n[i].exports}for(var u="function"==typeof require&&require,i=0;i<t.length;i++)o(t[i]);return o}return r})()({1:[function(_dereq_,module,exports){
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

},{}],2:[function(_dereq_,module,exports){
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

},{"d3-array":1}],3:[function(_dereq_,module,exports){
(function (global){
!function(e){if("object"==typeof exports&&"undefined"!=typeof module)module.exports=e();else if("function"==typeof define&&define.amd)define([],e);else{("undefined"!=typeof window?window:"undefined"!=typeof global?global:"undefined"!=typeof self?self:this).fitAspect=e()}}(function(){return function a(o,s,d){function c(t,e){if(!s[t]){if(!o[t]){var i="function"==typeof _dereq_&&_dereq_;if(!e&&i)return i(t,!0);if(l)return l(t,!0);var n=new Error("Cannot find module '"+t+"'");throw n.code="MODULE_NOT_FOUND",n}var r=s[t]={exports:{}};o[t][0].call(r.exports,function(e){return c(o[t][1][e]||e)},r,r.exports,a,o,s,d)}return s[t].exports}for(var l="function"==typeof _dereq_&&_dereq_,e=0;e<d.length;e++)c(d[e]);return c}({1:[function(e,t,i){"use strict";var n=[{names:["square","1:1","instagram"],description:"Square",decimal:1,orientation:"landscape"},{names:["4:3","fullscreen","four three","1.33:1","ipad","pythagorean"],description:"Traditional TVs",decimal:1.333333,orientation:"landscape"},{names:["a4","√2:1","paper","lichtenberg","1:1.41"],description:"A4 paper",decimal:1.41},{names:["imax","1.43:1"],description:"IMAX film",decimal:1.43,orientation:"landscape"},{names:["3:2","35mm","photo","1.5:1","1.5"],description:"35mm photos",decimal:1.5,orientation:"landscape"},{names:["business card","bank card","1.58:1"],description:"Bank Cards",decimal:1.58577,orientation:"landscape"},{names:["golden","kepler","1.618","1.6:1"],description:"Golden ratio",decimal:1.61803,orientation:"landscape"},{names:["16:9","hd","hdtv","fhd","tv","computer","iphone","4k","8k","1.78:1"],description:"HD video",decimal:1.77777,orientation:"landscape"},{names:["widescreen","1.85:1"],description:"Movie-theatres",decimal:1.85,orientation:"landscape"},{names:["2:1","univisium","mobile","18:9"],description:"2:1",decimal:2,orientation:"landscape"},{names:["cinemascope","widescreen","wide","2.35:1","2.39:1"],description:"Widescreen",decimal:2.35,orientation:"landscape"},{names:["silver","1 + √2","2.41:1"],description:"Silver ratio",decimal:2.41,orientation:"landscape"}],r=n.map(function(e){return(e=Object.assign({},e)).decimal=1/e.decimal,e.orientation="portrait",e}),a={};n.forEach(function(t){t.names.forEach(function(e){a[e]=t})}),t.exports={lookup:a,portraits:r,list:n}},{}],2:[function(e,t,i){"use strict";var n=e("./aspects");t.exports=function(e,t){var i=e/t;return(i=parseInt(100*i,10)/100)<1?function(e,t){for(var i=0;i<t.length;i+=1)if(e>t[i].decimal){if(t[i-1]){var n=Math.abs(e-t[i].decimal);if(Math.abs(e-t[i-1].decimal)<n)return t[i-1]}return t[i]}return t[t.length-1]}(i,n.portraits):function(e,t){for(var i=0;i<t.length;i+=1)if(e<=t[i].decimal){if(t[i-1]){var n=Math.abs(e-t[i].decimal);if(Math.abs(e-t[i-1].decimal)<n)return t[i-1]}return t[i]}return t[t.length-1]}(i,n.list)}},{"./aspects":1}],3:[function(e,t,i){"use strict";var n=function(e,t){var i=1/t.decimal,n=e.orientation||"landscape";"portrait"===n&&(i=1/i);var r=e.width*i;return r=Math.round(r),{closest:t,width:e.width,height:r,orientation:n,original:e}},r=function(e,t){var i=t.decimal,n=e.orientation||"landscape";"portrait"===n&&(i=1/i);var r=e.height*i;return{closest:t,width:r=Math.round(r),height:e.height,orientation:n,original:e}};t.exports={both:function(e,t){var i=r(e,t);return i.width>e.width?n(e,t):i},width:r,height:n}},{}],4:[function(i,n,e){(function(e){"use strict";var o=i("./find-best-ratio"),s=i("./parse-ratio"),d=i("./fit"),t=function(){var e=0<arguments.length&&void 0!==arguments[0]?arguments[0]:{};if(!e.aspect&&!e.ratio){var t=o(e.width,e.height),i=1/t.decimal,n=e.width*i,r=(n-e.height)/e.height;return r=parseInt(1e3*r,10)/10,n=Math.round(n),{closest:t,percent_change:r,width:e.width,height:n}}var a=s(e.aspect||e.ratio||"");return null===a?(console.error("find-aspect-ratio error: Could not find a given aspect ratio."),e):"number"==typeof e.width&&"number"==typeof e.height?d.both(e,a):"number"==typeof e.width?d.height(e,a):"number"==typeof e.height?d.width(e,a):(console.error("find-aspect-ratio error: Please supply a height, width, or ratio value."),e)};"undefined"!=typeof self?self.nlp=t:"undefined"!=typeof window?window.nlp=t:void 0!==e&&(e.nlp=t),void 0!==n&&(n.exports=t)}).call(this,"undefined"!=typeof global?global:"undefined"!=typeof self?self:"undefined"!=typeof window?window:{})},{"./find-best-ratio":2,"./fit":3,"./parse-ratio":5}],5:[function(e,t,i){"use strict";var n=e("./aspects"),r=/^[0-9\.]+:[0-9\.]+$/;t.exports=function(e){if(e=(e=(e=(e=e.toLowerCase()).trim()).replace(" ratio","")).replace("-"," "),!0===n.lookup.hasOwnProperty(e))return n.lookup[e];if(!0!==r.test(e))return null;var t=e.split(":");return{description:"custom",decimal:parseFloat(t[0])/parseFloat(t[1])}}},{"./aspects":1}]},{},[4])(4)});

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{}],4:[function(_dereq_,module,exports){
!function(){var n=function(t,e,r,u){for(var o=1;o<e.length;o++){var f=e[o++],s="number"==typeof f?r[f]:f;1===e[o]?u[0]=s:2===e[o]?(u[1]=u[1]||{})[e[++o]]=s:3===e[o]?u[1]=Object.assign(u[1]||{},s):u.push(e[o]?t.apply(null,n(t,s,r,["",null])):s)}return u},t=function(n){for(var t,e,r=1,u="",o="",f=[0],s=function(n){1===r&&(n||(u=u.replace(/^\s*\n\s*|\s*\n\s*$/g,"")))?f.push(n||u,0):3===r&&(n||u)?(f.push(n||u,1),r=2):2===r&&"..."===u&&n?f.push(n,3):2===r&&u&&!n?f.push(!0,2,u):4===r&&e&&(f.push(n||u,2,e),e=""),u=""},p=0;p<n.length;p++){p&&(1===r&&s(),s(p));for(var h=0;h<n[p].length;h++)t=n[p][h],1===r?"<"===t?(s(),f=[f],r=3):u+=t:o?t===o?o="":u+=t:'"'===t||"'"===t?o=t:">"===t?(s(),r=1):r&&("="===t?(r=4,e=u,u=""):"/"===t?(s(),3===r&&(f=f[0]),r=f,(f=f[0]).push(r,4),r=0):" "===t||"\t"===t||"\n"===t||"\r"===t?(s(),r=2):u+=t)}return s(),f},e="function"==typeof Map,r=e?new Map:{},u=e?function(n){var e=r.get(n);return e||r.set(n,e=t(n)),e}:function(n){for(var e="",u=0;u<n.length;u++)e+=n[u].length+"-"+n[u];return r[e]||(r[e]=t(n))},o=function(t){var e=n(this,u(t),arguments,[]);return e.length>1?e:e[0]};"undefined"!=typeof module?module.exports=o:self.htm=o}();

},{}],5:[function(_dereq_,module,exports){
(function (global){
!function(e){if("object"==typeof exports&&"undefined"!=typeof module)module.exports=e();else if("function"==typeof define&&define.amd)define([],e);else{("undefined"!=typeof window?window:"undefined"!=typeof global?global:"undefined"!=typeof self?self:this).spencerColor=e()}}(function(){return function u(i,a,c){function f(r,e){if(!a[r]){if(!i[r]){var o="function"==typeof _dereq_&&_dereq_;if(!e&&o)return o(r,!0);if(d)return d(r,!0);var n=new Error("Cannot find module '"+r+"'");throw n.code="MODULE_NOT_FOUND",n}var t=a[r]={exports:{}};i[r][0].call(t.exports,function(e){return f(i[r][1][e]||e)},t,t.exports,u,i,a,c)}return a[r].exports}for(var d="function"==typeof _dereq_&&_dereq_,e=0;e<c.length;e++)f(c[e]);return f}({1:[function(e,r,o){"use strict";r.exports={blue:"#6699cc",green:"#6accb2",yellow:"#e1e6b3",red:"#cc7066",pink:"#F2C0BB",brown:"#705E5C",orange:"#cc8a66",purple:"#d8b3e6",navy:"#335799",olive:"#7f9c6c",fuscia:"#735873",beige:"#e6d7b3",slate:"#8C8C88",suede:"#9c896c",burnt:"#603a39",sea:"#50617A",sky:"#2D85A8",night:"#303b50",rouge:"#914045",grey:"#838B91",mud:"#C4ABAB",royal:"#275291",cherry:"#cc6966",tulip:"#e6b3bc",rose:"#D68881",fire:"#AB5850",greyblue:"#72697D",greygreen:"#8BA3A2",greypurple:"#978BA3",burn:"#6D5685",slategrey:"#bfb0b3",light:"#a3a5a5",lighter:"#d7d5d2",fudge:"#4d4d4d",lightgrey:"#949a9e",white:"#fbfbfb",dimgrey:"#606c74",softblack:"#463D4F",dark:"#443d3d",black:"#333333"}},{}],2:[function(e,r,o){"use strict";var n=e("./colors"),t={juno:["blue","mud","navy","slate","pink","burn"],barrow:["rouge","red","orange","burnt","brown","greygreen"],roma:["#8a849a","#b5b0bf","rose","lighter","greygreen","mud"],palmer:["red","navy","olive","pink","suede","sky"],mark:["#848f9a","#9aa4ac","slate","#b0b8bf","mud","grey"],salmon:["sky","sea","fuscia","slate","mud","fudge"],dupont:["green","brown","orange","red","olive","blue"],bloor:["night","navy","beige","rouge","mud","grey"],yukon:["mud","slate","brown","sky","beige","red"],david:["blue","green","yellow","red","pink","light"],neste:["mud","cherry","royal","rouge","greygreen","greypurple"],ken:["red","sky","#c67a53","greygreen","#dfb59f","mud"]};Object.keys(t).forEach(function(e){t[e]=t[e].map(function(e){return n[e]||e})}),r.exports=t},{"./colors":1}],3:[function(e,r,o){"use strict";var n=e("./colors"),t=e("./combos"),u={colors:n,list:Object.keys(n).map(function(e){return n[e]}),combos:t};r.exports=u},{"./colors":1,"./combos":2}]},{},[3])(3)});

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{}],6:[function(_dereq_,module,exports){
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

},{}],7:[function(_dereq_,module,exports){
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


},{}],8:[function(_dereq_,module,exports){
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
        style: 'overflow:hidden; margin: 10px 20px 25px 25px;' // border:1px solid lightgrey;

      };
      return h(_templateObject(), attrs, elements);
    }
  }]);

  return World;
}();

module.exports = World;

},{"./shapes/Dot":17,"./shapes/Line":18,"./shapes/Shape":19,"./shapes/Text":20,"d3-geo":2,"fit-aspect-ratio":3,"htm":4,"vhtml":7}],9:[function(_dereq_,module,exports){
"use strict";

module.exports = {
  shapes: {
    'great-lakes': _dereq_('./shapes/great-lakes'),
    'north-america': _dereq_('./shapes/north-america')
  },
  points: _dereq_('./points/index')
};

},{"./points/index":11,"./shapes/great-lakes":14,"./shapes/north-america":15}],10:[function(_dereq_,module,exports){
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

},{}],11:[function(_dereq_,module,exports){
"use strict";

module.exports = Object.assign({}, _dereq_('./cities'), _dereq_('./ontario'), _dereq_('./manitoba'));

},{"./cities":10,"./manitoba":12,"./ontario":13}],12:[function(_dereq_,module,exports){
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

},{}],13:[function(_dereq_,module,exports){
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

},{}],14:[function(_dereq_,module,exports){
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

},{}],15:[function(_dereq_,module,exports){
module.exports={"type":"Topology","arcs":[[[52674,16904],[-17,-20],[-7,6],[0,30],[-23,66],[0,18],[54,-65],[1,-17],[-8,-18]],[[52949,17057],[-11,0],[5,15],[22,28],[8,-11],[-1,-14],[-23,-18]],[[53053,17157],[-15,-5],[21,43],[6,-4],[-12,-34]],[[63296,41297],[-227,-14],[-181,61],[-135,30],[-131,54],[-285,174],[-33,61],[-27,74],[-55,67],[-59,56],[-302,170],[-24,59],[61,39],[70,18],[61,-3],[202,-66],[255,-58],[110,-44],[125,-67],[123,-82],[272,-216],[49,-17],[123,-106],[44,-80],[22,-65],[-27,-33],[-31,-12]],[[63675,38155],[18,-8],[74,48],[38,-2],[-1,-34],[-63,-38],[-29,-29],[36,-26],[0,-18],[-44,-43],[-22,-46],[18,-45],[70,44],[29,0],[39,-9],[37,13],[21,22],[123,173],[6,22],[-132,-36],[-16,24],[87,104],[-7,55],[44,87],[41,52],[29,29],[43,26],[29,-41],[10,-74],[73,10],[71,-15],[52,-31],[9,-18],[0,-28],[-17,-52],[-29,-42],[57,-54],[-7,-23],[-95,-61],[-53,-62],[-50,-76],[-98,-88],[-156,-63],[-49,0],[-59,19],[-58,-4],[-57,-24],[-56,3],[-27,-13],[-26,2],[-22,24],[-46,71],[-22,47],[-25,224],[8,119],[39,110],[58,73],[33,59],[141,347],[27,79],[34,67],[60,66],[78,113],[24,23],[45,11],[44,-7],[-13,-40],[4,-40],[51,-155],[-1,-31],[-28,-123],[-53,-199],[-14,-110],[7,-33],[-22,-56],[-23,-42],[-91,-80],[-47,-18],[-44,-30],[-105,-99]],[[6558,52564],[50,-50],[26,-3],[84,15],[31,-10],[30,-21],[18,-33],[4,-57],[-14,-52],[4,-72],[-3,-31],[44,-41],[14,-55],[8,-60],[-97,-20],[-96,-3],[-84,-40],[-17,-30],[14,-44],[-23,-11],[-21,9],[-41,41],[-44,20],[-155,30],[-194,119],[-84,24],[-85,87],[-77,111],[50,18],[51,9],[226,-17],[28,80],[29,20],[72,22],[67,44],[30,-1],[31,-17],[63,26],[34,6],[27,-13]],[[8004,47172],[53,-165],[23,-32],[34,-18],[47,-17],[29,-25],[24,-37],[4,-17],[-150,67],[-94,-98],[-28,-12],[-267,-5],[-53,-18],[-34,-31],[-61,-89],[-31,-35],[-33,-20],[-69,-23],[-83,3],[-44,12],[-23,42],[-20,83],[0,23],[8,40],[75,55],[24,29],[96,188],[28,26],[29,7],[82,-13],[70,54],[151,84],[33,10],[109,2],[29,-13],[23,-23],[19,-34]],[[54889,10724],[66,-36],[69,-18],[29,-1],[28,-12],[63,-86],[50,-47],[191,-105],[63,-181],[12,-57],[-49,-33],[-62,-12],[-59,-2],[-55,35],[-24,27],[-57,12],[14,25],[-25,11],[-32,-2],[-23,-70],[-27,-55],[-49,5],[-20,47],[-26,-21],[-21,-35],[-25,-130],[-41,65],[-44,54],[-56,22],[-113,4],[-53,17],[-44,110],[-17,32],[-45,28],[-44,126],[-16,18],[-120,26],[-24,69],[7,63],[40,76],[19,22],[67,-3],[63,23],[28,32],[29,22],[230,-55],[52,-1],[51,-9]],[[60943,10712],[17,-20],[16,3],[-13,41],[12,0],[104,-25],[67,-42],[69,-21],[4,-138],[-53,-56],[-35,-58],[-28,-71],[-74,-83],[-90,-25],[-59,-2],[-22,3],[-22,14],[-45,-14],[-56,37],[-47,-9],[-95,8],[-35,-31],[-34,-7],[-34,6],[-28,14],[-70,-1],[-30,27],[12,158],[1,71],[-17,60],[-19,37],[-13,44],[27,29],[23,42],[7,63],[25,15],[29,8],[134,-30],[340,-16],[19,-6],[13,-25]],[[13757,50009],[4,-54],[22,6],[83,57],[44,16],[57,3],[45,-26],[9,-20],[-4,-23],[-36,-48],[2,-31],[38,-57],[95,-31],[12,-17],[-1,-21],[-65,-95],[-24,-22],[-17,-5],[-119,17],[-109,31],[-44,5],[-16,-10],[-30,-29],[22,-8],[95,-7],[34,-43],[14,-31],[8,-34],[-19,-15],[-41,-10],[-49,0],[-61,-38],[-33,-45],[-121,-11],[-92,-59],[-33,-29],[-11,-35],[-35,-25],[-78,-26],[47,-22],[7,-20],[1,-26],[-7,-23],[-62,-102],[-117,-84],[-30,4],[-12,10],[-11,16],[0,16],[151,169],[-6,8],[-41,6],[-65,46],[-46,-30],[-11,1],[13,40],[30,47],[-6,14],[-16,12],[-38,7],[-59,2],[-44,-11],[-29,-25],[-2,-11],[61,3],[16,-12],[17,-24],[10,-27],[3,-31],[-13,-41],[-29,-52],[-42,9],[-87,116],[-38,169],[-75,129],[-3,31],[21,79],[75,113],[81,31],[57,48],[56,14],[34,-1],[47,-20],[19,-44],[-11,-20],[4,-11],[33,-26],[36,-91],[41,-80],[28,-33],[36,-20],[-37,61],[-23,76],[-11,153],[-11,40],[21,10],[60,-5],[-2,22],[-63,51],[-38,43],[-15,33],[1,29],[35,43],[19,12],[20,6],[40,-9],[18,-14],[51,-96],[24,-30],[20,0],[18,16],[17,32],[18,19],[18,6],[57,-14],[19,5],[9,24],[0,44],[14,15],[3,32],[-31,47],[36,14],[118,-36],[51,-39],[-26,-73]],[[14018,50545],[20,-9],[20,60],[15,2],[64,-51],[40,12],[25,-62],[23,-6],[21,8],[13,-5],[-5,-65],[-46,-66],[-22,-17],[-31,17],[-11,6],[-21,30],[-15,37],[-9,0],[-36,-43],[1,-22],[15,-32],[-2,-20],[-39,-10],[-39,5],[-46,-28],[-12,17],[-7,50],[-14,-6],[-22,-60],[-24,-39],[-44,-32],[-10,-15],[-32,-2],[-48,-20],[-29,3],[-176,66],[-41,24],[144,151],[76,58],[44,-3],[44,-18],[23,3],[2,67],[-41,50],[2,21],[90,34],[34,-5],[37,-18],[36,-28],[33,-39]],[[25676,47679],[-18,-109],[-38,-112],[-57,-61],[-40,14],[-30,48],[-27,-2],[-30,10],[-16,40],[15,52],[-13,40],[-15,-35],[-26,-32],[-65,-42],[-44,-81],[-22,-52],[-26,57],[-17,135],[-2,58],[47,86],[61,83],[12,245],[194,124],[18,-7],[62,-92],[67,-129],[18,-58],[0,-101],[-8,-79]],[[24270,48528],[103,-22],[94,2],[32,-41],[22,-43],[12,-41],[3,-38],[-2,-27],[-12,-29],[3,-10],[184,-95],[87,-100],[35,-52],[20,-44],[36,-111],[77,-129],[40,-39],[23,-38],[-12,-2],[-55,29],[-117,86],[-9,-3],[-11,-47],[-18,-41],[-27,-30],[22,-8],[94,18],[79,-84],[31,-15],[30,-60],[1,-24],[-17,-44],[-14,-18],[5,-13],[22,-6],[88,12],[15,-22],[-13,-174],[12,-64],[0,-29],[-10,-39],[0,-33],[8,-33],[1,-30],[-23,-79],[-23,-13],[-38,-1],[-30,23],[-42,67],[-40,105],[-16,15],[-54,15],[-9,13],[-34,2],[-25,43],[3,57],[-21,57],[3,26],[-24,11],[-19,-16],[10,-57],[-12,-44],[-42,19],[-71,139],[-80,113],[-32,26],[8,33],[40,17],[32,-1],[7,19],[-67,108],[2,31],[24,55],[-30,23],[-84,-18],[-30,12],[-24,44],[-14,38],[-73,7],[-27,-5],[-48,58],[-21,36],[8,18],[44,33],[25,-5],[49,-34],[20,1],[48,46],[8,42],[36,34],[-6,36],[-20,62],[-44,17],[-91,-37],[-80,-56],[-31,22],[-7,35],[85,94],[37,52],[-7,30],[-28,40],[-2,100],[18,23]],[[24378,49191],[36,-32],[20,33],[37,-1],[68,-30],[41,-43],[22,-50],[2,-30],[-7,-68],[5,-69],[-2,-36],[-9,-30],[-16,-23],[-17,-3],[-53,62],[-60,112],[-46,35],[-2,-12],[12,-32],[38,-61],[7,-36],[26,-44],[12,-34],[7,-44],[0,-39],[-8,-33],[-12,-21],[-17,-9],[-93,9],[-55,-22],[-64,12],[-16,20],[-10,33],[-5,79],[-16,114],[4,87],[-41,79],[-36,48],[-51,42],[-34,43],[9,34],[53,25],[86,-6],[185,-59]],[[23506,49538],[46,-110],[34,-85],[30,-102],[50,-211],[23,-80],[7,-44],[6,-115],[-8,-24],[-15,-23],[-3,-33],[14,-87],[1,-133],[-13,-75],[-15,-11],[-37,24],[-31,41],[-23,42],[-55,133],[-17,62],[-1,45],[9,32],[18,20],[32,54],[-5,9],[-24,-12],[-49,-7],[-43,43],[-34,22],[7,77],[-9,22],[-66,-24],[-25,21],[-5,29],[1,43],[12,38],[63,95],[-6,18],[-31,4],[-40,33],[-18,106],[-44,61],[-25,-6],[-58,-172],[-29,-37],[-82,-24],[17,48],[7,42],[-29,130],[-1,50],[20,37],[58,15],[30,22],[24,35],[6,35],[45,92],[21,17],[56,-1],[117,-101],[35,-15],[52,-65]],[[23664,50346],[138,-22],[101,5],[92,-149],[57,-121],[34,-84],[19,-82],[25,-78],[-2,-12],[-54,54],[-38,107],[-20,42],[-19,19],[-20,40],[-40,102],[-1,29],[-18,27],[-20,11],[-23,-4],[-8,-11],[3,-70],[18,-79],[100,-172],[67,-98],[13,-32],[9,-90],[-29,-40],[35,-83],[-1,-16],[-8,-16],[-95,-36],[-87,-153],[-95,-90],[-44,-14],[-21,15],[-20,34],[-12,46],[-3,57],[23,36],[48,188],[1,62],[-59,86],[-35,70],[-19,98],[-32,258],[-15,83],[-21,69],[-27,55],[-21,60],[-14,66],[5,26],[47,-34],[57,-95],[29,-64]],[[23093,50429],[77,-97],[1,-23],[-16,-67],[-42,-18],[12,-27],[31,-20],[23,17],[82,93],[25,20],[15,2],[100,-28],[88,-45],[25,-36],[15,-62],[-23,-135],[-72,-24],[-34,3],[-36,19],[-58,-47],[48,-36],[147,-8],[45,-76],[13,-59],[-32,-108],[-83,30],[-74,63],[-151,88],[-36,4],[-24,-15],[-7,-54],[2,-116],[-40,-61],[-119,27],[-47,88],[-44,139],[-164,165],[-44,33],[-59,99],[23,79],[7,45],[32,12],[46,35],[26,76],[41,-62],[56,-59],[1,56],[26,44],[54,-2],[26,9],[35,42],[52,21],[31,-24]],[[26996,45135],[25,-318],[-2,-102],[-36,-67],[-26,-111],[-30,-50],[-28,66],[-2,111],[-8,88],[-9,43],[11,164],[-13,-12],[-33,-72],[-38,-5],[-67,82],[-33,66],[-6,69],[-44,73],[-5,26],[4,27],[36,73],[15,49],[13,101],[15,39],[34,-6],[61,-45],[64,-50],[58,-66],[44,-173]],[[24765,46321],[50,-59],[119,37],[22,-7],[23,-22],[25,-48],[27,-73],[6,-80],[-11,-29],[-21,-32],[-191,-126],[-5,-13],[4,-11],[18,-12],[38,1],[151,32],[8,22],[10,100],[21,53],[2,40],[-14,95],[1,39],[105,7],[66,36],[68,65],[15,-2],[-10,-118],[-9,-36],[-64,-145],[-37,-127],[-18,-126],[-4,-208],[-16,-71],[-29,-43],[-183,-76],[-94,5],[-84,69],[-39,47],[30,58],[20,2],[59,-11],[46,-21],[20,-1],[-3,13],[-130,108],[-95,49],[-29,55],[-1,43],[-7,23],[-76,149],[-15,59],[-10,83],[0,85],[19,138],[8,15],[31,-1],[54,-17],[129,-13]],[[25255,45391],[55,-93],[17,-83],[-7,-97],[-87,-37],[-46,29],[-20,-5],[-30,-30],[36,-13],[52,-48],[45,-62],[63,-11],[85,-43],[-64,-78],[-10,-45],[80,-125],[7,-31],[26,-7],[59,10],[8,-9],[0,-26],[-37,-73],[4,-14],[33,-11],[64,-1],[14,-72],[-57,-65],[-109,84],[-49,84],[-28,78],[-33,44],[-101,97],[-154,211],[-40,30],[-39,83],[-12,41],[1,26],[15,15],[46,9],[1,42],[-177,72],[-19,15],[-23,51],[12,7],[97,-9],[105,25],[64,18],[25,23],[52,29],[22,-1],[54,-34]],[[62203,38678],[15,-14],[26,25],[30,82],[80,-22],[42,-36],[24,8],[24,-4],[45,-48],[85,-38],[90,6],[137,22],[16,9],[141,19],[141,9],[49,-21],[18,-20],[9,-24],[-80,-66],[-81,-77],[-112,-76],[-14,-37],[7,-67],[-2,-70],[22,-5],[13,-23],[-29,-23],[-115,-10],[-33,6],[-40,28],[-14,67],[-50,-10],[-14,8],[69,57],[-32,72],[-34,-6],[-22,34],[2,46],[31,22],[9,25],[-42,-21],[-33,-43],[-42,-16],[-42,-37],[64,-11],[-33,-29],[-34,-6],[-159,56],[-39,21],[-49,58],[-38,78],[20,3],[7,14],[-4,14],[-55,10],[-87,-3],[-50,20],[3,137],[-16,37],[-54,32],[-83,9],[-9,51],[27,78],[40,66],[31,65],[36,54],[89,106],[-2,-79],[9,-69],[-58,-137],[100,-136],[13,-30],[9,-37],[-7,-34],[-16,-29],[39,-15],[12,-25]],[[57473,33206],[-39,-64],[35,-7],[30,19],[29,38],[66,52],[56,23],[18,5],[27,-37],[53,29],[55,17],[-237,-166],[-49,-19],[-69,-50],[-65,-35],[-48,-12],[-234,-124],[-19,-3],[-20,12],[-193,-63],[-80,-7],[-72,-22],[54,51],[1,19],[-13,15],[-29,-4],[-29,-53],[-47,-18],[-9,58],[16,44],[21,42],[46,66],[67,42],[33,36],[24,-31],[5,43],[18,25],[19,13],[47,0],[26,7],[18,14],[19,3],[51,-19],[50,5],[41,27],[42,9],[111,6],[111,20],[45,35],[93,98],[53,28],[-83,-114],[-45,-53]],[[52400,15420],[143,-46],[115,13],[55,27],[-5,-28],[51,-69],[18,-5],[75,35],[194,13],[20,-19],[34,-67],[50,-42],[51,-31],[54,-8],[53,14],[51,-7],[62,-65],[20,-8],[56,18],[-16,-60],[94,-84],[70,-165],[50,-68],[54,-61],[45,-41],[50,-19],[153,8],[36,-5],[32,-24],[31,-9],[18,9],[295,-257],[94,-137],[58,-72],[124,-103],[50,-22],[26,13],[-5,23],[-37,57],[-5,21],[47,-18],[84,-116],[23,-43],[42,-39],[43,-29],[-20,-46],[-35,-4],[-66,19],[52,-75],[10,-54],[24,-5],[36,60],[23,50],[93,-129],[50,-60],[-13,-34],[-4,-35],[56,32],[21,-3],[20,-19],[23,-56],[52,-12],[52,2],[107,-47],[101,-93],[95,-19],[96,-4],[48,-49],[21,-67],[-24,-47],[-13,-49],[36,-60],[-78,-26],[-11,-36],[4,-40],[16,-21],[44,19],[64,-17],[102,-15],[68,12],[138,-41],[42,-22],[81,-77],[38,-51],[81,-138],[71,-54],[60,-13],[21,9],[20,-15],[16,-19],[17,-60],[-9,-63],[-35,-51],[-19,-38],[-87,-4],[-121,-17],[-117,-56],[-57,-44],[-26,-30],[-62,-27],[-4,23],[1,30],[-16,54],[-14,-49],[-23,-36],[-38,-30],[-142,-2],[-58,41],[-58,28],[-214,29],[-52,-2],[-142,-31],[-144,-16],[-60,-19],[-60,-28],[-115,1],[-137,-33],[-137,-6],[88,227],[185,217],[35,47],[25,60],[6,46],[-9,38],[-44,68],[-9,51],[-13,33],[-64,29],[-65,17],[-68,0],[-144,24],[-76,2],[-65,46],[-107,165],[-51,46],[-26,38],[-20,42],[-25,243],[-21,117],[-33,101],[-49,77],[-52,26],[-199,-66],[-47,10],[-45,22],[-301,158],[-124,86],[-50,43],[-43,61],[-45,100],[-50,90],[0,-37],[-8,-23],[-251,-11],[-41,21],[-25,24],[-19,36],[-13,73],[-24,61],[-8,-65],[-12,-60],[-34,-33],[-38,-6],[-47,80],[-204,16],[-18,14],[-67,77],[-57,96],[57,34],[117,45],[26,30],[14,38],[-10,57],[-24,41],[-24,24],[-26,16],[-35,6],[-454,10],[-27,-31],[-40,-63],[-81,-81],[-53,-84],[-20,-43],[-25,-31],[-56,-52],[-47,-80],[-58,-36],[-32,22],[-31,0],[-23,-20],[-23,-9],[-117,-10],[-17,-20],[-17,-58],[-19,-111],[-17,-37],[-59,-14],[-55,-31],[-114,-106],[-29,-16],[6,78],[-5,76],[-32,3],[-37,-13],[-30,-21],[-56,-57],[-28,-14],[-27,29],[5,37],[188,137],[21,10],[33,-10],[32,4],[26,39],[-31,181],[12,123],[43,95],[87,144],[42,47],[428,301],[44,15],[278,61],[42,21],[129,89],[136,36],[143,-27]],[[27733,42841],[270,-125],[270,-61],[199,-73],[121,-23],[44,-16],[29,-25],[33,-62],[58,-148],[44,-94],[91,-164],[72,-116],[16,-46],[-15,-15],[1,-27],[54,-113],[102,-102],[80,-48],[169,-79],[103,-77],[32,-53],[45,-51],[19,-36],[37,-132],[68,-127],[70,-241],[13,20],[9,72],[8,16],[15,7],[14,-27],[12,-64],[45,-151],[-14,-44],[-13,-5],[-61,21],[-21,-27],[-28,-55],[-20,-22],[-12,11],[-175,53],[-108,49],[-142,79],[-170,82],[-97,56],[-81,58],[-57,49],[-10,42],[2,19],[109,134],[47,72],[17,55],[10,58],[-7,71],[-5,-6],[-9,-68],[-16,-60],[-20,-47],[-13,-17],[-130,-23],[-106,7],[-53,-57],[-16,-7],[-29,19],[-64,76],[-92,62],[9,16],[60,32],[32,46],[-6,7],[-21,-2],[-19,9],[-37,60],[-21,17],[-45,-27],[-19,-2],[-17,40],[25,92],[1,21],[-46,-33],[-15,11],[-15,30],[-13,12],[-38,-6],[-40,27],[-15,-10],[-5,-40],[-13,-10],[-62,67],[-15,2],[-30,-51],[-10,-3],[-16,22],[-8,124],[3,35],[8,12],[54,29],[155,30],[13,23],[-116,-12],[-30,17],[-33,42],[-33,0],[-18,14],[-19,30],[-49,112],[-34,29],[-57,18],[-29,21],[-12,-10],[-12,-32],[-16,-19],[-39,-12],[-36,9],[-28,30],[-16,39],[-7,43],[16,58],[0,23],[-7,26],[-13,22],[-18,16],[-11,-8],[-1,-34],[-10,-25],[-33,-19],[-26,33],[-17,46],[-21,33],[-113,0],[-52,-42],[-25,-4],[-25,10],[-5,22],[24,62],[-6,81],[-6,21],[-52,12],[-9,20],[31,99],[17,19],[23,7],[104,8],[34,-14],[50,-61],[-2,23],[-19,68],[-2,41],[34,47],[-33,13],[-122,11],[1,-30],[10,-42],[-73,-37],[-54,-6],[-51,6],[-42,22],[-72,89],[-45,87],[2,47],[25,50],[32,34],[76,30],[100,2],[112,-39],[280,-180]],[[57775,11986],[65,-25],[30,21],[24,17],[16,60],[21,54],[28,28],[32,18],[63,-1],[87,-47],[25,1],[83,42],[70,24],[65,-27],[27,-36],[54,-58],[27,-17],[85,1],[23,-6],[72,-95],[60,-38],[35,-2],[63,37],[31,-1],[36,-82],[7,-116],[30,-106],[46,-68],[225,29],[50,-56],[-17,-46],[-32,-25],[-107,11],[-47,-5],[-9,-46],[-1,-43],[63,-10],[62,-21],[62,-35],[64,-23],[72,-15],[70,-25],[118,-83],[130,-190],[35,-44],[23,-59],[-11,-73],[-46,-120],[-27,-39],[-38,-24],[-26,-49],[-26,-84],[-15,-7],[-19,4],[-31,47],[-22,73],[-63,69],[-75,-9],[-110,41],[-66,-20],[-67,-4],[-68,20],[-68,7],[-69,-25],[-66,-44],[-25,-28],[-42,-69],[-23,-25],[-161,-34],[-47,50],[-43,68],[-62,10],[-90,-53],[-56,-20],[-23,-22],[-7,-26],[0,-96],[-13,-58],[-87,-220],[-50,-155],[-20,-48],[-24,-11],[-43,89],[-27,33],[-35,16],[-14,47],[1,68],[-9,65],[-21,51],[-31,34],[-46,80],[-51,66],[-30,26],[-31,17],[-242,-9],[-27,-11],[-21,-22],[-22,-10],[-67,-20],[-66,-5],[-154,54],[-61,28],[-61,17],[-71,-5],[-70,-17],[-56,-38],[-42,-69],[-8,-63],[-25,-16],[-57,101],[-52,71],[-59,54],[-122,77],[-23,47],[-9,57],[49,174],[56,32],[31,6],[69,-21],[68,-40],[61,-26],[96,-10],[52,-43],[366,-66],[70,-21],[27,7],[24,26],[19,47],[23,35],[109,8],[23,16],[16,49],[-1,51],[-64,69],[-100,150],[-88,177],[38,60],[-14,109],[14,101],[21,99],[-87,85],[-103,84],[-143,27],[-44,21],[-23,63],[21,85],[46,48],[53,29],[54,20],[132,24],[130,-27],[112,-88],[115,-68],[144,-23]],[[61641,10597],[-41,-9],[-56,47],[44,17],[29,-11],[24,-44]],[[31480,26259],[-29,-2],[-44,11],[15,23],[25,17],[8,-13],[25,-36]],[[31953,25453],[-25,-1],[-33,9],[-17,53],[27,4],[25,-7],[20,-42],[3,-16]],[[57313,11044],[-10,-71],[-139,84],[-113,106],[5,57],[58,13],[55,-35],[81,-71],[63,-83]],[[64391,36154],[-73,-36],[-63,3],[-41,32],[-2,14],[99,-13],[38,7],[75,56],[-33,-63]],[[44147,20149],[-12,-3],[27,82],[4,33],[44,103],[22,16],[10,-44],[-44,-73],[-51,-114]],[[48566,21956],[-30,-52],[1,20],[22,52],[16,20],[-9,-40]],[[48746,22458],[-7,-11],[-48,21],[-30,19],[-5,20],[81,-35],[9,-14]],[[49010,22495],[-48,-22],[-70,2],[-15,8],[29,14],[84,19],[20,-21]],[[62638,10313],[-6,-26],[-37,22],[-18,23],[5,26],[28,25],[21,-2],[8,-9],[-1,-59]],[[61325,10373],[-42,-6],[-28,8],[-9,30],[51,27],[61,-4],[34,-16],[5,-11],[-72,-28]],[[62899,9391],[-27,-21],[-24,29],[6,70],[22,1],[22,-30],[1,-49]],[[61886,10731],[-14,-7],[-10,2],[-3,14],[10,40],[55,4],[-38,-53]],[[61742,10621],[-36,-26],[-24,4],[-10,9],[20,31],[50,-18]],[[61778,10670],[-43,-4],[-12,13],[24,30],[44,4],[13,-8],[-26,-35]],[[61684,10063],[46,-45],[55,0],[-58,-44],[-110,-4],[3,70],[64,23]],[[54842,17267],[-61,-21],[-44,21],[-11,16],[18,28],[42,23],[66,2],[29,-27],[4,-12],[-43,-30]],[[52832,16971],[-25,-24],[-26,17],[28,24],[85,25],[-32,-31],[-30,-11]],[[55278,26888],[-12,-2],[-22,11],[-29,21],[-7,15],[29,-5],[41,-40]],[[30096,39412],[-6,-20],[-7,2],[-17,39],[-2,29],[14,20],[20,-58],[-2,-12]],[[49908,8649],[-87,-77],[-27,1],[40,60],[64,52],[55,25],[45,-11],[-90,-50]],[[56898,37680],[-42,-20],[-71,19],[-79,-26],[-22,-1],[59,75],[90,45],[89,140],[25,3],[-34,-158],[-7,-56],[-8,-21]],[[56828,37796],[-65,-20],[-23,8],[72,98],[83,23],[-67,-109]],[[63753,37701],[-22,-8],[-6,14],[-34,41],[-3,21],[30,19],[62,-10],[-22,-52],[-5,-25]],[[63235,39492],[19,-19],[35,2],[23,-8],[-33,-37],[-64,-4],[-31,16],[45,190],[53,44],[109,124],[43,38],[40,15],[43,-9],[-43,-73],[-58,-4],[-56,-61],[-36,-68],[-44,-39],[-30,-47],[-15,-60]],[[61824,40093],[-14,-73],[-47,-62],[-24,-4],[-11,6],[14,40],[-2,69],[39,9],[15,-6],[30,21]],[[61842,40165],[-63,-52],[28,78],[11,20],[11,8],[10,-6],[3,-48]],[[58280,39081],[-49,-8],[11,35],[68,61],[49,35],[30,-1],[-48,-75],[-61,-47]],[[60864,36506],[-27,-35],[6,34],[34,87],[22,13],[-35,-99]],[[59587,36410],[-21,-20],[-22,6],[-1,50],[7,16],[7,8],[12,-14],[18,-46]],[[60598,36894],[-73,-52],[29,134],[23,42],[31,-14],[-5,-82],[-5,-28]],[[59824,36546],[-32,-20],[-35,9],[1,-53],[-4,-19],[-38,27],[-14,18],[1,70],[34,65],[26,27],[33,-19],[26,-73],[2,-32]],[[4369,49371],[-56,-47],[-51,17],[-16,34],[0,15],[147,38],[-24,-57]],[[9715,47528],[-7,-78],[-75,50],[-20,24],[13,19],[71,5],[18,-20]],[[9522,47505],[8,0],[17,8],[35,56],[12,2],[0,-17],[-17,-56],[29,-74],[26,-36],[-3,-12],[-66,-26],[-50,19],[-27,-7],[-24,-26],[-16,29],[-12,137],[4,24],[27,48],[36,21],[15,-6],[13,-19],[3,-20],[-10,-45]],[[8645,47039],[-13,-5],[-38,31],[-13,22],[-10,36],[76,51],[16,1],[14,-29],[3,-23],[-21,-64],[-14,-20]],[[8505,46594],[-47,-22],[-50,23],[-43,42],[-4,50],[95,-33],[20,-15],[29,-45]],[[10241,47164],[-17,-5],[-15,10],[-20,57],[-2,24],[39,-18],[15,-41],[0,-27]],[[9964,47320],[-34,-22],[-10,-28],[-25,-12],[-21,-22],[-72,-103],[-31,-19],[35,88],[5,28],[1,18],[-11,66],[21,-3],[17,14],[34,58],[31,6],[34,69],[18,6],[9,-10],[-15,-42],[32,-39],[-8,-37],[-10,-16]],[[10158,47343],[-3,-80],[-8,-12],[-14,21],[-31,-24],[-16,18],[7,27],[-3,22],[23,2],[5,39],[-3,17],[11,36],[16,8],[16,-74]],[[4624,45043],[-9,-40],[-8,-15],[-84,22],[-56,-8],[-6,24],[5,21],[88,32],[36,0],[23,-17],[11,-19]],[[4057,44778],[-34,-32],[-10,12],[-6,39],[19,30],[59,67],[41,-13],[12,-17],[-1,-25],[-15,-33],[-19,-17],[-23,0],[-23,-11]],[[6718,46264],[-21,-18],[-16,-4],[-13,10],[-56,-12],[-11,8],[-25,59],[-2,31],[10,24],[25,22],[41,20],[40,-4],[70,-55],[32,-32],[7,-20],[-25,-19],[-56,-10]],[[6870,46330],[-24,-8],[-5,11],[-3,44],[-18,69],[34,25],[22,6],[9,-10],[25,-52],[25,-14],[19,-11],[-33,-16],[-51,-44]],[[5563,45540],[-166,-107],[-54,-78],[-42,-75],[-32,-42],[-24,-7],[-26,-20],[-55,-52],[-23,-7],[-177,-122],[-12,-2],[8,32],[54,45],[36,42],[39,68],[22,25],[7,35],[3,69],[10,26],[38,52],[27,29],[37,9],[73,-8],[31,26],[10,20],[-19,19],[-4,31],[4,55],[22,49],[37,43],[52,33],[65,23],[48,1],[85,-50],[13,-23],[-21,-50],[-12,-47],[-54,-42]],[[6297,46094],[23,-47],[41,30],[30,41],[23,54],[15,20],[20,-28],[57,-38],[-48,-59],[-92,-88],[-31,-59],[-2,-25],[89,19],[25,-2],[16,-21],[-24,-23],[-49,-22],[-42,-42],[-98,-73],[-37,-61],[-44,-23],[-60,-5],[-105,-40],[-64,-37],[-16,-20],[-20,-9],[-24,1],[-25,-17],[-28,-33],[-23,-16],[-38,-3],[-20,-14],[-22,1],[-61,40],[-15,24],[54,47],[38,15],[59,7],[58,44],[120,57],[37,32],[24,109],[27,19],[14,43],[66,-1],[31,-49],[10,-7],[6,5],[3,39],[33,28],[-19,20],[-61,25],[-46,12],[-29,-1],[-24,14],[-19,30],[-8,30],[1,29],[16,33],[28,37],[33,20],[70,15],[62,25],[33,3],[25,-10],[7,-95]],[[6517,45917],[-7,-3],[-14,25],[-1,23],[9,16],[25,38],[19,14],[22,7],[6,-10],[-20,-46],[-24,-30],[-15,-34]],[[3116,44469],[-40,-15],[-44,16],[20,52],[21,29],[40,34],[47,-15],[38,-44],[-82,-57]],[[2524,44333],[107,-40],[132,5],[49,-10],[0,-12],[-84,-16],[-29,5],[-74,-22],[-51,-4],[-115,21],[-89,-14],[-23,5],[-28,18],[-32,32],[-2,20],[29,8],[79,-27],[8,14],[67,25],[56,-8]],[[1913,44232],[-292,-41],[-44,29],[44,16],[52,8],[110,47],[135,41],[105,49],[92,32],[25,53],[-80,27],[-16,21],[38,25],[32,36],[75,42],[67,-53],[15,-35],[-6,-42],[-14,-44],[-59,-23],[-9,-23],[32,-65],[-121,-57],[-181,-43]],[[1189,44010],[-46,-22],[-61,45],[6,51],[66,-43],[35,-31]],[[1037,43990],[-34,-59],[-25,27],[-10,80],[20,21],[53,-58],[-4,-11]],[[871,44064],[3,-33],[62,4],[20,-17],[0,-66],[-9,-18],[-9,-5],[-22,15],[-26,-34],[-116,-82],[-36,46],[-68,-72],[48,186],[55,29],[20,21],[-5,54],[26,92],[55,-5],[25,-37],[-1,-25],[-22,-53]],[[1182,44200],[-13,-29],[-53,31],[-19,25],[-4,26],[16,44],[42,0],[25,-17],[24,-33],[7,-21],[-25,-26]],[[569,43915],[-16,-13],[-28,-10],[-83,11],[-51,-4],[-55,-7],[-43,-17],[-8,24],[2,20],[181,55],[42,29],[26,36],[23,68],[19,20],[12,-1],[25,-26],[-9,-36],[-22,-30],[-8,-29],[-7,-90]],[[172,43848],[-13,-33],[-13,1],[-72,55],[-11,18],[43,27],[12,20],[-5,26],[-32,37],[-59,47],[-22,34],[14,21],[28,13],[89,2],[48,-58],[36,-20],[84,-14],[-44,-24],[-24,-24],[-31,-92],[-28,-36]],[[28334,41557],[15,-90],[-66,14],[-24,18],[-1,43],[13,41],[50,-15],[13,-11]],[[30345,39602],[-2,-22],[-21,-18],[-11,4],[-1,27],[-5,4],[-23,-32],[2,64],[11,67],[10,1],[14,-44],[26,-51]],[[30289,39801],[-3,-19],[-30,23],[-10,21],[1,47],[8,31],[5,6],[18,-13],[5,-8],[6,-88]],[[30134,40877],[7,-22],[-22,-21],[-15,-3],[-25,34],[-11,4],[9,-52],[-3,-18],[-51,32],[-9,26],[14,26],[32,28],[11,4],[63,-38]],[[30113,40636],[-9,-10],[-29,14],[-17,22],[-6,28],[11,53],[15,13],[9,-3],[4,-46],[25,-50],[-3,-21]],[[30008,40706],[15,-33],[-59,21],[-24,19],[-8,18],[-9,60],[4,20],[27,7],[48,-75],[6,-37]],[[30248,40362],[27,-131],[11,55],[74,-95],[0,-46],[-9,-16],[-15,-5],[-15,13],[-13,33],[-16,17],[-36,11],[-18,37],[-7,26],[-1,73],[-9,23],[-19,5],[-18,17],[-28,52],[-4,13],[13,42],[31,70],[22,33],[13,-3],[17,-21],[20,-39],[-4,-27],[-81,-53],[-3,-12],[40,-15],[14,-13],[14,-44]],[[29388,41734],[8,-21],[-121,78],[-52,46],[-20,33],[-13,19],[-61,49],[-11,22],[15,17],[41,-10],[69,-37],[61,-60],[84,-136]],[[29813,41090],[-7,-11],[-85,71],[-56,92],[-24,57],[112,-141],[57,-46],[3,-22]],[[29779,40959],[-23,-26],[-12,3],[-10,18],[-35,176],[15,-4],[47,-55],[-11,-21],[36,-53],[8,-37],[-15,-1]],[[28940,42231],[-13,-9],[-13,113],[16,41],[5,21],[-2,21],[29,-51],[11,-34],[5,-48],[0,-12],[-38,-42]],[[28035,41808],[-21,-4],[-34,12],[-39,28],[-67,77],[-7,17],[6,15],[19,12],[4,20],[-15,55],[54,34],[48,-28],[22,-35],[26,-63],[12,-71],[1,-48],[-9,-21]],[[28828,42299],[-6,-53],[-35,85],[-54,181],[7,43],[24,60],[22,4],[35,-28],[31,-50],[7,-20],[19,-52],[10,-47],[-21,-57],[-39,-66]],[[32544,25065],[-32,-9],[-35,20],[-30,96],[-33,76],[18,21],[27,-72],[67,-110],[18,-22]],[[32603,21187],[-24,-38],[-62,208],[-1,50],[19,25],[30,-5],[0,-52],[25,-44],[10,-43],[3,-101]],[[31623,26154],[-38,-14],[-29,13],[-46,96],[98,12],[42,-41],[6,-12],[-33,-54]],[[31711,26314],[111,-51],[59,25],[11,-25],[-7,-22],[-134,-38],[-42,27],[-3,37],[-14,35],[19,12]],[[32545,25622],[28,-74],[-40,9],[-41,-4],[-13,40],[-12,56],[-8,14],[-29,5],[-2,5],[-3,27],[8,13],[89,-62],[23,-29]],[[58095,33704],[-15,-18],[-5,48],[10,56],[11,1],[5,-30],[-6,-57]],[[56590,31904],[-64,-151],[-1,29],[79,187],[-14,-65]],[[56560,32744],[-26,-4],[26,96],[48,43],[17,-9],[1,-34],[-7,-28],[-32,-45],[-27,-19]],[[58561,33596],[-150,-49],[-24,31],[37,15],[47,75],[32,8],[49,-42],[9,-38]],[[58163,33711],[-27,-27],[-30,5],[15,36],[6,55],[15,59],[8,18],[18,16],[-5,-162]],[[58850,33485],[-42,-16],[-97,37],[80,31],[13,11],[10,46],[1,23],[31,-99],[4,-33]],[[55694,27422],[-99,-71],[-11,4],[64,51],[46,16]],[[55773,28087],[-8,-21],[-36,111],[37,-36],[6,-30],[1,-24]],[[55938,30115],[-25,-16],[83,200],[48,167],[21,58],[-20,-117],[-37,-108],[-70,-184]],[[55823,27472],[-73,-27],[-7,9],[84,57],[26,200],[4,93],[-13,162],[1,34],[14,-52],[12,-152],[-5,-116],[-24,-168],[-19,-40]],[[55301,26876],[-14,-11],[50,124],[99,158],[26,25],[-82,-135],[-79,-161]],[[44283,20399],[-20,-4],[25,54],[40,27],[88,104],[36,7],[19,36],[8,5],[-5,-44],[-71,-62],[-120,-123]],[[43962,19549],[-16,-58],[4,86],[44,194],[90,256],[38,42],[-104,-280],[-56,-240]],[[45221,21391],[-27,-10],[118,154],[25,51],[32,-3],[-53,-86],[-95,-106]],[[44062,18410],[-8,-46],[-45,216],[-73,490],[-3,279],[12,96],[19,-394],[81,-501],[17,-140]],[[50731,21886],[-55,-36],[-59,27],[37,4],[26,-10],[67,51],[36,39],[41,15],[-93,-90]],[[52628,23212],[-24,-243],[-11,86],[-1,84],[19,49],[17,24]],[[48384,22327],[1,-47],[-26,24],[-40,2],[18,15],[12,16],[6,17],[50,57],[-14,-42],[-7,-42]],[[46986,21745],[-20,-15],[-90,87],[-6,37],[45,34],[27,-3],[43,-44],[15,-12],[7,-18],[-4,-28],[-17,-38]],[[48599,22051],[-15,-32],[15,152],[-21,129],[22,-57],[7,-67],[-8,-125]],[[54673,16504],[1,-23],[-54,-63],[40,-44],[36,97],[30,-79],[16,-149],[-3,-26],[3,-22],[6,-29],[1,-41],[-30,-130],[-107,13],[-3,110],[-16,21],[-25,157],[-34,50],[-47,129],[27,33],[37,-11],[18,15],[50,13],[31,17],[23,-38]],[[54667,14212],[-22,-31],[-24,44],[-16,5],[-22,17],[-41,50],[-10,50],[34,4],[44,-9],[76,-29],[-6,-57],[-13,-44]],[[54473,14544],[-12,-17],[-29,38],[-43,16],[-25,58],[-24,22],[-1,21],[39,14],[27,-6],[31,-44],[18,-81],[19,-21]],[[54145,14811],[74,-22],[26,13],[25,4],[26,-9],[37,-83],[-32,-10],[-26,0],[-19,15],[-67,4],[-44,24],[-24,20],[-12,25],[36,19]],[[54908,18156],[-11,-9],[-47,100],[-38,29],[59,71],[26,60],[-1,133],[15,72],[-5,63],[13,64],[-18,72],[-49,57],[-99,227],[-155,56],[-80,2],[44,36],[41,-4],[63,-22],[75,-10],[46,-67],[44,-89],[39,-36],[16,-22],[-1,-26],[6,-24],[52,-42],[52,-67],[15,-196],[-70,-93],[-13,-284],[-19,-51]],[[56631,14981],[13,-17],[-35,-40],[-78,49],[-18,-3],[-15,52],[-6,39],[4,36],[46,-28],[25,-51],[64,-37]],[[56429,16323],[-43,-108],[-23,9],[13,136],[30,22],[12,-1],[11,-58]],[[57192,13454],[-18,-74],[-57,-139],[-129,-35],[-141,-7],[-11,38],[-3,34],[10,52],[-1,21],[-6,21],[52,22],[34,65],[54,11],[66,-46],[36,-1],[54,50],[43,108],[25,-14],[-8,-106]],[[56550,14473],[-38,-30],[8,52],[74,87],[40,77],[23,28],[8,21],[32,28],[16,50],[-4,40],[-34,67],[1,47],[11,34],[57,15],[-15,-50],[22,-142],[-75,-177],[-64,-55],[-62,-92]],[[55222,17740],[90,-113],[76,-42],[83,-141],[35,-50],[6,-46],[-14,-209],[-18,-125],[3,-111],[-20,33],[-19,72],[-32,42],[-11,21],[58,5],[5,114],[28,89],[-4,93],[-68,103],[-46,90],[-71,29],[-66,90],[-39,11],[-48,-16],[18,54],[12,71],[8,13],[34,-77]],[[55757,15707],[-22,-6],[-41,26],[-95,121],[-44,11],[14,68],[34,-24],[77,-104],[29,-53],[48,-39]],[[56205,15152],[-3,-26],[-69,200],[-86,48],[-50,48],[11,28],[34,12],[6,63],[-15,68],[-45,138],[-26,93],[-12,22],[-3,78],[54,-121],[22,-108],[37,-105],[24,-183],[69,-62],[50,-88],[2,-105]],[[55951,16455],[4,-51],[-37,10],[-54,-19],[-19,-1],[12,35],[37,47],[2,45],[-46,64],[-54,161],[-25,38],[-12,60],[-45,65],[10,34],[9,8],[30,-16],[69,-231],[5,-22],[114,-227]],[[54553,14387],[-18,-3],[-54,42],[-16,34],[19,47],[4,51],[7,3],[9,-61],[41,-26],[3,-14],[25,-52],[-20,-21]],[[52292,18704],[-20,-26],[-39,19],[-21,34],[-9,67],[34,-71],[12,-16],[43,-7]],[[54627,16961],[-2,-121],[6,-90],[-5,-32],[-59,-61],[-15,-33],[-56,-34],[-33,-48],[-17,77],[-33,48],[-5,81],[-25,-28],[-36,17],[-59,61],[-37,84],[52,14],[11,-52],[41,64],[-9,33],[-8,5],[-13,62],[62,163],[14,105],[-29,168],[27,11],[70,-59],[32,-58],[1,-80],[29,-62],[43,-148],[53,-87]],[[52548,16854],[-36,-10],[29,39],[9,60],[17,-47],[-1,-28],[-18,-14]],[[53298,19527],[9,-73],[-50,170],[-62,267],[-33,207],[22,-56],[22,-116],[92,-399]],[[52007,13832],[-51,-53],[-108,-74],[-58,-3],[-58,28],[-41,62],[-23,62],[2,29],[37,-49],[31,-24],[26,16],[19,27],[-60,198],[3,42],[47,109],[128,-33],[22,-19],[18,-69],[29,-54],[33,-145],[4,-50]],[[52266,18802],[0,-58],[-28,98],[-18,108],[26,-35],[20,-113]],[[53753,14922],[1,-26],[-95,73],[-41,76],[-16,18],[27,2],[107,-126],[17,-17]],[[54219,18979],[66,-31],[35,4],[22,20],[97,-8],[79,30],[11,-53],[-2,-28],[-167,-26],[-153,-78],[-85,-52],[-39,-6],[-30,28],[-102,160],[27,-16],[75,-91],[47,17],[43,60],[6,45],[-6,22],[20,72],[56,-69]],[[53192,17395],[-108,-187],[12,47],[42,100],[14,47],[28,30],[27,54],[2,64],[38,44],[13,6],[-68,-205]],[[49076,10193],[-26,-19],[21,58],[4,37],[29,153],[22,-1],[6,-13],[-56,-215]],[[49129,9692],[-42,-140],[-2,40],[17,103],[24,37],[14,38],[4,45],[20,-22],[-5,-45],[-30,-56]],[[47046,10944],[-61,-23],[-11,21],[123,102],[21,-4],[8,-14],[-64,-48],[-16,-34]],[[52185,1666],[-6,-46],[-42,84],[25,13],[9,-1],[14,-50]],[[52655,11614],[17,-20],[23,12],[6,21],[84,-15],[13,-42],[-64,-1],[-29,-27],[-14,-5],[-55,7],[-8,96],[15,10],[12,-36]],[[49626,12566],[-28,-31],[-16,110],[23,107],[27,62],[54,7],[36,20],[4,-27],[-29,-83],[-71,-165]],[[62845,9509],[-14,-15],[-24,61],[-40,17],[-35,36],[-1,8],[1,18],[7,21],[17,17],[44,-50],[21,-62],[19,-28],[5,-23]],[[63789,6098],[-30,-104],[-60,66],[-7,82],[6,50],[36,95],[28,62],[20,20],[12,-82],[-5,-189]],[[63827,6769],[-6,-56],[-14,-12],[-20,48],[-90,-7],[-13,43],[-1,20],[43,72],[-51,19],[-20,31],[-42,152],[3,44],[18,23],[29,3],[54,-49],[41,-69],[14,-1],[4,-20],[-8,-49],[24,-42],[11,-31],[24,-119]],[[63108,9011],[-4,-59],[-36,18],[-1,52],[17,53],[9,5],[15,-69]],[[62644,10489],[-87,-50],[4,29],[69,69],[26,-5],[-12,-43]],[[63326,9844],[-8,-27],[-45,48],[-13,89],[1,19],[8,9],[17,-17],[24,-6],[14,-29],[2,-86]],[[63607,8162],[-31,-4],[-13,9],[-4,60],[23,41],[14,10],[21,-46],[5,-39],[-15,-31]],[[63343,9307],[-18,-40],[-60,16],[-13,50],[-2,35],[38,70],[42,-30],[17,-33],[12,-7],[0,-29],[-5,-20],[-11,-12]],[[63579,7522],[-51,-21],[-22,172],[-36,125],[7,78],[6,30],[75,-48],[23,-58],[15,-154],[-17,-124]],[[63554,8502],[-64,-11],[-42,8],[-10,72],[22,60],[-15,74],[9,44],[22,28],[35,-38],[6,-55],[22,-50],[99,-107],[-84,-25]],[[63412,8279],[-44,-45],[-22,14],[-27,86],[-19,238],[15,40],[10,14],[58,-29],[24,-34],[25,-21],[-13,-43],[6,-179],[-13,-41]],[[63594,5466],[-13,21],[5,78],[24,43],[23,25],[24,3],[8,-64],[-6,-91],[-2,-5],[-63,-10]],[[57184,14688],[34,-15],[18,1],[62,-30],[37,-41],[9,-17],[-20,-36],[-57,69],[-50,9],[-71,-2],[-28,13],[19,75],[47,-26]],[[57934,14025],[-2,-13],[-31,38],[-59,1],[-9,52],[23,9],[75,-19],[17,-46],[-14,-22]],[[57252,13767],[-72,-49],[-8,58],[37,46],[43,-55]],[[57815,14100],[-9,-10],[-32,34],[-4,28],[-11,1],[-19,25],[4,32],[43,1],[18,-89],[10,-22]],[[57569,14111],[62,-55],[38,8],[3,-11],[-21,-13],[-5,-10],[-59,-15],[-19,3],[-4,37],[5,56]],[[9395,50761],[-40,-16],[-43,8],[-33,99],[26,3],[53,65],[118,53],[30,6],[-111,-218]],[[13697,49312],[-69,-32],[-12,1],[-43,-65],[-32,-28],[-44,51],[11,80],[38,53],[190,-18],[15,-15],[0,-12],[-14,-11],[-40,-4]],[[14234,52542],[-27,-4],[36,50],[24,96],[33,-14],[6,-17],[-54,-99],[-18,-12]],[[13981,50669],[-16,-6],[-40,30],[-26,33],[17,24],[77,52],[37,1],[15,-8],[6,-16],[-3,-24],[-17,-30],[-50,-56]],[[13570,50035],[-15,-27],[-14,7],[-31,32],[-62,48],[-28,32],[-2,14],[22,15],[72,-38],[31,-35],[27,-48]],[[12306,48011],[-21,-32],[-42,2],[-23,11],[-7,27],[63,83],[15,12],[11,-4],[5,-34],[-1,-65]],[[12786,48624],[-37,-23],[-13,8],[-1,19],[9,31],[17,31],[57,59],[58,40],[29,-3],[10,-25],[-36,-51],[-93,-86]],[[13044,48703],[-27,-2],[-40,26],[5,32],[57,38],[58,-6],[5,-21],[-3,-24],[-5,-14],[-18,-14],[-32,-15]],[[17294,52630],[12,-28],[104,7],[30,-5],[12,-13],[-14,-18],[-40,-25],[-118,-43],[-96,-56],[-12,5],[-18,62],[-18,25],[-10,35],[0,13],[17,24],[35,35],[25,13],[91,-31]],[[17987,52517],[-17,-24],[-48,9],[-25,15],[85,64],[14,-13],[-9,-51]],[[18288,52000],[-26,-6],[39,65],[53,73],[49,45],[64,19],[-7,-33],[-84,-61],[-88,-102]],[[16564,51995],[-60,-15],[-14,30],[31,73],[25,42],[19,10],[69,82],[77,61],[70,88],[71,126],[14,46],[33,5],[54,-31],[34,-43],[-16,-34],[-180,-179],[-15,-23],[-17,-61],[-15,-21],[-23,-10],[-18,-26],[-10,-44],[-23,-22],[-33,-2],[-24,-12],[-13,-22],[-36,-18]],[[16408,52246],[-28,-30],[-108,18],[23,60],[82,38],[90,-59],[-59,-27]],[[16606,52631],[0,-27],[-17,-25],[16,-46],[-28,-80],[-11,-50],[-15,-32],[-15,-12],[-16,7],[-3,17],[9,28],[-36,-1],[-11,70],[20,22],[7,29],[3,20],[23,89],[8,6],[3,-21],[8,-6],[12,8],[19,40],[8,4],[16,-40]],[[16458,53005],[-69,-8],[-32,13],[-4,14],[13,50],[-1,21],[35,8],[40,-24],[11,-24],[7,-50]],[[24652,47086],[14,-14],[13,10],[22,36],[39,-8],[28,-11],[18,-15],[-10,-52],[-7,-84],[-16,-30],[-16,-42],[-55,25],[-45,53],[-64,92],[-37,67],[-2,28],[-23,20],[-45,115],[-25,90],[-40,10],[-51,26],[-19,49],[13,43],[72,22],[108,-111],[17,-48],[39,-55],[7,-78],[20,-30],[45,-108]],[[24039,49032],[36,-69],[51,6],[28,-52],[22,-79],[-16,-50],[-23,11],[-25,-29],[-15,-96],[7,-97],[-8,-96],[-29,-99],[-6,-65],[-11,-20],[-13,-7],[-16,18],[-24,14],[-30,-56],[-36,0],[-31,126],[27,210],[61,43],[-36,57],[-76,66],[6,37],[-57,107],[-4,25],[9,89],[54,79],[72,15],[49,-34],[28,-29],[6,-25]],[[24412,47734],[11,-28],[1,-18],[-78,-66],[-2,-14],[-18,-41],[-17,-15],[-29,-44],[-56,-48],[8,144],[-55,83],[56,42],[36,-12],[61,-4],[60,36],[22,-15]],[[25480,47271],[56,-130],[3,-46],[-53,-16],[-42,7],[-21,15],[-6,21],[14,66],[-28,38],[-32,14],[-28,-23],[-2,65],[22,46],[-12,64],[0,48],[7,15],[28,-2],[59,-49],[35,-133]],[[24715,48714],[-6,-15],[-69,2],[-25,12],[-9,42],[6,40],[16,30],[20,59],[15,98],[102,-109],[31,-49],[16,-60],[-36,-22],[-43,-11],[-18,-17]],[[25060,48298],[-11,-165],[-22,9],[-20,0],[-42,-23],[-44,10],[-21,18],[-7,22],[7,49],[-24,27],[-83,10],[-31,12],[-17,52],[-4,68],[13,25],[42,18],[32,83],[18,11],[69,163],[34,-11],[60,-100],[76,-143],[-25,-135]],[[24697,48436],[-28,-3],[-33,15],[-78,81],[-1,24],[11,27],[45,50],[18,12],[106,-5],[34,-14],[8,-23],[-2,-23],[-11,-24],[-2,-25],[6,-26],[-13,-26],[-60,-40]],[[23864,50414],[-4,-25],[-74,2],[-75,37],[-37,47],[7,22],[70,20],[66,-45],[47,-58]],[[26787,44706],[-17,-45],[-73,110],[-26,31],[-54,116],[-10,50],[2,29],[10,10],[17,-7],[15,-13],[105,-129],[29,-61],[2,-91]],[[25705,46671],[-13,-1],[-5,21],[3,43],[18,73],[8,17],[70,-13],[11,-5],[2,-14],[-5,-23],[-23,-33],[-66,-65]],[[25649,44159],[-10,-2],[-18,21],[-12,33],[-8,87],[6,36],[5,14],[47,-55],[-10,-134]],[[26080,46152],[-16,-36],[-39,-56],[-25,-22],[-13,11],[-34,6],[-36,56],[-29,23],[-19,1],[-9,-20],[-3,-27],[5,-34],[-3,-16],[-11,3],[-11,16],[-10,29],[-3,29],[7,30],[26,40],[83,82],[25,15],[28,-3],[44,-39],[10,-11],[33,-77]],[[26291,45363],[-11,-3],[-35,12],[-118,168],[-84,62],[-57,83],[-59,54],[36,86],[31,-10],[109,-71],[87,-68],[49,-45],[103,-191],[-7,-28],[-44,-49]],[[27338,43672],[-10,-16],[-21,0],[-35,17],[-25,37],[-31,115],[3,20],[11,20],[49,41],[18,-4],[7,-30],[29,-68],[9,-20],[0,-79],[-4,-33]],[[27096,44598],[-42,-14],[14,54],[4,31],[-7,30],[-6,58],[-2,135],[41,86],[64,2],[-1,-43],[-27,-192],[-14,-89],[-10,-33],[-14,-25]],[[26582,45188],[-8,-8],[-23,13],[-21,26],[-37,74],[-12,31],[-8,50],[7,9],[16,-5],[12,-8],[58,-122],[16,-60]],[[26662,45313],[-3,-7],[-57,0],[-16,10],[-9,21],[-4,32],[9,38],[33,73],[1,31],[7,14],[24,-37],[9,-34],[6,-141]],[[21291,61738],[-45,-37],[-71,39],[-19,19],[83,52],[36,-2],[77,-31],[28,-27],[-89,-13]],[[38987,13871],[-15,-82],[-42,33],[-14,51],[-3,85],[23,14],[33,-35],[7,-24],[11,-42]],[[36777,17257],[15,-112],[-28,17],[-36,60],[-25,78],[2,35],[5,6],[52,-46],[15,-38]],[[37191,16406],[-12,-3],[-27,53],[-13,130],[6,13],[52,-161],[-1,-19],[-5,-13]],[[35369,21297],[-4,-17],[-55,62],[-126,210],[-46,105],[-7,50],[3,110],[44,-13],[50,-74],[21,-68],[1,-79],[94,-37],[13,-169],[12,-80]],[[34532,23945],[-18,-4],[-33,46],[2,42],[8,4],[33,-37],[12,-32],[-4,-19]],[[34273,20316],[-7,-32],[-92,67],[51,116],[-8,122],[22,26],[20,-41],[26,-155],[-12,-103]],[[36161,16648],[-7,-47],[-164,187],[40,17],[45,-13],[86,-144]],[[36487,18272],[7,-36],[-4,-11],[-22,25],[-38,-149],[-11,-14],[23,205],[24,29],[26,6],[-5,-55]],[[35887,21250],[-41,-235],[-42,4],[-86,74],[-9,46],[33,273],[25,36],[76,37],[11,-34],[8,-80],[25,-121]],[[35966,16800],[-10,-11],[-47,115],[-7,79],[-19,34],[-47,26],[41,162],[34,333],[15,-61],[-36,-338],[2,-44],[17,-42],[18,-70],[2,-75],[32,-70],[5,-38]],[[49816,63169],[22,-22],[124,28],[106,32],[164,77],[98,26],[299,0],[51,-15],[-23,-45],[-12,-13],[10,-20],[32,-27],[64,-30],[25,27],[19,64],[44,265],[19,80],[9,76],[-2,71],[-21,45],[-76,27],[-105,-4],[-54,7],[-65,14],[-48,22],[-31,30],[-63,89],[-46,50],[-118,90],[-53,30],[27,35],[107,41],[65,38],[76,114],[46,18],[164,-15],[223,-89],[140,-76],[37,-8],[1,13],[-36,36],[-160,95],[-73,69],[-35,50],[16,21],[91,22],[12,25],[-124,30],[-62,-1],[-50,-21],[-54,-2],[-101,39],[-27,22],[-58,67],[-31,58],[-33,36],[-12,28],[-7,89],[3,52],[14,45],[24,37],[65,68],[37,20],[68,9],[148,-34],[399,-123],[-10,40],[-446,166],[-157,42],[-39,60],[238,230],[218,54],[109,67],[178,3],[167,-43],[3,12],[-75,79],[6,20],[94,48],[174,55],[212,45],[42,23],[54,16],[100,14],[249,7],[139,-7],[186,-33],[108,-62],[34,-36],[57,-118],[47,-166],[69,-68],[111,-38],[76,-41],[42,-45],[12,-56],[-20,-68],[15,-69],[49,-71],[38,-40],[84,-46],[1,-25],[-26,-28],[-55,-38],[-136,-120],[-176,-133],[-126,-114],[-6,-34],[261,179],[81,-6],[4,-25],[-53,-87],[-65,-78],[-65,-50],[12,-19],[124,-88],[-23,-14],[-60,7],[-24,-8],[-18,-16],[-11,-24],[-1,-34],[10,-43],[-1,-33],[-12,-21],[13,-9],[36,4],[30,17],[53,59],[174,161],[111,60],[36,5],[102,-39],[24,2],[-113,123],[-10,32],[23,46],[14,16],[63,33],[51,19],[29,-8],[46,-63],[22,-43],[38,-18],[85,23],[57,53],[70,-35],[105,-83],[-9,-84],[0,-84],[5,-61],[126,-112],[88,-50],[16,0],[-2,17],[-18,37],[-48,38],[-44,57],[-39,71],[23,164],[66,86],[64,-22],[83,-50],[66,-4],[104,5],[212,-100],[114,-2],[-10,40],[-87,20],[-126,55],[-196,66],[-90,75],[-17,36],[3,38],[11,33],[20,29],[39,29],[190,87],[134,37],[102,12],[169,0],[197,-16],[106,-25],[123,-62],[154,-61],[55,-11],[65,2],[74,13],[70,-5],[223,-90],[59,-47],[35,-56],[27,-55],[17,-53],[-7,-43],[-186,-186],[-80,-32],[-54,-71],[-79,-133],[-68,-72],[-6,-14],[14,-4],[41,33],[70,92],[50,80],[94,66],[152,78],[133,37],[114,-3],[95,-11],[77,-19],[46,-16],[15,-13],[31,-59],[-2,-40],[-19,-45],[-37,-51],[-166,-56],[-92,-44],[-56,-17],[-170,-16],[8,-18],[126,-24],[140,8],[-3,-28],[-66,-76],[-22,-66],[19,-54],[-4,-44],[-49,-92],[-56,-84],[21,-12],[129,120],[34,131],[52,115],[61,63],[46,24],[144,10],[80,67],[68,22],[29,0],[58,-25],[-4,-26],[-84,-121],[-179,-194],[73,22],[49,46],[67,46],[75,69],[49,-62],[76,-47],[46,-105],[74,-51],[44,-40],[-6,67],[-64,134],[17,54],[50,27],[155,113],[108,-38],[66,-33],[34,9],[84,-5],[135,-18],[131,-32],[128,-45],[98,-52],[69,-60],[42,-42],[14,-23],[24,-60],[-18,-39],[-97,-92],[-53,-42],[-54,-19],[-143,19],[-44,-11],[-47,-29],[-149,-126],[-82,-54],[-81,-35],[-19,-19],[174,2],[48,38],[41,70],[76,73],[145,33],[202,-72],[101,3],[76,72],[86,49],[34,10],[18,-6],[65,-51],[20,-45],[-2,-104],[-9,-32],[-57,-78],[-142,-118],[-92,-44],[-103,-24],[-112,-40],[-39,-32],[-39,-46],[-38,-31],[-49,-25],[64,-37],[25,1],[23,22],[65,88],[48,38],[27,8],[28,-3],[27,-18],[27,-30],[-2,-75],[-82,-297],[14,1],[50,80],[145,309],[36,62],[70,63],[158,94],[121,49],[138,42],[73,16],[84,-11],[54,-48],[74,-9],[90,12],[57,-6],[66,-19],[56,-36],[95,-41],[215,-77],[27,-16],[24,-29],[24,-42],[-3,-41],[-29,-42],[-36,-25],[-43,-9],[-44,-22],[-82,-59],[-27,-10],[-129,-25],[-119,-12],[-74,-24],[-143,-64],[-198,-118],[2,-28],[79,-14],[64,18],[88,82],[82,31],[129,25],[178,22],[76,-4],[14,-4],[10,-19],[6,-34],[-29,-45],[-34,-21],[-92,-101],[61,-26],[83,-11],[47,27],[43,62],[48,34],[54,7],[47,15],[40,25],[11,16],[-59,33],[-5,19],[24,48],[44,53],[45,33],[33,3],[111,-36],[76,-61],[191,-185],[25,-36],[67,-138],[12,-61],[-11,-42],[-15,-26],[-21,-10],[-42,0],[-255,56],[-117,-7],[-51,-16],[-41,-24],[-31,-29],[-23,-36],[-45,-21],[-162,0],[-91,-20],[-156,-49],[-56,-27],[-13,-36],[96,6],[157,46],[148,13],[248,-101],[81,-16],[46,15],[54,5],[198,-7],[67,-13],[101,-38],[153,-85],[29,-24],[17,-27],[5,-27],[-1,-67],[-16,-23],[-52,-15],[-220,18],[-67,14],[-83,-18],[-67,6],[-86,27],[-94,48],[-141,-45],[-114,29],[-115,-26],[-229,-108],[25,-18],[314,92],[61,-6],[99,-33],[157,-67],[44,-27],[1,-105],[-24,-70],[-48,-79],[-72,10],[-168,50],[-68,6],[-51,-8],[-67,-31],[-32,-1],[-268,62],[-61,3],[-7,-6],[13,-12],[244,-97],[179,-11],[113,-17],[67,-29],[32,-22],[2,-65],[60,-65],[54,-27],[34,-1],[60,24],[60,4],[48,-17],[61,-36],[73,-10],[64,-22],[50,-4],[140,10],[60,-14],[16,-12],[-27,-21],[-127,-50],[-20,-48],[72,-63],[38,-47],[-2,-36],[-38,-81],[-10,-33],[13,-3],[93,66],[14,-8],[10,-91],[12,5],[31,75],[-14,101],[54,39],[174,30],[-30,-157],[-4,-82],[-76,-136],[-63,-44],[2,-9],[45,-17],[28,-2],[27,21],[64,105],[130,110],[23,2],[0,-39],[-17,-74],[60,-35],[57,35],[31,29],[72,-4],[33,-14],[10,-33],[-33,-136],[6,-33],[76,-91],[7,5],[-14,44],[-16,108],[15,47],[63,60],[128,87],[47,17],[30,-12],[47,-42],[-14,-24],[-51,-26],[-38,-47],[-25,-68],[27,-37],[104,-3],[105,56],[60,-28],[71,-72],[131,-117],[74,32],[92,-89],[-123,-70],[37,-148],[-160,6],[-91,-11],[-60,13],[-65,-5],[60,-35],[116,-14],[11,-45],[91,1],[68,9],[124,-2],[6,51],[81,30],[46,32],[38,-19],[112,-22],[150,-101],[-66,-61],[-18,-57],[-23,-48],[-12,-44],[-26,-30],[-215,-172],[36,-1],[90,41],[177,62],[99,25],[70,-17],[36,0],[31,22],[58,-26],[122,-23],[139,141],[84,-28],[79,-87],[168,-153],[88,-89],[29,-40],[-4,-40],[-79,-42],[-41,-9],[-107,80],[-98,40],[-59,-4],[-60,-31],[19,-17],[237,-122],[42,-91],[3,-39],[-159,-60],[-51,-4],[-110,29],[-64,53],[-53,20],[-74,6],[-23,-10],[80,-91],[-8,-27],[-41,-18],[-21,-44],[159,-79],[118,-80],[18,-32],[-80,-24],[-58,-6],[-121,12],[-67,17],[-18,-18],[69,-42],[26,-29],[21,-40],[12,-37],[4,-35],[-57,-29],[-68,-80],[-26,-84],[-61,-8],[-25,16],[-83,-25],[-108,35],[-39,38],[-119,157],[-3,-18],[30,-79],[-6,-47],[-126,-34],[0,-14],[78,-25],[93,-19],[-15,-73],[1,-313],[-21,-111],[-46,-97],[-65,-93],[-71,61],[-29,62],[-24,32],[-34,26],[-43,12],[-47,0],[-49,-55],[-54,48],[-51,58],[19,151],[22,77],[-9,-1],[-29,-36],[-71,-111],[-46,-136],[-60,52],[-54,65],[-45,66],[-72,75],[-70,88],[-37,105],[-17,21],[-40,87],[-17,25],[-14,8],[-35,54],[13,58],[55,68],[50,50],[83,48],[97,27],[45,63],[54,114],[59,79],[64,45],[-32,8],[-82,-38],[-57,-56],[-69,-93],[-64,-60],[-163,-69],[-60,-14],[-70,-7],[-153,9],[-36,24],[19,66],[109,118],[-17,8],[-39,-42],[-53,-29],[-45,-14],[-68,5],[-78,74],[-38,22],[-77,26],[-31,25],[-128,180],[-26,48],[-15,47],[-41,40],[-67,32],[-16,-6],[24,-40],[1,-34],[-59,-22],[-61,7],[-64,37],[-6,-49],[69,-88],[1,-110],[-19,-12],[-47,-6],[-31,13],[-104,83],[-99,58],[-69,32],[-8,-23],[45,-100],[52,-98],[86,-82],[136,-96],[62,-56],[-49,-79],[-42,-26],[-27,-8],[-82,0],[-151,44],[-71,48],[-103,116],[-170,119],[-37,0],[-120,-49],[19,-8],[78,-3],[56,-16],[136,-94],[11,-40],[-34,-44],[2,-56],[38,-67],[39,-43],[80,-31],[40,-4],[15,-19],[-48,-151],[-4,-41],[14,-17],[17,-1],[102,62],[43,15],[37,3],[44,-17],[49,-38],[29,-39],[9,-40],[14,-26],[101,-43],[-9,-20],[-104,-63],[-6,-11],[21,-4],[66,-38],[60,-60],[37,-71],[7,-35],[0,-33],[8,-20],[32,-3],[13,12],[15,-2],[16,-17],[17,-55],[37,-159],[19,-45],[11,-1],[6,159],[16,27],[64,-28],[94,-63],[66,-55],[8,-26],[-50,-50],[11,-23],[36,-32],[34,12],[25,56],[42,56],[49,39],[93,-33],[77,-82],[12,-29],[50,-35],[44,20],[85,-95],[-40,-43],[-89,-62],[-9,-22],[21,5],[170,0],[45,-25],[11,-48],[-75,-133],[-69,12],[-91,3],[-47,-7],[7,-17],[127,-61],[36,-51],[48,-52],[24,-42],[-1,-20],[-20,-29],[9,-10],[87,-20],[54,18],[67,6],[60,-4],[4,-19],[-8,-48],[-61,-45],[16,-11],[71,13],[33,-21],[43,-107],[47,-83],[-39,-20],[-43,-7],[6,-107],[29,-109],[0,-105],[-9,-94],[-39,-21],[-43,3],[-17,23],[-103,279],[-26,51],[-31,44],[-109,120],[4,-19],[27,-56],[24,-83],[32,-165],[14,-107],[-6,-40],[-23,-10],[-6,-19],[11,-29],[84,-109],[41,-65],[28,-68],[27,-46],[25,-24],[-6,-20],[-37,-15],[-64,-7],[-29,9],[-114,63],[-17,-20],[64,-230],[-2,-55],[-32,-20],[-39,23],[-47,65],[-71,73],[-95,81],[-92,63],[-21,-2],[-14,-19],[-15,-3],[-17,13],[-31,47],[-31,32],[-134,107],[-13,1],[12,-32],[14,-70],[-15,-15],[-35,5],[-66,31],[-45,70],[-56,121],[-30,46],[-3,-30],[16,-115],[-3,-39],[-33,-11],[-14,11],[-14,31],[-13,51],[-32,38],[-50,27],[-28,28],[-14,49],[-9,12],[-89,-12],[-44,36],[-127,140],[-116,152],[-74,81],[-27,19],[40,-99],[42,-145],[11,-66],[-19,-3],[-43,29],[-221,187],[-136,89],[-76,15],[-123,11],[-28,-48],[66,-108],[65,-82],[63,-54],[98,-107],[90,-137],[37,-42],[122,-59],[65,-15],[66,-5],[6,-21],[-32,-39],[-8,-24],[147,-61],[55,-34],[53,-56],[31,-14],[126,-143],[32,-23],[113,-46],[37,-29],[63,-92],[39,-47],[55,-111],[41,-49],[101,-56],[44,-16],[19,-22],[-13,-50],[-13,-21],[-57,-35],[9,-48],[32,-87],[-1,-53],[-35,-21],[-73,-25],[-36,2],[-55,21],[-69,35],[-137,86],[-205,61],[-77,32],[-25,29],[-39,17],[-510,83],[-86,21],[-53,26],[-49,37],[-195,86],[-24,19],[-130,148],[-98,172],[-32,23],[-106,23],[-89,-15],[-59,-19],[-90,7],[-58,28],[-124,77],[-125,41],[-109,68],[-57,24],[6,17],[81,101],[-25,-1],[-142,-77],[-51,24],[-84,60],[-63,61],[-129,167],[-75,61],[11,14],[83,5],[67,-5],[45,14],[86,68],[37,43],[4,25],[-72,6],[-16,12],[-13,28],[-33,36],[-53,43],[-61,19],[-210,-16],[-37,19],[1,29],[41,82],[22,33],[8,18],[-9,3],[-28,-2],[-122,-73],[-26,6],[-48,77],[-28,88],[-22,31],[-28,10],[-99,87],[-142,164],[-53,51],[-58,46],[-42,20],[6,24],[90,137],[4,22],[-77,-8],[-117,28],[-55,-34],[-35,-2],[-41,19],[-23,-7],[-21,-112],[-17,-28],[-24,-16],[-22,2],[-18,19],[0,27],[-17,137],[-41,21],[-115,5],[-24,11],[-28,25],[-24,48],[-20,68],[-22,38],[-26,7],[-21,-6],[-15,-18],[-36,-11],[-56,-4],[-1,-26],[52,-49],[51,-69],[48,-90],[-29,-61],[-109,-31],[-94,-9],[-80,13],[-62,21],[-87,50],[-123,-16],[-28,-132],[-27,-7],[-117,3],[-47,-12],[-157,-73],[-48,-10],[-36,9],[-36,-18],[-53,-42],[-72,-4],[-91,33],[-77,14],[-64,-5],[-65,20],[-67,44],[-55,19],[-72,-3],[-17,7],[-105,94],[-33,38],[-71,118],[-12,47],[-2,50],[6,37],[25,56],[26,131],[22,44],[33,39],[65,50],[233,89],[47,35],[-2,23],[-53,108],[1,28],[18,16],[37,63],[17,18],[42,10],[85,-32],[74,-13],[97,-4],[162,-44],[225,-85],[130,-58],[98,-86],[70,-85],[10,-42],[-32,-66],[-17,-20],[1,-22],[18,-26],[57,-37],[13,14],[-4,45],[12,37],[27,31],[4,39],[-22,49],[-27,42],[-32,35],[-146,122],[-14,41],[49,18],[214,-41],[81,9],[31,47],[34,33],[36,17],[71,8],[101,-22],[49,-3],[43,8],[57,25],[84,86],[53,20],[81,14],[61,1],[110,-34],[68,1],[-6,57],[-45,110],[-55,114],[-44,38],[-112,72],[-133,135],[-68,84],[-17,42],[9,28],[23,41],[240,150],[190,148],[83,76],[40,54],[42,38],[43,24],[91,29],[26,37],[6,63],[15,55],[86,147],[66,40],[99,27],[66,36],[79,120],[-8,30],[-37,23],[-28,35],[-121,316],[-81,152],[-97,133],[-87,162],[-144,160],[-2,42],[26,48],[-13,10],[-148,-69],[-35,-4],[-57,31],[-39,38],[-32,66],[3,35],[22,32],[28,81],[1,41],[-10,40],[-13,27],[-16,15],[-45,10],[-75,4],[-25,-14],[83,-122],[-13,-30],[-105,-13],[-47,6],[-44,14],[-39,24],[-123,127],[-26,48],[8,35],[-10,18],[-25,-12],[-34,0],[-46,12],[-10,15],[86,69],[5,21],[-39,23],[-60,4],[-15,21],[20,21],[80,38],[29,25],[-48,18],[-27,2],[-54,-41],[-82,-84],[-59,-31],[-82,39],[-51,13],[-36,-9],[-54,-65],[-119,-47],[-214,-112],[-91,-36],[-99,7],[-19,22],[2,39],[7,31],[12,25],[3,30],[-7,128],[17,35],[34,21],[62,22],[158,-27],[74,5],[52,30],[52,42],[51,56],[11,53],[-55,87],[-20,19],[-141,68],[-78,24],[-69,11],[-50,20],[-31,28],[-30,47],[-3,32],[5,42],[29,29],[126,33],[-1,9],[-104,25],[-48,-3],[-42,-28],[-52,-65],[-31,-18],[-94,38],[-57,6],[-38,18],[-21,18],[13,18],[47,18],[82,56],[5,30],[-56,49],[-30,12],[-117,18],[-143,-18],[-53,9],[-23,55],[-15,65],[-7,76],[-25,129],[-29,67],[-37,9],[-171,-28],[-41,0],[-28,10],[-113,86],[-45,31],[-26,7],[-82,92],[-32,18],[-37,45],[-44,72],[-47,23],[-50,-29],[-51,-40],[-51,-53],[-28,-44],[-4,-37],[32,-28],[179,-47],[45,-32],[39,-51],[29,-64],[19,-74],[-2,-56],[-22,-37],[-39,-33],[-110,-52],[-114,-31],[-116,-7],[-54,7],[-298,101],[-53,1],[-69,14],[-154,41],[-84,5],[-149,34],[-250,20],[-51,-16],[67,-47],[58,-23],[51,0],[72,-42],[94,-83],[54,-49],[44,-59],[2,-20],[-45,-40],[-349,211],[-214,-74],[-98,-27],[-85,-4],[-105,29],[-238,102],[-90,35],[-32,5],[-207,-44],[-179,-2],[-360,42],[-133,29],[-35,29],[-43,15],[-78,0],[-207,33],[-189,-74],[-227,68],[-67,40],[-22,28],[-66,115],[-9,62],[19,56],[18,38],[19,21],[-124,-64],[-43,-11],[-57,-2],[-171,23],[-27,-12],[9,-22],[45,-33],[6,-19],[-95,-16],[-144,16],[-63,-7],[-29,-9],[-64,-51],[-27,-12],[-34,6],[-151,116],[-122,74],[-141,28],[-66,24],[-35,28],[-196,236],[-27,51],[-62,185],[-21,40],[-25,26],[50,5],[185,-22],[179,1],[97,-15],[113,-46],[147,-33],[105,-7],[169,12],[192,31],[23,24],[-125,41],[-110,55],[-102,70],[-62,30],[-102,19],[-287,13],[-267,48],[-184,64],[-150,71],[-61,39],[-22,29],[-23,94],[-24,156],[-24,106],[-23,53],[-3,47],[51,100],[145,109],[4,17],[-29,5],[-61,28],[-20,40],[-7,64],[-1,55],[7,44],[24,56],[63,100],[90,122],[97,114],[16,37],[9,101],[13,74],[13,52],[21,39],[60,74],[75,70],[117,60],[10,22],[2,31],[7,22],[12,16],[291,191],[132,78],[113,49],[134,37],[383,74],[197,21],[248,-5],[456,-41],[56,-30],[14,-15],[20,-42],[-16,-27],[-125,-91],[-157,-75],[-102,-66],[-174,-149],[-46,-52],[-215,-300],[-51,-49],[-29,-40],[-21,-108],[6,-38],[33,-63],[117,-136],[31,-64],[0,-59],[-14,-139],[-1,-71],[5,-68],[24,-97],[43,-127],[99,-127],[156,-130],[115,-86],[115,-63],[135,-93],[30,-45],[-62,-50],[-146,-77],[-191,-32],[-103,-31],[-127,-67],[-160,-53],[-63,-32]],[[45619,63922],[74,-7],[46,13],[150,-5],[35,-24],[-2,-26],[-17,-42],[10,-36],[102,-69],[91,-48],[83,-60],[122,-125],[27,-35],[18,-38],[32,-140],[4,-53],[-12,-153],[-10,-28],[-32,-35],[11,-14],[97,-40],[77,-79],[40,-26],[98,-43],[17,-16],[22,-26],[55,-115],[90,-103],[7,-22],[-20,-49],[14,-15],[35,-18],[31,9],[28,37],[30,10],[32,-15],[25,-27],[33,-66],[49,-55],[-4,-16],[-24,-14],[-132,-15],[-74,10],[-69,27],[-48,27],[-44,43],[-16,-5],[-24,-34],[-49,-51],[-31,-46],[35,-21],[174,2],[37,-13],[45,-33],[-51,-55],[-117,-91],[-253,-178],[-75,-46],[18,-13],[28,-4],[87,8],[82,24],[98,-10],[44,-20],[-15,-19],[27,-31],[162,-72],[102,15],[104,69],[81,34],[99,-4],[28,-8],[-11,-17],[-73,-39],[-66,-43],[-8,-12],[83,17],[184,-27],[89,-7],[65,7],[61,-11],[56,-29],[20,-19],[-54,-11],[-50,0],[-42,-18],[-36,-35],[-25,-46],[-15,-56],[-38,-23],[-61,10],[-24,16],[13,22],[-18,3],[-49,-16],[-38,0],[-10,-17],[267,-179],[86,-159],[59,-65],[6,-18],[-39,-44],[-2,-33],[18,-97],[-8,-78],[-26,-136],[24,-42],[57,-38],[35,-48],[23,-16],[17,-39],[21,-24],[23,-10],[17,13],[9,35],[22,33],[61,60],[58,92],[10,31],[-9,73],[6,31],[55,110],[18,76],[17,119],[29,83],[62,70],[109,142],[39,28],[44,14],[77,-3],[57,-46],[75,-84],[96,-77],[176,-104],[50,-39],[99,-104],[42,-102],[28,-144],[25,-87],[21,-29],[9,-44],[-1,-59],[-8,-45],[-14,-33],[-22,-19],[-53,-6],[-65,9],[-19,15],[-35,69],[-14,4],[-61,-51],[-7,-28],[22,-94],[-3,-177],[6,-38],[65,-184],[108,-140],[271,-270],[15,-31],[29,-110],[14,-23],[18,-14],[21,-4],[29,11],[100,84],[86,88],[61,46],[34,4],[37,16],[39,28],[27,32],[13,34],[18,133],[15,63],[42,87],[16,25],[209,220],[18,27],[89,255],[31,117],[5,70],[-12,62],[6,52],[23,41],[26,29],[44,29],[24,41],[15,4],[35,0],[48,-28],[34,-5],[222,32],[0,17],[-130,54],[1,27],[11,36],[41,42],[50,13],[12,26],[1,32],[17,53],[-16,19],[-121,73],[-70,-3],[-18,9],[-61,60],[-22,86],[-2,35],[8,57],[8,16],[-4,26],[-15,36],[-1,31],[12,28],[-7,33],[-29,38],[-11,33],[30,96],[1,30],[-28,41],[-20,16],[15,10],[50,4],[60,-13],[69,-31],[86,-1],[101,31],[103,14],[177,-6],[42,-8],[176,-90],[137,-45],[62,4],[305,-18],[134,9],[69,-4],[133,-49],[-9,-41],[-58,-69],[-74,-14],[-66,-23],[61,-36],[181,-48],[42,-78],[13,-35],[-21,-32],[10,-16],[42,0],[108,27],[120,-18],[174,-60],[18,-12],[31,-47],[-4,-19],[-154,-117],[-80,-46],[-105,-47],[-3,-26],[147,-4],[115,-14],[52,-15],[27,-22],[37,-48],[6,-37],[-4,-51],[-12,-35],[-134,-101],[-61,-31],[-104,-38],[-46,-28],[-50,4],[-54,36],[-56,8],[-102,-29],[-55,1],[-27,-10],[-3,-22],[49,-64],[27,-25],[11,-19],[-20,-32],[4,-10],[17,-11],[92,-140],[20,-10],[19,5],[40,40],[24,16],[11,-2],[0,-21],[-43,-121],[-6,-33],[1,-29],[21,-61],[47,-66],[58,-60],[87,-80],[118,-85],[44,-41],[64,-98],[13,-38],[-17,-102],[-47,-168],[-30,-96],[-13,-23],[-89,-68],[-50,-15],[-83,1],[-27,-13],[-44,-54],[-59,-95],[-46,-60],[-34,-26],[-62,-29],[-97,-88],[-48,-34],[-167,-37],[-136,-119],[-54,-39],[-59,-21],[-63,-2],[-37,19],[-21,72],[-12,24],[-48,50],[-98,144],[-43,50],[-28,11],[-57,-8],[-29,4],[-63,46],[-23,29],[3,11],[48,17],[-20,24],[-85,64],[-35,34],[-4,11],[-84,44],[-84,12],[-105,-72],[-40,-49],[1,-16],[51,-19],[22,9],[42,45],[23,14],[65,-7],[54,-31],[20,-27],[7,-19],[147,-144],[52,-29],[22,-35],[16,-57],[32,-64],[71,-106],[75,-130],[15,-51],[-38,-25],[-20,-2],[-58,19],[-153,61],[-17,-1],[-40,-32],[-33,-73],[-11,-7],[-81,29],[-153,63],[-102,54],[-52,44],[-62,71],[-74,97],[-88,31],[-102,-36],[-148,-12],[-310,11],[-40,-8],[-16,-13],[26,-53],[-6,-15],[-21,-9],[-4,-16],[33,-57],[55,-38],[154,-52],[102,-45],[61,-38],[20,-32],[3,-36],[-29,-70],[-17,-26],[-358,-351],[-137,-142],[-69,-87],[-61,-58],[-54,-28],[-87,-16],[-123,-4],[-159,14],[-82,45],[-149,122],[-104,71],[-47,24],[-38,67],[-36,13],[-74,11],[-77,36],[-181,122],[-94,48],[-85,27],[-77,4],[-29,-7],[52,-61],[-23,-4],[-62,14],[-61,0],[-108,44],[-108,-6],[-77,9],[-93,24],[-99,12],[-161,-1],[-58,-5],[-9,-11],[78,-54],[132,-64],[-17,55],[4,15],[45,19],[210,-34],[238,-72],[61,-7],[67,-26],[74,-43],[102,-89],[195,-201],[62,-51],[83,-47],[423,-68],[145,0],[293,-19],[154,-36],[44,-27],[13,-89],[-14,-45],[-84,-136],[-53,-101],[-328,-425],[-43,-97],[-19,-59],[-58,-60],[-149,-94],[-149,-80],[-90,-17],[-79,19],[-52,23],[-77,79],[-5,-8],[56,-124],[-13,-13],[-45,16],[-103,56],[-33,-11],[-20,-15],[-27,0],[-36,16],[-63,46],[-17,23],[-15,69],[-11,12],[-125,-40],[-21,-13],[51,-27],[18,-20],[50,-101],[3,-22],[-36,-14],[-120,39],[-14,-4],[59,-101],[23,-48],[2,-25],[-77,-114],[-49,-48],[-67,-17],[-42,12],[-48,29],[-34,-4],[-21,-38],[-39,-27],[-57,-16],[-73,7],[-88,30],[-234,108],[-74,16],[-137,15],[-16,15],[1,15],[19,14],[-5,12],[-29,9],[-31,-9],[-31,-28],[-54,-9],[-78,12],[-115,42],[-229,107],[-250,91],[-145,119],[55,-107],[-5,-36],[-28,-31],[-4,-32],[57,-76],[78,-27],[79,3],[2,12],[-33,19],[-29,28],[-14,41],[15,7],[70,-21],[46,-25],[346,-141],[103,-26],[78,-27],[22,-16],[-27,-35],[-139,-87],[-2,-14],[95,8],[116,76],[65,36],[63,22],[85,-38],[107,-97],[86,-55],[124,-31],[72,-33],[122,-91],[19,-48],[11,-193],[-4,-47],[-16,-46],[-28,-46],[-50,-25],[-75,-6],[-58,-19],[-127,-102],[-54,-15],[-229,31],[-90,28],[-41,-1],[-23,-23],[-24,-10],[-90,-9],[-14,-19],[5,-28],[18,-39],[22,-22],[34,-28],[51,-18],[105,-21],[11,-52],[-5,-17],[-34,-35],[-40,4],[-69,41],[-34,2],[-29,-23],[-42,-7],[-53,7],[-29,-18],[-6,-43],[-18,-33],[-61,-52],[-33,-38],[1,-30],[35,-22],[41,-49],[45,-76],[9,-34],[-28,9],[-38,30],[-47,52],[-71,46],[-157,61],[-28,-2],[15,-15],[104,-64],[40,-41],[4,-29],[-87,-66],[-2,-22],[23,-19],[6,-16],[-31,-32],[-51,-27],[-99,-3],[-9,-16],[38,-33],[12,-20],[-32,-28],[-21,-4],[-114,13],[30,-70],[17,-25],[35,-35],[63,-32],[1,-12],[-22,-28],[-37,-34],[-157,-103],[-110,-121],[-16,-37],[28,-80],[1,-20],[-29,-36],[-65,10],[-12,-14],[15,-38],[3,-54],[-10,-73],[-47,-113],[-84,-154],[-64,-140],[-45,-126],[-32,-62],[-60,-4],[-45,-40],[31,-21],[19,-24],[13,-37],[-13,-115],[-39,-193],[-24,-154],[5,-474],[-6,-210],[-17,-116],[-28,-63],[-47,-19],[60,-19],[38,-29],[18,-46],[16,-71],[22,-34],[27,6],[24,-7],[22,-20],[65,-103],[72,-29],[4,-57],[-29,-320],[1,-41],[32,80],[35,248],[46,110],[37,24],[150,12],[160,-28],[60,-3],[53,15],[54,-31],[13,-31],[15,-130],[15,-74],[95,-265],[45,-148],[56,-232],[20,-66],[116,-308],[21,-84],[10,-64],[-4,-45],[-21,-70],[-38,-94],[-36,-72],[-33,-49],[-34,-37],[-35,-24],[2,-7],[39,12],[41,25],[76,63],[30,15],[83,9],[3,-23],[-40,-48],[8,-4],[58,38],[126,55],[493,177],[116,16],[166,-33],[135,-75],[147,-97],[154,-70],[244,-66],[71,-32],[143,-32],[67,-37],[79,-89],[127,-116],[96,-73],[106,-68],[107,-126],[173,-285],[43,-35],[106,-46],[200,-60],[295,-141],[129,-55],[84,-22],[84,-39],[84,-56],[63,-61],[44,-65],[39,-43],[68,-48],[35,-34],[3,-52],[-81,-201],[-2,-17],[82,145],[47,42],[36,19],[77,-2],[117,-24],[102,0],[88,23],[74,11],[62,-3],[45,7],[28,17],[33,0],[132,-50],[54,-2],[193,-48],[124,17],[21,-9],[43,-61],[37,-5],[61,9],[60,-16],[100,-81],[46,-70],[45,-142],[4,-42],[-83,-329],[-25,-128],[-4,-111],[17,-64],[71,-107],[13,-28],[42,-159],[11,-68],[-5,-78],[-22,-128],[6,-98],[19,-148],[-8,-100],[-32,-53],[-22,-51],[-18,-94],[0,-35],[17,-69],[32,-40],[51,-44],[48,-66],[88,-159],[63,-88],[78,-131],[15,-65],[-22,-43],[-27,-30],[-62,-40],[-28,-30],[7,-7],[90,22],[54,-2],[45,-32],[37,-63],[62,-52],[86,-44],[86,-73],[143,-174],[26,-40],[38,-91],[50,-141],[26,-93],[3,-44],[-29,-43],[-96,-80],[-97,-140],[30,7],[64,59],[109,116],[60,23],[55,-14],[88,-33],[78,-48],[68,-61],[100,-166],[106,-132],[59,-112],[-18,71],[-39,83],[-99,131],[-43,69],[-9,33],[-4,35],[8,66],[19,94],[25,67],[31,39],[21,41],[10,44],[17,32],[86,58],[23,-5],[18,-71],[20,-14],[40,-11],[34,-25],[27,-38],[20,-38],[11,-39],[26,-127],[18,-59],[3,72],[22,112],[16,46],[53,67],[-3,29],[-22,39],[-109,168],[-3,41],[30,24],[20,46],[10,68],[24,50],[69,69],[59,107],[30,73],[25,39],[24,9],[-36,30],[-7,20],[-1,88],[-16,92],[-23,43],[-67,93],[-11,27],[-12,105],[6,52],[18,44],[-10,43],[-62,77],[-24,67],[-27,162],[-24,204],[-28,150],[-31,96],[-7,60],[17,25],[23,76],[20,18],[31,-5],[1,9],[-48,41],[-23,51],[1,19],[42,51],[-8,21],[-31,28],[-95,47],[34,18],[22,40],[-4,12],[-38,17],[-43,31],[-33,46],[-41,73],[-24,58],[-25,97],[-42,110],[-18,27],[-21,18],[-24,8],[1,17],[26,25],[411,184],[34,26],[202,103],[93,58],[95,82],[129,85],[63,55],[40,53],[205,211],[87,107],[51,93],[73,111],[96,129],[60,110],[25,92],[32,159],[9,141],[5,207],[-3,184],[-26,289],[-15,91],[-30,108],[-71,217],[-12,60],[-45,100],[-144,260],[-181,176],[-34,44],[-72,52],[-109,59],[-70,48],[-184,182],[-60,21],[-24,47],[-5,33],[7,87],[11,59],[13,45],[15,28],[101,133],[57,106],[39,59],[44,45],[79,59],[44,75],[-10,31],[-35,33],[-8,32],[60,82],[8,23],[-7,77],[11,17],[74,5],[105,-111],[26,10],[-33,30],[-41,75],[6,31],[78,81],[2,37],[-22,51],[-3,40],[46,97],[-12,20],[-128,19],[-22,27],[8,13],[60,32],[5,12],[-107,217],[-17,64],[45,80],[51,36],[-6,20],[-68,4],[-41,11],[-41,61],[16,40],[15,17],[40,93],[38,19],[-7,16],[-146,-39],[-69,32],[-67,-8],[-32,9],[11,33],[123,151],[57,81],[35,73],[19,49],[2,24],[-13,162],[8,44],[47,42],[73,77],[-100,70],[-62,68],[-42,34],[-31,33],[-40,70],[-31,91],[-32,184],[-6,101],[8,74],[13,35],[22,37],[92,71],[161,105],[126,41],[90,-22],[178,-25],[143,-59],[434,-150],[77,-66],[-72,-56],[10,-14],[164,107],[43,20],[37,5],[124,-41],[49,-6],[63,-35],[150,-114],[10,10],[-42,58],[24,28],[118,60],[122,50],[86,48],[92,64],[62,35],[32,4],[41,-19],[109,-85],[71,-44],[56,-46],[80,-80],[30,-18],[61,-55],[79,6],[27,-7],[8,-11],[14,-37],[7,-24],[0,-25],[-18,-72],[-57,-114],[24,-2],[36,24],[47,44],[37,16],[79,-35],[73,-56],[26,-30],[29,-48],[24,-25],[22,-48],[-1,-16],[-21,-23],[-90,-39],[17,-15],[105,25],[32,24],[22,40],[31,12],[124,-71],[18,-25],[-8,-19],[-20,-22],[-54,-24],[-46,-61],[-8,-26],[33,-19],[81,-8],[-1,-14],[-46,-27],[-6,-35],[104,-123],[70,-55],[41,-9],[94,-3],[76,-21],[169,-70],[100,-13],[85,22],[57,4],[50,-28],[16,-20],[7,-39],[-1,-59],[28,-42],[56,-23],[45,3],[59,47],[50,8],[17,35],[15,64],[14,35],[37,10],[29,-19],[17,-29],[31,-89],[9,-39],[-3,-36],[-16,-32],[-31,-35],[-46,-37],[-36,-53],[-44,-122],[-16,-80],[-4,-47],[2,-51],[8,-56],[17,-47],[41,-65],[3,-22],[3,-53],[-4,-23],[-24,-46],[-66,-46],[-91,-12],[-297,-3],[-79,11],[19,-42],[83,-13],[76,0],[283,-26],[39,-27],[33,-47],[24,-49],[23,-99],[5,-47],[-13,-53],[-29,-57],[-20,-77],[-10,-95],[16,-51],[153,-4],[31,-35],[-4,-25],[-54,-97],[-5,-28],[24,-65],[-3,-19],[-15,-19],[-15,-49],[-13,-79],[-17,-51],[-42,-40],[-22,-9],[-16,11],[-41,108],[-17,16],[-16,-10],[-8,-17],[0,-23],[-7,-25],[-14,-25],[-61,-39],[-98,-27],[3,-29],[66,-15],[84,-48],[48,-8],[76,39],[146,122],[60,33],[53,13],[60,1],[66,-9],[133,16],[33,-12],[40,-27],[48,-42],[33,-41],[18,-39],[30,-148],[40,-38],[10,-29],[3,-43],[-3,-86],[-44,-171],[-22,-64],[-61,-87],[-70,-39],[-125,-37],[-65,-31],[-49,-43],[-3,-23],[142,73],[155,36],[45,40],[33,39],[34,86],[62,225],[35,70],[49,12],[22,-24],[50,-131],[0,-34],[-12,-28],[-83,-128],[30,13],[83,120],[18,38],[9,55],[27,39],[10,-19],[26,-142],[1,-103],[4,-33],[-9,-98],[10,-18],[25,85],[9,64],[11,46],[14,27],[103,86],[119,71],[78,62],[64,30],[97,28],[62,58],[28,87],[23,61],[19,33],[64,61],[35,3],[33,-22],[38,-45],[42,-68],[25,-54],[8,-38],[7,-136],[8,0],[38,104],[5,37],[-3,39],[-12,39],[-39,84],[-15,53],[4,33],[41,19],[61,7],[10,14],[-45,35],[-1,19],[43,62],[26,4],[49,-10],[-10,32],[1,21],[14,8],[82,-20],[10,24],[70,2],[7,22],[-61,31],[-60,20],[-18,16],[-14,25],[-19,61],[5,16],[17,0],[29,-16],[16,31],[17,75],[18,31],[55,-36],[2,16],[-44,117],[8,22],[68,9],[41,-15],[108,-86],[20,10],[-17,24],[-55,50],[-50,33],[-45,16],[-34,28],[-37,78],[-6,31],[3,42],[26,86],[15,19],[26,13],[38,9],[41,-8],[86,-57],[15,18],[-47,31],[-25,28],[-12,37],[5,41],[37,85],[17,69],[77,189],[23,35],[24,22],[15,22],[61,5],[112,-67],[34,-40],[10,-57],[-59,-77],[-101,-57],[-30,-25],[19,-14],[95,45],[82,22],[68,-1],[54,-91],[8,-126],[-31,-105],[41,52],[51,29],[42,-69],[5,-56],[22,-52],[48,-71],[50,-61],[-55,-65],[-65,-39],[13,-29],[90,-30],[12,-31],[-9,-41],[13,0],[62,64],[53,-9],[68,-138],[-50,-79],[-76,-36],[-60,-16],[-84,2],[-33,-11],[17,-27],[80,0],[123,20],[92,32],[39,2],[42,-13],[15,-11],[-45,-23],[-3,-8],[17,-23],[34,-75],[-3,-17],[-33,-43],[53,-10],[73,20],[23,-22],[45,-90],[28,-93],[-124,-126],[-63,-26],[-94,-67],[-26,-54],[-54,-69],[35,1],[96,110],[47,26],[35,-6],[15,-19],[-7,-31],[30,4],[130,64],[54,13],[71,4],[7,-22],[-43,-154],[-75,-119],[-137,-73],[-48,-43],[-60,-69],[23,-13],[130,91],[89,36],[125,29],[55,-4],[98,-182],[57,-17],[45,9],[87,-52],[32,-50],[-8,-37],[-29,-21],[-15,-34],[35,-101],[-20,-56],[-63,-50],[-45,-25],[-48,-8],[-46,-44],[-21,-7],[-64,11],[22,-27],[32,-14],[50,-7],[60,14],[57,-2],[91,-32],[38,-39],[1,-11],[-20,-23],[-28,-73],[-21,-26],[17,-20],[46,-31],[34,-11],[45,10],[47,-13],[161,-172],[-7,-89],[-24,-69],[9,-77],[1,-94],[-87,-27],[-289,45],[-165,68],[-8,20],[46,45],[-41,4],[-48,-18],[-20,-17],[55,-71],[152,-63],[68,-76],[74,-7],[23,-13],[41,-45],[-12,-15],[-76,-5],[-60,-52],[38,-30],[135,-27],[96,-6],[49,-30],[-40,-33],[-113,-39],[-4,-57],[84,-23],[75,14],[31,-6],[22,-140],[12,-29],[-80,-24],[0,-27],[53,-22],[89,-18],[29,-25],[6,-42],[19,-22],[50,-4],[57,52],[33,43],[48,-16],[3,-55],[58,-61],[21,-11],[16,-87],[47,78],[35,-17],[39,-4],[-14,-75],[-23,-60],[31,-37],[23,-55],[63,-76],[-17,-36],[-74,-78],[-40,-123],[-9,-43],[-38,-71],[-53,-69],[33,8],[118,126],[69,42],[154,23],[37,36],[56,14],[35,-40],[3,-73],[46,-23],[47,24],[44,-21],[-26,-46],[-140,-187],[-40,-75],[-12,-53],[48,73],[175,168],[18,25],[38,72],[36,47],[94,-17],[48,-34],[23,-94],[38,-102],[57,-114],[153,-55],[55,-9],[95,38],[15,53],[75,17],[52,-7],[18,-102],[56,-54],[55,-45],[54,-24],[78,-10],[42,-49],[0,-20],[-44,-53],[-42,-77],[-74,-54],[-103,-2],[-143,-34],[-5,-30],[-32,-34],[-76,-33],[-41,-25],[-67,-125],[-41,-53],[-47,-10],[-66,5],[-43,-12],[-32,-22],[-18,-34],[-14,-13],[-89,-34],[-160,-95],[-85,-3],[-53,11],[-40,-8],[-27,-26],[-77,-47],[-23,-28],[-13,-32],[-11,-67],[-9,-24],[-15,-14],[-64,14],[-72,43],[14,-45],[114,-78],[32,-43],[-30,-37],[-73,-58],[-8,-31],[29,-17],[-11,-26],[-40,-28],[4,-12],[5,-11],[98,40],[88,86],[57,87],[29,25],[113,32],[67,38],[96,69],[104,100],[115,131],[145,102],[177,74],[130,38],[82,1],[5,13],[-75,23],[-61,4],[-76,-16],[-24,40],[3,17],[25,30],[64,26],[314,-37],[108,-28],[118,-238],[27,-77],[8,-55],[-12,-35],[-47,-42],[-134,-82],[-19,-21],[-1,-12],[58,-17],[18,-22],[30,-91],[60,59],[114,144],[93,66],[78,19],[94,8],[32,-1],[12,-47],[49,-93],[45,-25],[87,-12],[79,-116],[29,-81],[29,-46],[-2,-32],[4,-26],[21,-40],[10,-34],[-5,-77],[-46,-134],[34,-122],[-16,-55],[-7,-87],[29,-59],[8,-34],[-24,-19],[-174,-49],[-68,-1],[-17,-29],[52,-9],[96,2],[115,-30],[51,-33],[22,-46],[-5,-37],[-34,-27],[-65,5],[-62,25],[4,-24],[93,-60],[27,-30],[50,-38],[10,-52],[-12,-52],[-175,-208],[-144,-132],[-145,-116],[-233,-223],[-23,-11],[-42,-4],[-110,36],[-89,-9],[-167,-43],[-46,-27],[-92,-77],[-36,-11],[-98,-16],[-94,11],[-37,-11],[-45,-38],[-12,-20],[-12,-65],[-227,-292],[-60,-99],[-116,-104],[-127,-183],[-111,-74],[-38,-102],[-106,-61],[-194,-16],[-93,-18],[-108,29],[-81,-44],[-122,-14],[-59,10],[-237,-98],[-60,93],[-46,36],[-134,6],[-106,38],[-98,8],[-95,17],[-62,-1],[-65,-10],[-101,3],[-56,-51],[-190,15],[-79,47],[-66,9],[-88,-10],[-84,-35],[-184,40],[-195,-33],[-171,22],[-47,22],[-269,-61],[-105,35],[-92,-93],[-64,20],[-69,-14],[-22,17],[-46,-12],[-30,-51],[-39,-5],[-65,-90],[-109,-72],[-159,-391],[-15,-150],[-60,-103],[-53,-13],[-43,-3],[-276,-75],[-123,-60],[33,-47],[-40,-35],[-65,-15],[-70,-43],[-46,-49],[-22,-68],[-142,-110],[-164,-255],[-78,-187],[-96,-135],[-67,-51],[-48,-8],[-49,16],[-81,63],[-59,7],[-149,88],[-345,89],[52,-33],[46,-55],[91,-14],[93,1],[193,-110],[94,-38],[58,-33],[49,-74],[-35,-146],[-36,-120],[-48,-92],[-166,-236],[-81,-80],[-140,-283],[-145,-133],[-78,-81],[-83,-129],[-194,-97],[-72,-25],[-66,13],[-81,-79],[-96,-48],[-28,-74],[-231,-197],[-88,-25],[-75,-53],[-23,-89],[-67,-54],[-18,-41],[-57,-125],[-105,-161],[-128,-27],[-47,-56],[-54,-91],[-76,-62],[-151,29],[37,-38],[135,-60],[14,-88],[-68,-21],[-141,-117],[-191,-202],[78,37],[161,147],[119,53],[156,154],[112,29],[21,34],[19,127],[10,46],[53,125],[63,106],[51,146],[92,93],[139,78],[129,171],[71,52],[69,37],[28,69],[43,40],[113,81],[125,21],[126,67],[97,36],[59,62],[87,33],[257,180],[72,85],[92,172],[81,88],[28,93],[117,152],[121,200],[59,143],[90,80],[174,227],[93,91],[38,10],[105,81],[66,84],[105,85],[190,104],[178,125],[241,108],[283,162],[228,86],[161,13],[195,40],[69,-4],[305,-70],[146,-87],[166,-182],[25,-48],[4,-68],[-88,33],[-78,2],[54,-37],[92,-113],[-4,-140],[-52,-127],[-155,-63],[-39,-49],[-32,-83],[-31,-31],[-75,-37],[-42,-53],[-121,-85],[-55,-10],[-63,20],[-152,80],[-93,77],[-47,-42],[-38,-44],[-90,15],[-41,-20],[-68,22],[-139,-97],[40,-11],[110,56],[37,-7],[82,-72],[196,-77],[51,-51],[48,-163],[33,-27],[67,17],[76,81],[63,43],[123,36],[-24,-54],[93,4],[93,-72],[-34,-51],[-47,-103],[-32,-201],[-95,-135],[-127,-132],[32,-32],[37,-20],[82,39],[54,-2],[61,-26],[-19,-102],[-22,-70],[13,-65],[36,-124],[49,-27],[20,-159],[26,-86],[-4,-70],[50,-44],[8,-71],[180,-20],[36,-28],[124,-27],[24,-19],[22,-39],[-122,-86],[99,-62],[93,-101],[74,20],[32,-3],[82,-63],[23,-32],[12,-28],[42,6],[59,25],[107,-6],[114,-36],[-9,-54],[-18,-38],[90,12],[56,-39],[19,19],[14,24],[111,66],[143,137],[17,-17],[6,-52],[19,-84],[55,-59],[65,-13],[89,45],[36,-39],[42,-75],[39,-97],[-2,-35],[-51,-30],[-47,-44],[193,-18],[20,-19],[21,-38],[-20,-39],[-18,-19],[-35,23],[-64,-21],[-56,-50],[-61,-28],[-38,-4],[-43,-23],[-39,-36],[-41,-10],[-126,-90],[-130,-57],[-135,-93],[-138,-59],[-144,-70],[-31,-6],[-36,3],[-82,-69],[-41,10],[-41,-12],[-48,15],[-32,28],[25,-73],[7,-66],[-12,-30],[-23,-34],[-82,6],[-33,25],[-38,35],[-18,58],[-41,41],[-25,-57],[0,-43],[-30,-58],[-36,99],[-65,-36],[-28,-105],[14,-30],[20,-81],[-32,-42],[-24,12],[-49,-118],[-60,-43],[-61,-121],[-73,-91],[-20,-62],[-122,-140],[-47,4],[-34,-5],[-51,-58],[-8,-118],[-23,15],[-23,-4],[-12,-37],[-17,-6],[-45,35],[-53,-19],[-42,27],[-52,173],[-28,61],[-50,19],[-13,-36],[-19,-36],[-48,71],[-37,266],[0,64],[51,223],[126,201],[-40,6],[-112,-139],[12,34],[19,35],[37,57],[57,53],[76,31],[52,5],[36,29],[52,52],[10,28],[-46,-32],[-77,-31],[20,41],[19,22],[410,359],[82,60],[165,75],[22,50],[-22,32],[64,-28],[-5,-41],[-10,-30],[-4,-51],[6,-49],[65,-24],[54,-91],[-26,124],[49,70],[188,93],[156,10],[50,44],[-134,29],[-158,-16],[-99,33],[-135,-21],[-144,20],[-44,-27],[-35,-58],[-47,25],[-22,5],[-22,20],[47,100],[145,150],[89,130],[25,28],[20,52],[-49,-9],[-43,-20],[-29,60],[-52,80],[-5,-34],[26,-99],[-101,-175],[-65,-12],[-86,-82],[-123,-71],[-143,-135],[-184,-115],[-39,0],[-84,94],[24,43],[22,58],[-21,-17],[-14,-25],[-50,-41],[41,-78],[-20,-29],[-59,-38],[-54,-56],[-48,-37],[-39,47],[-107,-60],[-89,-16],[-20,30],[-6,48],[-31,12],[-58,-13],[-23,25],[-3,-30],[16,-51],[11,-99],[-18,-45],[4,-59],[51,-17],[12,-18],[2,-22],[-111,-152],[-94,21],[-51,-40],[-53,-12],[-24,-67],[-29,-15],[-40,4],[-35,19],[-26,-9],[-37,-122],[-30,10],[-12,-44],[-16,-19],[-23,-16],[-21,54],[-12,52],[-19,11],[-25,13],[-26,0],[-17,-8],[-22,-33],[-31,-29],[-23,24],[-19,39],[-15,-62],[-23,-66],[4,-76],[-10,-45],[-22,12],[-21,40],[-61,32],[-48,-3],[10,42],[45,61],[-14,12],[-22,-9],[-10,9],[16,55],[2,61],[-21,-22],[-25,-64],[-62,-51],[2,-86],[-58,-175],[-3,-75],[-37,-59],[-48,-51],[-65,14],[-49,-45],[-25,-51],[-22,-7],[-11,65],[-8,20],[-18,-96],[-19,-6],[-7,68],[-8,45],[-26,-39],[-16,-103],[-18,9],[-5,38],[-13,12],[-4,-44],[6,-61],[-9,-33],[-17,17],[-18,30],[-29,-23],[-26,-9],[0,30],[5,37],[-53,-20],[-63,-68],[-50,-94],[17,-16],[20,-30],[-86,-146],[-87,-131],[-67,-214],[-26,-25],[-23,-39],[-24,-129],[-28,-115],[16,-51],[10,-53],[25,-52],[21,-5],[23,10],[16,-2],[11,-22],[-5,-27],[-26,-6],[-49,-47],[-43,-17],[-22,-56],[-32,-65],[-63,-100],[27,-31],[97,-35],[43,-36],[66,-189],[-15,-18],[-6,-35],[58,-48],[19,-135],[48,-46],[71,-28],[87,40],[73,57],[-3,46],[-45,107],[-11,50],[-34,33],[-13,-28],[-22,36],[-2,20],[20,10],[24,-4],[28,-19],[71,-116],[20,-154],[4,-97],[-8,-33],[-21,7],[-40,-7],[-187,-50],[-42,-45],[-96,-48],[-6,24],[7,50],[-6,102],[-18,5],[-149,-166],[-57,-11],[-49,-48],[-10,27],[-9,124],[30,104],[-16,-1],[-50,-63],[-22,39],[-11,42],[-15,24],[-17,9],[14,-93],[-34,-69],[-9,-179],[-43,-75],[-134,-48],[-88,11],[-78,-15],[-104,-35],[-58,21],[-59,-37],[-200,-9],[-42,19],[-54,-69],[-86,-40],[-218,-154],[-48,-56],[-58,-87],[-40,-47],[-32,-15],[-20,-39],[-21,-26],[21,87],[22,74],[19,142],[-5,115],[-24,48],[-24,31],[28,-113],[5,-140],[-10,-81],[-53,-158],[-23,-37],[-27,-32],[-20,-14],[-18,-25],[-22,-40],[-20,-79],[12,-72],[104,-27],[28,22],[15,-51],[8,-72],[-8,-78],[-18,-79],[-13,-98],[-11,-150],[-17,-134],[-2,41],[10,163],[-17,-17],[-11,-37],[-32,-212],[-44,-112],[-40,-78],[-42,13],[10,-62],[-12,-32],[-10,-67],[-24,-45],[-24,5],[-33,-31],[-13,-24],[-1,-45],[-23,-39],[-81,-206],[-70,-60],[-16,8],[18,97],[13,99],[-43,42],[-41,23],[-46,-3],[-52,76],[-66,55],[-93,151],[2,41],[-2,70],[28,110],[27,77],[38,40],[108,41],[27,61],[16,52],[-53,-89],[-81,-30],[-43,-33],[-35,-50],[-20,-64],[-47,-76],[3,-51],[8,-37],[-3,-76],[29,-74],[58,-121],[11,-188],[45,-126],[68,-147],[53,-42],[2,-54],[-24,-90],[-32,-42],[41,9],[21,-21],[20,-75],[-1,-77],[-7,-43],[-13,-18],[1,45],[-9,15],[-15,-19],[-9,-22],[-4,-86],[-10,-43],[-36,-13],[-36,-113],[-34,-64],[-132,-433],[5,-72],[-24,-23],[-36,-19],[-37,-44],[-25,-47],[-23,-129],[-43,-144],[-28,60],[-7,52],[12,134],[48,220],[52,137],[40,65],[32,132],[-41,20],[-63,-2],[12,61],[18,54],[-33,54],[-19,6],[-20,22],[23,45],[12,47],[-7,58],[10,43],[-17,-7],[-26,-46],[-16,-18],[-10,41],[-12,-10],[-7,-27],[-17,-16],[-35,38],[-52,44],[-29,75],[-16,58],[16,105],[36,19],[46,-17],[61,0],[-8,23],[-22,-4],[-64,86],[-21,51],[-35,14],[-17,-49],[-18,-13],[23,108],[29,4],[42,30],[-12,63],[-27,28],[-49,-35],[1,44],[9,57],[37,0],[32,-19],[27,91],[2,41],[-46,-59],[-10,128],[45,123],[42,53],[54,-1],[54,9],[-34,22],[-35,12],[27,49],[22,9],[22,42],[-53,-6],[6,80],[-26,-16],[-30,-8],[-12,-34],[2,-56],[-9,-37],[-24,-30],[-40,-24],[-4,41],[-14,18],[-5,-86],[-10,-30],[-30,81],[-9,-16],[1,-23],[-7,-40],[-26,-20],[2,-51],[-10,-28],[-81,44],[-2,-15],[46,-95],[33,-33],[4,-52],[-28,-43],[-40,37],[-7,-3],[22,-63],[13,-56],[-14,-47],[3,-58],[-3,-52],[-9,-45],[19,-211],[24,-57],[22,-54],[13,-51],[-25,-8],[-38,42],[-34,32],[-41,103],[-7,41],[-9,32],[4,-74],[15,-83],[127,-186],[23,-71],[18,-56],[-5,-54],[-33,38],[-28,49],[-76,54],[-95,35],[-54,127],[0,-53],[-12,-45],[-33,55],[-21,46],[-7,52],[-41,-4],[-43,-44],[-41,10],[-5,87],[11,46],[47,109],[44,56],[19,72],[-6,111],[-9,-113],[-25,-57],[-39,-42],[-53,-77],[-12,-70],[-16,-133],[22,-45],[22,-11],[66,30],[35,-14],[76,-159],[142,-63],[52,-39],[42,-83],[64,-48],[49,-70],[2,-45],[-18,-54],[-6,-72],[-21,-46],[-51,-5],[-30,11],[-163,256],[-20,23],[-60,134],[-71,71],[-22,-1],[101,-133],[41,-92],[73,-130],[52,-55],[38,-86],[36,-40],[97,-57],[-34,-41],[54,-35],[8,-65],[-5,-73],[-75,29],[-3,-54],[7,-32],[-33,-27],[-46,36],[-119,196],[1,-26],[10,-31],[69,-126],[61,-75],[53,-34],[40,-64],[14,-38],[10,-58],[-30,-39],[-34,-22],[-33,39],[-25,42],[-52,70],[-15,79],[-40,-4],[-165,99],[-132,12],[13,-20],[16,-13],[106,-25],[42,-46],[87,-41],[51,-11],[20,-125],[70,-86],[10,-64],[48,-7],[84,62],[55,-22],[78,-18],[18,-50],[14,-97],[27,-108],[73,-427],[108,-349],[13,-60],[-25,53],[-80,231],[-45,167],[-45,295],[-13,66],[-16,26],[-10,-21],[-5,-38],[8,-29],[-18,-97],[8,-44],[29,-46],[31,-115],[26,-155],[-34,63],[-37,33],[-57,26],[-50,44],[3,-64],[-5,-69],[-39,21],[-26,23],[23,-74],[-51,22],[-34,-4],[-22,-66],[-29,-39],[-44,-13],[-65,60],[-21,72],[-9,81],[-4,-95],[12,-100],[-4,-76],[63,-14],[58,13],[79,-3],[52,14],[31,24],[74,-21],[5,-92],[-8,-90],[-5,-97],[21,0],[24,31],[12,174],[68,64],[23,-5],[22,-55],[7,-57],[8,-77],[-16,-119],[-105,-138],[-74,-128],[-39,-26],[-55,15],[-62,32],[-31,7],[-23,-11],[-14,39],[-10,72],[-24,24],[-18,-3],[-13,-76],[-58,-22],[-80,32],[-83,64],[36,-69],[206,-128],[23,-24],[22,-36],[-29,-55],[-22,-62],[-4,-48],[-8,-31],[-82,-83],[-44,15],[-114,149],[52,-129],[41,-55],[84,-29],[157,48],[51,-53],[-42,-93],[-42,-66],[-55,-7],[-49,-18],[-14,-45],[-34,-3],[-54,-2],[-84,-4],[-46,10],[-64,-92],[-24,-13],[-34,18],[-14,74],[-15,36],[0,-138],[5,-38],[13,-28],[-75,-75],[-72,-94],[-26,-25],[-29,-47],[-60,-135],[-15,-99],[-21,-110],[-3,49],[4,84],[-15,95],[-10,-175],[-23,-81],[-213,5],[-94,-44],[-143,-149],[-43,-65],[-119,-252],[-30,-162],[-24,68],[6,51],[0,42],[-30,-89],[29,-130],[-26,-50],[-78,-93],[-43,-15],[-48,-26],[-15,-92],[-65,-84],[-38,-37],[-70,22],[21,-81],[-25,-61],[-44,-48],[-55,-30],[-31,3],[-27,-16],[-21,-39],[-52,-36],[-53,20],[-61,12],[-33,-22],[56,-36],[31,-52],[-7,-71],[-15,-27],[-35,-37],[-16,5],[-10,33],[-11,70],[-17,-15],[-3,-32],[-14,-12],[-51,111],[3,-85],[17,-64],[18,-33],[17,-20],[4,-30],[-35,-73],[-17,-17],[-32,-12],[-18,-45],[5,-39],[-28,-84],[-66,-53],[-20,2],[-17,-15],[10,-38],[16,-27],[0,-26],[-18,-34],[-34,-10],[-20,-39],[6,-38],[12,-20],[-3,-36],[-39,-36],[-9,-35],[19,-11],[15,11],[11,-8],[-23,-59],[-21,-36],[-21,-65],[-46,-18],[1,-20],[27,-19],[22,-50],[-42,-92],[-26,8],[-15,20],[-11,-72],[4,-39],[-10,-79],[-15,-95],[-11,-39],[2,-73],[7,-70],[25,-91],[39,-370],[27,-128],[47,-347],[79,-336],[111,-406],[183,-493],[22,-70],[-24,-59],[-7,-62],[-2,-93],[6,-90],[22,-111],[41,-169],[-23,34],[-60,242],[-7,143],[9,202],[-15,-5],[-11,-66],[-6,-77],[-15,-30],[-21,118],[1,53],[22,62],[-6,23],[-36,32],[-7,50],[4,49],[-20,26],[-16,-1],[11,-122],[17,-74],[21,-180],[33,-108],[20,-91],[231,-972],[54,-124],[20,-89],[21,-186],[5,-238],[-37,-436],[-9,-298],[-5,9],[-4,31],[-9,5],[-32,-137],[-45,-122],[-15,-192],[-21,-95],[-64,-101],[-40,2],[-97,-76],[-68,20],[-82,-43],[-53,5],[-31,90],[5,41],[12,40],[21,10],[72,-95],[13,40],[-21,47],[-42,27],[-31,29],[-62,215],[-64,149],[-11,99],[-110,60],[-80,91],[-52,163],[-30,288],[-36,33],[-15,22],[35,107],[36,90],[-29,-22],[-21,-34],[-27,-79],[-19,-12],[-19,12],[-20,152],[6,187],[29,70],[-45,2],[-46,-27],[6,-62],[-6,-34],[-34,8],[-26,22],[-35,65],[-47,124],[-97,340],[-19,48],[-33,51],[16,15],[27,10],[63,153],[49,93],[16,64],[-3,27],[-21,40],[-29,-35],[-12,10],[-32,81],[-31,22],[-21,-17],[22,-66],[21,-24],[-8,-96],[-8,-31],[-19,-28],[-30,15],[-15,-24],[-18,25],[-17,42],[-20,69],[52,390],[48,249],[5,283],[4,42],[-4,76],[-64,163],[-284,400],[-219,473],[-191,178],[-144,-39],[-25,-36],[-11,-47],[9,-53],[-13,-21],[-39,2],[-52,-12],[-136,-125],[-48,5],[-44,-32],[-32,-24],[-86,-14],[-72,-27],[-31,15],[-20,72],[0,75],[16,-58],[26,-44],[11,18],[5,39],[-26,78],[-82,100],[-93,146],[29,-5],[7,31],[-29,41],[12,47],[20,50],[-39,-7],[-35,-35],[-1,-43],[-7,-34],[-19,5],[-36,42],[-173,118],[-152,66],[116,30],[63,-23],[-7,36],[-15,22],[-50,29],[-64,-11],[-40,14],[-41,-29],[-45,-42],[-40,-22],[-156,-30],[-127,-33],[20,34],[22,23],[75,34],[11,71],[-18,68],[-19,-16],[-21,-54],[-26,39],[-28,0],[-7,-85],[-37,-57],[-16,-57],[-106,-45],[-13,15],[31,54],[-3,30],[-35,-26],[-59,-103],[-209,-34],[11,24],[44,4],[62,33],[-12,55],[-24,61],[-22,6],[-15,36],[1,111],[-14,66],[-34,67],[-11,-13],[-25,-115],[-21,-151],[-11,-48],[-61,-4],[-55,11],[-186,-18],[-70,51],[-29,9],[-17,-1],[-81,-47],[-92,-35],[-22,11],[-31,2],[-67,-122],[-79,-57],[-199,102],[-49,83],[-44,17],[-54,10],[-58,-101],[-44,-137],[70,-75],[59,-36],[99,30],[54,67],[45,-3],[21,14],[19,35],[38,-28],[2,-27],[-27,-39],[-34,-32],[-21,-39],[39,-77],[61,-26],[23,11],[14,87],[37,56],[51,-12],[-7,-35],[7,-33],[24,-57],[-3,-81],[5,-19],[-55,-36],[-41,-12],[-33,-47],[17,-27],[-33,-24],[-23,9],[-11,-9],[-3,-28],[-18,-27],[25,-80],[52,-53],[36,-66],[146,-86],[35,2],[35,-87],[28,-30],[27,-16],[-3,-60],[-48,-44],[-13,-52],[-12,-29],[-22,37],[-22,27],[-51,-82],[-25,-18],[12,89],[-19,35],[-30,89],[-42,55],[-31,18],[-23,35],[-28,14],[-25,-4],[-41,21],[-2,47],[-12,35],[-32,42],[-153,79],[-1,-33],[11,-24],[22,-17],[26,-31],[0,-95],[-12,-40],[-4,-57],[-11,-58],[-18,-45],[-42,-31],[-19,26],[-30,124],[-42,40],[-67,4],[-45,-28],[-50,-121],[-40,-19],[-137,62],[-156,95],[4,32],[25,10],[47,-13],[-3,33],[-48,106],[-9,48],[6,59],[-15,-1],[-29,-49],[-100,41],[-27,50],[-59,141],[-83,4],[-37,85],[-68,-35],[-34,-40],[-30,-61],[12,-32],[30,-50],[-14,-24],[-96,-36],[-223,40],[-65,37],[-88,80],[-121,64],[-59,11],[-57,-13],[-167,-7],[-38,-17],[-33,-27],[-22,30],[-10,54],[20,9],[21,32],[20,63],[2,38],[-14,25],[-26,3],[-57,-165],[33,-92],[-2,-33],[-114,-19],[-258,-185],[-101,-101],[5,34],[122,130],[-43,20],[-69,-33],[-25,13],[29,107],[-9,94],[-49,3],[-32,-75],[-21,3],[-29,32],[-22,-10],[16,-171],[31,-70],[26,-90],[-71,-111],[-65,-92],[-7,-88],[-66,-115],[-62,-65],[-146,-154],[-42,-33],[-66,-71],[-90,-53],[-88,-85],[-29,-13],[56,72],[66,71],[-57,-10],[-87,33],[-53,2],[-1,-26],[-40,-37],[-42,54],[-19,36],[-8,31],[-18,8],[-17,-15],[62,-219],[27,-10],[30,-22],[-37,-51],[-40,-39],[-62,-25],[-53,80],[-11,-101],[-7,-100],[-18,-26],[-28,-37],[-16,28],[-7,39],[-18,-35],[-27,-26],[-43,-5],[-33,-14],[0,-41],[8,-42],[58,33],[-21,-108],[-53,-106],[-44,-25],[-67,15],[-16,-10],[-15,-22],[78,-167],[-50,-250],[-32,-91],[-22,-12],[-24,-3],[-86,81],[-47,63],[41,-170],[113,-50],[6,-64],[-1,-55],[-22,-65],[-21,-85],[15,-60],[18,-148],[15,-67],[17,-206],[18,-89],[102,-328],[35,-3],[6,-35],[-4,-68],[-10,-206],[-33,-169],[-108,-352],[-45,-218],[-87,-623],[-27,-409],[-6,-192],[-8,-28],[8,-28],[-21,-425],[12,-363],[-8,-56],[-31,-110],[-23,-151],[9,-68],[0,-46],[32,-231],[11,-173],[97,-297],[54,-103],[68,-90],[25,-50],[-12,-126],[-28,-64],[-12,-101],[-13,83],[5,109],[22,58],[0,43],[-40,46],[-72,149],[-89,264],[63,-422],[22,-68],[17,-28],[28,-29],[7,-46],[1,-34],[78,-294],[88,-303],[4,-83],[36,-102],[224,-425],[138,-318],[47,-302],[29,-94],[14,-129],[90,-145],[28,-93],[48,-52],[39,-156],[73,-91],[-17,-2],[-64,46],[4,-30],[54,-65],[106,-64],[26,0],[-42,33],[-35,45],[12,6],[74,-56],[207,-18],[90,-130],[118,-56],[64,-165],[73,-174],[47,-8],[37,-1],[110,30],[172,109],[59,53],[116,72],[176,14],[55,-21],[132,45],[63,56],[21,50],[11,37],[122,53],[24,11],[124,9],[60,20],[70,11],[51,-78],[0,-38],[-34,-36],[16,-35],[60,-58],[111,-23],[35,9],[51,85],[90,83],[-2,96],[-16,52],[-26,4],[-6,29],[18,71],[-5,23],[-56,-67],[-12,0],[6,32],[11,25],[164,147],[42,61],[56,53],[117,200],[26,377],[23,66],[78,116],[8,34],[5,79],[-2,198],[4,156],[-3,176],[14,157],[13,44],[44,251],[93,111],[161,132],[37,22],[511,139],[72,34],[88,90],[64,31],[118,-3],[36,12],[7,13],[0,12],[22,12],[67,-12],[127,-54],[46,-14],[114,-63],[125,-25],[17,11],[17,19],[13,37],[-13,32],[-13,-3],[-21,-17],[-25,-2],[-49,27],[10,22],[50,-2],[32,11],[48,39],[51,-29],[67,-130],[48,-41],[3,-187],[8,-34],[17,-49],[-24,-145],[-26,-120],[-34,-99],[-72,-154],[-88,-124],[-109,-275],[-25,-129],[1,-104],[18,-100],[-5,-37],[-14,-37],[-22,3],[-43,-47],[-55,-143],[-1,-43],[24,-40],[32,19],[41,2],[23,12],[24,-3],[-5,-82],[-27,-57],[-16,-19],[-29,-9],[-33,-33],[-17,-31],[1,-94],[19,-7],[38,70],[23,-3],[4,-30],[-50,-241],[-32,-247],[-44,-144],[-15,-208],[-23,-89],[-27,-88],[-16,5],[-42,167],[-43,43],[-10,40],[25,202],[-11,112],[-23,-5],[-29,-61],[-38,-53],[1,-77],[-44,-128],[-11,-42],[-29,-113],[29,-15],[27,11],[63,-4],[25,-125],[-7,-104],[-60,-275],[-7,-94],[-28,-142],[37,-93],[-34,-123],[-12,-80],[-3,-120],[18,-229],[-28,-330],[-50,-144],[-31,-54],[-55,-143],[-72,-43],[-100,-231],[-18,-60],[10,-65],[29,-22],[23,-6],[49,-56],[57,-42],[37,85],[-20,51],[-14,26],[2,23],[199,-221],[53,-27],[41,63],[24,22],[57,76],[17,16],[94,32],[46,-1],[40,-77],[31,-43],[60,37],[49,8],[206,-72],[82,32],[150,6],[68,-17],[95,101],[61,21],[72,47],[-9,49],[-17,22],[109,-21],[164,-103],[174,18],[62,56],[41,16],[178,-106],[47,-82],[37,-8],[28,19],[8,17],[-35,17],[-16,26],[140,-50],[264,-385],[6,-31],[-113,113],[-60,-8],[-15,-19],[3,-62],[6,-29],[25,-3],[19,16],[47,-20],[31,-41],[37,-63],[23,-69],[24,-1],[24,41],[45,5],[29,-46],[20,3],[-28,71],[-68,72],[16,3],[150,-129],[43,-160],[36,-37],[37,-49],[-16,-37],[-16,-24],[-35,-119],[-12,-11],[-2,89],[-21,11],[-24,-31],[-14,-45],[21,-59],[19,-1],[22,-17],[61,-408],[-13,-72],[-38,-114],[-35,-96],[-36,-60],[-44,-258],[-40,-417],[28,-376],[-14,-347],[13,-82],[3,-102],[-29,-18],[-17,3],[-17,63],[2,55],[18,64],[7,88],[-8,46],[-18,-100],[-30,-45],[-20,-15],[-20,-51],[21,-95],[26,-69],[9,-50],[-10,-60],[-6,-203],[-9,6],[-10,27],[-28,2],[-3,-81],[2,-46],[-24,-35],[-9,-36],[20,-24],[22,-15],[26,3],[22,-100],[7,-82],[-51,-75],[-17,-63],[-29,-76],[-16,-74],[-5,-53],[20,-170],[35,-120],[29,-76],[39,-16],[14,-40],[15,-62],[7,-81],[69,-268],[56,-150],[120,-273],[53,-50],[88,-220],[30,-37],[18,-65],[91,-53],[26,-39],[34,-53],[70,-95],[5,-47],[-7,-44],[19,-128],[37,-18],[37,24],[9,-23],[-7,-24],[-19,-26],[-5,-110],[61,-52],[30,-45],[100,22],[36,-12],[26,13],[-28,88],[-38,65],[3,29],[29,-21],[21,-45],[50,-55],[90,-191],[104,-47],[82,7],[77,25],[122,75],[88,134],[70,60],[229,128],[81,133],[34,18],[32,17],[72,101],[39,78],[41,40],[120,-29],[78,-37],[54,5],[52,-26],[23,-58],[24,-24],[36,2],[-172,-419],[-36,9],[-106,-14],[-35,-36],[-25,-45],[-11,-21],[-62,-52],[-24,-76],[-15,-63],[-31,-73],[35,-43],[-204,-246],[-40,-35],[-92,-25],[-21,-27],[-28,-48],[-4,-74],[5,-63],[26,-48],[24,-31],[57,-146],[65,-119],[-8,-21],[-486,98],[-8,90],[-16,-23],[-35,-22],[-8,-31],[-227,45],[-47,147],[-10,56],[-18,66],[-72,27],[-61,50],[-68,8],[-34,-28],[-35,35],[-6,81],[-70,-36],[-91,12],[-81,35],[-55,-19],[-46,-57],[6,-146],[-13,-29],[-37,111],[-51,106],[-45,65],[-3,152],[-18,82],[-67,76],[-58,54],[-42,-11],[26,-88],[67,-112],[5,-43],[-1,-57],[-47,8],[-41,24],[-49,7],[-34,35],[-70,134],[50,114],[15,75],[-1,155],[-11,76],[-55,114],[-86,126],[-121,103],[-57,83],[-141,63],[-54,43],[-42,78],[-6,56],[15,86],[-39,110],[-169,216],[-94,79],[-21,47],[-15,14],[15,-149],[41,-89],[108,-84],[29,-49],[12,-63],[-62,-122],[-32,-31],[-9,-66],[-21,-20],[-21,38],[-87,191],[-169,91],[-31,56],[-63,174],[-28,159],[10,106],[69,165],[22,71],[-5,45],[3,65],[-26,45],[-64,59],[-41,48],[11,24],[73,63],[5,58],[0,19],[-1,27],[-45,110],[-72,132],[-276,406],[-102,243],[-54,174],[-52,91],[-149,186],[-35,74],[-147,249],[-113,146],[-1,62],[46,77],[22,-3],[25,-56],[40,-63],[19,0],[28,29],[1,29],[2,106],[-44,42],[-24,88],[-22,59],[7,36],[-2,42],[-62,33],[-57,-26],[-33,17],[-25,22],[-13,-13],[10,-100],[-32,-61],[-28,-44],[-51,-12],[-85,-4],[-129,49],[-95,68],[-51,0],[16,-22],[41,-14],[53,-48],[-15,-13],[-193,99],[-224,194],[-133,31],[-153,51],[-90,123],[-68,53],[-209,164],[-70,28],[-293,-3],[-126,64],[-143,125],[-97,113],[-226,316],[-16,23],[-145,271],[-151,298],[-60,98],[-57,74],[-78,137],[-203,302],[-105,138],[-99,164],[-89,92],[-87,60],[-39,34],[-34,45],[-19,3],[-10,-61],[31,-31],[37,-26],[29,-2],[30,-22],[90,-83],[15,-43],[-256,167],[-104,15],[-11,27],[52,87],[-16,32],[-19,4],[-55,-61],[-21,-4],[-3,40],[2,38],[-36,55],[-22,-3],[-19,-38],[-48,-73],[1,-28],[95,-31],[33,-17],[-7,-20],[-82,0],[-101,-33],[-179,-201],[-167,-87],[-241,-194],[-106,-10],[-55,-31],[-162,74],[-205,182],[-310,58],[-209,239],[-207,98],[-132,230],[-79,9],[-51,37],[-188,83],[-186,55],[-182,200],[-119,64],[-103,80],[-226,136],[-84,75],[-79,117],[-129,121],[-55,101],[-62,36],[-88,190],[-46,80],[-39,38],[-42,13],[-120,-16],[-180,84],[-83,21],[-173,124],[-230,138],[-76,159],[-64,148],[-116,195],[-73,84],[-126,99],[-69,80],[-108,61],[-182,156],[-58,134],[-34,119],[-97,144],[-107,269],[-27,99],[-21,152],[-25,88],[-29,69],[15,51],[54,62],[90,13],[64,67],[8,55],[-4,35],[-41,84],[-50,22],[-39,2],[-11,32],[30,35],[34,82],[50,100],[35,92],[7,131],[-5,130],[14,110],[-121,128],[-14,54],[-39,145],[-66,170],[2,338],[-79,300],[-83,149],[-42,52],[-116,231],[-91,135],[-90,253],[-88,161],[-112,270],[-82,134],[-369,455],[21,0],[108,-111],[18,9],[3,54],[-13,65],[-20,16],[-29,-15],[-39,14],[-19,21],[-58,14],[-74,76],[-31,78],[-3,90],[-105,191],[-40,107],[21,-8],[27,-43],[29,-12],[33,-1],[23,18],[-8,31],[-23,27],[-152,100],[-51,71],[-125,117],[-30,42],[-19,113],[-31,5],[-27,-31],[-74,-31],[-18,36],[-3,35],[52,37],[48,107],[-1,37],[-27,-43],[-40,-48],[-41,-27],[-61,-22],[-31,16],[-27,24],[-44,94],[-22,305],[38,105],[46,105],[38,62],[23,-47],[23,-6],[-17,53],[-37,50],[-14,49],[-1,45],[-18,84],[-113,176],[-105,-13],[-42,6],[-40,67],[-34,113],[-17,94],[-2,51],[-10,50],[-182,83],[-54,71],[-55,89],[-22,73],[-23,54],[-17,94],[-13,110],[21,141],[26,69],[-125,50],[-48,3],[-40,-29],[-35,37],[-74,41],[-87,148],[-103,268],[-113,87],[-37,93],[-46,84],[-41,104],[-6,45],[-12,26],[-57,72],[-63,123],[-17,99],[-17,151],[-42,53],[-42,25],[-8,72],[2,40],[-14,76],[-86,183],[-44,150],[-24,47],[-22,68],[-11,140],[-36,174],[-69,210],[-58,141],[-28,142],[12,144],[-9,89],[-8,21],[6,29],[19,-16],[16,26],[-2,92],[-20,28],[-55,28],[-25,20],[-136,38],[-77,52],[-5,122],[-37,55],[-32,34],[-103,72],[-16,-37],[-14,-67],[-42,-15],[-37,-3],[-63,47],[-154,179],[-33,29],[-48,15],[-24,29],[-105,94],[21,-50],[30,-52],[27,-151],[-31,-109],[-19,-381],[21,-76],[45,-121],[31,-193],[10,-144],[29,-114],[-9,-268],[10,-82],[44,-134],[80,-125],[16,-66],[105,-96],[64,-125],[127,-169],[40,-72],[114,-265],[4,-79],[21,-96],[64,20],[29,-74],[-3,-34],[7,-25],[34,5],[29,-20],[61,-286],[32,-39],[42,-17],[47,-32],[3,-73],[-2,-59],[40,-84],[-7,-115],[32,-97],[-5,-94],[13,-75],[99,-168],[122,-134],[25,-175],[49,-161],[51,-40],[51,-67],[-6,-69],[4,-43],[68,-126],[11,-162],[59,-105],[17,-9],[13,17],[-43,106],[-21,69],[-3,109],[15,14],[123,-172],[13,-128],[41,-72],[3,-98],[25,-59],[7,-84],[40,-140],[1,-193],[21,-142],[77,-216],[63,-46],[11,-106],[65,-276],[75,-149],[38,-127],[5,-79],[-28,-117],[-3,-81],[41,-248],[61,-127],[68,-31],[12,-18],[-6,-34],[23,-30],[26,39],[13,55],[-13,64],[-3,47],[12,33],[22,5],[131,-169],[22,-66],[48,-76],[45,-94],[19,-73],[36,-61],[18,-142],[91,-64],[48,-117],[3,-75],[-23,-190],[-21,-55],[-73,-81],[-53,-97],[-52,-59],[-54,-36],[-45,8],[-44,111],[-51,335],[-35,71],[-24,105],[-40,87],[-145,132],[-74,140],[-71,93],[-76,134],[-209,224],[-86,113],[-57,113],[-37,-2],[-28,-11],[-11,31],[-1,58],[-13,38],[-122,170],[-25,93],[-5,109],[27,281],[12,164],[-7,84],[-13,11],[-12,46],[-3,135],[-29,147],[-111,301],[-81,59],[-71,43],[-198,266],[-53,133],[-13,76],[-7,153],[-27,-89],[-36,-65],[-84,4],[-94,-74],[-56,70],[-30,79],[-46,96],[-51,18],[-33,3],[-62,118],[-50,37],[-71,15],[-61,60],[-19,65],[-10,93],[-23,54],[-95,109],[-78,119],[-73,77],[-22,62],[-2,43],[116,-12],[138,-46],[66,13],[41,45],[39,31],[7,-34],[-7,-63],[37,-57],[52,-47],[37,4],[-36,51],[-21,104],[9,37],[0,52],[-51,-12],[-8,27],[44,78],[51,208],[24,204],[-53,179],[-89,125],[-193,363],[-115,187],[-34,70],[-30,33],[-94,42],[-79,105],[-138,148],[-59,76],[-41,179],[-32,24],[10,124],[-14,219],[-24,56],[-74,55],[-18,148],[-4,142],[-14,98],[-128,166],[-7,81],[0,76],[-13,75],[-68,158],[-82,138],[-29,66],[-3,133],[-29,36],[11,9],[27,-3],[15,17],[2,93],[-124,146],[-36,201],[-65,106],[-16,39],[-35,189],[-1,6],[-4,110],[-25,38],[-33,-24],[-15,142],[8,67],[-4,66],[-30,161],[-81,195],[-174,242],[-89,81],[-70,102],[-44,29],[-56,8],[-16,-46],[-64,31],[10,114],[-62,159],[-50,18],[-127,-11],[-169,87],[-50,52],[-18,94],[-79,81],[-105,79],[-58,-18],[-76,12],[-109,57],[-63,7],[-124,-17],[-46,12],[-42,72],[-47,36],[10,89],[-6,80],[7,63],[-21,137],[16,127],[-13,46],[-26,35],[-82,52],[-15,65],[13,90],[-21,60],[-67,55],[-63,126],[-79,70],[-33,115],[-49,72],[-17,63],[-108,226],[-116,177],[-18,102],[-4,139],[45,85],[25,74],[-3,69],[-7,50],[-40,88],[-154,51],[-125,217],[-8,165],[-49,169],[-1,110],[-7,119],[37,26],[34,-10],[-4,-47],[11,-85],[40,-64],[37,-28],[34,-62],[26,-19],[26,-4],[-14,40],[-15,25],[-19,83],[-34,105],[-41,58],[-20,106],[-17,25],[-11,39],[39,47],[52,33],[71,9],[201,-16],[43,28],[36,-9],[26,3],[-55,28],[-31,-9],[-36,6],[-72,-6],[-29,12],[-32,33],[-21,4],[-66,-58],[-30,7],[-70,63],[-30,9],[-50,-36],[-6,-155],[16,-115],[-30,-12],[-34,48],[-53,28],[-43,43],[-61,80],[-32,30],[-37,-67],[-1,31],[18,77],[-5,130],[55,-104],[-17,74],[-43,80],[-33,28],[-40,144],[-92,86],[-73,140],[-151,231],[-10,203],[-54,257],[23,146],[-3,104],[-27,156],[-28,85],[-122,234],[-117,157],[-18,119],[-8,119],[25,107],[23,112],[16,30],[6,-12],[-4,-24],[16,-8],[6,50],[10,25],[-17,4],[2,15],[10,32],[36,147],[-3,186],[39,228],[-2,75],[-25,162],[-25,97],[-44,69],[19,100],[-1,96],[-78,138],[-30,181],[-6,76],[8,203],[-21,86],[-52,142],[23,124],[24,75],[58,329],[14,27],[25,-1],[43,56],[-20,13],[-30,-27],[27,130],[30,111],[19,40],[10,364],[17,277],[28,92],[-10,94],[12,128],[-8,129],[60,622],[-8,75],[18,101],[-17,265],[8,297],[-16,38],[-7,41],[14,6],[28,-43],[129,0],[83,40],[30,-13],[35,-54],[44,-11],[55,9],[-17,14],[-26,3],[-57,50],[-33,50],[-101,-3],[-21,32],[-112,-32],[-35,33],[-62,-21],[15,93],[-3,117],[4,115],[15,-84],[38,-88],[18,100],[13,127],[-37,48],[-62,36],[-22,118],[146,100],[-78,21],[-30,45],[-38,6],[-3,-35],[-12,-45],[-13,60],[-4,71],[-15,122],[-60,196],[-37,253],[-45,125],[-88,120],[-23,70],[-21,177],[12,134],[-16,94],[42,-5],[110,-74],[138,-58],[42,-43],[67,-31],[370,-49],[25,5],[48,30],[20,-4],[54,-69],[28,-8],[35,3],[26,14],[45,47],[6,-17],[-1,-44],[16,-63],[33,-81],[12,-51],[-66,-142],[-13,-3],[-2,48],[-8,9],[-125,-240],[-44,-114],[-4,-52],[1,-30],[18,-8],[40,12],[59,48],[3,10],[-55,-17],[-27,-1],[4,54],[6,25],[36,80],[38,48],[54,51],[31,42],[21,61],[60,74],[11,20],[-3,60],[4,12],[29,-8],[12,-104],[-7,-46],[-52,-57],[-6,-20],[9,-77],[-8,-7],[-20,9],[-6,-5],[49,-84],[15,-65],[2,-58],[-13,-112],[-14,-18],[-24,7],[-32,35],[-8,-12],[-25,-86],[-9,7],[-15,103],[-9,8],[-50,-47],[-20,-45],[-17,-72],[-21,-34],[61,-7],[56,15],[45,-35],[15,-1],[40,34],[13,22],[34,109],[16,19],[26,1],[24,16],[36,60],[1,23],[-13,133],[4,75],[-7,24],[-16,25],[2,24],[13,40],[1,36],[-11,31],[5,37],[34,78],[6,34],[42,77],[-11,32],[-30,38],[-20,33],[-19,53],[-15,17],[-5,-8],[21,-86],[-4,-5],[-54,46],[-13,29],[-6,40],[4,29],[29,30],[35,10],[-3,25],[-43,80],[-29,37],[-22,18],[-30,4],[-14,13],[-3,19],[6,25],[16,8],[46,-10],[25,17],[-2,32],[-8,18],[2,114],[-18,92],[-9,16],[-10,2],[-11,-13],[-29,-3],[-18,30],[-19,59],[-37,140],[-20,35],[-53,46],[-21,0],[-22,-14],[-13,-22],[-12,-45],[-8,-16],[-7,3],[-5,13],[-17,63],[4,28],[18,34],[-2,12],[-38,-10],[-17,9],[-8,18],[3,72],[-21,41],[25,17],[63,14],[65,1],[19,30],[18,76],[-46,-70],[-28,-7],[-86,26],[-56,-4],[-7,15],[2,16],[12,15],[9,53],[13,147],[18,54],[6,29],[-4,7],[-75,-103],[-6,-32],[7,-28],[-7,-58],[-34,-17],[-21,9],[-39,-49],[-12,-5],[-178,86],[-18,11],[-31,40],[-44,68],[-13,59],[16,50],[17,24],[18,0],[20,-18],[41,-81],[16,-50],[43,7],[68,64],[17,24],[-69,-25],[-29,2],[-30,26],[-31,52],[-16,58],[0,186],[10,36],[34,26],[21,44],[-2,19],[-20,38],[-30,29],[-29,15],[-7,-5],[45,-81],[-1,-30],[-50,-80],[-8,-23],[0,-77],[-7,-17],[-40,-21],[-45,-61],[-76,-21],[-72,6],[-38,30],[-119,149],[-44,63],[-1,52],[-81,185],[0,40],[-28,57],[-32,9],[-7,54],[65,119],[41,100],[4,31],[-2,48],[-9,109],[8,46],[-40,-61],[-8,-46],[9,-46],[-6,-53],[-23,-74],[-40,-78],[-82,-37],[-145,21],[-17,10],[-10,27],[-7,114],[-9,-14],[-17,-61],[-13,-87],[-17,-20],[-30,-2],[-24,14],[-17,30],[-37,2],[-68,-23],[-32,9],[-38,1],[-77,26],[-92,6],[-25,20],[2,38],[17,19],[96,17],[94,41],[92,20],[-4,20],[-41,7],[-212,-45],[-67,5],[-12,8],[-1,45],[24,42],[41,40],[12,30],[-23,13],[-38,-8],[-19,24],[21,94],[-18,96],[-25,-91],[-37,-50],[-179,-22],[-29,-26],[-24,0],[-114,49],[-49,29],[-46,44],[-81,98],[-64,63],[-3,118],[14,74],[30,84],[116,180],[40,35],[36,11],[170,16],[126,23],[24,11],[-186,13],[-166,-10],[-56,-28],[-73,-116],[-18,-48],[-19,-36],[-13,-1],[-25,13],[-8,15],[-10,38],[-39,60],[-17,69],[-10,102],[3,46],[18,58],[56,114],[-72,-3],[9,96],[26,104],[66,62],[67,44],[61,59],[107,39],[36,-80],[91,-24],[26,-37],[33,-66],[40,-63],[48,-60],[13,-4],[-21,52],[-78,112],[-3,40],[-20,43],[-103,63],[-18,24],[-18,62],[-8,41],[12,39],[104,120],[24,63],[0,30],[-8,34],[-23,57],[-6,0],[7,-87],[-4,-35],[-12,-38],[-17,-29],[-22,-20],[-224,-273],[-22,-16],[-84,-25],[-42,-29],[-23,-38],[-37,-100],[-51,-200],[-58,-162],[-50,210],[-89,160],[174,159],[5,24],[-17,85],[4,26],[17,37],[44,55],[-2,5],[-51,-19],[-80,-123],[-31,-40],[-16,-5],[-2,72],[43,187],[35,183],[12,52],[31,52],[-29,-4],[-141,-80],[-47,50],[-40,264],[-69,103],[-118,84],[-115,39],[-27,75],[-23,91],[31,107],[51,51],[47,23],[44,-11],[2,-39],[-30,-105],[40,-10],[158,-129],[34,-10],[66,49],[36,-2],[86,-40],[29,-48],[83,-94],[-12,55],[-90,115],[-50,38],[-88,7],[-55,-20],[-23,6],[-48,29],[-41,49],[-41,106],[-9,49],[2,37],[10,33],[17,30],[33,21],[51,14],[15,13],[-64,60],[-30,0],[-101,-89],[-20,-6],[-9,17],[-8,1],[-30,-44],[-23,-19],[-82,-136],[-13,-65],[-5,-98],[-10,-60],[-14,-25],[-97,-47],[-55,-95],[-67,82],[-73,79],[-49,139],[-89,24],[-102,78],[-39,70],[55,143],[80,108],[12,130],[11,28],[138,32],[89,64],[-92,6],[-56,-10],[-101,-45],[-112,89],[-58,81],[-18,69],[20,59],[3,60],[10,80],[10,35],[24,45],[48,30],[42,92],[18,65],[87,193],[32,84],[59,115],[120,182],[-38,-10],[-20,-16],[-19,2],[-18,19],[-17,40],[-16,62],[-9,-27],[-2,-115],[-12,-98],[-21,-69],[-59,-138],[-34,-55],[-24,51],[18,87],[34,69],[6,94],[-32,113],[-20,91],[-7,69],[-1,63],[6,56],[13,62],[22,67],[-3,8],[-27,-52],[-20,-57],[-14,-60],[-5,-65],[4,-69],[10,-60],[34,-112],[11,-69],[1,-45],[-74,-161],[-25,-77],[2,-34],[6,-79],[-79,-132],[-101,-65],[-22,20],[-22,22],[-127,16],[-46,143],[-25,110],[-37,97],[1,21],[34,65],[127,53],[1,21],[-47,14],[-12,23],[-13,104],[4,91],[-3,61],[-22,124],[-31,76],[-82,147],[-7,38],[36,47],[22,43],[-138,-76],[-189,-79],[-81,-56],[-17,-22],[-6,-19],[15,-53],[-3,-16],[-16,-31],[-20,-88],[-41,-93],[-21,-19],[-73,35],[-20,30],[-37,121],[8,32],[27,27],[37,59],[46,90],[88,230],[58,1],[101,46],[-159,22],[-24,12],[-21,32],[-19,50],[-33,57],[-60,19],[-27,21],[-40,68],[-27,31],[-14,38],[-3,45],[-11,23],[-42,8],[-23,16],[-6,116],[-82,30],[-34,26],[-55,73],[-15,35],[-5,29],[13,81],[-5,15],[-48,-8],[-301,124],[16,164],[-56,217],[-59,88],[12,34],[12,18],[27,1],[115,-64],[109,-76],[14,11],[-174,161],[-43,48],[-11,57],[0,31],[13,16],[163,-14],[10,12],[-165,47],[-34,0],[-36,-68],[-17,-16],[-35,4],[-12,10],[-42,82],[-40,56],[-74,79],[-14,56],[-4,83],[10,78],[61,178],[24,31],[7,19],[-19,-3],[-18,-17],[-50,-82],[-52,-136],[-42,-46],[-26,11],[-41,55],[-83,68],[-97,17],[-61,69],[-91,192],[-12,96],[-12,23],[-49,31],[-30,46],[-46,233],[-62,163],[-15,85],[5,85],[-8,10],[-21,-67],[-6,-34],[-38,-10],[36,-67],[9,-34],[-18,3],[-37,-8],[63,-115],[28,-177],[42,-132],[27,-106],[12,-81],[18,-78],[49,-171],[7,-34],[-6,-28],[-16,-33],[-28,-12],[-88,22],[-33,43],[-47,77],[-66,36],[-164,-18],[-13,6],[0,64],[19,113],[-15,45],[-85,167],[2,32],[117,76],[-57,6],[-46,-30],[-18,19],[-28,107],[-18,39],[-9,9],[-5,-101],[20,-53],[3,-31],[-4,-43],[-12,-31],[-21,-18],[-22,-5],[-40,21],[-44,41],[-38,19],[-15,16],[-18,44],[-30,34],[-144,43],[-86,50],[-7,-13],[27,-54],[2,-32],[-21,-10],[-39,-51],[11,-7],[41,17],[46,-2],[76,-31],[68,-41],[25,-23],[10,-34],[9,-12],[67,-40],[3,-21],[-43,-62],[89,7],[52,-22],[66,-96],[23,-53],[3,-69],[-14,-19],[-26,-14],[-181,-23],[-66,-82],[-13,-1],[-50,22],[-90,66],[-113,62],[-257,186],[-7,9],[-4,36],[-17,18],[-35,16],[-49,47],[-62,78],[-38,61],[-14,44],[-36,51],[-117,105],[-60,40],[-54,23],[-47,5],[-12,14],[21,23],[3,14],[-103,21],[-98,49],[-248,139],[-128,87],[-75,41],[-32,23],[-14,19],[18,20],[50,21],[34,23],[54,88],[4,28],[-29,64],[-13,58],[1,32],[6,31],[8,21],[23,21],[16,10],[20,-7],[63,-80],[8,-29],[-2,-109],[18,-127],[5,9],[6,42],[4,81],[7,38],[13,38],[23,20],[71,-11],[33,6],[-138,58],[-87,108],[-16,11],[-48,5],[-50,-44],[-130,-142],[-36,-25],[-163,-79],[-111,-16],[-124,13],[-106,25],[-266,125],[-41,29],[62,77],[2,24],[-21,79],[-18,22],[-25,12],[-7,-8],[-1,-24],[7,-43],[-20,-23],[-45,-24],[-77,-25],[-235,63],[-242,53],[-216,11],[-305,-42],[-163,-42],[-94,-4],[-92,7],[-7,30],[42,17],[-3,21],[-52,66],[-80,41],[-108,14],[-61,19],[-16,25],[-38,23],[-60,22],[-27,40],[21,124],[21,74],[21,51],[52,85],[-18,-7],[-76,-61],[-66,-64],[-61,-83],[-36,-39],[-46,-35],[-73,9],[-99,52],[-84,27],[-70,1],[-28,9],[48,47],[27,37],[38,60],[9,29],[-258,9],[-10,32],[0,23],[-8,19],[-38,14],[-52,-10],[-85,-38],[-36,29],[13,16],[27,11],[57,54],[-76,28],[-39,32],[-19,26],[1,95],[20,61],[171,59],[-54,23],[-108,-8],[-72,-50],[-86,-72],[-58,-27],[-30,18],[-39,6],[-48,-5],[-33,-19],[-16,-31],[-20,-21],[-21,-11],[-16,3],[-24,32],[-49,20],[-23,25],[-14,-16],[-17,-46],[-18,-23],[-82,-24],[-46,4],[-54,58],[-8,20],[19,51],[119,199],[-12,-1],[-38,-31],[-77,-80],[-35,-24],[-59,-3],[-27,9],[-34,-7],[-39,-21],[-25,-23],[-12,-26],[8,-4],[59,29],[33,8],[10,-14],[-46,-90],[-28,-86],[-27,-22],[-42,4],[-46,-9],[-1,-24],[86,-69],[32,-9],[40,-25],[6,-24],[-15,-66],[-12,-26],[-18,-13],[-70,2],[-23,-7],[-47,-41],[-24,-35],[9,-3],[41,29],[59,15],[78,2],[58,15],[38,28],[38,-8],[36,-44],[11,-38],[-15,-33],[-30,-24],[-45,-14],[-29,-21],[-11,-28],[-7,-42],[-2,-55],[12,-100],[-9,-13],[-17,-8],[-25,-1],[-23,-23],[-53,-133],[-19,-14],[-22,14],[-20,-2],[-17,-17],[-38,-13],[-58,-9],[-50,3],[-88,29],[-36,20],[-28,33],[-79,-35],[-20,16],[-50,91],[-10,-5],[-10,-99],[-15,-35],[-48,-72],[-27,-123],[-8,-4],[-9,18],[-30,110],[-16,25],[-44,-64],[-5,-23],[12,-82],[-10,-13],[-89,45],[-22,2],[-6,-8],[30,-63],[-3,-23],[-126,-124],[-33,5],[-21,12],[-22,-2],[-80,-46],[-22,2],[-31,27],[-14,-1],[-7,-28],[-1,-55],[-30,-53],[-95,-85],[-25,-39],[-20,-53],[-14,-6],[-56,35],[-65,22],[-9,-11],[20,-33],[-4,-20],[-28,-7],[-36,3],[-42,13],[-61,-15],[-77,-43],[-65,1],[-90,71],[-25,6],[-7,20],[17,57],[26,44],[19,20],[85,55],[98,21],[61,33],[75,69],[40,52],[78,134],[-6,11],[-18,7],[-171,-127],[-25,-12],[-34,1],[-137,49],[-28,20],[-20,61],[38,139],[26,67],[67,104],[87,110],[30,72],[46,191],[-3,87],[-20,106],[-1,63],[19,20],[200,98],[95,74],[184,108],[50,-1],[37,-37],[42,-30],[49,-22],[63,2],[77,27],[121,-10],[250,-71],[54,-4],[2,9],[-39,50],[-172,29],[-73,29],[-204,127],[-46,50],[19,23],[49,19],[18,18],[7,32],[28,44],[51,55],[76,54],[145,80],[-56,4],[-105,-15],[-38,-15],[-70,-58],[-27,-40],[-39,-78],[-16,-15],[-73,-12],[-197,-8],[-33,41],[-19,6],[-24,-6],[-182,-102],[-65,-53],[-46,-59],[-72,-44],[-96,-28],[-73,-34],[-76,-69],[-26,-53],[-2,-25],[19,-78],[-19,-14],[-44,-6],[-71,-52],[-149,-154],[-20,-56],[1,-19],[24,-43],[-17,-29],[-42,-44],[-93,-71],[-62,-27],[-40,-1],[-38,10],[-68,45],[-55,3],[-4,-6],[75,-49],[77,-63],[47,-52],[19,-41],[1,-43],[-17,-44],[-54,-76],[-53,-23],[-136,-23],[-43,-18],[-14,-14],[93,-32],[9,-17],[-13,-63],[-25,-21],[-77,-38],[-70,-11],[-10,7],[13,50],[-4,12],[-26,11],[-37,-20],[-93,-73],[-10,-12],[34,-20],[-7,-17],[-50,-53],[-21,-35],[-34,-36],[-149,-110],[12,-27],[-39,-96],[-22,-85],[27,-35],[125,-42],[61,-10],[71,-29],[130,-79],[43,-51],[6,-24],[-4,-27],[-15,-35],[-41,-68],[-98,-100],[-44,-28],[-67,-22],[-22,-16],[-85,-95],[-24,-51],[4,-45],[-16,-31],[-111,-61],[4,-11],[40,-5],[-15,-54],[-6,-75],[-19,-12],[-69,0],[-88,-29],[-6,-8],[-2,-54],[-229,-40],[-50,-102],[-27,-32],[-90,-74],[-55,-30],[-62,-19],[-33,-25],[-3,-32],[-18,-28],[-55,-47],[-27,-58],[-19,-9],[-101,-14],[-21,-18],[-9,-79],[-19,-3],[-36,19],[-47,-15],[-105,-89],[-23,-32],[2,-17],[17,-17],[24,-53],[-1,-35],[-41,-100],[-14,-15],[-49,-25],[-20,-55],[-46,6],[-36,-10],[-24,-37],[-26,-21],[-28,-6],[-36,-29],[-42,-52],[-40,-34],[-36,-14],[-35,-4],[-35,7],[-30,-6],[-28,-19],[-26,-31],[-22,-86],[-27,-39],[-17,-6],[-35,5],[-52,19],[-54,-7],[-86,-53],[-28,-40],[55,-9],[27,-11],[-1,-11],[-28,-11],[-48,1],[-30,-10],[-35,-23],[-89,-24],[-33,-18],[-67,-100],[-8,-23],[8,-5],[38,11],[44,-17],[23,-21],[15,-26],[14,-51],[8,-7],[-85,-85],[-24,-35],[-15,-14],[-11,11],[-10,94],[-6,16],[-20,1],[-20,-29],[-42,-112],[-46,-56],[-348,-144],[-51,-32],[-10,-62],[-14,-53],[-23,-42],[-27,-27],[-6,19],[2,150],[-7,29],[-35,19],[-15,-2],[-21,-9],[-35,-32],[-22,-8],[-26,3],[-45,-32],[-108,-103],[-70,-25],[-19,-21],[-30,-56],[-20,-21],[-30,-1],[-39,17],[-31,-12],[-24,-41],[-24,-16],[-68,30],[-30,-21],[-40,-52],[-40,-35],[-43,-17],[-111,-17],[-45,11],[-9,15],[2,67],[18,48],[17,23],[22,20],[32,3],[61,-15],[-7,16],[-22,19],[-56,33],[-55,17],[-31,-11],[-45,-25],[-30,-30],[-16,-33],[-20,-109],[-12,-29],[-129,-193],[-51,-59],[-51,5],[-24,-23],[-35,-48],[-31,-23],[-29,3],[-23,9],[-16,16],[3,15],[21,14],[-7,38],[-38,63],[-25,34],[-48,4],[-7,-28],[16,-146],[-3,-33],[-30,-42],[-79,-47],[-25,5],[-71,92],[-67,18],[-5,-30],[16,-61],[-17,-57],[-49,-53],[-37,-26],[-25,2],[-2,37],[23,73],[6,60],[-11,49],[1,37],[13,26],[90,72],[37,11],[21,-18],[25,-3],[30,11],[19,23],[8,35],[38,44],[69,53],[80,99],[89,145],[104,124],[120,105],[130,83],[262,114],[21,-7],[-24,-37],[16,-23],[26,-3],[96,18],[38,24],[11,-23],[-13,-29],[-58,-30],[2,-24],[83,-115],[27,-18],[22,2],[9,15],[-7,83],[28,16],[58,4],[38,-11],[18,-26],[33,-21],[49,-16],[30,5],[12,27],[-21,32],[-94,71],[-25,29],[-7,41],[14,53],[29,78],[45,104],[41,73],[83,82],[56,40],[142,125],[273,126],[68,82],[92,89],[39,23],[0,-35],[12,-31],[62,-21],[40,-7],[18,6],[5,33],[-8,60],[-2,56],[4,54],[9,41],[41,75],[60,85],[84,98],[52,45],[49,24],[48,43],[83,101],[26,17],[59,20],[22,-9],[12,-25],[16,-16],[60,-14],[40,22],[-7,12],[-32,8],[-21,15],[-20,60],[-39,37],[-9,41],[7,65],[34,151],[6,155],[30,89],[61,32],[136,22],[-80,40],[-29,0],[-52,19],[-19,97],[0,71],[34,81],[126,138],[139,95],[-19,8],[-17,28],[64,191],[62,170],[-84,-145],[-97,-111],[-285,-129],[-194,-108],[-92,-26],[-60,28],[-48,103],[-27,37],[-35,68],[15,88],[28,60],[60,10],[68,-29],[59,-3],[-76,60],[-110,53],[-50,-17],[-38,-85],[-51,-58],[-45,20],[-26,24],[18,-71],[-33,-109],[-13,-75],[48,-198],[-9,-79],[-88,-36],[-72,65],[-150,251],[-52,71],[-117,118],[-39,-16],[-49,-59],[-48,-16],[-127,86],[-58,66],[-56,79],[-85,-44],[-75,-52],[-87,-83],[-58,1],[-159,-72],[-17,-1],[-22,-39],[-22,-17],[-18,-74],[-214,-57],[-212,32],[74,41],[83,32],[72,77],[-31,103],[-5,52],[1,66],[78,93],[-81,0],[-53,-33],[-49,70],[-23,137],[56,82],[26,62],[22,86],[2,74],[-43,126],[-125,265],[-57,198],[-97,105],[72,173],[81,157],[105,70],[-9,11],[-57,-1],[-38,-9],[-34,-51],[-35,-39],[-111,-200],[-71,-98],[-47,-28],[75,-37],[11,-32],[15,-73],[-19,-88],[-20,-48],[-88,4],[-79,-71],[-185,-77],[-251,-44],[-123,5],[-128,90],[0,52],[6,45],[-185,155],[-104,154],[-75,4],[-65,41],[-77,64],[7,51],[12,37],[-47,25],[-61,-3],[-70,18],[184,198],[63,133],[51,19],[67,-20],[92,-53],[78,-23],[28,-24],[29,-47],[-31,-78],[-27,-54],[34,14],[97,85],[71,74],[35,-7],[22,-13],[40,-77],[50,-78],[109,74],[59,93],[-49,40],[-61,24],[-154,32],[38,27],[99,-3],[37,25],[-39,35],[-49,32],[-134,-105],[-243,5],[-170,61],[-169,-10],[-27,12],[-33,33],[95,77],[68,43],[4,25],[-40,4],[-74,-21],[-32,36],[5,62],[-12,-6],[-29,-34],[-42,17],[-35,28],[18,30],[37,41],[-16,6],[-33,-9],[-32,-53],[7,-44],[-1,-62],[-54,-11],[-46,7],[-34,63],[-35,134],[-93,36],[-23,68],[59,87],[-26,45],[-63,15],[-73,-44],[-31,39],[-6,43],[-3,61],[20,7],[17,-12],[145,34],[13,17],[-114,52],[-32,54],[47,31],[86,3],[120,32],[-50,58],[-11,32],[-10,53],[20,88],[141,203],[138,169],[43,39],[63,21],[58,-16],[61,-36],[12,16],[-21,14],[-26,70],[85,27],[50,78],[4,23],[-54,-33],[-57,-53],[-14,53],[-15,124],[25,117],[20,52],[47,50],[135,20],[24,-10],[5,24],[-81,73],[33,58],[30,29],[164,47],[89,-15],[113,-54],[65,-67],[-10,-35],[-16,-20],[-34,-23],[-12,-17],[6,-14],[48,40],[79,49],[44,-21],[35,-39],[38,1],[123,33],[62,35],[77,92],[101,59],[142,186],[42,77],[49,12],[44,-7],[30,-63],[45,-18],[255,15],[130,29],[90,61],[94,102],[55,69],[26,89],[-34,116],[-34,96],[-46,219],[-126,145],[-90,44],[-57,-6],[41,92],[121,-10],[78,18],[64,45],[20,33],[32,69],[-11,74],[-17,40],[-44,43],[-52,65],[-36,21],[-31,-1],[-152,-129],[-91,-2],[-68,23],[-60,-73],[-165,-65],[-88,-65],[-164,-161],[-41,-74],[-52,-3],[-38,142],[-178,135],[-54,-46],[30,-42],[40,-30],[67,-14],[-29,-41],[-21,-54],[-67,51],[-119,74],[-124,39],[-321,-5],[-211,-76],[-19,16],[-21,6],[-35,-18],[-15,-31],[-23,-20],[-43,-7],[-87,12],[-167,47],[-379,70],[-99,43],[-85,102],[2,70],[37,29],[-3,99],[-74,27],[-150,143],[-55,60],[12,7],[27,-16],[51,-13],[126,20],[43,92],[93,27],[87,-13],[-20,25],[-22,20],[-224,47],[-30,-15],[-402,84],[-317,145],[-26,28],[-29,62],[43,61],[43,29],[2,-34],[7,-33],[181,77],[95,101],[180,18],[42,28],[56,54],[80,92],[113,49],[77,44],[100,25],[85,-42],[27,-6],[155,-9],[51,18],[22,14],[16,22],[-153,78],[16,43],[19,31],[178,91],[137,30],[73,-3],[212,117],[116,34],[219,22],[179,6],[49,-42],[-97,9],[-42,-8],[30,-14],[34,-30],[-10,-39],[-59,-114],[5,-91],[-39,-30],[-37,-41],[184,-132],[285,-8],[155,24],[89,-40],[74,-9],[202,20],[153,-28],[64,11],[141,197],[55,30],[60,-34],[78,-28],[50,21],[41,-51],[-19,106],[-28,39],[-230,73],[-155,-36],[-48,41],[16,81],[-165,199],[-69,41],[-81,2],[-42,69],[-34,89],[70,36],[63,17],[59,-29],[66,-117],[62,-17],[-18,-117],[77,-107],[173,-100],[139,37],[98,-1],[58,-21],[144,-90],[73,-11],[227,47],[3,88],[-19,64],[-54,40],[-154,-8],[-119,66],[-102,-18],[-189,-101],[-94,40],[-60,54],[-95,54],[-12,104],[80,118],[59,57],[-53,41],[-133,29],[-232,-30],[-11,41],[1,43],[-94,-85],[-97,18],[-131,-9],[-288,75],[-103,93],[-43,75],[-78,206],[-99,129],[-686,438],[-311,110],[-151,122],[-94,30],[-90,12],[-115,39],[77,49],[54,16],[-56,-51],[42,-12],[67,29],[37,35],[53,147],[55,224],[-15,88],[380,-18],[253,15],[84,20],[320,34],[83,25],[153,75],[181,133],[155,175],[24,47],[10,-12],[14,7],[17,67],[20,156],[77,147],[327,335],[152,133],[51,60],[53,44],[37,-41],[18,-13],[10,-20],[-31,-9],[-51,-43],[-71,-28],[-17,-15],[41,3],[125,31],[70,38],[350,70],[189,116],[8,26],[281,144],[39,-5],[45,-18],[-78,-95],[54,-25],[-48,-114],[102,-2],[23,-52],[5,45],[-1,65],[8,63],[15,44],[72,-20],[161,48],[-196,6],[-117,103],[-65,1],[218,151],[199,92],[45,-2],[22,-17],[5,-28],[-43,-18],[-43,-32],[20,-29],[29,-4],[95,24],[43,29],[204,-2],[60,21],[15,21],[264,4],[48,15],[166,81],[152,98],[71,53],[120,137],[104,88],[170,89],[41,-11],[-55,-17],[-39,-38],[53,-50],[358,-103],[90,-5],[36,-62],[-30,-59],[-92,-66],[-186,-68],[57,-25],[37,-60],[55,-8],[89,23],[70,37],[145,120],[46,67],[34,17],[121,-16],[69,-34],[78,-62],[-29,-59],[-32,-33],[102,-46],[112,-10],[107,-37],[151,76],[118,16],[111,-3],[144,42],[243,-57],[62,15],[98,-10],[104,-34],[36,-36],[-111,-77],[-18,-79],[39,-34],[70,-5],[9,-47],[44,-11],[220,3],[-17,-22],[-11,-26],[-68,-60],[392,-33],[52,33],[81,13],[172,45],[65,-20],[76,-46],[71,-9],[66,9],[153,66],[178,3],[72,-21],[77,9],[231,-75],[85,-9],[114,-98],[58,-3],[67,41],[58,-1],[56,-40],[92,-12],[43,-63],[46,-23],[350,-47],[173,22],[252,-6],[122,-30],[127,4],[208,-109],[111,-17],[21,-25],[315,-27],[110,57],[192,15],[172,48],[98,0],[114,-12],[44,5],[32,21],[277,-82],[156,-94],[68,-70],[325,-99],[94,-55],[64,-62],[38,-6],[26,18],[114,-6],[43,-8],[77,-16],[247,-32],[233,19],[433,-106],[267,-199],[217,-97],[89,-67],[140,-59],[332,-128],[103,-14],[192,-61],[119,8],[204,-15],[140,-50],[274,-136],[56,-12],[15,10],[-95,135],[-15,13],[-111,50],[-130,25],[-10,9],[-24,48],[8,18],[28,8],[97,0],[57,8],[8,18],[-41,5],[-50,18],[-58,32],[-33,30],[119,199],[41,-20],[63,46],[113,-29],[20,16],[14,101],[17,24],[31,18],[157,18],[195,-18],[20,10],[-18,68],[-3,26],[12,61],[12,32],[23,17],[91,-13],[28,-30],[31,-52],[31,-29],[96,-30],[11,-20],[-37,-78],[-38,-41],[-80,-108],[-5,-27],[123,48],[138,67],[119,37],[99,7],[71,21],[43,37],[30,38],[62,121],[40,21],[171,-8],[40,4],[27,12],[-5,15],[-36,19],[-49,5],[-1,9],[18,19],[27,12],[84,14],[55,-44],[37,-3],[125,48],[192,129],[76,35],[67,7],[56,-24],[43,6],[57,73],[22,38],[35,34],[142,75],[91,16],[55,-14],[66,-31],[55,-12],[71,9],[54,-4],[25,15],[92,86],[29,1],[29,-25],[46,-63],[0,-31],[-60,-76],[-441,-217],[-135,-94],[-68,-35],[-69,-19],[-135,-16],[-54,-19],[-90,-17],[-212,-30],[-41,-15],[-28,-17],[-76,-115],[-37,-38],[-73,-56],[-81,-35],[-112,-13],[-71,-54],[-82,-104],[-66,-73],[-76,-61],[-82,-81],[-21,-42],[24,-56],[14,-18],[82,-30],[32,6],[-29,30],[-69,43],[-10,16],[18,12],[325,-32],[70,32],[25,28],[-6,15],[-88,5],[-19,27],[-14,50],[-3,39],[8,30],[20,37],[95,62],[102,26],[78,35],[43,32],[117,56],[47,46],[25,35],[3,17],[-21,13],[17,30],[86,27],[37,3],[120,-26],[21,-20],[-11,-54],[16,2],[46,69],[26,22],[27,5],[27,-9],[27,-22],[15,-66],[2,-109],[6,-44],[31,76],[21,35],[117,161],[78,88],[89,86],[128,64],[297,106],[167,29],[84,26],[42,23],[26,29],[47,33],[8,-3],[-18,-69],[-12,-19],[-109,-43],[-10,-32],[12,-50],[18,-33],[26,-16],[45,11],[65,39],[80,59],[173,150],[15,27],[44,120],[99,53],[180,61],[44,38],[-157,33],[-33,22],[-5,13],[30,35],[-73,34],[-26,21],[1,61],[22,44],[48,43],[26,7],[71,-25],[59,-32],[204,-148],[82,-72],[48,-57],[115,-177],[51,-103],[41,-105],[40,-76],[39,-47],[197,-185],[101,-78],[86,-48],[97,-38],[111,-29],[75,-2],[117,78],[2,52],[-51,86],[-52,60],[6,36],[69,70],[-5,25],[15,71],[47,-13],[19,2],[26,27],[34,50],[43,41],[52,33],[14,21],[-51,16],[-32,0],[-23,7],[-15,14],[21,15],[113,38],[21,37],[36,24],[46,10],[28,-10],[32,-30],[2,-49],[-14,-79],[-3,-64],[35,-151],[32,-33],[122,-44],[-8,-37],[-141,-159],[-30,-39],[-15,-29],[5,-25],[25,-21],[48,-15],[123,-6],[34,14],[239,5],[44,12],[37,30],[54,77],[61,23],[19,22],[38,91],[19,105],[18,44],[28,28],[37,8],[93,-10],[44,9],[173,-9],[172,8],[179,-19],[114,-21],[106,-35],[204,-81],[80,-43],[284,-196],[83,-40],[156,-38],[535,-85],[68,-23],[141,-89],[97,-52],[115,-50],[144,-43],[282,-65],[46,-21],[52,-6],[58,8],[258,-37],[68,2],[50,-7],[60,-27],[89,-9],[-3,18],[-102,102],[5,16],[42,2],[125,-17],[29,29],[42,-1],[95,-14],[103,-32],[110,-49],[134,-41],[203,-104],[112,-86],[106,-108],[59,-74],[10,-42],[22,-21],[34,1],[13,-16],[-30,-93],[-18,-23],[-23,-16],[-97,-19],[-267,22],[-47,-75],[-150,-63],[-27,-27],[-5,-21],[10,-65],[-19,-20],[-122,-75],[-4,-21],[79,-31],[85,-52],[66,-13],[84,7],[105,-18],[127,-44],[89,-20],[49,4],[68,-8],[86,-20],[115,-7],[254,2],[76,-15],[106,-7],[205,2],[37,3],[66,35],[42,11],[73,1],[213,25],[73,0],[68,19],[87,39],[54,7],[20,-24],[37,-11],[53,4],[100,41],[236,122],[85,0],[62,38],[15,0],[17,-15],[58,-89],[17,-15],[40,-7],[39,-46],[40,-68],[30,-19],[221,-3],[78,-19],[23,-20],[24,-54],[15,-104],[9,-39],[32,-55],[22,-16],[20,15],[54,145],[19,23],[36,-8],[11,-7],[56,-107],[78,-80],[195,-147],[32,-53],[11,-40],[-11,-36],[-33,-32],[-53,-25],[-72,-19],[-67,7],[-63,33],[-20,2],[22,-30],[129,-121],[33,-49],[31,-32],[28,-16],[26,-26],[24,-37],[107,-97],[30,-47],[122,-145],[58,-57],[45,-32],[18,-4],[-11,27],[-155,193],[-80,121],[-11,29],[-5,45],[-3,147],[11,23],[54,19],[69,-67],[26,-10],[18,4],[9,18],[40,-16],[70,-50],[24,-1],[-53,95],[-38,46],[-14,32],[36,49],[-20,24],[-89,71],[-46,74],[-42,111],[-3,44],[6,46],[-6,37],[-56,75],[-61,52],[-48,61],[-10,32],[7,86],[37,38],[70,50],[18,52],[-32,54],[-5,24],[20,-4],[137,26],[34,-9],[52,11],[69,33],[54,15],[71,-3],[39,8],[47,15],[24,16],[44,62],[23,9],[73,-7],[41,-13],[19,5],[-3,87],[14,31],[72,64],[76,6],[50,18],[58,36],[41,32],[41,51],[17,65],[-13,19],[-86,26],[-51,-14],[-115,-46],[-120,-60],[-46,-56],[-13,-71],[-22,-32],[-94,30],[-40,-1],[-50,-12],[-53,-27],[-56,-43],[-82,-8],[-110,27],[-65,8],[-67,-45],[4,-34],[31,-49],[-31,-29],[-159,-11],[-42,6],[-85,-19],[-34,4],[-24,24],[-174,99],[-17,20],[42,81],[161,219],[17,13],[298,38],[180,40],[329,120],[63,10],[212,80],[87,20],[81,-13],[118,-42],[61,-36],[44,-46],[36,-63],[45,-141],[15,-119],[28,-45],[99,-83],[51,-32],[32,-10],[27,19],[18,3],[13,-7],[13,-53],[18,-5],[60,7],[63,-22],[9,-16],[-13,-65],[17,-27],[77,-59],[74,-21],[86,-10],[159,9],[132,28],[100,46],[82,-51],[164,-123],[98,-88],[81,-41],[165,-49],[38,-26],[60,-3],[83,20],[94,-8],[106,-36],[73,-16],[249,69],[38,4],[93,33],[60,10],[70,-1],[53,10],[34,21],[133,-1],[239,-23],[163,-31],[98,-39],[79,-21],[63,-5],[60,7],[60,18],[62,36],[133,18],[22,9],[-3,19],[-27,30],[-76,55],[-53,54],[-10,35],[1,41],[17,24],[32,7],[49,-26],[68,-58],[192,-218],[46,-31],[26,-28],[175,-80],[83,-14],[99,49],[43,31],[21,29],[-1,29],[10,41],[-10,26],[-27,33],[-70,46],[-114,59],[-105,18],[-95,-25],[-107,-50],[-45,20],[-133,141],[-34,53],[0,14],[62,-16],[3,17],[-37,68],[-23,23],[-77,107],[-10,32],[48,7],[22,13],[29,-1],[136,-65],[70,29],[161,41],[-64,62],[-15,61],[8,13],[52,9],[103,-51],[50,-7],[36,21],[39,1],[40,-17],[38,-27],[71,-74],[35,-45],[39,-68],[12,-10],[190,-5],[107,60],[-2,-20],[-25,-47],[-133,-181],[2,-23],[71,10],[33,14],[20,22],[18,51],[12,15],[198,85],[57,13],[-36,-91],[-74,-326],[-15,-113],[-16,-39],[-77,-125],[1,-44],[85,-105],[15,-29],[9,-86],[15,-17],[70,-1],[72,27],[87,19],[13,-18],[-47,-105],[3,-9],[82,27],[38,3],[15,-6],[61,-53],[6,-40],[-1,-60],[-6,-42],[-21,-24],[-25,-10],[-31,-7],[-28,3],[-86,-9],[-50,12],[-50,32],[-36,8],[-41,-25],[-66,4],[-73,72],[-29,-6],[-10,-12],[1,-15],[33,-47],[258,-248],[39,-51],[9,-73],[5,0],[24,74],[-15,35],[-133,144],[-16,53],[5,14],[35,15],[189,-36],[73,8],[49,25],[25,31],[18,172],[34,109],[-20,99],[-51,156],[-41,92],[-92,94],[-9,33],[103,284],[18,24],[24,10],[81,4],[59,24],[93,-33],[51,-9],[63,29],[141,119],[56,37],[70,70],[85,102],[93,74],[150,69],[91,56],[19,20],[-85,5],[-21,9],[-18,53],[9,98],[-1,54],[-11,49],[-18,44],[-27,38],[-25,23],[-23,8],[-15,-4],[-8,-14],[-23,-93],[-29,-69],[-40,-35],[-83,-25],[-142,-17],[-59,32],[-7,28],[20,108],[48,47],[130,91],[83,73],[1,13],[-76,0],[-19,15],[-16,90],[5,34],[12,38],[53,28],[164,35],[128,40],[4,-14],[-103,-121],[-10,-29],[39,-26],[98,71],[64,58],[11,20],[-58,7],[-3,24],[11,44],[-5,30],[-64,38],[-79,-21],[-66,-40],[-54,-11],[-80,-1],[-59,10],[-37,19],[-44,42],[-51,67],[-65,66],[-24,7],[-19,-8],[-42,-63],[-18,-8],[-255,89],[-109,50],[-52,38],[-65,23],[-78,7],[-62,17],[-47,29],[-37,40],[-29,53],[-53,67],[-121,133],[-32,85],[-5,33],[8,84],[114,142],[21,42],[38,30],[57,19],[40,6],[93,-18],[-56,44],[-4,24],[55,77],[-10,4],[-153,-60],[-39,4],[-55,37],[-103,129],[-1,81],[32,113],[12,67],[-30,57],[11,16],[32,16],[14,17],[-14,64],[22,33],[73,66],[72,57],[42,17],[37,-3],[38,-18],[40,-31],[68,-33],[50,-10],[38,18],[63,123],[23,32],[-22,14],[-124,-2],[-54,10],[-31,12],[-23,47],[19,25],[120,87],[57,91],[169,127],[171,60],[83,19],[67,4],[29,-8],[36,-63],[7,-67],[93,-83]],[[53002,44943],[-39,-13],[-113,27],[-47,19],[-139,72],[-264,105],[-91,53],[-18,39],[48,82],[27,34],[29,20],[278,38],[109,-18],[128,-169],[73,-113],[30,-92],[1,-44],[-12,-40]],[[56868,59949],[-265,-2],[-144,7],[-58,16],[-50,23],[-57,77],[-37,79],[-2,34],[10,30],[13,18],[178,27],[146,-33],[125,-39],[162,-6],[49,-15],[19,-11],[13,-20],[18,-91],[2,-48],[-5,-37],[-117,-9]],[[54554,55644],[46,-42],[48,3],[27,-35],[63,-109],[6,-18],[-3,-35],[-33,-45],[-35,-24],[-72,-35],[-82,-15],[-45,24],[-126,101],[-116,118],[-37,66],[16,27],[49,19],[100,19],[163,-10],[31,-9]],[[52312,55129],[21,-28],[7,-42],[-9,-56],[-14,-51],[-19,-47],[-48,-77],[-149,-133],[-56,-73],[-42,-43],[-243,-193],[-31,-9],[-31,3],[-67,29],[-68,5],[-175,-77],[-8,13],[-8,83],[-18,47],[-77,98],[-5,23],[1,29],[6,22],[87,92],[198,336],[47,16],[97,-37],[46,-12],[33,2],[141,70],[134,-8],[122,41],[58,0],[44,-7],[26,-16]],[[53647,54588],[43,-27],[71,-91],[27,-46],[8,-62],[-19,-82],[-10,-77],[-26,-58],[-49,-74],[-43,-86],[-38,-98],[-31,-65],[-25,-32],[-27,-16],[-28,-2],[-44,36],[-59,72],[-47,44],[-62,31],[-33,40],[-6,41],[-2,130],[3,65],[8,55],[14,43],[31,60],[85,130],[51,50],[32,11],[85,-9],[34,3],[28,15],[29,-1]],[[39647,65074],[-27,-5],[-52,23],[-75,51],[-124,104],[-147,103],[-23,62],[-36,46],[-188,109],[-123,44],[-93,22],[-15,30],[65,89],[74,71],[44,25],[137,23],[462,47],[105,2],[111,-23],[152,-97],[63,-11],[38,-22],[34,-36],[17,-37],[2,-75],[-17,-112],[-21,-42],[-92,-144],[-97,-78],[-18,-51],[-39,-41],[-69,-50],[-48,-27]],[[42522,66098],[459,-214],[65,18],[138,11],[145,32],[199,26],[122,47],[52,14],[88,8],[48,0],[139,-26],[56,-18],[29,-19],[32,-34],[55,-84],[8,-31],[-3,-9],[-51,-51],[-33,-23],[-70,-21],[-60,-7],[-52,-37],[-58,10],[-16,-35],[7,-21],[15,-10],[30,3],[33,14],[65,-7],[35,-22],[31,-38],[-23,-34],[-115,-47],[-170,-54],[-207,-169],[-108,-71],[-23,-22],[-11,-23],[4,-42],[5,-17],[30,-7],[101,59],[65,28],[66,15],[117,1],[48,-9],[87,-36],[80,-54],[18,-20],[-8,-20],[-33,-20],[-4,-13],[77,-27],[84,-75],[5,-45],[-37,-44],[-10,-30],[17,-16],[42,10],[99,50],[108,26],[43,-3],[27,-12],[29,-68],[23,-78],[3,-65],[-18,-52],[-25,-42],[-64,-51],[-59,-19],[-30,0],[3,-9],[69,-42],[29,-34],[12,-33],[-4,-31],[-10,-27],[-80,-100],[4,-14],[23,-7],[50,-57],[6,-134],[-181,-42],[-43,-31],[-50,-49],[-57,-38],[-130,-39],[-66,-4],[-326,32],[-33,19],[-22,34],[-13,49],[-3,39],[5,28],[-1,16],[-10,5],[-36,-28],[-38,-51],[21,-57],[103,-157],[20,-68],[3,-28],[-5,-22],[-116,-92],[-67,-31],[-70,-15],[-65,11],[-63,39],[-47,17],[-99,-2],[-30,20],[-29,37],[-69,132],[-97,94],[-83,106],[-212,154],[-110,92],[-146,148],[-61,34],[-51,13],[-100,11],[-23,18],[-36,51],[-63,37],[-23,5],[-37,-9],[-97,-35],[-123,37],[-28,26],[-16,44],[-15,25],[-42,21],[-35,55],[-233,109],[-139,124],[-28,45],[-2,18],[14,60],[32,68],[42,67],[26,28],[91,58],[72,13],[100,-7],[53,-13],[47,-39],[21,-45],[23,-30],[75,-33],[40,-27],[59,-62],[46,-74],[42,-25],[101,-8],[104,13],[224,46],[9,6],[14,28],[25,170],[15,1],[76,-79],[22,-8],[33,16],[19,38],[-2,18],[-48,92],[-28,40],[-25,25],[-28,8],[-61,-8],[-47,18],[-10,25],[8,35],[25,36],[28,20],[53,11],[63,-11],[86,-43],[55,-10],[76,13],[-98,25],[-137,101],[-59,18],[-71,-43],[-49,-16],[-91,-21],[-73,-3],[-306,155],[-17,15],[-22,40],[3,19],[30,29],[76,38],[113,24],[76,4],[66,-29],[96,-77],[85,-45],[7,17],[-15,44],[-39,65],[-27,18],[-68,20],[-63,45],[-30,33],[-15,36],[-2,38],[12,26],[26,14],[234,38],[164,-42],[104,-6],[43,51],[-15,10],[-54,-13],[-62,0],[-39,29],[-1,16],[48,40],[74,17]],[[50724,57431],[19,-12],[23,7],[39,49],[86,142],[25,13],[36,-2],[128,-91],[47,-51],[25,-70],[26,-28],[101,-37],[96,-12],[126,-37],[46,-28],[100,-135],[12,-9],[114,-55],[176,-124],[44,-18],[171,-41],[62,-30],[59,-46],[66,-85],[77,-133],[60,-213],[5,-42],[-7,-25],[-22,-27],[-99,-87],[8,-15],[93,6],[207,53],[126,-38],[44,-6],[10,2],[46,67],[50,-11],[73,-65],[47,-52],[21,-40],[-5,-23],[-49,-8],[119,-38],[103,-61],[-23,-40],[-109,-88],[-113,-77],[-132,-116],[-33,-18],[-17,-1],[-73,22],[-104,55],[-322,126],[-99,26],[-127,16],[-18,30],[-30,190],[-57,33],[-195,40],[-56,23],[-3,38],[12,65],[-26,32],[-66,-1],[-64,-14],[-104,-44],[-48,-41],[-18,-45],[-12,-95],[-13,-45],[-36,-60],[-160,-153],[-65,-46],[-64,-13],[-26,-15],[-43,-56],[-65,-139],[-26,-40],[-43,-41],[-88,-62],[-90,-49],[-151,-58],[-84,-20],[-56,20],[-38,131],[-81,386],[-13,26],[-17,17],[-19,6],[-270,-50],[-149,6],[-148,-87],[-37,-6],[-75,2],[-54,14],[-13,10],[-9,37],[2,40],[19,42],[67,115],[54,72],[25,22],[252,127],[62,42],[31,43],[0,45],[-12,56],[-44,138],[-11,127],[0,62],[17,97],[63,233],[21,114],[41,406],[21,116],[31,107],[31,63],[80,128],[62,51],[79,35],[17,-4],[15,-15],[29,-54],[109,-50],[36,-48],[25,-54],[13,-70],[-13,-31],[-52,-46],[-9,-19],[1,-16],[100,-73],[74,-176]],[[55751,60486],[284,-88],[28,-32],[13,-29],[8,-32],[1,-65],[-6,-27],[-28,-63],[-2,-20],[23,-214],[-3,-116],[-19,-97],[-41,-78],[-61,-59],[-47,-34],[-206,-83],[-147,-21],[-154,-4],[-196,-22],[-90,4],[-47,10],[-33,17],[-38,52],[-44,89],[-37,100],[-43,176],[-1,21],[42,144],[56,97],[98,143],[111,140],[29,23],[50,25],[126,39],[104,-9],[46,5],[57,18],[64,5],[103,-15]],[[53651,65807],[93,-13],[587,25],[122,-19],[370,-110],[96,-37],[47,-49],[43,-78],[19,-17],[134,-46],[55,-55],[20,-29],[28,-66],[61,-37],[69,-22],[22,-19],[-10,-83],[29,-39],[65,-46],[25,-31],[-51,-39],[-118,-22],[-333,23],[-447,53],[-261,-15],[-130,-24],[-316,-86],[-100,-14],[-99,-1],[-174,70],[-63,37],[-21,29],[-41,84],[-35,101],[-17,83],[-20,63],[-60,21],[-176,25],[-60,35],[-27,29],[-26,46],[1,48],[14,43],[11,10],[22,1],[-49,52],[-17,56],[-1,79],[6,51],[14,22],[33,14],[77,10],[113,-2],[159,-58],[126,-6],[191,-47]],[[46237,66313],[214,-48],[104,-30],[51,-21],[97,-70],[50,-20],[189,41],[133,15],[295,-19],[250,-57],[92,-43],[57,-40],[-15,-44],[-46,-71],[-54,-67],[-108,-106],[-92,-52],[-23,-26],[-14,-35],[-36,-51],[-99,-112],[-26,-19],[-139,-48],[47,-22],[22,-19],[-18,-50],[-88,-119],[-91,-110],[-64,-66],[-115,-96],[-64,-27],[-85,-8],[-516,83],[-130,-1],[-344,-43],[33,-22],[126,-32],[81,-35],[108,-109],[14,-28],[7,-31],[-4,-63],[-8,-16],[-171,-168],[-56,-122],[-35,-101],[-58,-28],[-192,43],[-62,-2],[-216,-29],[-101,15],[15,152],[-14,164],[-32,157],[-161,279],[-18,50],[-12,53],[-6,57],[1,57],[11,116],[1,59],[-8,153],[-22,229],[-2,81],[2,32],[6,25],[35,33],[66,23],[34,4],[212,-70],[97,-23],[65,1],[3,8],[-57,15],[-55,30],[-88,88],[-41,77],[-8,25],[-2,26],[5,27],[13,25],[44,39],[36,17],[133,45],[134,27],[296,18],[83,-13],[128,49],[76,11],[130,-17]],[[43916,61804],[17,-12],[31,10],[25,32],[14,7],[23,-7],[76,-58],[59,-61],[62,-44],[98,-39],[214,-126],[64,-86],[66,-132],[60,-102],[53,-70],[56,-56],[90,-62],[80,45],[35,12],[30,-19],[28,-47],[-15,-21],[-35,-27],[-58,-31],[-81,-2],[-38,-9],[-64,-50],[-50,-59],[-70,-20],[-133,-99],[-73,-37],[-107,-10],[-223,78],[-139,-11],[-113,16],[-126,82],[-99,46],[-190,67],[-12,10],[-9,23],[-3,35],[-9,23],[-13,12],[-30,-1],[-31,-23],[-58,-20],[-90,5],[-39,13],[-29,22],[-16,26],[-3,29],[-8,23],[-14,16],[-32,0],[-51,-15],[-20,-19],[11,-22],[-9,-13],[-89,0],[-35,13],[-67,41],[-29,42],[-38,74],[4,21],[24,45],[31,31],[202,19],[94,18],[102,51],[120,89],[26,26],[3,20],[-9,21],[-37,51],[-13,35],[12,16],[48,2],[-25,20],[-22,28],[-7,17],[1,28],[38,6],[47,-14],[91,-80],[36,-16],[62,-12],[-66,55],[-68,117],[-9,40],[2,23],[18,62],[16,26],[21,17],[65,36],[104,25],[54,3],[54,-20],[47,-39],[105,-64],[16,-26],[-2,-12],[-40,-17],[-5,-16],[16,-24]],[[34626,64748],[35,-13],[63,11],[91,33],[118,27],[144,22],[38,-26],[24,5],[43,42],[7,28],[-5,31],[3,69],[22,41],[86,86],[46,31],[73,14],[174,-9],[163,-49],[220,-49],[323,-122],[101,-52],[10,-45],[-57,-96],[-139,-136],[-111,-49],[-43,-30],[72,-21],[47,-35],[72,50],[52,57],[74,47],[25,4],[7,-10],[-13,-24],[-6,-24],[2,-23],[9,-14],[45,-7],[24,9],[99,66],[96,102],[146,66],[40,33],[128,28],[-2,20],[6,77],[-45,34],[-149,69],[-74,83],[16,63],[82,-9],[226,-7],[47,-8],[216,-108],[77,-67],[61,-33],[128,-49],[43,-39],[30,-16],[10,-17],[-9,-17],[-5,-41],[24,-13],[83,-16],[23,-16],[31,-52],[38,-87],[34,-94],[52,-178],[106,-238],[35,-149],[12,-28],[24,-18],[68,-28],[51,-38],[62,-13],[15,4],[15,33],[38,53],[185,103],[10,16],[-22,23],[-7,17],[3,10],[38,9],[-128,130],[-84,124],[-53,153],[-8,45],[-8,95],[-18,25],[-29,23],[-12,29],[4,37],[-5,31],[-35,64],[-131,453],[1,44],[18,33],[46,19],[76,4],[24,11],[-28,16],[-49,48],[-7,22],[33,45],[168,-20],[122,-40],[208,-96],[21,5],[23,48],[44,31],[67,-11],[188,-69],[218,-124],[146,-61],[103,-83],[70,-79],[44,-60],[1,-24],[-10,-24],[11,-32],[30,-39],[18,-34],[13,-75],[28,-96],[7,-49],[193,-434],[37,-76],[24,-36],[135,-168],[72,-122],[7,-82],[10,-23],[3,-38],[-3,-52],[-16,-44],[-28,-37],[-28,-52],[-41,-116],[-4,-28],[30,-40],[189,-136],[116,-165],[55,-28],[146,-103],[158,-58],[53,-26],[50,-37],[15,-1],[30,7],[9,10],[2,15],[-44,78],[-3,31],[21,5],[163,-133],[87,-53],[120,-55],[206,-128],[29,-11],[111,12],[31,-8],[19,-13],[8,-17],[4,-76],[31,-37],[175,16],[50,-3],[31,-12],[26,-25],[38,-81],[35,-160],[2,-58],[-16,-96],[-26,-35],[-33,-11],[-94,11],[-66,30],[-34,39],[-31,85],[-15,16],[-13,-17],[-31,-78],[-20,-34],[-25,-23],[-46,6],[-69,34],[-130,86],[-45,22],[-29,-4],[-62,-29],[-95,-55],[-39,-40],[16,-26],[11,-32],[6,-39],[-3,-29],[-13,-18],[-31,-22],[-67,-4],[-95,17],[-76,30],[-133,79],[-30,11],[-42,-17],[-16,-24],[26,-32],[67,-43],[82,-71],[23,-14],[21,1],[7,-13],[10,-38],[-5,-66],[-39,-131],[-4,-31],[16,7],[112,129],[58,36],[126,57],[54,42],[160,11],[58,-23],[37,-39],[1,-18],[-42,-47],[-8,-23],[-2,-30],[4,-26],[10,-24],[29,-21],[51,9],[14,-5],[27,-23],[18,-35],[1,-50],[-38,-108],[-67,-35],[-205,-66],[-71,-34],[-136,-23],[-52,-32],[-33,-10],[-145,5],[-167,-20],[-191,40],[-135,17],[-155,63],[-58,-16],[-61,-40],[-290,48],[-35,35],[12,23],[69,74],[4,15],[-3,13],[-132,13],[-148,40],[-147,19],[-111,-5],[-72,14],[-70,33],[-39,29],[-7,27],[-1,29],[6,59],[-9,42],[-32,32],[-65,29],[-65,-3],[-55,-32],[-52,-60],[-97,-166],[-48,-29],[-126,-120],[-47,-30],[-230,-47],[-273,-20],[-102,-37],[-96,-70],[-118,-67],[-286,-83],[-264,-47],[-277,-20],[-207,-31],[-59,15],[-93,-5],[-100,-47],[-112,-9],[-428,-16],[-197,-31],[-107,-8],[-85,4],[-59,13],[-56,41],[-58,63],[-118,168],[-34,71],[14,121],[-9,70],[-39,152],[-8,13],[-206,56],[-136,18],[-204,4],[-250,-8],[-250,18],[-132,20],[-131,34],[-224,87],[-13,9],[-17,29],[-22,49],[-56,64],[-153,143],[-60,85],[-10,22],[-14,62],[-20,102],[-6,63],[18,39],[15,7],[316,75],[557,82],[510,55],[231,-5],[136,-26],[137,-12],[247,-5],[312,-39],[62,3],[140,26],[41,21],[221,-2],[43,12],[39,23],[-50,44],[-212,94],[-560,169],[-137,36],[-196,39],[-114,5],[-144,-22],[-54,1],[-142,-33],[-135,-21],[-256,-20],[-370,-15],[-51,5],[-77,25],[-55,8],[-361,-20],[-324,26],[-368,258],[-61,79],[13,32],[45,35],[183,100],[65,23],[272,54],[271,65],[214,60],[105,22],[101,2],[82,20],[-16,19],[-67,22],[0,32],[35,15],[134,15],[143,-19],[71,7],[20,21],[-19,17],[-135,35],[-649,-102],[-303,-9],[-210,-44],[-115,1],[-138,44],[-18,13],[-2,18],[42,59],[147,35],[74,98],[-79,2],[-263,-21],[-115,9],[-156,38],[-45,44],[-19,33],[-4,39],[5,109],[14,59],[8,14],[193,181],[120,38],[85,57],[3,24],[-21,25],[-79,58],[-31,30],[-18,28],[13,45],[45,61],[131,99],[316,198],[161,83],[155,44],[217,96],[555,158],[497,159],[183,-42],[52,-33],[23,-28],[20,-39],[17,-50],[24,-109],[3,-56],[-4,-57],[-12,-51],[-18,-46],[-38,-55],[-55,-66],[-119,-113],[-13,-33]],[[31790,66264],[5,-4],[140,93],[85,5],[59,-8],[19,-12],[12,-18],[5,-41],[4,-99],[8,-12],[19,5],[31,24],[153,147],[65,40],[44,12],[187,22],[127,0],[141,-14],[105,-21],[171,-60],[135,-70],[124,-74],[416,-279],[176,-82],[67,-46],[30,-34],[27,-43],[8,-40],[-11,-37],[-19,-26],[-42,-22],[-254,-93],[-134,-29],[-133,-41],[-317,-146],[-217,-69],[-282,-136],[-532,-217],[-64,-44],[-29,-31],[-150,-248],[-57,-58],[-138,-59],[-176,-14],[-49,-17],[-8,-85],[-63,-143],[-30,-96],[-42,-256],[-10,-26],[-32,-48],[-54,-51],[-169,-59],[-125,-32],[-170,-26],[-40,18],[-41,41],[-43,3],[-26,-6],[-222,-178],[-213,-73],[-93,-65],[-65,-31],[-53,-9],[-86,6],[-63,29],[-57,46],[-42,50],[-110,204],[-47,69],[-40,35],[-108,124],[-29,25],[-409,158],[-200,88],[-48,31],[-44,18],[-256,-12],[-34,5],[-9,13],[31,46],[11,29],[4,29],[-4,46],[3,8],[98,46],[-16,9],[-12,18],[-7,26],[11,18],[28,11],[34,42],[41,74],[30,42],[42,27],[74,72],[54,29],[45,36],[1,16],[-18,13],[-6,28],[8,86],[0,45],[8,37],[16,30],[22,20],[189,62],[10,18],[2,21],[-5,24],[-10,16],[-30,14],[-51,4],[-44,36],[-10,17],[18,49],[85,79],[29,39],[92,175],[169,108],[45,117],[127,124],[0,17],[-41,42],[-117,29],[-56,45],[-38,51],[-171,295],[-29,21],[-9,35],[-35,22],[7,21],[668,88],[460,28],[476,77],[132,4],[103,-15],[101,-40],[135,-69],[178,-68],[332,-99],[207,-21],[-83,-79],[-12,-23],[0,-18]],[[66369,38971],[-11,-10],[-39,14],[18,31],[14,9],[17,4],[9,-9],[-2,-23],[-6,-16]],[[53734,44150],[-22,-8],[-52,9],[-41,25],[-26,32],[168,87],[35,-11],[1,-16],[-26,-47],[-6,-31],[-13,-24],[-18,-16]],[[54038,48334],[-28,-14],[-20,2],[4,33],[27,66],[16,57],[3,50],[13,44],[20,37],[20,19],[31,-1],[6,-121],[-9,-57],[-20,-47],[-28,-39],[-35,-29]],[[53978,48455],[-44,-101],[-36,-97],[-50,-182],[-28,-7],[-25,44],[71,213],[3,24],[-3,21],[-22,31],[-21,-36],[-101,-241],[-26,-38],[-22,-22],[-17,-4],[-43,5],[-86,-69],[145,286],[1,22],[-27,13],[-11,-7],[-117,-180],[-67,-69],[-46,21],[-11,18],[4,21],[115,182],[105,130],[44,82],[18,77],[8,57],[-1,62],[6,17],[6,-3],[7,-24],[2,-65],[-25,-135],[-19,-66],[-23,-54],[10,-12],[42,31],[36,65],[29,100],[18,86],[18,137],[6,-4],[9,-28],[19,-19],[30,-11],[17,-18],[13,-46],[12,-20],[45,-17],[17,-15],[13,-50],[-1,-27],[6,-17],[13,-7],[-16,-54]],[[53472,48962],[-15,-17],[-27,19],[-1,50],[26,39],[20,-2],[19,-20],[-5,-27],[-17,-42]],[[53554,49702],[-9,-8],[-23,7],[-9,-66],[-9,-6],[-16,40],[12,35],[-2,24],[5,17],[23,40],[14,9],[8,-4],[7,-49],[-1,-39]],[[53412,48396],[-28,-8],[-33,15],[17,73],[29,30],[71,31],[12,19],[23,9],[33,0],[37,27],[41,55],[14,8],[-29,-80],[-30,-60],[-157,-119]],[[53245,51806],[-18,-3],[-4,12],[14,41],[23,5],[25,45],[24,-16],[-9,-25],[-31,-38],[-24,-21]],[[53433,45499],[0,-30],[-36,6],[-19,17],[-15,27],[-4,20],[14,20],[40,-12],[20,-48]],[[53365,51952],[-56,-7],[24,60],[21,28],[23,18],[47,7],[31,-24],[-28,-43],[-62,-39]],[[53662,48845],[-20,-13],[-13,1],[15,88],[-17,31],[-1,17],[7,14],[10,3],[22,-27],[11,-29],[5,-28],[-1,-28],[-7,-19],[-11,-10]],[[44273,65091],[-88,-10],[-81,70],[-2,65],[5,36],[10,32],[30,27],[83,32],[37,-24],[14,-28],[13,-11],[53,-24],[26,-28],[-3,-32],[-17,-49],[-19,-32],[-21,-14],[-40,-10]],[[43881,66647],[-83,-30],[-36,23],[-5,7],[124,0]],[[43464,66021],[-157,-21],[-72,9],[-38,-28],[-30,-11],[-86,-6],[-176,50],[-47,18],[-18,16],[8,15],[31,15],[135,22],[51,17],[20,23],[34,19],[47,14],[127,13],[285,70],[142,10],[55,-6],[17,-18],[4,-19],[-7,-18],[-45,-47],[-57,-37],[-156,-80],[-67,-20]],[[47853,61581],[12,-62],[-49,-84],[-15,-15],[-19,-6],[-18,10],[-55,66],[-15,41],[20,21],[42,26],[30,12],[38,-11],[13,18],[16,-16]],[[47694,61383],[-45,-12],[-28,42],[-22,8],[-10,28],[-47,5],[3,44],[13,21],[43,17],[34,-7],[33,-43],[16,-34],[14,-44],[-4,-25]],[[56662,54794],[-29,-9],[-108,12],[-134,47],[-69,44],[3,14],[31,7],[34,-7],[58,-31],[155,-15],[50,-18],[16,-26],[-7,-18]],[[56184,60512],[-43,-6],[-62,62],[-129,70],[-49,51],[-2,24],[4,39],[14,47],[45,52],[48,8],[68,-11],[50,-37],[53,-103],[36,-50],[11,-37],[-17,-17],[1,-17],[8,-10],[-3,-19],[-14,-28],[-19,-18]],[[54198,52908],[-74,-12],[-1,15],[31,41],[116,35],[86,11],[-20,-35],[-51,-27],[-87,-28]],[[53953,60357],[-44,-11],[-60,53],[-1,29],[12,71],[109,18],[46,-43],[22,-42],[-84,-75]],[[55206,55568],[-57,-10],[-75,22],[-74,44],[-167,138],[125,93],[202,-108],[60,-70],[-14,-109]],[[61653,54734],[104,-10],[63,3],[28,-15],[25,-49],[-32,-69],[-37,-27],[-61,-7],[-98,22],[-34,15],[-30,37],[14,27],[48,9],[8,11],[-12,20],[0,18],[14,15]],[[62818,59223],[-68,-28],[-36,34],[26,9],[37,41],[50,35],[21,29],[85,13],[29,-1],[11,-11],[-48,-44],[-107,-77]],[[59798,52421],[-49,-8],[-22,21],[-2,61],[16,46],[57,95],[50,106],[30,26],[59,-18],[35,-30],[36,-51],[16,-39],[-14,-58],[-42,-52],[-49,-35],[-121,-64]],[[58655,54724],[-38,-3],[-74,7],[-79,21],[-43,23],[-39,52],[-7,56],[-73,83],[-82,28],[-46,58],[47,4],[66,-13],[97,-25],[87,-33],[126,-73],[41,-68],[41,-49],[13,-38],[-11,-17],[-26,-13]],[[61648,53544],[-13,-11],[-13,2],[-40,53],[-55,22],[-20,21],[-164,110],[-18,48],[-3,38],[55,19],[109,17],[96,0],[89,-23],[18,-25],[48,-44],[-12,-54],[-3,-68],[-20,-33],[-30,-25],[-24,-47]],[[59972,61702],[-14,-6],[-143,46],[-10,36],[70,40],[56,22],[44,3],[43,-7],[41,-44],[-49,-39],[-38,-51]],[[53709,61949],[22,-58],[14,-18],[-20,-27],[-82,-54],[-179,-22],[-90,25],[41,-77],[9,-33],[-12,-14],[-38,4],[-62,22],[-36,27],[-9,30],[-13,7],[-17,-17],[-17,4],[-37,46],[-27,17],[-180,27],[-9,12],[10,21],[28,30],[40,10],[101,-13],[9,8],[6,37],[8,15],[70,-3],[44,8],[25,-20],[24,-45],[34,8],[50,-6],[55,16],[84,40],[65,15],[89,-22]],[[54471,61876],[28,-50],[5,-26],[-38,-31],[-146,-56],[-88,-49],[-45,-11],[-60,11],[-72,-23],[-29,3],[32,40],[115,116],[96,12],[31,24],[27,-8],[15,20],[3,30],[33,22],[30,0],[63,-24]],[[54682,56165],[-39,-47],[-117,17],[-16,14],[-4,17],[19,21],[120,21],[50,2],[26,-7],[4,-8],[-43,-30]],[[51701,58451],[54,-12],[41,1],[10,-15],[-44,-48],[-26,-10],[-48,35],[-36,43],[-10,26],[-3,28],[8,5],[54,-53]],[[53829,61008],[-38,-7],[-44,19],[-16,33],[-8,33],[8,16],[20,16],[26,37],[34,57],[54,38],[116,35],[17,12],[52,100],[18,17],[59,10],[7,13],[-21,24],[0,26],[21,26],[29,19],[75,19],[68,-3],[18,-8],[15,-17],[21,-47],[3,-10],[-32,-42],[-82,-63],[-52,-53],[-10,-17],[-4,-23],[-20,-27],[-59,-65],[-40,-63],[-40,-35],[-110,-33],[-85,-37]],[[55034,61306],[-69,-12],[-51,6],[-33,24],[-25,32],[-31,80],[11,37],[4,67],[6,25],[12,13],[71,24],[43,-3],[62,-26],[135,-7],[34,-24],[8,-14],[0,-18],[-10,-21],[-67,-60],[-32,-42],[-23,-50],[-45,-31]],[[61540,54056],[12,-9],[15,11],[11,-9],[8,-29],[12,-19],[32,-26],[11,-18],[-1,-18],[-27,-28],[-17,-1],[-129,65],[-35,66],[-3,34],[13,30],[20,15],[27,-1],[31,-16],[20,-47]],[[49640,62273],[62,-8],[59,10],[43,-9],[26,-29],[18,-29],[9,-31],[-24,-21],[-96,-19],[-65,9],[-70,25],[-33,-11],[-80,19],[-40,25],[-32,36],[0,22],[85,25],[33,20],[105,-34]],[[51373,57966],[70,-39],[70,-23],[112,-12],[16,-8],[0,-18],[-16,-29],[-38,-37],[-26,-1],[-62,26],[-23,13],[-25,30],[-13,3],[-16,-10],[-4,-13],[7,-16],[-10,-5],[-75,11],[-12,10],[6,31],[53,46],[-47,14],[-14,17],[-68,-29],[-38,-6],[-60,20],[-6,105],[-7,39],[-28,26],[-16,27],[-26,21],[-54,21],[-44,53],[-9,24],[6,18],[27,24],[155,-53],[93,-50],[89,-62],[47,-45],[4,-30],[-12,-29],[-27,-30],[21,-34]],[[49813,59901],[-23,-1],[-37,15],[-85,60],[-17,26],[-8,31],[0,35],[7,36],[26,72],[-49,57],[-12,32],[6,19],[27,44],[8,28],[27,39],[72,75],[72,-17],[64,-63],[17,-42],[-5,-44],[5,-65],[16,-85],[5,-61],[-8,-38],[-27,-72],[-24,-34],[-30,-31],[-27,-16]],[[50857,57744],[-28,-11],[-30,6],[-26,29],[-21,52],[-34,38],[-77,48],[-13,19],[-22,65],[-4,64],[-16,58],[-1,29],[14,43],[65,10],[50,-17],[10,-11],[17,-25],[11,-31],[61,-82],[35,-66],[49,-135],[0,-26],[-13,-28],[-27,-29]],[[46306,54021],[-22,-2],[-50,51],[-11,26],[66,16],[45,-45],[-4,-22],[-24,-24]],[[66746,43735],[-41,-99],[-27,-49],[-26,-15],[-54,-15],[-115,-15],[-49,-15],[-7,-66],[8,-35],[16,-28],[22,-7],[47,16],[18,-3],[14,-14],[10,-25],[6,-34],[0,-45],[-8,-53],[-39,-126],[-49,-69],[-63,-57],[-16,-22],[-8,-25],[-8,-83],[-32,-66],[-102,-167],[-39,-37],[0,-30],[-16,-79],[-30,-64],[-84,-147],[-21,-52],[-9,-41],[1,-58],[-4,-25],[-19,-49],[-27,-46],[-5,-23],[10,-40],[11,-14],[1,-37],[-8,-60],[36,38],[79,136],[60,82],[40,27],[28,36],[29,80],[41,77],[37,24],[17,-14],[14,-38],[-3,-47],[-20,-58],[2,-17],[46,42],[81,35],[29,-4],[59,-53],[49,6],[79,31],[15,-14],[-14,-46],[-30,-45],[-73,-63],[-176,-124],[-54,-85],[10,2],[38,37],[40,19],[42,3],[17,-11],[-7,-24],[-5,-65],[-106,-128],[25,5],[123,58],[74,-81],[103,28],[62,27],[-2,-16],[13,-36],[0,-56],[6,-8],[29,19],[6,21],[-2,100],[9,11],[20,-16],[12,-26],[4,-74],[-13,-73],[-19,-68],[-45,-99],[6,-42],[-12,-46],[10,-2],[45,43],[3,18],[-4,41],[5,20],[37,45],[62,52],[20,8],[8,-12],[-3,-33],[18,8],[40,48],[36,29],[34,10],[35,33],[36,56],[39,46],[42,37],[18,2],[-7,-60],[9,-69],[1,-58],[8,-13],[32,63],[18,23],[22,9],[25,-5],[170,23],[52,-16],[58,-41],[74,-63],[27,-58],[5,-72],[-7,-50],[-53,-64],[-48,-42],[-27,-42],[-9,-42],[-11,-26],[-32,-36],[-140,-101],[34,-3],[81,22],[53,4],[3,-14],[-22,-29],[-40,-29],[-5,-14],[2,-18],[44,-22],[55,11],[47,-16],[-4,-24],[-38,-79],[-10,-50],[-50,-42],[-99,-65],[-25,-25],[5,-6],[91,48],[46,13],[29,0],[33,47],[51,15],[50,-29],[77,80],[27,10],[47,-9],[30,14],[50,55],[39,26],[7,-2],[9,-22],[3,-63],[-9,-55],[-12,-36],[-41,-78],[-25,-28],[-24,-10],[-41,4],[-18,-12],[-39,-62],[-68,-62],[-43,-24],[27,-34],[10,-65],[-15,-20],[-73,-20],[-4,-12],[-26,-14],[-60,-23],[41,-10],[77,16],[8,-11],[-11,-46],[-20,-47],[-91,-122],[-1,-12],[14,-59],[18,-46],[22,-31],[50,-2],[37,14],[55,81],[120,254],[107,70],[90,79],[20,-16],[10,-19],[-4,-19],[-45,-64],[-23,-53],[-62,-164],[-23,-77],[-12,-82],[3,-140],[7,-24],[19,-33],[36,28],[61,69],[39,67],[30,108],[19,42],[20,-1],[19,-23],[4,-53],[17,-71],[10,-72],[-8,-79],[-8,-44],[-125,-322],[13,-57],[4,-35],[-4,-37],[-40,-154],[-38,-97],[-20,-41],[-25,-26],[-28,-9],[-26,14],[-22,37],[-20,19],[-17,2],[-32,-7],[-84,-79],[-17,-4],[-12,10],[-15,42],[12,208],[8,68],[-17,53],[17,89],[1,34],[-10,13],[-21,-9],[-32,-44],[-43,-81],[-47,-73],[-81,-99],[-37,-20],[-14,5],[-16,14],[-23,41],[1,37],[9,50],[34,119],[67,175],[55,125],[11,53],[-16,23],[-12,46],[-21,135],[-27,111],[-32,50],[-78,54],[-15,-13],[-8,-75],[-93,-215],[-16,-94],[-12,-34],[-18,-24],[-40,-29],[11,50],[43,111],[-6,11],[-55,-89],[-41,-50],[-51,-13],[-31,4],[-30,-14],[-127,-211],[-5,-70],[-22,-57],[-63,-104],[-33,-36],[-46,-7],[-43,18],[-28,-3],[-66,-32],[-74,-15],[-30,7],[-20,13],[-38,41],[-4,28],[2,17],[19,44],[44,54],[36,21],[88,27],[65,41],[49,61],[22,37],[92,190],[117,67],[58,55],[40,69],[6,24],[-58,-35],[-30,-9],[-48,13],[-22,24],[-66,-7],[-92,11],[-13,-19],[-12,-92],[-12,-49],[-14,-16],[-21,-10],[-42,-10],[-108,33],[-20,18],[-28,14],[-119,-30],[-25,3],[24,21],[117,68],[13,195],[-8,31],[-32,-27],[-56,-28],[-38,8],[-17,17],[-16,-14],[-38,-103],[-23,-13],[-34,-4],[-73,-38],[-143,-24],[-28,-27],[-97,9],[-283,57],[-101,-6],[-122,34],[-23,15],[-171,-6],[-51,8],[4,43],[-6,11],[-49,-47],[-44,-31],[-57,-26],[-178,-46],[-96,-10],[-55,32],[-22,32],[-34,102],[-22,129],[0,23],[11,45],[38,61],[169,164],[136,165],[58,86],[55,31],[91,71],[3,9],[-88,-9],[-62,20],[-63,7],[-121,-19],[-121,0],[0,37],[56,69],[121,119],[12,0],[-37,-55],[-10,-41],[16,-28],[18,-17],[69,-7],[16,24],[25,127],[52,147],[28,106],[49,81],[25,12],[21,-15],[72,-20],[75,-74],[24,-6],[8,7],[-27,21],[-23,35],[-9,33],[27,102],[32,30],[5,20],[-63,0],[-51,29],[-15,46],[3,82],[16,48],[41,64],[50,43],[30,-11],[57,-57],[35,16],[-4,18],[-54,92],[-17,68],[2,32],[116,323],[57,174],[77,264],[18,42],[39,78],[17,22],[50,0],[32,9],[-46,35],[-17,24],[-2,25],[13,26],[18,20],[60,42],[42,70],[26,82],[-4,28],[-13,28],[1,15],[32,17],[82,98],[11,19],[30,130],[37,58],[35,29],[55,37],[168,89],[99,80],[67,-5],[20,-55],[96,-37],[16,40],[-22,48],[19,19],[78,18],[14,-7],[24,-27],[-3,-26]],[[66703,42920],[-18,-11],[-17,0],[-15,12],[-3,19],[16,41],[42,20],[32,-5],[-2,-20],[-17,-34],[-18,-22]],[[67237,41791],[-84,-58],[-19,-23],[-23,-11],[-17,18],[-25,62],[5,20],[23,3],[13,-8],[3,-19],[10,-10],[17,0],[62,59],[33,10],[11,-12],[-9,-31]],[[66799,44087],[-26,0],[-6,11],[11,38],[29,44],[39,13],[-11,-65],[-36,-41]],[[67488,41947],[40,-65],[21,-18],[-140,-70],[-17,-5],[-9,9],[-1,65],[6,51],[10,7],[32,-30],[34,62],[24,-6]],[[66306,39047],[-47,-43],[-17,24],[4,28],[24,68],[0,20],[-29,132],[5,21],[7,10],[42,-28],[5,-36],[-20,-81],[14,-54],[18,-39],[-6,-22]],[[67415,39648],[-26,-34],[-27,1],[3,30],[33,59],[17,43],[1,25],[6,20],[25,22],[22,39],[-10,-73],[-44,-132]],[[42405,60969],[-17,-31],[-21,-9],[-43,-37],[-17,-5],[-24,23],[-21,39],[-9,6],[-13,-2],[-28,-24],[-13,0],[-11,17],[-5,33],[1,49],[13,76],[1,27],[-7,22],[6,19],[19,17],[24,8],[58,-7],[46,-30],[22,-40],[45,-31],[15,-23],[-21,-97]],[[42526,61176],[-12,-60],[-67,16],[-30,22],[-25,49],[-3,12],[5,19],[28,42],[18,16],[44,-18],[20,-25],[16,-40],[6,-33]],[[42356,62655],[-7,-8],[-118,37],[-45,22],[-15,16],[-10,33],[-7,50],[24,24],[53,-2],[56,-19],[86,-54],[-24,-17],[-1,-35],[9,-34],[-1,-13]],[[44963,61735],[72,-67],[0,-33],[-10,-54],[-20,-41],[-32,-28],[-45,-15],[-58,0],[-25,12],[19,42],[14,13],[2,36],[-9,62],[-11,37],[-33,22],[-23,1],[-6,-20],[10,-41],[-7,-53],[-25,-64],[-19,-31],[-34,16],[-17,24],[5,41],[-10,37],[9,39],[23,58],[33,39],[43,20],[49,-2],[57,-22],[48,-28]],[[41886,61559],[-44,-9],[-8,2],[3,27],[-3,14],[-12,10],[39,21],[6,17],[-13,13],[-53,25],[-15,22],[2,19],[21,17],[37,-2],[80,-33],[37,-45],[16,-34],[-26,-5],[-21,-13],[-21,-29],[-25,-17]],[[41520,60750],[-23,-2],[-31,18],[-175,61],[-23,18],[21,26],[63,32],[43,34],[33,51],[101,-26],[38,-25],[14,-21],[6,-29],[-6,-72],[-33,-15],[-28,-50]],[[40054,60569],[-30,-3],[-56,16],[-83,35],[-62,38],[-41,40],[-6,27],[32,15],[47,8],[112,-5],[54,-15],[70,-58],[15,-33],[3,-21],[-9,-18],[-46,-26]],[[38227,59567],[-27,-83],[-11,8],[-18,40],[-38,19],[-43,44],[1,95],[16,43],[-4,61],[43,36],[32,-43],[8,-72],[-7,-36],[31,-37],[14,-9],[8,-30],[-5,-36]],[[37539,60147],[61,-11],[45,8],[33,-40],[13,-41],[-4,-14],[-15,-6],[-95,46],[-36,27],[-12,23],[10,8]],[[38122,59172],[69,-8],[88,1],[-16,-77],[-33,-49],[-27,-14],[-11,27],[-52,61],[-18,59]],[[37453,60155],[-20,-3],[-75,59],[16,51],[69,-52],[10,-32],[0,-23]],[[63735,48228],[7,-24],[-84,18],[-30,15],[-2,16],[4,25],[16,29],[39,23],[21,-12],[44,-30],[6,-19],[-21,-41]],[[63328,49741],[45,-30],[12,-109],[-85,7],[-98,73],[-20,52],[2,10],[13,9],[25,-13],[24,19],[25,6],[57,-24]],[[59295,51223],[-33,-72],[-44,9],[-16,-15],[-12,-1],[20,68],[2,46],[-10,46],[15,24],[58,1],[1,-53],[7,-20],[14,-10],[-2,-23]],[[61879,52547],[-19,-69],[-63,26],[-98,52],[-38,34],[-14,38],[-1,53],[29,9],[74,4],[62,-73],[18,-11],[50,-63]],[[50207,8732],[-42,-57],[-7,26],[20,53],[25,31],[24,-3],[-6,-23],[-14,-27]],[[4589,48823],[72,-20],[40,13],[35,-10],[6,-24],[-61,-51],[-25,3],[-73,62],[6,27]],[[57389,12301],[22,-23],[-8,-29],[-55,18],[-57,32],[-19,-8],[-11,4],[-33,31],[29,23],[30,8],[33,-2],[69,-54]],[[6572,58396],[-21,-6],[1,15],[62,41],[114,50],[-4,-11],[-61,-39],[-91,-50]],[[49748,13500],[10,-48],[-9,6],[-13,36],[-9,46],[8,1],[6,-16],[7,-25]],[[52430,16799],[-14,-2],[-1,15],[23,19],[16,-1],[0,-21],[-24,-10]]],"transform":{"scale":[0.001838860850736404,0.0010020709717534302],"translate":[-178.19453125,7.710991655433217]},"objects":{"ne_50m_land":{"type":"GeometryCollection","geometries":[{"arcs":[[0]],"type":"Polygon"},{"arcs":[[1]],"type":"Polygon"},{"arcs":[[2]],"type":"Polygon"},{"arcs":[[3]],"type":"Polygon"},{"arcs":[[4]],"type":"Polygon"},{"arcs":[[5]],"type":"Polygon"},{"arcs":[[6]],"type":"Polygon"},{"arcs":[[7]],"type":"Polygon"},{"arcs":[[8]],"type":"Polygon"},{"arcs":[[9]],"type":"Polygon"},{"arcs":[[10]],"type":"Polygon"},{"arcs":[[11]],"type":"Polygon"},{"arcs":[[12]],"type":"Polygon"},{"arcs":[[13]],"type":"Polygon"},{"arcs":[[14]],"type":"Polygon"},{"arcs":[[15]],"type":"Polygon"},{"arcs":[[16]],"type":"Polygon"},{"arcs":[[17]],"type":"Polygon"},{"arcs":[[18]],"type":"Polygon"},{"arcs":[[19]],"type":"Polygon"},{"arcs":[[20]],"type":"Polygon"},{"arcs":[[21]],"type":"Polygon"},{"arcs":[[22]],"type":"Polygon"},{"arcs":[[23]],"type":"Polygon"},{"arcs":[[24]],"type":"Polygon"},{"arcs":[[25]],"type":"Polygon"},{"arcs":[[26]],"type":"Polygon"},{"arcs":[[27]],"type":"Polygon"},{"arcs":[[28]],"type":"Polygon"},{"arcs":[[29]],"type":"Polygon"},{"arcs":[[30]],"type":"Polygon"},{"arcs":[[31]],"type":"Polygon"},{"arcs":[[32]],"type":"Polygon"},{"arcs":[[33]],"type":"Polygon"},{"arcs":[[34]],"type":"Polygon"},{"arcs":[[35]],"type":"Polygon"},{"arcs":[[36]],"type":"Polygon"},{"arcs":[[37]],"type":"Polygon"},{"arcs":[[38]],"type":"Polygon"},{"arcs":[[39]],"type":"Polygon"},{"arcs":[[40]],"type":"Polygon"},{"arcs":[[41]],"type":"Polygon"},{"arcs":[[42]],"type":"Polygon"},{"arcs":[[43]],"type":"Polygon"},{"arcs":[[44]],"type":"Polygon"},{"arcs":[[45]],"type":"Polygon"},{"arcs":[[46]],"type":"Polygon"},{"arcs":[[47]],"type":"Polygon"},{"arcs":[[48]],"type":"Polygon"},{"arcs":[[49]],"type":"Polygon"},{"arcs":[[50]],"type":"Polygon"},{"arcs":[[51]],"type":"Polygon"},{"arcs":[[52]],"type":"Polygon"},{"arcs":[[53]],"type":"Polygon"},{"arcs":[[54]],"type":"Polygon"},{"arcs":[[55]],"type":"Polygon"},{"arcs":[[56]],"type":"Polygon"},{"arcs":[[57]],"type":"Polygon"},{"arcs":[[58]],"type":"Polygon"},{"arcs":[[59]],"type":"Polygon"},{"arcs":[[60]],"type":"Polygon"},{"arcs":[[61]],"type":"Polygon"},{"arcs":[[62]],"type":"Polygon"},{"arcs":[[63]],"type":"Polygon"},{"arcs":[[64]],"type":"Polygon"},{"arcs":[[65]],"type":"Polygon"},{"arcs":[[66]],"type":"Polygon"},{"arcs":[[67]],"type":"Polygon"},{"arcs":[[68]],"type":"Polygon"},{"arcs":[[69]],"type":"Polygon"},{"arcs":[[70]],"type":"Polygon"},{"arcs":[[71]],"type":"Polygon"},{"arcs":[[72]],"type":"Polygon"},{"arcs":[[73]],"type":"Polygon"},{"arcs":[[74]],"type":"Polygon"},{"arcs":[[75]],"type":"Polygon"},{"arcs":[[76]],"type":"Polygon"},{"arcs":[[77]],"type":"Polygon"},{"arcs":[[78]],"type":"Polygon"},{"arcs":[[79]],"type":"Polygon"},{"arcs":[[80]],"type":"Polygon"},{"arcs":[[81]],"type":"Polygon"},{"arcs":[[82]],"type":"Polygon"},{"arcs":[[83]],"type":"Polygon"},{"arcs":[[84]],"type":"Polygon"},{"arcs":[[85]],"type":"Polygon"},{"arcs":[[86]],"type":"Polygon"},{"arcs":[[87]],"type":"Polygon"},{"arcs":[[88]],"type":"Polygon"},{"arcs":[[89]],"type":"Polygon"},{"arcs":[[90]],"type":"Polygon"},{"arcs":[[91]],"type":"Polygon"},{"arcs":[[92]],"type":"Polygon"},{"arcs":[[93]],"type":"Polygon"},{"arcs":[[94]],"type":"Polygon"},{"arcs":[[95]],"type":"Polygon"},{"arcs":[[96]],"type":"Polygon"},{"arcs":[[97]],"type":"Polygon"},{"arcs":[[98]],"type":"Polygon"},{"arcs":[[99]],"type":"Polygon"},{"arcs":[[100]],"type":"Polygon"},{"arcs":[[101]],"type":"Polygon"},{"arcs":[[102]],"type":"Polygon"},{"arcs":[[103]],"type":"Polygon"},{"arcs":[[104]],"type":"Polygon"},{"arcs":[[105]],"type":"Polygon"},{"arcs":[[106]],"type":"Polygon"},{"arcs":[[107]],"type":"Polygon"},{"arcs":[[108]],"type":"Polygon"},{"arcs":[[109]],"type":"Polygon"},{"arcs":[[110]],"type":"Polygon"},{"arcs":[[111]],"type":"Polygon"},{"arcs":[[112]],"type":"Polygon"},{"arcs":[[113]],"type":"Polygon"},{"arcs":[[114]],"type":"Polygon"},{"arcs":[[115]],"type":"Polygon"},{"arcs":[[116]],"type":"Polygon"},{"arcs":[[117]],"type":"Polygon"},{"arcs":[[118]],"type":"Polygon"},{"arcs":[[119]],"type":"Polygon"},{"arcs":[[120]],"type":"Polygon"},{"arcs":[[121]],"type":"Polygon"},{"arcs":[[122]],"type":"Polygon"},{"arcs":[[123]],"type":"Polygon"},{"arcs":[[124]],"type":"Polygon"},{"arcs":[[125]],"type":"Polygon"},{"arcs":[[126]],"type":"Polygon"},{"arcs":[[127]],"type":"Polygon"},{"arcs":[[128]],"type":"Polygon"},{"arcs":[[129]],"type":"Polygon"},{"arcs":[[130]],"type":"Polygon"},{"arcs":[[131]],"type":"Polygon"},{"arcs":[[132]],"type":"Polygon"},{"arcs":[[133]],"type":"Polygon"},{"arcs":[[134]],"type":"Polygon"},{"arcs":[[135]],"type":"Polygon"},{"arcs":[[136]],"type":"Polygon"},{"arcs":[[137]],"type":"Polygon"},{"arcs":[[138]],"type":"Polygon"},{"arcs":[[139]],"type":"Polygon"},{"arcs":[[140]],"type":"Polygon"},{"arcs":[[141]],"type":"Polygon"},{"arcs":[[142]],"type":"Polygon"},{"arcs":[[143]],"type":"Polygon"},{"arcs":[[144]],"type":"Polygon"},{"arcs":[[145]],"type":"Polygon"},{"arcs":[[146]],"type":"Polygon"},{"arcs":[[147]],"type":"Polygon"},{"arcs":[[148]],"type":"Polygon"},{"arcs":[[149]],"type":"Polygon"},{"arcs":[[150]],"type":"Polygon"},{"arcs":[[151]],"type":"Polygon"},{"arcs":[[152]],"type":"Polygon"},{"arcs":[[153]],"type":"Polygon"},{"arcs":[[154]],"type":"Polygon"},{"arcs":[[155]],"type":"Polygon"},{"arcs":[[156]],"type":"Polygon"},{"arcs":[[157]],"type":"Polygon"},{"arcs":[[158]],"type":"Polygon"},{"arcs":[[159]],"type":"Polygon"},{"arcs":[[160]],"type":"Polygon"},{"arcs":[[161]],"type":"Polygon"},{"arcs":[[162]],"type":"Polygon"},{"arcs":[[163]],"type":"Polygon"},{"arcs":[[164]],"type":"Polygon"},{"arcs":[[165]],"type":"Polygon"},{"arcs":[[166]],"type":"Polygon"},{"arcs":[[167]],"type":"Polygon"},{"arcs":[[168]],"type":"Polygon"},{"arcs":[[169]],"type":"Polygon"},{"arcs":[[170]],"type":"Polygon"},{"arcs":[[171]],"type":"Polygon"},{"arcs":[[172]],"type":"Polygon"},{"arcs":[[173]],"type":"Polygon"},{"arcs":[[174]],"type":"Polygon"},{"arcs":[[175]],"type":"Polygon"},{"arcs":[[176]],"type":"Polygon"},{"arcs":[[177]],"type":"Polygon"},{"arcs":[[178]],"type":"Polygon"},{"arcs":[[179]],"type":"Polygon"},{"arcs":[[180]],"type":"Polygon"},{"arcs":[[181]],"type":"Polygon"},{"arcs":[[182]],"type":"Polygon"},{"arcs":[[183]],"type":"Polygon"},{"arcs":[[184]],"type":"Polygon"},{"arcs":[[185]],"type":"Polygon"},{"arcs":[[186]],"type":"Polygon"},{"arcs":[[187]],"type":"Polygon"},{"arcs":[[188]],"type":"Polygon"},{"arcs":[[189]],"type":"Polygon"},{"arcs":[[190]],"type":"Polygon"},{"arcs":[[191]],"type":"Polygon"},{"arcs":[[192]],"type":"Polygon"},{"arcs":[[193]],"type":"Polygon"},{"arcs":[[194]],"type":"Polygon"},{"arcs":[[195]],"type":"Polygon"},{"arcs":[[196]],"type":"Polygon"},{"arcs":[[197]],"type":"Polygon"},{"arcs":[[198]],"type":"Polygon"},{"arcs":[[199]],"type":"Polygon"},{"arcs":[[200]],"type":"Polygon"},{"arcs":[[201]],"type":"Polygon"},{"arcs":[[202]],"type":"Polygon"},{"arcs":[[203]],"type":"Polygon"},{"arcs":[[204]],"type":"Polygon"},{"arcs":[[205]],"type":"Polygon"},{"arcs":[[206]],"type":"Polygon"},{"arcs":[[207]],"type":"Polygon"},{"arcs":[[208]],"type":"Polygon"},{"arcs":[[209]],"type":"Polygon"},{"arcs":[[210]],"type":"Polygon"},{"arcs":[[211]],"type":"Polygon"},{"arcs":[[212]],"type":"Polygon"},{"arcs":[[213]],"type":"Polygon"},{"arcs":[[214]],"type":"Polygon"},{"arcs":[[215]],"type":"Polygon"},{"arcs":[[216]],"type":"Polygon"},{"arcs":[[217]],"type":"Polygon"},{"arcs":[[218]],"type":"Polygon"},{"arcs":[[219]],"type":"Polygon"},{"arcs":[[220]],"type":"Polygon"},{"arcs":[[221]],"type":"Polygon"},{"arcs":[[222]],"type":"Polygon"},{"arcs":[[223]],"type":"Polygon"},{"arcs":[[224]],"type":"Polygon"},{"arcs":[[225]],"type":"Polygon"},{"arcs":[[226]],"type":"Polygon"},{"arcs":[[227]],"type":"Polygon"},{"arcs":[[228]],"type":"Polygon"},{"arcs":[[229]],"type":"Polygon"},{"arcs":[[230]],"type":"Polygon"},{"arcs":[[231]],"type":"Polygon"},{"arcs":[[232]],"type":"Polygon"},{"arcs":[[233]],"type":"Polygon"},{"arcs":[[234]],"type":"Polygon"},{"arcs":[[235]],"type":"Polygon"},{"arcs":[[236]],"type":"Polygon"},{"arcs":[[237]],"type":"Polygon"},{"arcs":[[238]],"type":"Polygon"},{"arcs":[[239]],"type":"Polygon"},{"arcs":[[240]],"type":"Polygon"},{"arcs":[[241]],"type":"Polygon"},{"arcs":[[242]],"type":"Polygon"},{"arcs":[[243]],"type":"Polygon"},{"arcs":[[244]],"type":"Polygon"},{"arcs":[[245]],"type":"Polygon"},{"arcs":[[246]],"type":"Polygon"},{"arcs":[[247]],"type":"Polygon"},{"arcs":[[248]],"type":"Polygon"},{"arcs":[[249]],"type":"Polygon"},{"arcs":[[250]],"type":"Polygon"},{"arcs":[[251]],"type":"Polygon"},{"arcs":[[252]],"type":"Polygon"},{"arcs":[[253]],"type":"Polygon"},{"arcs":[[254]],"type":"Polygon"},{"arcs":[[255]],"type":"Polygon"},{"arcs":[[256]],"type":"Polygon"},{"arcs":[[257]],"type":"Polygon"},{"arcs":[[258]],"type":"Polygon"},{"arcs":[[259]],"type":"Polygon"},{"arcs":[[260]],"type":"Polygon"},{"arcs":[[261]],"type":"Polygon"},{"arcs":[[262]],"type":"Polygon"},{"arcs":[[263]],"type":"Polygon"},{"arcs":[[264]],"type":"Polygon"},{"arcs":[[265]],"type":"Polygon"},{"arcs":[[266]],"type":"Polygon"},{"arcs":[[267]],"type":"Polygon"},{"arcs":[[268]],"type":"Polygon"},{"arcs":[[269]],"type":"Polygon"},{"arcs":[[270]],"type":"Polygon"},{"arcs":[[271]],"type":"Polygon"},{"arcs":[[272]],"type":"Polygon"},{"arcs":[[273]],"type":"Polygon"},{"arcs":[[274]],"type":"Polygon"},{"arcs":[[275]],"type":"Polygon"},{"arcs":[[276]],"type":"Polygon"},{"arcs":[[277]],"type":"Polygon"},{"arcs":[[278]],"type":"Polygon"},{"arcs":[[279]],"type":"Polygon"},{"arcs":[[280]],"type":"Polygon"},{"arcs":[[281]],"type":"Polygon"},{"arcs":[[282]],"type":"Polygon"},{"arcs":[[283]],"type":"Polygon"},{"arcs":[[284]],"type":"Polygon"},{"arcs":[[285]],"type":"Polygon"},{"arcs":[[286]],"type":"Polygon"},{"arcs":[[287]],"type":"Polygon"},{"arcs":[[288]],"type":"Polygon"},{"arcs":[[289]],"type":"Polygon"},{"arcs":[[290]],"type":"Polygon"}]}}}
},{}],16:[function(_dereq_,module,exports){
"use strict";

var World = _dereq_('./World'); //


var somehow = function somehow(obj) {
  return new World(obj);
};

module.exports = somehow;

},{"./World":8}],17:[function(_dereq_,module,exports){
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

},{"./Shape":19,"spencer-color":5}],18:[function(_dereq_,module,exports){
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

},{"../data":9,"./Shape":19,"d3-geo":2,"spencer-color":5}],19:[function(_dereq_,module,exports){
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

},{"../data":9,"d3-geo":2,"spencer-color":5,"topojson-client":6}],20:[function(_dereq_,module,exports){
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

},{"./Shape":19}]},{},[16])(16)
});
