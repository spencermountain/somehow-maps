const data = require('./data')

exports.bounds = function(arr) {
  let min = null
  let max = null
  arr.forEach(a => {
    if (min === null || a < min) {
      min = a
    }
    if (max === null || a > max) {
      max = a
    }
  })
  return {
    min: min,
    max: max
  }
}

exports.parsePoint = function(input) {
  if (typeof input === 'string') {
    return data.points[input].reverse()
  } else {
    return input
  }
}

exports.parseBounds = function(input) {
  if (typeof input === 'string') {
    if (data.bounds.hasOwnProperty(input)) {
      return data.bounds[input]
    }
    return data.points[input]
  } else {
    return exports.bounds(input)
  }
}

const zoomLevels = {
  city: 100,
  country: 200,
  world: 500
}
exports.parseZoom = function(input) {
  if (typeof input === 'number') {
    return input
  }
  if (zoomLevels.hasOwnProperty(input)) {
    return zoomLevels[input]
  }
  return 500
}
