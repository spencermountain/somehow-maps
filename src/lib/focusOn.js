import * as d3Geo from 'd3-geo'

const focusOn = function (shape, projection, width, height) {
  let path = d3Geo.geoPath().projection(projection)

  var b = path.bounds(shape)
  let s = 0.95 / Math.max((b[1][0] - b[0][0]) / width, (b[1][1] - b[0][1]) / height)
  let t = [(width - s * (b[1][0] + b[0][0])) / 2, (height - s * (b[1][1] + b[0][1])) / 2]
  projection.scale(s).translate(t)
  return projection
}
export default focusOn
