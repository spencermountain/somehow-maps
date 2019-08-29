/* Draw GeoJSON

Iterates through the latitude and longitude values, converts the values to XYZ coordinates,
and draws the geoJSON geometries.

*/
const fns = require('./fns')

function interpolatePoints(interpolation_array) {
  //This function is recursive. It will continue to add midpoints to the
  //interpolation array until needsInterpolation() returns false.
  var temp_array = []
  var point1, point2

  for (let point_num = 0; point_num < interpolation_array.length - 1; point_num++) {
    point1 = interpolation_array[point_num]
    point2 = interpolation_array[point_num + 1]

    if (fns.needsInterpolation(point2, point1)) {
      temp_array.push(point1)
      temp_array.push(fns.getMidpoint(point1, point2))
    } else {
      temp_array.push(point1)
    }
  }
  return temp_array
}

function drawThreeGeo(json, radius, shape, materalOptions, container) {
  container = container || window.scene

  var x_values = []
  var y_values = []
  var z_values = []

  function createCoordinateArray(feature) {
    //Loop through the coordinates and figure out if the points need interpolation.
    var temp_array = []
    var interpolation_array = []

    for (let point_num = 0; point_num < feature.length; point_num++) {
      var point1 = feature[point_num]
      var point2 = feature[point_num - 1]

      if (point_num > 0) {
        if (fns.needsInterpolation(point2, point1)) {
          interpolation_array = [point2, point1]
          temp_array = interpolatePoints(interpolation_array)

          for (
            var inter_point_num = 0;
            inter_point_num < interpolation_array.length;
            inter_point_num++
          ) {
            temp_array.push(interpolation_array[inter_point_num])
          }
        } else {
          temp_array.push(point1)
        }
      } else {
        temp_array.push(point1)
      }
    }
    return temp_array
  }

  var json_geom = fns.createGeometryArray(json)
  //An array to hold the feature geometries.
  var convertCoordinates = fns.getConversionFunctionName(shape)

  //Whether you want to convert to spherical or planar coordinates.
  var coordinate_array = []
  //Re-usable array to hold coordinate values. This is necessary so that you can add
  //interpolated coordinates. Otherwise, lines go through the sphere instead of wrapping around.

  for (var geom_num = 0; geom_num < json_geom.length; geom_num++) {
    if (json_geom[geom_num].type === 'Point') {
      convertCoordinates(json_geom[geom_num].coordinates, radius, x_values, y_values, z_values)
      fns.drawParticle(x_values[0], y_values[0], z_values[0], materalOptions, container)
    } else if (json_geom[geom_num].type == 'MultiPoint') {
      for (let point_num = 0; point_num < json_geom[geom_num].coordinates.length; point_num++) {
        convertCoordinates(
          json_geom[geom_num].coordinates[point_num],
          radius,
          x_values,
          y_values,
          z_values
        )
        fns.drawParticle(x_values[0], y_values[0], z_values[0], materalOptions, container)
      }
    } else if (json_geom[geom_num].type === 'LineString') {
      coordinate_array = createCoordinateArray(json_geom[geom_num].coordinates)

      for (let point_num = 0; point_num < coordinate_array.length; point_num++) {
        convertCoordinates(coordinate_array[point_num], radius, x_values, y_values, z_values)
      }
      fns.drawLine(x_values, y_values, z_values, materalOptions, container)
    } else if (json_geom[geom_num].type == 'Polygon') {
      for (
        let segment_num = 0;
        segment_num < json_geom[geom_num].coordinates.length;
        segment_num++
      ) {
        coordinate_array = createCoordinateArray(json_geom[geom_num].coordinates[segment_num])
        for (let point_num = 0; point_num < coordinate_array.length; point_num++) {
          convertCoordinates(coordinate_array[point_num], radius, x_values, y_values, z_values)
        }
        fns.drawLine(x_values, y_values, z_values, materalOptions, container)
      }
    } else if (json_geom[geom_num].type === 'MultiLineString') {
      for (
        let segment_num = 0;
        segment_num < json_geom[geom_num].coordinates.length;
        segment_num++
      ) {
        coordinate_array = createCoordinateArray(json_geom[geom_num].coordinates[segment_num])

        for (let point_num = 0; point_num < coordinate_array.length; point_num++) {
          convertCoordinates(coordinate_array[point_num], radius, x_values, y_values, z_values)
        }
        fns.drawLine(x_values, y_values, z_values, materalOptions, container)
      }
    } else if (json_geom[geom_num].type === 'MultiPolygon') {
      for (
        var polygon_num = 0;
        polygon_num < json_geom[geom_num].coordinates.length;
        polygon_num++
      ) {
        for (
          let segment_num = 0;
          segment_num < json_geom[geom_num].coordinates[polygon_num].length;
          segment_num++
        ) {
          coordinate_array = createCoordinateArray(
            json_geom[geom_num].coordinates[polygon_num][segment_num]
          )

          for (let point_num = 0; point_num < coordinate_array.length; point_num++) {
            convertCoordinates(coordinate_array[point_num], radius, x_values, y_values, z_values)
          }
          fns.drawLine(x_values, y_values, z_values, materalOptions, container)
        }
      }
    } else {
      throw new Error('The geoJSON is not valid.')
    }
  }

  // temp_array.push(interpolation_array[interpolation_array.length - 1])

  // if (temp_array.length > interpolation_array.length) {
  //   temp_array = interpolatePoints(temp_array)
  // } else {
  //   return temp_array
  // }
  return []
}

module.exports = drawThreeGeo
