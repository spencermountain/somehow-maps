const THREE = require('../assets/three')

exports.createGeometryArray = function(json) {
  var geometry_array = []

  if (json.type === 'Feature') {
    geometry_array.push(json.geometry)
  } else if (json.type === 'FeatureCollection') {
    for (var feature_num = 0; feature_num < json.features.length; feature_num++) {
      geometry_array.push(json.features[feature_num].geometry)
    }
  } else if (json.type === 'GeometryCollection') {
    for (var geom_num = 0; geom_num < json.geometries.length; geom_num++) {
      geometry_array.push(json.geometries[geom_num])
    }
  } else {
    throw new Error('The geoJSON is not valid.')
  }
  //alert(geometry_array.length);
  return geometry_array
}

exports.needsInterpolation = function(point2, point1) {
  //If the distance between two latitude and longitude values is
  //greater than five degrees, return true.
  var lon1 = point1[0]
  var lat1 = point1[1]
  var lon2 = point2[0]
  var lat2 = point2[1]
  var lon_distance = Math.abs(lon1 - lon2)
  var lat_distance = Math.abs(lat1 - lat2)

  if (lon_distance > 5 || lat_distance > 5) {
    return true
  } else {
    return false
  }
}

exports.getMidpoint = function(point1, point2) {
  var midpoint_lon = (point1[0] + point2[0]) / 2
  var midpoint_lat = (point1[1] + point2[1]) / 2
  var midpoint = [midpoint_lon, midpoint_lat]

  return midpoint
}

exports.createVertexForEachPoint = function(
  object_geometry,
  values_axis1,
  values_axis2,
  values_axis3
) {
  for (var i = 0; i < values_axis1.length; i++) {
    object_geometry.vertices.push(
      new THREE.Vector3(values_axis1[i], values_axis2[i], values_axis3[i])
    )
  }
}

exports.drawParticle = function(x, y, z, options, container) {
  var particle_geom = new THREE.Geometry()
  particle_geom.vertices.push(new THREE.Vector3(x, y, z))

  var particle_material = new THREE.ParticleSystemMaterial(options)

  var particle = new THREE.ParticleSystem(particle_geom, particle_material)
  container.add(particle)

  x_values.length = 0
  y_values.length = 0
  z_values.length = 0
}

exports.drawLine = function(x_values, y_values, z_values, options, container) {
  var line_geom = new THREE.Geometry()
  exports.createVertexForEachPoint(line_geom, x_values, y_values, z_values)

  var line_material = new THREE.LineBasicMaterial(options)
  var line = new THREE.Line(line_geom, line_material)
  container.add(line)

  x_values.length = 0
  y_values.length = 0
  z_values.length = 0
}

function convertToSphereCoords(coordinates_array, sphere_radius, x_values, y_values, z_values) {
  var lon = coordinates_array[0]
  var lat = coordinates_array[1]

  x_values.push(Math.cos((lat * Math.PI) / 180) * Math.cos((lon * Math.PI) / 180) * sphere_radius)
  y_values.push(Math.cos((lat * Math.PI) / 180) * Math.sin((lon * Math.PI) / 180) * sphere_radius)
  z_values.push(Math.sin((lat * Math.PI) / 180) * sphere_radius)
}

function convertToPlaneCoords(coordinates_array, radius, x_values, y_values, z_values) {
  var lon = coordinates_array[0]
  var lat = coordinates_array[1]

  z_values.push((lat / 180) * radius)
  y_values.push((lon / 180) * radius)
}

exports.getConversionFunctionName = function(shape) {
  var conversionFunctionName

  if (shape === 'sphere') {
    conversionFunctionName = convertToSphereCoords
  } else if (shape === 'plane') {
    conversionFunctionName = convertToPlaneCoords
  } else {
    throw new Error('The shape that you specified is not valid.')
  }
  return conversionFunctionName
}
