<script>
  import * as d3Geo from 'd3-geo'
  // import countries from './data/more.geo.json'
  import landJson from './data/land-50m.js'
  const toronto = [-79.43, 43.65]
  const edmonton = [-113.5129, 53.5323]
  const dr = [-68.9906, 18.4096]
  // import north from './data/points/north-america.js'
  import { geoPath, geoOrthographic } from 'd3-geo'

  import { feature } from 'topojson-client'
  let width = 500
  let height = 500
  const projection = geoOrthographic().rotate([90, -20]).scale(228).translate([210, 200])
  const path = geoPath().projection(projection)
  path.prototype.toString = function () {
    console.log(this._.length)
    let str = this._.replace(/\d+\.\d+/g, (s) => parseFloat(s).toFixed(1))
    console.log(this._.length)
    return str
  }
  const graticule = d3Geo.geoGraticule10()
  const land = feature(landJson, landJson.objects.land)
  // data = path(land)

  let circle = d3Geo
    .geoCircle()
    .center(edmonton)
    .radius((120 / (6371 * Math.PI * 2)) * 360)
    .precision(1)
  let dot1 = path(circle())

  let circle2 = d3Geo
    .geoCircle()
    .center(dr)
    .radius((120 / (6371 * Math.PI * 2)) * 360)
    .precision(1)
  let dot2 = path(circle2())
</script>

<svg viewBox="0 0 {width} {height}" style="display: block;">
  <g clip-path="url(${new URL('#clip', location)})">
    <path d={path(graticule)} stroke="#9c9ea1" fill="none" opacity="0.3" strokeWidth="0.1px" />
    <path d={path(land)} fill="#7c8691" opacity="0.8" />
  </g>
  <use xlink:href="${new URL('#outline', location)}" fill="none" stroke="#000" />
  <path d={dot1} fill="#cc7066" stroke="#D68881" fill-opacity="0.4" />
  <path d={dot2} fill="#cc7066" stroke="#D68881" fill-opacity="0.4" />
</svg>
