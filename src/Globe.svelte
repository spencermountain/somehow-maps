<script>
  import { setContext } from 'svelte'
  import Shape from './shapes/Shape.svelte'
  import * as d3Geo from 'd3-geo'
  import * as d3Zoom from 'd3-zoom'

  import findPoint from './lib/findPoint'
  export let rotate = 0
  export let tilt = 0
  export let width = 500
  // export let zoom = 75
  export let height = 500
  export let flat = false
  export let show = ''
  show = findPoint(show) || show || ''
  // let projection = d3Geo.geoOrthographic().scale(180)

  // if (flat) {
  let projection = d3Geo.geoMercator() //.scale(zoom)
  projection.scale(1).translate([0, 0])
  // }
  // if (show) {
  //   projection.center(show)
  // }
  // if (rotate || tilt) {
  //   projection.rotate([rotate, tilt, 3])
  // }

  // projection.translate([width / 2, height / 2])
  // projection.translate([200, 200])

  import countries from './data/more.geo.json'
  let path = d3Geo.geoPath().projection(projection)
  let shape = countries.features[2]

  var b = path.bounds(shape),
    s = 0.95 / Math.max((b[1][0] - b[0][0]) / width, (b[1][1] - b[0][1]) / height),
    t = [(width - s * (b[1][0] + b[0][0])) / 2, (height - s * (b[1][1] + b[0][1])) / 2]
  console.log(b)
  projection.scale(s).translate(t)

  setContext('projection', projection)
</script>

<style>
  svg {
    margin: 10px 20px 25px 25px;
    border: 1px solid lightgrey;
  }
</style>

<div>scale: {JSON.stringify(s)}</div>
<div>translate: {JSON.stringify(t)}</div>
<svg
  viewBox="0,0,{width},{height}"
  preserveAspectRatio="xMidYMid meet"
  style="margin: 10px 20px 25px 25px;">
  <slot />
</svg>
