<script>
  import { setContext } from 'svelte'
  import Shape from './shapes/Shape.svelte'
  import * as d3Geo from 'd3-geo'

  import findPoint from './lib/findPoint'
  export let rotate = 0
  export let tilt = 0
  export let width = 500
  export let zoom = 75
  export let height = 500
  export let flat = false
  export let show = ''
  show = findPoint(show) || show || [0, 0]
  let projection = d3Geo.geoOrthographic().scale(180)

  if (flat) {
    projection = d3Geo.geoMercator().scale(zoom)
    projection.center(show)
    console.log(show)
  }
  if (rotate || tilt) {
    projection.rotate([rotate, tilt, 3])
  }

  projection.translate([width / 2, height / 2])
  // projection.translate([200, 200])
  setContext('projection', projection)
</script>

<style>
  svg {
    margin: 10px 20px 25px 25px;
    border: 1px solid lightgrey;
  }
</style>

<svg
  viewBox="0,0,{width},{height}"
  preserveAspectRatio="xMidYMid meet"
  style="margin: 10px 20px 25px 25px;">
  <slot />
</svg>
