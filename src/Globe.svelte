<script>
  import { setContext } from 'svelte'
  import Shape from './shapes/Shape.svelte'
  import * as d3Geo from 'd3-geo'
  import findPoint from './lib/findPoint'
  import focusOn from './lib/focusOn'
  import countries from './data/more.geo.json'
  // export let rotate = 0
  // export let tilt = 0
  export let focus = countries
  export let width = 500
  export let height = 500
  export let show = ''
  show = findPoint(show) || show || ''
  let projection = d3Geo.geoOrthographic()
  projection.scale(1).translate([0, 0])

  let path = d3Geo.geoPath().projection(projection)
  let shape = countries.features[2]

  focusOn(focus, projection, width, height)
  // console.log(projection.rotate([0.1, 0.1, 0]))
  // if (rotate || tilt) {
  // projection.rotate([rotate, tilt, 3])
  // projection.rotate([0, 0, 0])
  // }

  setContext('projection', projection)
</script>

<style>
  svg {
    margin: 10px 20px 25px 25px;
    /* border: 1px solid lightgrey; */
  }
</style>

<svg
  viewBox="0,0,{width},{height}"
  width="100%"
  height="100%"
  preserveAspectRatio="xMidYMid meet"
  style="margin: 10px 20px 25px 25px;">
  <slot />
</svg>
