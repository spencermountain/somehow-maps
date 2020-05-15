<script>
  import { setContext } from 'svelte'
  import Shape from './shapes/Shape.svelte'
  import * as d3Geo from 'd3-geo'
  import * as d3Zoom from 'd3-zoom'

  import findPoint from './lib/findPoint'
  import focusOn from './lib/focusOn'
  export let width = 500
  export let height = 500
  export let focus = []

  let projection = d3Geo.geoMercator()
  projection.scale(1).translate([0, 0])

  focusOn(focus, projection, width, height)

  // let path = d3Geo.geoPath().projection(projection)
  // console.log(focus)

  // var b = path.bounds(focus),
  //   s = 0.95 / Math.max((b[1][0] - b[0][0]) / width, (b[1][1] - b[0][1]) / height),
  //   t = [(width - s * (b[1][0] + b[0][0])) / 2, (height - s * (b[1][1] + b[0][1])) / 2]
  // console.log(b)
  // projection.scale(s).translate(t)

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
