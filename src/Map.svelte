<script>
  import { setContext } from 'svelte'
  import Shape from './shapes/Shape.svelte'
  import * as d3Geo from 'd3-geo'
  import * as d3Zoom from 'd3-zoom'

  import findPoint from './lib/findPoint'
  import focusOn from './lib/focusOn'
  export let width = 500
  export let height = 300
  export let focus = []

  let projection = d3Geo.geoMercator()
  projection.scale(1).translate([0, 0])

  focusOn(focus, projection, width, height)

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
