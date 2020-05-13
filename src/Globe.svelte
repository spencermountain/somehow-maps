<script>
  import { setContext } from 'svelte'
  import Shape from './Shape.svelte'
  import * as d3Geo from 'd3-geo'
  import countries from './shapes/countries'
  export let rotate = 0
  export let tilt = 0
  export let width = 500
  export let height = 500
  export let showCountries = true
  export let flat = false

  // import { rotate, tilt } from './stores.js'
  let projection = d3Geo.geoOrthographic().scale(180)
  projection.rotate([rotate, tilt, 3])

  if (flat) {
    projection = d3Geo.geoMercator().scale(75)
  }
  projection.translate([200, 200])
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
  {#if showCountries}
    <Shape data={countries} />
  {/if}
  <slot />
</svg>
