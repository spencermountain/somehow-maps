<script>
  import { getContext } from 'svelte'
  import * as d3Geo from 'd3-geo'
  let projection = getContext('projection')

  export let lat = 4
  export let lon = 4
  export let color = '#d7d5d2'
  const toPath = d3Geo.geoPath().projection(projection)
  let grat = d3Geo.geoGraticule().step([lon, lat])
  let lines = grat.lines()
</script>

<g class="graticule">
  {#each lines as line}
    <path d={toPath(line)} fill="none" stroke={color} class="line" />
  {/each}
</g>
