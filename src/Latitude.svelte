<script>
  import findPoint from './findPoint'
  import * as d3Geo from 'd3-geo'
  import { getContext } from 'svelte'
  import c from 'spencer-color'
  export let at = 0
  export let color = 'blue'
  export let label = ''
  color = c.colors[color] || color

  let projection = getContext('projection')

  const toPath = d3Geo.geoPath().projection(projection)

  let arr = []
  for (let lat = -180; lat <= 180; lat += 10) {
    arr.push([lat, at]) //create a ton of small line-segments
  }
  let geoJSON = {
    type: 'Feature',
    geometry: {
      type: 'LineString',
      coordinates: arr
    }
  }

  let d = toPath(geoJSON)
</script>

<path {d} fill="none" stroke={color}>
  <line>{label}</line>
</path>
