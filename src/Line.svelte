<script>
  import findPoint from './findPoint'
  import { getContext } from 'svelte'
  import * as d3Geo from 'd3-geo'
  import c from 'spencer-color'
  export let from = ''
  export let to = ''
  export let color = 'blue'
  color = c.colors[color] || color

  let projection = getContext('projection')
  from = findPoint(from).reverse()
  to = findPoint(to).reverse()
  console.log(from, to)

  const toPath = d3Geo.geoPath().projection(projection)
  let geoJSON = {
    type: 'Feature',
    geometry: {
      type: 'LineString',
      coordinates: [from, to]
    }
  }
  let d = toPath(geoJSON)
</script>

<path {d} stroke-width="4" fill="none" stroke={color} stroke-linecap="round" />
