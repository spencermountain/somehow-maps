<script>
  import findPoint from '../lib/findPoint'
  import * as d3Geo from 'd3-geo'
  import { getContext } from 'svelte'
  let c = {
    colors: {}
  }
  export let at = ''
  export let radius = 2
  export let opacity = 0.5
  export let color = 'blue'
  export let label = ''
  color = c.colors[color] || color

  let projection = getContext('projection')

  var path = d3Geo.geoPath().projection(projection)

  at = findPoint(at)
  at = at.reverse()

  let circle = d3Geo.geoCircle().radius(radius).center(at)
  let d = path(circle(at))
</script>

<path {d} fill={color} stroke={color} fill-opacity={opacity}>
  <line>{label}</line>
</path>
