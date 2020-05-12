<script>
  import findPoint from './findPoint'
  import * as d3Geo from 'd3-geo'
  import { getContext } from 'svelte'
  import c from 'spencer-color'
  export let at = ''
  export let radius = 4
  export let opacity = 0.5
  export let color = 'blue'
  export let label = ''
  color = c.colors[color] || color

  let projection = getContext('projection')

  var path = d3Geo.geoPath().projection(projection)

  at = findPoint(at)

  let circle = d3Geo.geoCircle().radius(15)

  let d = path(circle(at))
</script>

<path {d} fill={color} stroke="red" {opacity}>
  <line>{label}</line>
</path>
