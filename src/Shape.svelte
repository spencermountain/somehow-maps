<script>
  import { getContext } from 'svelte'
  import * as d3Geo from 'd3-geo'
  import * as topojson from 'topojson-client'
  import c from 'spencer-color'
  export let data = ''
  export let stroke = 'lightgrey'
  export let fill = 'white'
  fill = c.colors[fill] || fill
  stroke = c.colors[stroke] || stroke

  let projection = getContext('projection')
  const toPath = d3Geo.geoPath().projection(projection)

  let key = Object.keys(data.objects)[0]
  let geoJSON = topojson.feature(data, data.objects[key])
  let d = toPath(geoJSON)
</script>

<path {d} {fill} {stroke} shape-rendering="geometricPrecision" />
