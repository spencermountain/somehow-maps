<script>
  import { getContext } from 'svelte'
  // import c from 'spencer-color'
  export let shape = ''
  export let stroke = 'lightgrey'
  export let fill = 'white'
  let c = {
    colors: {}
  }
  fill = c.colors[fill] || fill
  stroke = c.colors[stroke] || stroke

  let projection = getContext('projection')
  const toPath = d3Geo.geoPath().projection(projection)

  // console.log(shape)
  // let key = Object.keys(shape.objects)[0]
  // let geoJSON = topojson.feature(shape, shape.objects[key])
  // console.log(geoJSON)
  let d = toPath(shape)
  // console.log(d)
</script>

<path {d} {fill} {stroke} shape-rendering="geometricPrecision" />
