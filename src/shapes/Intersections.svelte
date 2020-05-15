<script>
  import { ShapeInfo, Intersection } from 'kld-intersections'
  import { getContext } from 'svelte'
  import * as d3Geo from 'd3-geo'
  export let shape = []

  let projection = getContext('projection')
  const toPath = d3Geo.geoPath().projection(projection)

  let allPaths = shape.features.map(obj => {
    let tmp = {
      type: 'FeatureCollection',
      features: [obj]
    }
    let d = toPath(tmp)
    return ShapeInfo.path(d)
  })
  let allIntersections = []
  allPaths.forEach((one, i) => {
    allPaths.forEach((two, o) => {
      if (i == o) {
        return
      }
      let found = Intersection.intersect(one, two).points
      if (found.length) {
        allIntersections = allIntersections.concat(found)
      }
    })
  })
  console.log(allIntersections)

  let dufferin = {
    type: 'FeatureCollection',
    features: shape.features.filter(d => d.properties.name && d.properties.name.match(/dufferin/i))
  }
  let allD = toPath(dufferin)

  let start = [15, 75]
  let end = [355, 140]
  const path = ShapeInfo.path(allD)
  const line = ShapeInfo.line(start, end)
  const intersections = Intersection.intersect(path, line)

  let dots = allIntersections
</script>

<style>

</style>

<g>
  {#each dots as dot}
    <circle r="2" fill-opacity="0.5" fill="red" cx={dot.x} cy={dot.y} />
  {/each}
</g>
