<script>
  import { tweened } from 'svelte/motion'
  import { cubicOut } from 'svelte/easing'
  import { Globe, Countries, Line, Dot, Shape } from '../../src'
  import roads from './roads'
  let road = [[-79.455808, 43.6750582], [-79.4548503, 43.6752713], [-79.4547401, 43.6752928]]
  var geoJSON = {
    type: 'FeatureCollection',
    features: [
      {
        type: 'Feature',
        geometry: {
          type: 'Polygon',
          coordinates: road
        },
        properties: {}
      }
    ]
  }
</script>

<style>

</style>

<div>
  <Globe tilt={0} flat="true" show={[-79.4548503, 43.6752713]} zoom={182500}>
    <!-- <Countries /> -->
    <!-- <Latitude at={40} /> -->
    <Shape shape={roads} color="blue" />
  </Globe>
</div>
