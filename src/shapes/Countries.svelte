<script>
  import Shape from './Shape.svelte'
  import * as topojson from 'topojson-client'

  import countries from '../data/countries'
  let key = Object.keys(countries.objects)[0]
  let geoJSON = topojson.feature(countries, countries.objects[key])
  console.log(geoJSON)
</script>

<Shape shape={geoJSON} />
