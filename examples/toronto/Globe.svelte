<script>
  import { tweened } from 'svelte/motion'
  import { cubicOut } from 'svelte/easing'
  import { Map, Globe, Countries, Line, Dot, Shape } from '../../src'
  import roads from './roads'

  // import countries from '../../src/data/more.geo.json'
  // let focus = countries.features[3]
</script>

<style>

</style>

<div>
  <Map focus={roads}>
    <Countries />
    <Shape shape={roads} color="blue" />
  </Map>
</div>
