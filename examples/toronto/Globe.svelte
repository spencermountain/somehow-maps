<script>
  import { tweened } from 'svelte/motion'
  import { cubicOut } from 'svelte/easing'
  import { Globe, Countries, Line, Dot, Shape } from '../../src'
  import roads from './roads'
</script>

<style>

</style>

<div>
  <Globe tilt={0} flat="true">
    <Countries />
    <!-- <Latitude at={40} /> -->
    <!-- <Shape shape={roads} color="blue" /> -->
  </Globe>
</div>
