<script>
  import { tweened } from 'svelte/motion'
  import { cubicOut } from 'svelte/easing'
  import { Globe, Line, Graticule, Dot, Latitude, Longitude } from '../../src'
  const r = tweened(0, {
    duration: 1400,
    easing: cubicOut
  })
  r.set(180)
</script>

<style>

</style>

<div>
  <Globe rotate={$r} tilt={-10}>
    <Graticule />
    <Latitude at={40} />
    <!-- <Longitude at={20} /> -->
    <!-- <Line from="toronto" to={'jamaica'} /> -->
    <Dot at={[90, 0]} color="lightblue" radius={50} />
    <Dot at={[-90, 0]} color="lightblue" radius={50} />
  </Globe>
</div>
