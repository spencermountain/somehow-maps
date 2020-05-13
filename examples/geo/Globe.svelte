<script>
  import { Globe, Line, Graticule, Dot, Latitude, Longitude } from '../../src'
</script>

<style>

</style>

<div>
  <Globe rotate={90} tilt={-10}>
    <Graticule />
    <Latitude at={40} />
    <!-- <Longitude at={20} /> -->
    <!-- <Line from="toronto" to={'jamaica'} /> -->
    <Dot at={[90, 0]} color="lightblue" radius={50} />
    <Dot at={[-90, 0]} color="lightblue" radius={50} />
  </Globe>
</div>
