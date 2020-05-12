<script>
  import { Globe, Line, Graticule, Dot } from '../../src'
</script>

<style>

</style>

<div>
  <Globe rotate={0} tilt={40}>
    <Graticule />
    <Line from="toronto" to={'jamaica'} />
    <Dot at="vancouver" color="red" />
  </Globe>
</div>
