<div align="center">
  <img src="https://cloud.githubusercontent.com/assets/399657/23590290/ede73772-01aa-11e7-8915-181ef21027bc.png" />
  <div>make a map, without thinking</div>

  <a href="https://npmjs.org/package/somehow-maps">
    <img src="https://img.shields.io/npm/v/somehow-maps.svg?style=flat-square" />
  </a>
  <a href="https://unpkg.com/somehow-maps">
    <img src="https://badge-size.herokuapp.com/spencermountain/somehow/master/builds/somehow-maps.min.js" />
  </a>
</div>

<div align="center">
  <code>npm install somehow-maps</code>
</div>

<div align="center">
**work-in-progress**
</div>

```js
const somehowMaps = require('somehow-maps')

let w = somehowMaps({
  height: 300,
  aspect: 'widescreen'
})

w.shape({
  shape: 'great-lakes'
})

w.line()
  .from('toronto')
  .to('winnipeg')
  .color('red')

w.dot()
  .at('barrie')
  .color('blue')

w.text('Toronto')
  .at('toronto')
  .color('red')

document.querySelector('#stage').innerHTML = w.build()
```

![image](https://user-images.githubusercontent.com/399657/57199810-8cfdfb00-6f51-11e9-8ccd-39e5154e4ebd.png)

---

to generate ad-hoc topojson shapes (without using QGIS):

- grab shapefiles from [natural earth](https://www.naturalearthdata.com/downloads/)
- draw a new shape in [geojson.io](http://geojson.io)
- add them both to [mapshaper.org](https://mapshaper.org/)
- run `mapshaper -clip myshape.geojson`
- export as **topojson**

---

to generate ad-hoc topojson shapes from openstreetmap:

- download some pbf data from [geofabrik.de](http://download.geofabrik.de/)

- use James Halliday's nodejs [pbf parser](https://github.com/substack/osm-pbf-parser) to filter down [specific tags](https://wiki.openstreetmap.org/wiki/Key:amenity)

### See also

- [somehow](https://github.com/spencermountain/somehow)
- [thensome.how](https://thensome.how)

MIT
