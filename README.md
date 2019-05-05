# somehow-geo

make a map without thinking

---

to generate ad-hoc topojson shapes (without using QGIS):

- grab shapefiles from [natural earth](https://www.naturalearthdata.com/downloads/)
- draw a new shape in [geojson.io](http://geojson.io)
- add them both to [mapshaper.org](https://mapshaper.org/)
- run `mapshaper -clip myshape.geojson`
- export as **topojson**
