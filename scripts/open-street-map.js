var fs = require('fs')
var through = require('through2')
var parseOSM = require('osm-pbf-parser')
const path = '/Users/spencer/data/openstreetmap/ontario-latest.osm.pbf'

const tag = 'amenity:prison'
let [key, val] = tag.split(':')
console.log(key + '  :  ' + val)

let nodes = []
let ways = []
let want = {}

const doit = function(each, cb) {
  var osm = parseOSM()
  fs.createReadStream(path)
    .pipe(osm)
    .pipe(through.obj(each, cb))
}

//first-pass, get all nodes, ways
const firstPass = function(items, _, next) {
  items.forEach(item => {
    if (item.tags[key] === val) {
      delete item.info
      if (item.type === 'node') {
        nodes.push(item)
      } else if (item.type === 'way') {
        ways.push(item)
      } else {
        console.log('skipping ' + item.type)
      }
    }
  })
  next()
}

//get all the way's nodes
const secondPass = function(items, _, next) {
  items.forEach(item => {
    delete item.info
    if (want[item.id] === true) {
      want[item.id] = item
    }
  })
  next()
}

//cleanup the data
const postProcess = function() {
  ways = ways.map(way => {
    return {
      name: way.tags.name,
      path: way.refs.map(id => {
        return [want[id].lat, want[id].lon]
      })
    }
  })
  nodes = nodes.map(node => {
    return {
      name: node.tags.name,
      point: [node.lat, node.lon]
    }
  })
  nodes = nodes.filter(node => node.name)
}

doit(firstPass, () => {
  console.log(' - done first-pass - ')
  ways.forEach(way => {
    way.refs.forEach(id => (want[id] = true))
  })
  console.log('want ' + Object.keys(want).length + ' refs')

  doit(secondPass, () => {
    console.log(' - done second-pass - ')
    postProcess()
    let results = {
      ways: ways,
      nodes: nodes
    }
    fs.writeFileSync('./tmp.json', JSON.stringify(results, null, 2))
  })
})
