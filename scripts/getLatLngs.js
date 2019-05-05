const wtf = require('/Users/spencer/mountain/wtf_wikipedia/')

let arr = [
  'Barrie',
  'Belleville, Ontario',
  'Brampton',
  'County of Brant',
  'Brantford',
  'Brockville',
  'Burlington, Ontario',
  'Cambridge, Ontario',
  'Clarence-Rockland',
  'Cornwall, Ontario',
  'Dryden, Ontario',
  'Elliot Lake',
  'Greater Sudbury',
  'Guelph',
  'Haldimand County',
  'Hamilton, Ontario',
  'Kawartha Lakes',
  'Kenora',
  'Kingston, Ontario',
  'Kitchener, Ontario',
  'London, Ontario',
  'Markham, Ontario',
  'Mississauga',
  'Niagara Falls, Ontario',
  'Norfolk County, Ontario',
  'North Bay, Ontario',
  'Orillia',
  'Oshawa',
  'Ottawa',
  'Owen Sound',
  'Pembroke, Ontario',
  'Peterborough, Ontario',
  'Pickering, Ontario',
  'Port Colborne',
  'Prince Edward County, Ontario',
  'Quinte West',
  'Richmond Hill, Ontario',
  'Sarnia',
  'Sault Ste. Marie, Ontario',
  'St. Catharines',
  'St. Thomas, Ontario',
  'Stratford, Ontario',
  'Temiskaming Shores',
  'Thorold',
  'Thunder Bay',
  'Timmins',
  'Toronto',
  'Vaughan',
  'Waterloo, Ontario',
  'Welland',
  'Windsor, Ontario',
  'Woodstock, Ontario'
]

wtf.fetch(arr, (err, docs) => {
  let places = docs.reduce((h, doc) => {
    let coord = (doc.templates('coord') || [])[0] || {}
    console.log(coord)
    h[doc.title()] = [coord.lat, coord.lon]
    return h
  }, {})
  console.log(JSON.stringify(places, null, 2))
})
