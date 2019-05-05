const wtf = require('/Users/spencer/mountain/wtf_wikipedia/')

wtf.fetch('List of cities in Ontario', (err, doc) => {
  let list = []
  let tables = doc.tables()[0]
  tables.json().map(row => {
    let links = row.Name.links || [{}]
    list.push(links[0].page)
  })
  list = list.filter(l => l)
  console.log(JSON.stringify(list, null, 2))
})
