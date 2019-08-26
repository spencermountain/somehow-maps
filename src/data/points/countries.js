const points = [
  ['afghanistan', 'kabul', 34.28, 69.11],
  ['albania', 'tirane', 41.18, 19.49],
  ['algeria', 'algiers', 36.42, 3.08],
  ['american samoa', 'pago pago', -14.16, -170.43],
  ['andorra', 'andorra la vella', 42.31, 1.32],
  ['angola', 'luanda', -8.5, 13.15],
  ['antigua and barbuda', 'west indies', 17.2, -61.48],
  ['argentina', 'buenos aires', -36.3, -60.0],
  ['armenia', 'yerevan', 40.1, 44.31],
  ['aruba', 'oranjestad', 12.32, -70.02],
  ['australia', 'canberra', -35.15, 149.08],
  ['austria', 'vienna', 48.12, 16.22],
  ['azerbaijan', 'baku', 40.29, 49.56],
  ['bahamas', 'nassau', 25.05, -77.2],
  ['bahrain', 'manama', 26.1, 50.3],
  ['bangladesh', 'dhaka', 23.43, 90.26],
  ['barbados', 'bridgetown', 13.05, -59.3],
  ['belarus', 'minsk', 53.52, 27.3],
  ['belgium', 'brussels', 50.51, 4.21],
  ['belize', 'belmopan', 17.18, -88.3],
  ['benin', 'porto novo', 6.23, 2.42],
  ['bhutan', 'thimphu', 27.31, 89.45],
  ['bolivia', 'la paz', -16.2, -68.1],
  ['bosnia and herzegovina', 'sarajevo', 43.52, 18.26],
  ['botswana', 'gaborone', -24.45, 25.57],
  ['brazil', 'brasilia', -15.47, -47.55],
  ['british virgin islands', 'road town', 18.27, -64.37],
  ['brunei darussalam', 'bandar seri begawan', 4.52, 115.0],
  ['bulgaria', 'sofia', 42.45, 23.2],
  ['burkina faso', 'ouagadougou', 12.15, -1.3],
  ['burundi', 'bujumbura', -3.16, 29.18],
  ['cambodia', 'phnom penh', 11.33, 104.55],
  ['cameroon', 'yaounde', 3.5, 11.35],
  ['canada', 'ottawa', 45.27, -75.42],
  ['cape verde', 'praia', 15.02, -23.34],
  ['cayman islands', 'george town', 19.2, -81.24],
  ['central african republic', 'bangui', 4.23, 18.35],
  ['chad', "n'djamena", 12.1, 14.59],
  ['chile', 'santiago', -33.24, -70.4],
  ['china', 'beijing', 39.55, 116.2],
  ['colombia', 'bogota', 4.34, -74.0],
  ['comros', 'moroni', -11.4, 43.16],
  ['congo', 'brazzaville', -4.09, 15.12],
  ['costa rica', 'san jose', 9.55, -84.02],
  ["cote d'ivoire", 'yamoussoukro', 6.49, -5.17],
  ['croatia', 'zagreb', 45.5, 15.58],
  ['cuba', 'havana', 23.08, -82.22],
  ['cyprus', 'nicosia', 35.1, 33.25],
  ['czech republic', 'prague', 50.05, 14.22],
  ['democratic republic of the congo', 'kinshasa', -4.2, 15.15],
  ['denmark', 'copenhagen', 55.41, 12.34],
  ['djibouti', 'djibouti', 11.08, 42.2],
  ['dominica', 'roseau', 15.2, -61.24],
  ['dominica republic', 'santo domingo', 18.3, -69.59],
  ['east timor', 'dili', -8.29, 125.34],
  ['ecuador', 'quito', -0.15, -78.35],
  ['egypt', 'cairo', 30.01, 31.14],
  ['el salvador', 'san salvador', 13.4, -89.1],
  ['equatorial guinea', 'malabo', 3.45, 8.5],
  ['eritrea', 'asmara', 15.19, 38.55],
  ['estonia', 'tallinn', 59.22, 24.48],
  ['ethiopia', 'addis ababa', 9.02, 38.42],
  ['falkland islands', 'stanley', -51.4, -59.51],
  ['faroe islands', 'torshavn', 62.05, -6.56],
  ['fiji', 'suva', -18.06, 178.3],
  ['finland', 'helsinki', 60.15, 25.03],
  ['france', 'paris', 48.5, 2.2],
  ['french guiana', 'cayenne', 5.05, -52.18],
  ['french polynesia', 'papeete', -17.32, -149.34],
  ['gabon', 'libreville', 0.25, 9.26],
  ['gambia', 'banjul', 13.28, -16.4],
  ['georgia', 'tbilisi', 41.43, 44.5],
  ['germany', 'berlin', 52.3, 13.25],
  ['ghana', 'accra', 5.35, -0.06],
  ['greece', 'athens', 37.58, 23.46],
  ['greenland', 'nuuk', 64.1, -51.35],
  ['guadeloupe', 'basse-terre', 16.0, -61.44],
  ['guatemala', 'guatemala', 14.4, -90.22],
  ['guernsey', 'st. peter port', 49.26, -2.33],
  ['guinea', 'conakry', 9.29, -13.49],
  ['guinea-bissau', 'bissau', 11.45, -15.45],
  ['guyana', 'georgetown', 6.5, -58.12],
  ['haiti', 'port-au-prince', 18.4, -72.2],
  ['honduras', 'tegucigalpa', 14.05, -87.14],
  ['hungary', 'budapest', 47.29, 19.05],
  ['iceland', 'reykjavik', 64.1, -21.57],
  ['india', 'new delhi', 28.37, 77.13],
  ['indonesia', 'jakarta', -6.09, 106.49],
  ['iran', 'tehran', 35.44, 51.3],
  ['iraq', 'baghdad', 33.2, 44.3],
  ['ireland', 'dublin', 53.21, -6.15],
  ['israel', 'jerusalem', 31.71, -35.1],
  ['italy', 'rome', 41.54, 12.29],
  ['jamaica', 'kingston', 18.0, -76.5],
  ['jordan', 'amman', 31.57, 35.52],
  ['kazakhstan', 'astana', 51.1, 71.3],
  ['kenya', 'nairobi', -1.17, 36.48],
  ['kiribati', 'tarawa', 1.3, 173.0],
  ['kuwait', 'kuwait', 29.3, 48.0],
  ['kyrgyzstan', 'bishkek', 42.54, 74.46],
  ['laos', 'vientiane', 17.58, 102.36],
  ['latvia', 'riga', 56.53, 24.08],
  ['lebanon', 'beirut', 33.53, 35.31],
  ['lesotho', 'maseru', -29.18, 27.3],
  ['liberia', 'monrovia', 6.18, -10.47],
  ['libyan arab jamahiriya', 'tripoli', 32.49, 13.07],
  ['liechtenstein', 'vaduz', 47.08, 9.31],
  ['lithuania', 'vilnius', 54.38, 25.19],
  ['luxembourg', 'luxembourg', 49.37, 6.09],
  ['macao, china', 'macau', 22.12, 113.33],
  ['madagascar', 'antananarivo', -18.55, 47.31],
  ['macedonia', 'skopje', 42.01, 21.26],
  ['malawi', 'lilongwe', -14.0, 33.48],
  ['malaysia', 'kuala lumpur', 3.09, 101.41],
  ['maldives', 'male', 4.0, 73.28],
  ['mali', 'bamako', 12.34, -7.55],
  ['malta', 'valletta', 35.54, 14.31],
  ['martinique', 'fort-de-france', 14.36, -61.02],
  ['mauritania', 'nouakchott', -20.1, 57.3],
  ['mayotte', 'mamoudzou', -12.48, 45.14],
  ['mexico', 'mexico', 19.2, -99.1],
  ['micronesia', 'palikir', 6.55, 158.09],
  ['moldova, republic of', 'chisinau', 47.02, 28.5],
  ['mozambique', 'maputo', -25.58, 32.32],
  ['myanmar', 'yangon', 16.45, 96.2],
  ['namibia', 'windhoek', -22.35, 17.04],
  ['nepal', 'kathmandu', 27.45, 85.2],
  ['netherlands', 'amsterdam', 52.23, 4.54],
  ['netherlands antilles', 'willemstad', 12.05, -69.0],
  ['new caledonia', 'noumea', -22.17, 166.3],
  ['new zealand', 'wellington', -41.19, 174.46],
  ['nicaragua', 'managua', 12.06, -86.2],
  ['niger', 'niamey', 13.27, 2.06],
  ['nigeria', 'abuja', 9.05, 7.32],
  ['norfolk island', 'kingston', -45.2, 168.43],
  ['north korea', 'pyongyang', 39.09, 125.3],
  ['northern mariana islands', 'saipan', 15.12, 145.45],
  ['norway', 'oslo', 59.55, 10.45],
  ['oman', 'masqat', 23.37, 58.36],
  ['pakistan', 'islamabad', 33.4, 73.1],
  ['palau', 'koror', 7.2, 134.28],
  ['panama', 'panama', 9.0, -79.25],
  ['papua new guinea', 'port moresby', -9.24, 147.08],
  ['paraguay', 'asuncion', -25.1, -57.3],
  ['peru', 'lima', -12.0, -77.0],
  ['philippines', 'manila', 14.4, 121.03],
  ['poland', 'warsaw', 52.13, 21.0],
  ['portugal', 'lisbon', 38.42, -9.1],
  ['puerto rico', 'san juan', 18.28, -66.07],
  ['qatar', 'doha', 25.15, 51.35],
  ['republic of korea', 'seoul', 37.31, 126.58],
  ['romania', 'bucuresti', 44.27, 26.1],
  ['russia', 'moscow', 55.45, 37.35],
  ['rawanda', 'kigali', -1.59, 30.04],
  ['saint kitts and nevis', 'basseterre', 17.17, -62.43],
  ['saint lucia', 'castries', 14.02, -60.58],
  ['saint pierre and miquelon', 'saint-pierre', 46.46, -56.12],
  ['saint vincent and the greenadines', 'kingstown', 13.1, -61.1],
  ['samoa', 'apia', -13.5, -171.5],
  ['san marino', 'san marino', 43.55, 12.3],
  ['sao tome and principe', 'sao tome', 0.1, 6.39],
  ['saudi arabia', 'riyadh', 24.41, 46.42],
  ['senegal', 'dakar', 14.34, -17.29],
  ['sierra leone', 'freetown', 8.3, -13.17],
  ['slovakia', 'bratislava', 48.1, 17.07],
  ['slovenia', 'ljubljana', 46.04, 14.33],
  ['solomon islands', 'honiara', -9.27, 159.57],
  ['somalia', 'mogadishu', 2.02, 45.25],
  ['south africa', 'pretoria', -25.44, 28.12],
  ['spain', 'madrid', 40.25, -3.45],
  ['sudan', 'khartoum', 15.31, 32.35],
  ['suriname', 'paramaribo', 5.5, -55.1],
  ['swaziland', 'mbabane', -26.18, 31.06],
  ['sweden', 'stockholm', 59.2, 18.03],
  ['switzerland', 'bern', 46.57, 7.28],
  ['syria', 'damascus', 33.3, 36.18],
  ['tajikistan', 'dushanbe', 38.33, 68.48],
  ['thailand', 'bangkok', 13.45, 100.35],
  ['togo', 'lome', 6.09, 1.2],
  ['tonga', "nuku'alofa", -21.1, -174.0],
  ['tunisia', 'tunis', 36.5, 10.11],
  ['turkey', 'ankara', 39.57, 32.54],
  ['turkmenistan', 'ashgabat', 38.0, 57.5],
  ['tuvalu', 'funafuti', -8.31, 179.13],
  ['uganda', 'kampala', 0.2, 32.3],
  ['ukraine', 'kiev', 50.3, 30.28],
  ['united arab emirates', 'abu dhabi', 24.28, 54.22],
  ['united kingdom', 'london', 51.36, -0.05],
  ['united republic of tanzania', 'dodoma', -6.08, 35.45],
  ['united states of america', 'washington dc', 39.91, -77.02],
  ['united states of virgin islands', 'charlotte amalie', 18.21, -64.56],
  ['uruguay', 'montevideo', -34.5, -56.11],
  ['uzbekistan', 'tashkent', 41.2, 69.1],
  ['vanuatu', 'port-vila', -17.45, 168.18],
  ['venezuela', 'caracas', 10.3, -66.55],
  ['viet nam', 'hanoi', 21.05, 105.55],
  ['yugoslavia', 'belgrade', 44.5, 20.37],
  ['zambia', 'lusaka', -15.28, 28.16],
  ['zimbabwe', 'harare', -17.43, 31.02]
]

let obj = {}
points.forEach(a => {
  obj[a[0]] = [a[2], a[3]]
  obj[a[1]] = [a[2], a[3]]
})
module.exports = obj
