import points from '../data/points'
const findPoint = function (input) {
  if (points.hasOwnProperty(input)) {
    return points[input]
  }
  return input
}
export default findPoint
