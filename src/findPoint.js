import points from './points'
const findPoint = function (input) {
  if (points.hasOwnProperty(input)) {
    return points[input]
  }
  return input
}
export default findPoint
