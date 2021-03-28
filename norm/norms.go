package norm

import "math"

func L1Norm(X []float64) float64 {
	var sum float64
	sum = 0
	for i := range X {
		sum += math.Abs(X[i])
	}
	return sum
}

func EuclideanNorm(X []float64) float64 {
	var sum float64
	sum = 0
	for i := range X {
		sum += math.Pow(X[i], 2)
	}

	return math.Sqrt(sum)
}

func InfiniteNorm(X []float64) float64 {
	var maxSoFar float64
	var cur float64
	maxSoFar = 0
	for i := range X {
		cur = math.Abs(X[i])
		if maxSoFar <= cur {
			maxSoFar = cur
		}
	}
	return maxSoFar
}
