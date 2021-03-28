package derivative

import "math"

type numVal float64


// All derivative approximations are found using Richardson Exprapolation
func First_Derivative(f []numVal, x, h numVal, n int, d [][]float64) {
	size := len(f) - 1

	for i := 0; i < n; i++ {
		fxplus := nested_Mutliplication(f[:], x+h, size)
		fxminus := nested_Mutliplication(f[:], x-h, size)
		
		d[i][0] = float64((fxplus - fxminus) / (2.0 * h))
		for j := 1; j < i; j++ {
			d[i][j] = d[i][j-1] + (d[i][j-1]-d[i-1][j-1])/(math.Pow(4, float64(j))-1.0)
		}
		h = h / 2.0
	}
}

func Second_Derivative(f []numVal, x, h numVal, n int, d [][]float64) {
	size := len(f) - 1

	for i := 0; i < n; i++ {
		fx := nested_Mutliplication(f[:], x, size)
		fxplus := nested_Mutliplication(f[:], x+h, size)
		fxminus := nested_Mutliplication(f[:], x-h, size)
		d[i][0] = float64((fxplus - 2*fx + fxminus) / (h * h))

		for j := 1; j < i; j++ {
			d[i][j] = d[i][j-1] + (d[i][j-1]-d[i-1][j-1])/(math.Pow(4, float64(j))-1.0)
		}
		h = h / 2
	}

}

// =========== helper function: ===========

//returns p(x) using nested multiplication in O(n) time
func nested_Mutliplication(coef []numVal, x numVal, N int) numVal {
	p := coef[N]
	for i := N - 1; i >= 0; i-- {
		p = coef[i] + x*p
	}
	return p
}
