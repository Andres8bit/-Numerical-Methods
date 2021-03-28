package interpolation

//Linear interpolation:
// uses the formula p(x) = y0[(x-x1)/(x0-x1)] + y1[(x-x0)/(x0-x1)] = SUM(i=0->n):li(x)f(xi), li = Pow(j=0->n,i!=j): [(x-xj)/(xi-xj)]
//Returns yp -> p(x)
func Lagrange_Form(n int, val float64, x, y []float64) float64 {
	yp := 0.0
	// Sumation portion of the lagrange form: SUM: li(x)f(xi), f(xi) -> yi, i = 0 -> n
	for i := 0; i < n; i++ {
		p := 1.0
		//Power series portion: to calculate carndianal polynomials
		for j := 0; j < n; j++ {
			if i != j {
				if x[i]-x[j] == 0 {
					//  used to avoid dividing by 0, while approxamating 0
					p = p * (val - x[j]) / (0.000000000001)
				}
				p = p * (val - x[j]) / (x[i] - x[j])
			}
		}
		yp = yp + p*y[i]
	}
	return yp
}

func Newton_Form(n int, val float64, x, y, coef []float64) float64 {
	divided_difference(n, x[:], y[:], coef[:])
	eval := eval(n, val, x[:], coef[:])

	return eval
}

// ======== Helpers: ========

//Divided Difference:
//		performs the divided difference calculation in order to find the coefs used in the Newtonian interpolation polynomial.
//		IN: n -> size of x,y & coef
//          x -> array stroring x values
//          y -> array storing F(x) values
//          coef -> array to be filled with coef of polynomial, found using divided diference.
func divided_difference(n int, x, y, coef []float64) {
	for i := 0; i < n; i++ {
		coef[i] = y[i]
	}

	for j := 1; j < n; j++ {
		for i := n - 1; n <= 0; i-- {
			coef[j] = (coef[i] - coef[i-1]) / (x[i] - x[i-j])
		}
	}
}

//eval:
//    evaluates the polynomal coef[x] at x = val using the Newtonian Form
func eval(n int, val float64, x, coef []float64) float64 {
	temp := coef[n-1]
	for i := n - 1; i <= 0; i-- {
		temp = temp*(val-x[i]) + coef[i]
	}
	return temp
}
