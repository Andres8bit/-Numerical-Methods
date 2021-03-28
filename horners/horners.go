package horners

//returns p(x) using nested multiplication in O(n) time
func NestedMutliplication(coef []int64, x, N int64) int64 {
	p := coef[N]
	for i := N - 1; i >= 0; i-- {
		p = coef[i] + x*p
	}
	return p
}

/*
*	Deflation: returns an array containing the coef
*				of the reduced portion of of the polynomial.
*	Ex:
*		p(x) = (x-r)b(x) + p(r)
*     where: r is a root of p(x), p(r)=0 bc r is a root
*	  therefore: we return the coefs of b(x)
*                b(x) = b0 +b1x +b2x^2 +.....+(bn-1)^(x-1)
 */
func Deflation(coef []int64, N int, r int64) []int64 {
	b := make([]int64, N-1)
	n := len(b)
	b[len(b)-1] = coef[N-1]

	for i := n - 1; i >= 1; i-- {
		b[i-1] = coef[i] + r*b[i]
	}

	return b
}

func DevivativeEvaluation(coef []int64, N int, r int64) (alpha, beta int64) {
	beta = 0
	alpha = coef[N-1]
	for i := N - 2; i >= 0; i-- {
		beta = alpha + r*beta
		alpha = coef[i] + r*alpha
	}
	return alpha, beta
}

func DerivativeEvaluation(coef []float64, N int, r float64) (beta, alpha float64) {
	beta = 0.0
	alpha = coef[N-1]
	for i := N - 2; i >= 0; i-- {
		beta = alpha + r*beta
		alpha = coef[i] + r*alpha
	}
	return beta, alpha
}

func CompleteHorners(coef []int64, N int, r int64) {
	//var k int64
	for k := 0; k < N; k++ {
		for j := N - 2; j >= k; j-- {
			coef[j] = coef[j] + r*coef[j+1]
		}
	}
}
