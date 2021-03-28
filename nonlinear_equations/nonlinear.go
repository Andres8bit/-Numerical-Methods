//  1. Currently using float to hold error thershold, a better option would be to get an interger (err) representing a power,
// then taking 1^-err as our threshold
//  2. Add hybird scheme
// expand to work on complex numbers

package nonlinear

import (
	"fmt"
	"math"

	"github.com/andres8bit/numerical_methods/src/horners"
)

//
//====== Root finders	======

//=========== Bisection: ======
//takes advantage of the Intermediate-Value Theorem:
// 		if f(x) is a continous function on the interval [a,b]
//		and if sign(f(a)) != sign(f(b))
//				then at some point f crosses the y-axis on [a,b]
//				a region in which f is continous
//				therefore it must cross the x-axis at which point we find our root.
// Bisection works bu choosing first a and b
// then splitting the interval into two at each iteration
// until a root is found, think binary search
func Bisection(f []float64, a, b, eps float64, nMax int) float64 {
	size := len(f) - 1
	fa := nestedMultiplication(f[:], a, size)
	fb := nestedMultiplication(f[:], b, size)
	c := 0.0
	if fa*fb > 0 {
		fmt.Println("error: function has the same sign at a & b")
	} else {
		error := float64(b - a)
		for i := 0; i < nMax; i++ {
			error = error / 2
			c = a + error
			fc := nestedMultiplication(f[:], c, size)

			//check for convergence before prep for next iteration
			if math.Abs(error) < eps {
				fmt.Println("convergence at:", c)
				return c
			}

			if fa*fc < 0 {
				b = c
				fb = fc
			} else {
				a = c
				fa = fc
			}

		}
	}

	return c

}

// False Point Bisection:
//	works much like normal Bisection, except rather than using the midpoint
// False Point uses the point where the secant lines of (a,f(a)) & (b,f(b))
func FalsePointBisection(f []float64, a, b, err float64, nMax int) float64 {
	size := len(f) - 1
	fa := nestedMultiplication(f[:], a, size)
	fb := nestedMultiplication(f[:], b, size)
	c := 0.0

	if fa*fb > 0 {
		fmt.Println("error: function has the same sign at a & b")
	} else {
		error := float64(b - a)
		for i := 0; i < nMax; i++ {
			error = error / 2
			afb := a * fb
			bfa := b * fa
			c = (afb - bfa) / (fb - fa) //point where secant line crosses x-axis

			fc := nestedMultiplication(f[:], c, size)

			//check for convergence before prep for next iteration
			if math.Abs(error) < err {
				fmt.Println("convergence at:", c)
				return c
			}

			if fa*fc < 0 { //diff sign
				b = c
				fb = fc
			} else {
				a = c
				fa = fc
			}

		}
	}
	return c
}

// Modified False Point:
// 			Unlike traditional false point, the modified version avoid the pitfall of false point repeatedly taking the same point.
//			if f(a)f(b) < 0 -> signs are diff: regular false point.
//			else: multiply first top and bottom terms by 2
//The goal is to change the slope of the secant line in order to both get closer to the root, and avoiding falling into local minima traps
func ModifiedFalsePoint(f []float64, a, b, err float64, nMax int) float64 {
	size := len(f) - 1
	fa := nestedMultiplication(f[:], a, size)
	fb := nestedMultiplication(f[:], b, size)
	c := 0.0

	if fa*fb > 0 {
		fmt.Println("error: function has the same sign at a & b")
	} else {
		error := float64(b - a)
		for i := 0; i < nMax; i++ {
			error = error / 2
			afb := a * fb
			bfa := b * fa
			c = (afb - bfa) / (fb - fa) //point where secant line crosses x-axis

			fc := nestedMultiplication(f[:], c, size)

			//check for convergence before prep for next iteration
			if math.Abs(error) < err {
				fmt.Println("convergence at:", c)
				return c
			}

			if fa*fc < 0 { //diff sign
				b = c
				fb = fc
			} else {
				afb = 2 * a * fb
				bfa = b * fa
				c = (afb - bfa) / (2*fb - fa)
				a = c
				fa = fc
			}

		}
	}
	return c
}

//Newton's Method:
func NewtonRalph(f []float64, x, err, eps float64, nMax int) float64 {
	for i := 0; i <= nMax; i++ {
		fp, fx := horners.DerivativeEvaluation(f[:], len(f), x)

		if math.Abs(fp) < eps {
			fmt.Println("error small derivative at:", x)
			return x
		}

		delta := fx / fp
		x = x - delta
		if math.Abs(delta) < err {
			fmt.Println("convergence at :", x)
			return x
		}

	}

	fmt.Println("failed to converge after ", nMax, "iterattions, end at x:", x)
	return x
}

// Steffensen's method:
// similar to Newton's but avoids to derivative calclation, instead it uses:
//											 xn+1 = xn - f(g)/g(x)
//											where g(x) = [f(x+f(x))-f(x)]/f(x)
func Steffensen(f []float64, x, err, eps float64, nMax int) float64 {
	n := len(f) - 1
	for i := 0; i <= nMax; i++ {
		fx := nestedMultiplication(f[:], x, n)
		xg := nestedMultiplication(f[:], x+fx, n)
		gx := ((xg - fx) - fx) / fx
		x = x - (fx / gx)

		if math.Abs(fx) < err {
			fmt.Println("convergence at :", x)
			return x
		}

	}

	fmt.Println("failed to converge after ", nMax, "iterattions, end at x:", x)
	return x
}

// Secant Method:
// The Secant Method uses a similar approach to both Bisection and Netwon's Method.
// However unlike Newton's the Secant Method using the Reiman Sum:  dir(f) = lim(h-> infinity) [f(x+h)-f(x)]/h
// Whith a small enough h used itereative to approxamte the Rieman Sum.
// Recall that Newton uses x(n+1) = xn - f(x)/dir(f(x)) to find roots.
// Rewritting our approxamationo for the remain sum using x = xn , h = x(n-1) - xn
// We have dir(f(x)) = [f(xn-1) - f(xn)]/(x(n-1) - xn)
// Rewriting Newton's method for x(n+1) we get:
// 			x(n+1) = xn - [xn-x(n-1)/(f(x) - f(xn-1))]f(x) <- Secant Method for finding roots
func Secant(f []float64, a float64, b float64, sig float64, maxItr int) float64 {
	n := len(f) - 1
	fa := nestedMultiplication(f[:], a, n)
	fb := nestedMultiplication(f[:], b, n)
	if math.Abs(fa) > math.Abs(fb) {
		temp := a
		a = b
		b = temp
		temp = fa
		fa = fb
		fb = temp
	}

	for i := 1; i < maxItr; i++ {
		if math.Abs(fa) > math.Abs(fb) {
			temp := a
			a = b
			b = temp
			temp = fa
			fa = fb
			fb = temp
		}

		delta := (b - a) / (fb - fa)
		b = a
		fb = fa
		delta = delta * fa

		if math.Abs(delta) < sig {
			fmt.Println("converaged at ", a)
			return a
		}
		a = a - delta
		fa = nestedMultiplication(f[:], a, n)
	}

	fmt.Println("after reaching:", maxItr, "failed to fully conveger ended at ", a)
	return a
}

// Illinous Algorithm :
// Is a hybrid scheme for root finding hoping to combine the speed of Secant with the saftey of Bisection.
func Illinois() {

}

//=========== helper funcitons: ================

func nestedMultiplication(f []float64, x float64, N int) float64 {
	p := f[N]
	for i := N - 1; i >= 0; i-- {
		//fmt.Println(" i:", i, " x:", x, "p:", p, "f[i]:", f[i], "f(x):", p)
		p = f[i] + x*p
	}
	return p
}
