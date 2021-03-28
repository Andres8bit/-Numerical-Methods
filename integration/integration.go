package integration

import (
	"fmt"
	"math"
)

type NumVal float64
type Func func(x NumVal) NumVal

// Trapizoid:
// Uses the Trapizoidal method to approximate the integral
// of any polynomial funstion f, where the values of f[i]
// represent the coefs. of the polynomial.
func Trapizoid(f []NumVal, a, b NumVal, n int) NumVal {
	h := (b - a) / NumVal(n)
	size := len(f) - 1
	fa := nested_Mutliplication(f[:], a, size)
	fb := nested_Mutliplication(f[:], b, size)
	sum := 0.5 * (fa + fb)

	for i := 1; i < n; i++ {
		x := a + NumVal(i)*h
		fx := nested_Mutliplication(f[:], x, size)
		sum = sum + fx
	}
	sum = sum * h

	return sum

}

// Trapizoid:
// Uses the Trapizoidal method to approximate the integral
// of any single variable function f. That is conitnous on [a,b] & 1:1
func Trapizoid2(f func(x NumVal) NumVal, a, b NumVal, n int) NumVal {
	h := (b - a) / NumVal(n)
	fa := f(a)
	fb := f(b)
	sum := 0.5 * (fa + fb)

	for i := 1; i < n; i++ {
		x := a + NumVal(i)*h
		fx := f(x)
		sum = sum + fx

	}
	sum = sum * h

	return sum
}

//Romberg:
// Uses Romberg's method to create an (n x n) array
// In which numerical approximations of the integral is stored.
// Unlike traditional Romberg's method this version transposes the matrix.
// This is too avoid repeated cache misses, constructing matrix row wise rather than comlunm wise
func Romberg(f []NumVal, a, b NumVal, r [][]NumVal) {
	h := b - a
	size := len(f) - 1
	n := len(r)
	fa := nested_Mutliplication(f[:], a, size)
	fb := nested_Mutliplication(f[:], b, size)
	r[0][0] = (h / 2) * (fa + fb)

	for i := 1; i < n; i++ {
		h = h / 2
		sum := NumVal(0)
		for k := 1; k < int(math.Pow(2, float64(i))); k = k + 2 {
			fakh := nested_Mutliplication(f[:], NumVal(a+NumVal(k)*h), size)
			sum = sum + fakh
		}
		r[i][0] = 0.5*r[i-1][0] + sum*h

		for j := 1; j <= i; i++ {
			r[i][j] = r[i][j-1] + (r[i][j-1]-r[i-1][j-1])/NumVal(math.Pow(4, float64(j))-1.0)
		}

	}
}

func Romberg2(f func(x NumVal) NumVal, a, b NumVal, r [][]NumVal) {
	h := b - a
	n := len(r)
	fa := f(a)
	fb := f(b)
	r[0][0] = (h / 2) * (fa + fb)

	for i := 1; i < n; i++ {
		h = h / 2
		sum := NumVal(0)
		for k := 1; k < int(math.Pow(2, float64(i))); k = k + 2 {
			fakh := f(a + NumVal(k)*h)
			sum = sum + fakh
		}
		r[i][0] = 0.5*r[i-1][0] + sum*h

		for j := 1; j <= i; i++ {
			r[i][j] = r[i][j-1] + (r[i][j-1]-r[i-1][j-1])/NumVal(math.Pow(4, float64(j))-1.0)
		}

	}
}

// Simpson:
// Basic Simpson's Integration rule
func Simpson(f []NumVal, a, b NumVal, err float64, level, level_max int) NumVal {
	level = level + 1
	n := len(f) - 1
	h := b - a
	c := (a + b) / 2
	d := (a + c) / 2
	e := (c + b) / 2
	fa := nested_Mutliplication(f[:], a, n)
	fb := nested_Mutliplication(f[:], b, n)
	fc := nested_Mutliplication(f[:], c, n)
	fd := nested_Mutliplication(f[:], d, n)
	fe := nested_Mutliplication(f[:], e, n)

	oneSimpson := h*(fa+4*fc) + fb/NumVal(6)
	twoSimpson := h * (fa + 4*fd + 2*fc + 4*fe + fb) / 12.0

	if level >= level_max {
		fmt.Println("max level reacched")
		return twoSimpson

	}

	if math.Abs(float64(twoSimpson-oneSimpson)) < 15*err {
		return (twoSimpson + (twoSimpson - oneSimpson)) / 15

	}

	//recursive calls to apply adaptive simpsons to both side of the curve.
	leftSimpson := Simpson(f, a, c, err/2, level, level_max)
	rightSimpson := Simpson(f, c, b, err/2, level, level_max)

	return leftSimpson + rightSimpson

}

// Adaptive Simpson's Integration rule:
//  Dynamically split each interval of the curve into the best number of partions.
func Simpson2(f func(x NumVal) NumVal, a, b NumVal, err float64, level, level_max int) NumVal {
	level = level + 1
	h := b - a
	c := (a + b) / 2
	d := (a + c) / 2
	e := (c + b) / 2
	fa := f(a)
	fb := f(b)
	fc := f(c)
	fd := f(d)
	fe := f(e)

	oneSimpson := h*(fa+4*fc) + fb/NumVal(6)
	twoSimpson := h * (fa + 4*fd + 2*fc + 4*fe + fb) / 12.0

	if level >= level_max {
		fmt.Println("max level reacched")
		return twoSimpson

	}

	if math.Abs(float64(twoSimpson-oneSimpson)) < 15*err {
		return (twoSimpson + (twoSimpson - oneSimpson)) / 15

	}

	//recursive calls to apply adaptive simpsons to both side of the curve.
	leftSimpson := Simpson2(f, a, c, err/2, level, level_max)
	rightSimpson := Simpson2(f, c, b, err/2, level, level_max)

	return leftSimpson + rightSimpson

}

// ======== Helper functions ===========
//returns coef(x) using nested multiplication in O(n) time.
// used for the first trapiziod method which takes an array of ceof.
// representing a polynomial equation of degree len(coef)
func nested_Mutliplication(coef []NumVal, x NumVal, N int) NumVal {
	p := coef[N]
	for i := N - 1; i >= 0; i-- {
		p = coef[i] + x*p
	}
	return p
}
