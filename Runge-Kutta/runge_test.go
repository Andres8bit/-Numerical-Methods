package runge_test

import (
	"fmt"
	"math"
	"testing"

	"github.com/andres8bit/numerical_methods/src/runge-kutta"
)

func TestBaseRunge(t *testing.T) {
	ans := 3.192937699
	n := 73.0
	a := 1.0
	b := 1.5625
	x := 2.0
	h := (b - a) / n
	ti := a

	f := func(y, x float64) float64 {
		tmp := x - y - 1
		tmp = tmp * tmp
		return 2 + tmp
	}
	val := runge.RK4(f, ti, x, h, int(n))
	err := math.Abs(ans-float64(val)) / ans
	fmt.Println("val:", val)

	if err*100 >= 0.5 {
		fmt.Println("err:", err*100, "%")
		t.Errorf("invalid solution error rate:")

	}
}

func TestRkf45(t *testing.T) {
	ans := 3.192937699
	n := 72
	a := 1.0
	b := 1.5625
	x := 2.0
	h := (b - a) / float64(n)
	ti := a
	errPercent := 0.0

	f := func(y, x float64) float64 {
		tmp := x - y - 1
		tmp = tmp * tmp
		return 2 + tmp
	}

	// moves from 1 to 1.5625 in 72 steps of equal h size
	for i := 0; i < n; i++ {
		xVal, tVal, _ := runge.RKF45(f, ti, x, h)
		ti = tVal
		x = xVal
	}

	errPercent = math.Abs(ans-x) / ans

	if errPercent >= 0.5 {
		fmt.Println("err:", errPercent*100, "%")
		t.Errorf("invalid solution error rate:")

	}
}

func TestAdaptiveRKF45(t *testing.T) {
	// f := func(y, x float64) float64 {
	// 	b := 5 * math.Sin(y)
	// 	c := 0.2 * x

	// 	return 3 + b + c
	// }

	f := func(y, x float64) float64 {
		tmp := x - y - 1
		tmp = tmp * tmp
		return 2 + tmp
	}
	n := 1000
	e_max := 0.0000001
	e_min := 0.00000001
	h_min := 0.000001
	h_max := 1.0
	ti := 0.0
	tb := 1.0
	x := 2.0
	h := (tb - ti) / float64(n)

	val := runge.RKFAdaptive(f, ti, x, h, tb, e_max, e_min, h_min, h_max, n)

	fmt.Println("val:", val)

}
