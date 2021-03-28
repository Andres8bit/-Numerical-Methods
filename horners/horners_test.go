package horners_test

import (
	"math/rand"
	"testing"
	"time"

	"github.com/andres8bit/numerical_methods/src/horners"
)

func TestNestedFilledCoef(t *testing.T) {
	a := [...]int64{-2, -5, 7, -4, 1} // coef of  p : p(x) = -2 -5x +7x^2 -4x^3 +x^4
	var size int64 = int64(len(a) - 1)
	y := horners.NestedMutliplication(a[:], 3, size)
	if y != int64(19) {
		t.Errorf("incorrect output:")
	}
}

func TestNestedCoefSingle(t *testing.T) {
	s := rand.NewSource(time.Now().UnixNano())
	r := rand.New(s)
	n := r.Int63n(100)
	a := [...]int64{n}
	var size int64 = int64(len(a) - 1)
	y := horners.NestedMutliplication(a[:], 3, size)

	if y != n {
		t.Errorf("incorrect output:")
	}
}

func TestDeflation(t *testing.T) {
	var r int64 = 2
	// coef of  p : p(x) = -2 -5x +7x^2 -4x^3 +x^4
	a := [...]int64{-2, -5, 7, -4, 1}
	size := len(a)
	b := horners.Deflation(a[:], size, r)
	n := len(b)

	match := [...]int64{1, 3, -2, 1}
	for i := 0; i < n; i++ {
		if b[i] != match[i] {
			t.Errorf("incorrect output:")
		}
	}
}

func TestDerivativeAt(t *testing.T) {
	// coef of  p : p(x) = -2 -5x +7x^2 -4x^3 +x^4
	a := [...]int64{-2, -5, 7, -4, 1}
	r := int64(3)
	size := len(a)
	y, z := horners.DevivativeEvaluation(a[:], size, r)
	var x int64 = 37 //dev(p(r)) = 37
	var c int64 = 19 //p(r) = 19
	if y != x && z != c {
		t.Errorf("incorrect output:")
	}
}

func TestCompleteHorners(t *testing.T) {
	a := [...]int64{2, -5, 7, -4, 1}
	r := int64(3)
	tester := [...]int64{23, 37, 25, 8, 1}
	size := len(a)
	horners.CompleteHorners(a[:], size, r)
	for i := 0; i < len(a); i++ {
		if a[i] != tester[i] {
			t.Errorf("incorrect output:")
		}
	}
}
