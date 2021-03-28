package nonlinear_test

import (
	"fmt"
	"testing"

	nonlinear "github.com/andres8bit/numerical_methods/src/nonlinear_equations"
)

func TestBisection(t *testing.T) {
	f := [...]float64{1, -3, 0, 1}
	nonlinear.Bisection(f[:], 0, 1, 0.0001, 100)
	fmt.Println()

}

func TestFalsePointPosition(t *testing.T) {
	f := [...]float64{1, -3, 0, 1}
	nonlinear.FalsePointBisection(f[:], 0, 1, 0.0001, 100)
	fmt.Println()
}

func TestModifiedFalsePosition(t *testing.T) {
	f := [...]float64{1, -3, 0, 1}
	nonlinear.ModifiedFalsePoint(f[:], 0, 1, 0.0001, 100)
	fmt.Println()
}

func TestNewtonRaph(t *testing.T) {
	f := [...]float64{-8, -8, 17, -1, -9, 1}
	nonlinear.NewtonRalph(f[:], 0, 0.000015, 0.000015, 10)
	fmt.Println()
}

func TestSteffensen(t *testing.T) {
	f := [...]float64{-20, 10, 2, 1}
	nonlinear.Steffensen(f[:], 0, 0.000015, 0.000015, 100)
	fmt.Println()
}

func TestSecantMethod(t *testing.T) {
	f := []float64{3, 0, 0, 1, 0, 1}
	nonlinear.Secant(f[:], -1.0, 1.0, 0.00000001, 10)
}
