package integration_test

import (
	"math"
	"testing"

	integral "github.com/andres8bit/numerical_methods/src/integration"
)

func testFunc(x integral.NumVal) integral.NumVal {
	tmp := float64(x)
	sineVal := integral.NumVal(math.Sin(tmp))
	if x == 0 {
		return 1.0
	}
	return sineVal / x
}
func TestTrapizoid(t *testing.T) {
	val := integral.Trapizoid2(testFunc, 0, 1, 6)
	err := math.Abs(0.94508-float64(val)) / 0.94508

	// If our percent error is greater than .1 percent fail test.
	if err >= 0.001 {
		t.Errorf("did not return expected val")
	}

}
