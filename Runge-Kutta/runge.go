package runge

import (
	"fmt"
	"math"
)

type Num float64

type F func(float64, float64) float64
type FV func(t float64, X []float64) []float64

// implements 4th order Runge-Kutta method:

func RK4(f F, t, x, h float64, n int) float64 {
	ta := t
	tx := x
	N := float64(n)
	for i := 0.0; i < N; i++ {
		k1 := h * f(t, tx)
		k2 := h * f(t+(0.5*h), tx+(0.5*k1))
		k3 := h * f(t+(0.5*h), tx+(0.5*k2))
		k4 := h * f(t+h, tx+k3)

		k := k1 + 2*k2 + 2*k3 + k4
		tx = tx + 0.1666666666666666666667*k
		t = ta + i*h

	}
	return tx
}

//implements Runge-Kutta-Fehlberg method:
//		Uses 5th order runge-kutta to approximate IVP
func RKF45(f F, t, x, h float64) (ans, tb, err float64) {
	//Constants:
	c20 := 0.25
	c21 := 0.25
	c30 := 0.3750000000000
	c31 := 0.0937500000000
	c32 := 0.2812500000000000
	c40 := 0.9230769230769230769230769
	c41 := 0.8793809740555302685480200
	c42 := -3.2771961766044606281292671825216
	c43 := 3.3208921256258534365043240782886
	c51 := 2.0324074074074074074074074074074
	c52 := -8.0000000000000
	c53 := 7.1734892787524366471734892787524
	c54 := -0.2058966861598
	c60 := 0.500000000000000
	c61 := -0.2962962962962962962962962962963
	c62 := 2.00000000000000
	c63 := -1.3816764132553606237816764132554
	c64 := 0.45297270955165
	c65 := -0.27500000000000
	a1 := 0.1157407407407407407
	a3 := 0.5489278752436647173489
	a4 := 0.53533138401559454191033
	a5 := -0.200000000000
	b1 := 0.118518518518518518518
	b3 := 0.5189863547758284600389863547
	b4 := 0.50613149034201665780613149034202
	b5 := -0.18000000000000
	b6 := 0.0363636363636363636363636364

	k1 := h * f(t, x)
	k2 := h * f(t+c20*h, x+c21*k1)
	k3 := h * f(t+c30*h, x+c31*k1+c32*k2)
	k4 := h * f(t+c40*h, x+c41*k1+c42*k2+c43*k3)
	k5 := h * f(t+h, x+c51*k1+c52*k2+c53*k3+c54*k4)
	k6 := h * f(t+c60*h, x+c61*k1+c62*k2+c63*k3+c64*k4+c65*k5)

	x4 := x + a1*k1 + a3*k3 + a4*k4 + a5*k5

	ans = x + b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6
	tb = t + h
	err = math.Abs(float64(ans - x4))
	return ans, tb, err
}

func RKFAdaptive(f F, t, x, h, tb, eps_max, eps_min, h_min, h_max float64, itr_max int) float64 {
	delta := 0.000005
	iflag := 1 // stop

	for k := 0; k <= itr_max; {
		k = k + 1
		if math.Abs(h) < h_min {
			if h < 0 {
				h = -1 * h_min
			} else {
				h = h_min
			}
		}

		if math.Abs(h) > h_max {
			if h < 0 {
				h = -1 * h_max
			} else {
				h = h_max
			}
		}

		d := math.Abs(tb - t)

		if d <= math.Abs(h) {
			iflag = 0 // sucessful march from ta-tb

			if d <= delta*math.Max(tb, t) {
				if h < 0 {
					h = -1 * d
				} else {
					h = d
				}
				break
			}
		}

		xSave := x
		tSave := t

		xTemp, tTemp, err := RKF45(f, t, x, h)
		x = xTemp
		t = tTemp

		if iflag == 0 {
			break
		}

		if err < eps_min {
			h = 2 * h
		}

		if err > eps_max {
			h = h / 2
			x = xSave
			t = tSave
			k = k - 1
		}
		fmt.Println("k:", k, " x:", x, " t:", t, "tb:", tb, " err:", err, " h:", h, " d:", d, "iflag: ", iflag)
	}

	return x

}

// 4th order Runge-Kutta on a System on n equations, with a limit on steps
// n is the size of the system.
// h is the step size
// t is time
// xi is a vector of initial values X
// Fv is F(t,X)
// maxIterations is the maximum number of iteration steps
// ODE Sytem of form:
// {
// | X' = F(t,X)
// | X(a) = S
// {
// Where X = [x1,...xn]^T, X' = [x'1,x'2,x'3]^T, F(t,x) = [f(t,x)1,...fn(t,x)]^T, S = [s1,...sn]^T
// user must supply FV  F(t,x) vector to evaluate for each K(Runge-Kutta) vector (k1-k4)
// returns solutions and places them in xi
func RK4System(h, t float64, xi []float64, system FV, n, maxItr int) {
	y := make([]float64, n)
	for j := 0; j < maxItr; j++ {

		k1 := system(t, xi[:])

		for i := 0; i < n; i++ {
			y[i] = xi[i] + 0.5*h*k1[i]

		}

		k2 := system(t+0.5*h, y[:])

		for i := 0; i < n; i++ {
			y[i] = xi[i] + 0.5*h*k2[i]
		}

		k3 := system(t+0.5*h, y[:])

		for i := 0; i < n; i++ {
			y[i] = xi[i] + h*k3[i]
		}

		k4 := system(t+h, y[:])
		for i := 0; i < n; i++ {
			xi[i] = xi[i] + 0.16666666666666666666666666666667*h*(k1[i]+2*k2[i]+2*k3[i]+k4[i])
		}
		t = t + h
	}
}
