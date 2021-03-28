package eulers

type Num float64

type F func(Num, Num) Num

// Eulers method to find approximate values of the solutions to the intial-value promblem:
//				x' = f(t,x(t))
//				x(a) = xa
func Eulers(n int, a, b, x Num, f F) Num {
	h := b - a/Num(n)
	t := a

	for k := 0; k < n; k++ {
		x = x + h*f(t, x)
		t = t + h
	}

	return x
}
