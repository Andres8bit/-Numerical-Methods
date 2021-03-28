package bspline

// ===== Package Types ====
type Num float64

type F func(x Num) Num

type BSpline struct {
	a []Num
	h []Num
}

type Schoenberg struct {
	d []Num
}

// ===== End Types =====

// ===== Methods/Functions: =====
// ===== B_Spline =====
func newBSpline(n int, t, y []Num) *BSpline {
	s := new(BSpline)
	s.a = make([]Num, n+1)
	s.h = make([]Num, n+1)
	s.a, s.h = coef(n, t, y)

	return s
}

func (b *BSpline) eval(t []Num, x Num, n int) Num {
	i := n - 2
	for ; i >= 0; i-- {
		if x-t[i] >= 0 {
			break
		}
	}
	i = i + 1
	d := (b.a[i+1]*(x-t[i-1]) + b.a[i]*(t[i]-x-b.h[i+1])) / (b.h[i] + b.h[i+1])
	e := (b.a[i]*(x-t[i-1]+b.h[i-1]) + b.a[i-1]*(t[i-1]-x+b.h[i])) / (b.h[i-1] + b.h[i])

	return (d*(x-t[i-1]) + e*(t[i]-x)) / b.h[i]
}

// ======= End B_Spline =====

// ===== Schoenberg: =====
func newSchoenberg(f F, a, b Num, n int) *Schoenberg {
	s := new(Schoenberg)
	s.d = make([]Num, n+3)
	s.d = schoenberg_coef(f, a, b, n)

	return s

}

func (s *Schoenberg) eval(a, b, x Num, n int) Num {
	h := int(b-a) / n
	k := int(x-a) / (h + (5 / 2))
	p := int(x-a) - (k-(5/2))*h
	c := (s.d[k+1]*Num(p) + s.d[k]*Num(2*h-p)) / Num(2*h)
	e := (s.d[k]*Num(p+h) + s.d[k-1]*Num(h-p)) / Num(2*h)

	return (c*Num(p) + e*Num(h-p)) / Num(h)
}

// ===== End B_Spline =====

// ====== helpers: =====
func coef(n int, t, y []Num) ([]Num, []Num) {
	size := n + 1
	a := make([]Num, size)
	h := make([]Num, size)

	for i := 1; i <= n; i++ {
		h[i] = t[i] - t[i-1]
	}

	h[0] = h[1]
	h[size] = h[n]

	del := Num(-1)
	gam := 2 * y[0]
	p := del * gam
	q := Num(2)

	for i := 0; i < n; i++ {
		r := h[i+1] / h[i]
		del = Num(-1) * r * del
		gam = Num(-1)*gam + (r+1)*y[i]
		p = p + gam*del
		q = q + del*del
	}

	a[0] = (Num(-1) * p) / q

	for i := 1; i < size; i++ {
		a[i] = ((h[i-1]+h[i])*y[i-1] - h[i]*a[i-1]) / h[i-1]

	}

	return a[:], h[:]
}

func schoenberg_coef(f F, a, b Num, n int) []Num {
	d := make([]Num, n+3)
	h := (b - a) / Num(n)

	for i := 1; i < n+2; i++ {
		d[i] = f(a + Num(i-2)*h)
	}
	d[0] = 2*d[1] - d[2]
	d[n+2] = 2*d[n+1] - d[n]
	return d
}
