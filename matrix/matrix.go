package matrix

import (
	"errors"
	"fmt"
	"math"
	"math/cmplx"
)

type F func(complex128) complex128
type Matrix struct {
	m []complex128 //matix made of Numtype
	r int          //row count of matrix
	c int          //col count of matrix
}

//============= Methods: ======================

//============ Constructors: ============
func (m *Matrix) NewMatrix(row, col int) {
	m = NewMatrix(row, col)

}
func NewMatrix(row, col int) *Matrix {
	matrix := new(Matrix)
	matrix.r = row
	matrix.c = col
	if row == 0 || col == 0 {
		if row != 0 {
			matrix.m = make([]complex128, row)
		} else {
			matrix.m = make([]complex128, col)
		}

	} else {
		matrix.m = make([]complex128, row*col)
	}
	return matrix
}

//creates a whole new matrix copying source
func (m *Matrix) Copy(source *Matrix) {
	m.r = source.r
	m.c = source.c
	m.m = make([]complex128, m.r*m.c)
	copy(m.m, source.m)
}

//  ========== End of constructors ==========

// Printer method to print matrix in r x c format:
//need to add formating tools to make output nice for any size matrix

func (m *Matrix) Print() {
	r := m.r
	c := m.c

	for i := 0; i < r; i++ {
		fmt.Print("")
		for j := 0; j < c; j++ {
			fmt.Print(m.Get(i, j), " ")
		}
		fmt.Println(" ")
	}
}

func (m *Matrix) Equals(a *Matrix) bool {
	r := m.r
	c := m.c

	if r != a.r || c != a.c {
		return false
	}

	for i := 0; i < r*c; i++ {
		if m.m[i] != a.m[i] {
			return false
		}

	}
	return true
}

// =========== Getters: ===========
func (m *Matrix) ColCount() int { return m.c }

func (m *Matrix) RowCount() int { return m.r }

func (m *Matrix) Get(row, col int) complex128 {
	if (row < m.r || row == 0) && (col < m.c || col == 0) {
		return m.m[m.c*row+col]
	}
	fmt.Println("Get: out of bounds")

	return 0
}

func (m *Matrix) Dims() (int, int) {
	return m.r, m.c
}
func (m *Matrix) GetRow(row int) []complex128 {
	n := m.c
	r := make([]complex128, n)
	if row >= m.r {
		fmt.Println("GetRow: out of bounds")
		return r
	}

	start := row * n

	for i := 0; i < n; i++ {
		r[i] = m.m[start+i]
	}

	return r
}

func (m *Matrix) GetCol(col int) []complex128 {
	rowC := m.r
	colC := m.c
	c := make([]complex128, rowC)

	if col >= colC {
		fmt.Println("GetCol: out of bounds")
		return c
	}

	for i := 0; i < rowC; i++ {
		c[i] = m.m[col]
		col += colC

	}

	return c
}

// ======== End of Getters: ========

// ======== Setters: ========
func (m *Matrix) Set(row, col int, val complex128) error {
	if (row < m.r || row == 0) && (col < m.c || col == 0) {
		m.m[row*m.c+col] = val

		return nil
	}
	return errors.New("Set:out of bounds")
}

func (m *Matrix) SetRow(vals []complex128, row int) error {

	if row >= m.r && row != 0 {
		return errors.New("SetRow:out of bounds")
	}

	for i := 0; i < m.c; i++ {
		m.m[(row*m.c)+i] = vals[i]
	}
	return nil
}

func (m *Matrix) SetCol(vals []complex128, col int) error {
	if col >= m.c {

		return errors.New("SetCol:out of bounds")
	}

	for i := 0; i < m.r; i++ {
		m.m[i*m.c+col] = vals[i]
	}

	return nil
}

// ============ End of Setters: ======================

//================== End of Methods =====================

//================== Functions ========================
func Naive(source, x, b *Matrix) {
	forwardElimination(source, b)
	backSubstitution(source, b, x)
}

func ScaledPartialPivoting(source, x, b *Matrix) {
	n := source.RowCount()
	var order []int       //pivot order of elimination
	var scaledOrder []int //holds the scaled values of each matrix row
	var scaledVal, index int
	var s []complex128
	var rmax, r float64
	var sourcelik, alik, alkk, xMult, alij, alkj complex128

	//sets the natural order of the matrix
	// & places largest values into scaledvector in order
	for i := 0; i < n; i++ {
		order = append(order, i)
		s = source.GetRow(i)

		scaledVal = max(s[:])
		scaledOrder = append(scaledOrder, scaledVal)
	}
	//main step loop
	for k := 0; k < n; k++ {
		rmax = 0.0
		index = 0
		//sets up pivot order of matrix
		for i := k; i < n; i++ {
			sourcelik = source.Get(order[i], k)
			r = math.Abs(real(sourcelik) / float64(scaledOrder[order[i]]))
			if r > rmax {
				rmax = r
				index = i
			}
		}
		swap(order[:], k, index)
		//forward elimination phase
		for i := k + 1; i < n; i++ {
			alik = source.Get(order[i], k)
			alkk = source.Get(order[k], k)
			xMult = alik / alkk
			//storing xMult in the "zero" locations of matrix
			source.Set(order[i], k, xMult)
			for j := k + 1; j < n; j++ {
				alij = source.Get(order[i], j)
				alkj = source.Get(order[k], j)

				source.Set(order[i], j, alij-xMult*alkj)
			}
		}
	}
	//Proforms backwards propagation phase
	//updating both the b matrix and x matrix placing
	// solutions into the x matrix
	solve(n, order[:], source, b, x)
}

//=============== Tri: =====================
// solves a tridiaganol matrix:
// the matrix is broken up into three  1D matrixes:
// a := lowest diagonal
// d := mid diagonal
// c := top diagonal
// b := b vector
// x := stores solutions
// WARNING: this funciton is not stable
//			that is to say it will overwrite
//			the matrices given
//O(n) = n
func tri(n int, a, d, c, x, b *Matrix) {
	var dCur, bCur, cCur, xMult complex128
	var bn, dn, aPrev, dPrev, bPrev, cPrev, xPrev complex128
	//forward Elimination updates b vector
	for i := 1; i < n; i++ {
		aPrev = a.Get(0, i-1)
		dPrev = d.Get(0, i-1)
		bPrev = b.Get(0, i-1)
		cPrev = c.Get(0, i-1)
		dCur = d.Get(0, i)
		bCur = b.Get(0, i)

		xMult = aPrev / dPrev
		d.Set(0, i, dCur-xMult*cPrev)
		b.Set(0, i, bCur-xMult*bPrev)

	}

	//backSubstitution phase:
	//places answers into x using backSubstitution
	bn = b.Get(0, n-1)
	dn = d.Get(0, n-1)
	x.Set(0, n-1, bn/dn)
	for i := n - 2; i <= 0; i-- {
		bCur = b.Get(0, i)
		cCur = c.Get(0, i)
		dCur = d.Get(0, i)
		xPrev = x.Get(0, i+1)
		x.Set(0, i, (bCur-cCur*xPrev)/dCur)
	}

}

//=============== Penta: =====================
// Solves a 5 diagonal matrix system:
// Accepts matix as 5 1D matrices, and 2 1D matrices x,b:
// e := lowest diagonal
// a := 2nd lowest diagonal
// d := mid diagonal
// c := 2nd topmost diagoanl
// f := topmost diagonal
// x := variable, where answers will be stored in backSubPhase
// b := hold RHS of linear system in order.
// WARNING: this funciton is not stable
//			that is to say it will overwrite
//			the matrices given
func Penta(n int, e, a, d, c, f, x, b *Matrix) {
	var r, s, t, xMult complex128
	var xn, dn, dn1, cn, cn1, bn, bn1 complex128
	var dPrev, cPrev, fPrev, bPrev, xPrev complex128
	var cCur, dCur, bCur, fCur complex128
	var bNext, dNext complex128

	r = a.Get(0, 1)
	s = a.Get(0, 2)
	t = e.Get(0, 1)

	for i := 1; i < n-2; i++ {
		dPrev = d.Get(0, i-1)
		dNext = d.Get(0, i+1)
		dCur = d.Get(0, i)
		cPrev = c.Get(0, i-1)
		cCur = c.Get(0, i)
		fPrev = f.Get(0, i-1)
		bCur = b.Get(0, i)
		bPrev = b.Get(0, i-1)
		bNext = b.Get(0, i+1)
		xMult = r / dPrev

		d.Set(0, i, dCur-xMult*cPrev)
		c.Set(0, i, cCur-xMult*fPrev)
		b.Set(0, i, bCur-xMult*bPrev)

		xMult = t / dPrev
		r = s - xMult*cPrev
		d.Set(0, i+1, dNext-xMult*fPrev)
		b.Set(0, i+1, bNext-xMult*bPrev)
		s = a.Get(0, i+1)
		t = e.Get(0, i)
	}
	dn = d.Get(0, n-1)
	dn1 = d.Get(0, n-2)
	cn = c.Get(0, n-1)
	cn1 = c.Get(0, n-2)
	bn = b.Get(0, n-1)
	bn1 = b.Get(0, n-2)
	xn = x.Get(0, n-1)

	xMult = r / dn1
	d.Set(0, n-1, dn-xMult*cn1)
	x.Set(0, n-1, bn-xMult*bn1/dn)
	x.Set(0, n-2, (bn1-cn*xn)/dn1)

	for i := n - 3; i <= 0; i-- {
		bCur = b.Get(0, i)
		fCur = f.Get(0, i)
		xPrev = x.Get(0, i+1)
		dCur = d.Get(0, i-1)
		cCur = c.Get(0, i)

		x.Set(0, i, (bCur-fCur*xPrev-cCur*xPrev)/dCur)
	}

}

//=============== Doolittle Factorization: =====================
func Doolittle(a *Matrix) (lower, upper *Matrix) {
	n := a.ColCount()
	var sum, lowerVal complex128
	lower = NewMatrix(n, n)
	upper = NewMatrix(n, n)
	for k := 0; k < n; k++ {
		lower.Set(k, k, 1)
		for j := k; k < n; k++ {
			sum = complex128(0)
			for s := 0; s < k; s++ {
				sum += lower.Get(k, s) * upper.Get(s, j)
			}
			upper.Set(k, j, a.Get(k, j)-sum)
		}

		for i := k + 1; i < n; i++ {
			sum = complex128(0)

			for s := 0; s < k; s++ {
				sum += lower.Get(i, s) * upper.Get(s, k)

			}
			lowerVal = (a.Get(i, k) - sum) / upper.Get(k, k)

			lower.Set(i, k, lowerVal)
		}
	}
	return lower, upper
}

// Foward Substitution phase:
// works with Dooolittle Factorization, taking the lower unit matrix and using it to sovle for z
// in Ax = b -> LUx = b -> Lz = b -> Ux = z
// INPUT : The user must pass both the lower matrix (l) and the b vector matrix
// OUTPUT : zi the z vector being solved for.

func LUForwardSub(l, b *Matrix) (z *Matrix) {
	var sum complex128
	n := b.ColCount()
	z = NewMatrix(1, n)
	z.Set(0, 0, b.Get(0, 0))

	for i := 1; i < n; i++ {
		sum = complex128(0)
		for j := 0; j < i; j++ {
			sum += l.Get(i, j) + z.Get(0, j)
		}

		z.Set(0, i, b.Get(0, i)-sum)
	}

	return z
}

func LUBackSub(u, z *Matrix) (x *Matrix) {
	var sum complex128
	n := u.ColCount()
	x = NewMatrix(0, n)
	x.Set(0, n-1, z.Get(0, n-1)/u.Get(n-1, n-1))

	for i := n - 2; i <= 0; i-- {
		sum = complex128(0)
		for j := i + 1; j < n; j++ {
			sum += u.Get(i, j) * x.Get(0, j)
		}
		sum = sum / u.Get(i, i)
		x.Set(0, i, z.Get(0, i)-sum)

	}

	return x
}

func LUUpdateRHSB(a, b *Matrix) {
	var val complex128
	n := a.ColCount()
	for k := 0; k < n; k++ {
		for i := k + 1; i < n; i++ {
			val = b.Get(0, i) - a.Get(i, k)*b.Get(0, k)
			b.Set(0, i, val)
		}
	}

}

//LDLTranspose Factorization:
// INPUT: a and its diagonal D matrix.
// OUTPUT: L lower triangular Matrix
func LDLTFactorization(a *Matrix) (l, d *Matrix) {
	var sum, ljv, dia complex128
	n := a.ColCount()
	l = NewMatrix(n, n)
	d = NewMatrix(0, n)

	for j := 0; j < n; j++ {
		l.Set(j, j, 1)
		sum = complex128(0)

		for v := 0; v < j; v++ { //does not iterate on first loop
			ljv = l.Get(j, v) // this leaves d(0,0) = a(0,0)
			ljv = ljv * ljv
			dia = d.Get(0, v)
			sum += dia * ljv
		}

		d.Set(0, j, a.Get(j, j)-sum)

		for i := j + 1; i < n; i++ {
			l.Set(j, i, 0)
			sum := complex128(0)

			for v := 1; v < j; v++ {
				sum += l.Get(i, v) * d.Get(0, v) * l.Get(j, v)
			}
			sum = (a.Get(i, j) - sum) / d.Get(0, j)
			l.Set(i, j, sum)
		}
	}

	return l, d
}

// Cholesky Factorization:
// Given: a Symmetric, Positive Definite (SPD) Matrix
// Returns: l a Lower Triangular Matrix
// that is to say given A returns L in
// A = LL^T
func CholeskyFactorization(a *Matrix) (l *Matrix) {
	n := a.ColCount()
	l = NewMatrix(n, n)

	for k := 0; k < n; k++ {
		sum := complex128(0)
		for s := 0; s < k; s++ {
			lval := l.Get(k, s)
			sum += lval * lval
		}

		sum = a.Get(k, k) - sum
		sum = cmplx.Sqrt(sum)
		l.Set(k, k, sum)

		for i := k + 1; k < n; i++ {
			sum := complex128(0)
			for s := 0; s < k; s++ {
				sum += l.Get(i, s) * l.Get(k, s)
			}
			sum = a.Get(i, k) - sum
			sum = sum / l.Get(k, k)
			l.Set(i, k, sum)
		}
	}

	return l
}

func Power(a, x, y *Matrix, varphi F, kMax int) (complex128, []complex128) {
	n := a.ColCount()
	eig := complex128(0)
	eigVect := make([]complex128, n)

	for i := 0; i < kMax; i++ {
		//y =

	}

	return eig, eigVect
}

// ungar's/3m method"
// Given a and b are matrices, such that a = A1 +iA2, b = B1 +iB2
// define C as C = A*B = C1 +iC2
// then we can construct C using the following:
// T1 = A1*B1  , T2 = A2*B2
// C1 = T1-T2, C2 = [A1+A2][B1+B2] -T1-T2 = [A1+A2][B1+B2] - C1
func Mult(a, b *Matrix) (result *Matrix) {

	return result

}

//================ End of Functions ======================

//================ helper funcitons ======================

//================ Norms ============================
func L1Norm(X []float64) float64 {
	var sum float64
	sum = 0
	for i := range X {
		sum += math.Abs(X[i])
	}
	return sum
}

func EuclideanNorm(X []float64) float64 {
	var sum float64
	sum = 0
	for i := range X {
		sum += math.Pow(X[i], 2)
	}

	return math.Sqrt(sum)
}

func InfiniteNorm(X []float64) float64 {
	var maxSoFar float64
	var cur float64
	maxSoFar = 0
	for i := range X {
		cur = math.Abs(X[i])
		if maxSoFar <= cur {
			maxSoFar = cur
		}
	}
	return maxSoFar
}

//================ End of Norms ======================
//	backwards propagation of partial pivoting
func solve(n int, l []int, a, b, x *Matrix) {

	//update b vector values
	for k := 0; k < n; k++ {
		for i := k + 1; i < n; i++ {
			tmp := b.Get(0, l[i])
			tmp2 := b.Get(0, l[k])
			coef := a.Get(l[i], k)
			b.Set(0, l[i], tmp-coef*tmp2)
		}
	}

	bln := b.Get(0, l[n-1])
	alnn := a.Get(l[n-1], l[n-1])
	x.Set(0, n-1, bln/alnn)

	//backprobagation phase
	for i := n - 1; i >= 0; i-- {
		sum := b.Get(0, l[i])
		for j := i + 1; j < n; j++ {
			xVal := x.Get(0, j)
			aVal := a.Get(l[i], j)
			sum = sum - xVal*aVal
		}
		aVal := a.Get(l[i], i)
		x.Set(0, i, sum/aVal)
	}

}

//helper function to return the largest coeficient:
func max(source []complex128) int {
	maxVal := 0
	n := len(source)
	for i := 0; i < n; i++ {
		tmp := real(source[i])
		val := int(math.Abs(tmp))

		if val > maxVal {
			maxVal = val
		}
	}

	return maxVal
}

//helper funcction to swap the values of a slice
//at indecies source and dest
func swap(a []int, source, dest int) {
	tmp := a[source]
	a[source] = a[dest]
	a[dest] = tmp

}

//helper funciton for gaussian elimination, pivots inorder
func forwardElimination(source, b *Matrix) {
	var xMult complex128
	n := int(source.ColCount())

	for k := 0; k < n; k++ {
		for i := k + 1; i < n; i++ {
			coef := source.Get(i, k)
			div := source.Get(k, k)
			xMult = coef / div
			source.Set(i, k, xMult)

			for j := k + 1; j < n; j++ {
				coef := source.Get(i, j)
				mult := source.Get(k, j)
				val := coef - (xMult)*mult
				source.Set(i, j, val)
			}

			finVal := b.Get(0, i)
			finMult := b.Get(0, k)
			b.Set(0, i, finVal-(xMult*finMult))
		}
	}
}

//helper function for gaussina elimination, finds solutions to source, plavins them in x
func backSubstitution(source, b, x *Matrix) {
	n := source.ColCount()
	bn := b.Get(0, n-1)
	ann := source.Get(n-1, n-1)
	x.Set(0, n-1, bn/ann)

	for i := n - 1; i >= 0; i-- {
		sum := b.Get(0, i)
		for j := i + 1; j < n; j++ {
			ann = source.Get(i, j)
			xj := x.Get(0, j)
			sum = sum - ann*xj
		}
		ann = source.Get(i, i)

		x.Set(0, i, sum/ann)
	}
}

// Parallel Matrix multiplication to be used Mult.
// making complex multiplication aseries of normal floating point matrix multiplication
func Mult(a, b []float64, n, blockSize int) []float64 {

}

//naive multiplication method used in Paralllel matrix multiplication
// treats the underlying math as rows of A by all the cols of B
func naiveMult(a, b *Matrix) *Matrix {
	r := a.RowCount()
	c := b.ColCount()
	if r != c {
		return nil
	}

	result := NewMatrix(r, c)
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			val := dot(a.GetRow(i), b.GetCol(j))
			result.Set(i, j, complex(val, 0))
		}
	}

	fmt.Println("new matrix:", result)
	return result
}

//returns dot product but only on real portion of complex number,
//returns a floating point number
// meaning: result := x + 0i, where x:= Sum (ai*bi) i->len(a)
func dot(a, b []complex128) float64 {
	result := 0.0
	n := len(a)
	for i := 0; i < n; i++ {
		result += real(a[i]) * real(b[i])
	}

	return result
}

// n: row
// m: col
func blocking(a, b []float64, n, m, bSize int) []float64 {
	sum := 0.0
	chunk := bSize * (n / bSize) // chunk size to split marticies into
	c := make([]float64, n*m)

	for i := 0; i < chunk; i += bSize {
		for j := 0; j < chunk; j += bSize {
			for k := 0; k < n; k++ {
				for l := j; l < j+bSize; j++ {
					sum = c[m*k+l]
					for kk := i; kk < i+bSize; kk++ {
						sum += a[kk]

					}
				}

			}
		}

	}

	return c

}

//================ End of Helpers ======================
