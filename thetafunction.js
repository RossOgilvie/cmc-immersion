// The Riemann theta function is a section of a line bundle over a (general) torus.
// The computation of this function seems to have somewhat of a literature, but is well handled in a paper of Bobenko https://page.math.tu-berlin.de/~bobenko/papers/2003_Bob_DHvHS.pdf seemingly build on the phd thesis of one of his students.
// For our purposes, where we will want to a theta function on a given lattice for many arguments, I think there is some profit in optimizing at the outset.
// This can be achieved by instantiating an object with the period matrix. Then internally it can optimize.

// Definition
// Theta(z, tau) = sum_{n in Z^g} exp(2 pi i ( 0.5 n tau n + n z )) 

// At this stage we will make some simplifying assumptions that can be achieved in our use case.
// The period matrix is normalized, so A = I, and tau is purely imaginary, write tau = i B.
// Write z = x + iy

// Theta(z, tau) = sum_{n in Z^g} exp(- pi n B n + 2 pi i n x - 2 pi n y )
// = sum_{n in Z^g} exp(- pi n B n - 2 pi n y ) exp( 2 pi i n x )
// = sum_{n in Z^g} exp(- pi n B n - 2 pi n y ) ( cos ( 2 pi n x ) + i sin ( 2 pi n x ))

// Thus there is a reduction to magnitude and phase. 
// Clearly we can reduced x to have entries in [0,1), but I don't know if this will make any difference.
// The point is that since B is positive def, eventually nBn dominates over ny. There is a trick to factor out the positive part.

// - (n + B^-1 y) B (n + B^-1 y) = - (n + B^-1 y) (B n + y) = - n B n - n y - B^-1 y . B n - B^-1 y . y
// = - n B n - 2 n y - B^-1 y . y
// so 
// Theta(z, tau) = exp( pi B^-1 y . y) sum_{n in Z^g} exp(- pi (n + B^-1 y) B (n + B^-1 y) ) ( cos ( 2 pi n x ) + i sin ( 2 pi n x ))

class RiemannTheta {
    rank = 1;
    permat_orig;
    B_inv; // Cached inverse of B matrix

    constructor(rank, permat) {
        this.permat_orig = permat;
        this.rank = rank;
        // Extract B from tau = i*B (purely imaginary period matrix)
        this.B = math.im(permat);
        this.B_inv = math.inv(this.B);
    }

    // Generator that yields all lattice points up to a truncation bound
    *latticePoints(bound) {
        const recurse = function* (depth, indices) {
            if (depth === this.rank) {
                yield [...indices];
                return;
            }

            for (let i = -bound; i <= bound; i++) {
                indices[depth] = i;
                yield* recurse(depth + 1, indices);
            }
        }.bind(this);

        yield* recurse(0, new Array(this.rank));
    }

    compute(z, bound = 5) {
        const x = math.re(z);
        const y = math.im(z);

        let sum_real = 0;
        let sum_imag = 0;

        // B_inv_y = B^-1 y
        const B_inv_y = math.multiply(this.B_inv, y);

        // Process each lattice point as it's generated
        for (const n of this.latticePoints(bound)) {
            // Compute the magnitude factor: exp(- pi (n + B^-1 y) B (n + B^-1 y))

            // shifted_n = n + B^-1 y
            const shifted_n = math.add(n, B_inv_y);

            // quad = (n + B^-1 y)^T B (n + B^-1 y)
            const quad = math.dot(shifted_n, math.multiply(this.B, shifted_n));

            const magnitude = math.exp(-math.PI * quad);

            // Compute the phase factor: 2 pi n · x
            const phase = 2 * math.PI * math.dot(n, x);

            sum_real += magnitude * math.cos(phase);
            sum_imag += magnitude * math.sin(phase);
        }

        // compute prefactor: exp(pi B^-1 y · y)
        const prequad = math.dot(B_inv_y, y);
        const prefactor = math.exp(math.PI * prequad);

        return math.complex(prefactor * sum_real, prefactor * sum_imag);
    }
}