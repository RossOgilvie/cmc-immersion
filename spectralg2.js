var props = {
    r1: 0.5,
    r2: 0.3,
    phi: math.PI / 4,
};

function dist_c(z, w) {
    return math.abs(math.subtract(z, w))
}

function make_pt(x, y) {
    return { x, y };
}

function rho(pt) {
    let cji = pt.x.conjugate().inverse();
    return make_pt(cji, math.multiply(math.pow(cji, 3), pt.y.conjugate()));
}


// the polynomial function in x, with the given coefficients, beginning from the constant term.
function poly_coeff(coeff, x) {
    let total = math.complex(0, 0);
    for (let i = coeff.length - 1; i >= 0; i--) {
        total = math.chain(total).multiply(x).add(coeff[i]).done();
    }
    return total;
}


class Hyperelliptic_Curve {
    genus;
    alpha1 = math.complex(0.5, 0.5);
    alpha2 = math.complex(0.5, -0.3);
    roots;


    update() {
        this.genus = 2;
        this.alpha1 = math.complex({ r: props.r1, phi: props.phi });
        this.alpha2 = math.complex({ r: props.r2, phi: -props.phi });
        const alpha1cji = this.alpha1.conjugate().inverse();
        const alpha2cji = this.alpha2.conjugate().inverse();
        this.roots = [this.alpha1, this.alpha2, alpha1cji, alpha2cji];
    }

    // The polynomial that defines the spectral curve, y^2 = x*a(x)
    a_poly(x) {
        let a_poly_val = math.complex(1, 0);
        for (let i = 0; i < this.roots.length; i++) {
            const factor = math.subtract(x, this.roots[i]);
            a_poly_val = math.multiply(a_poly_val, factor);
        }
        return a_poly_val;
    }
    // The polynomial x*a(x) = y^2, including the x factor that we usually handle separately in our working.
    xa_poly(x) {
        return math.multiply(x, this.a_poly(x));
    }

    // The polynomial a(x), with the drop_alpha_idx-th root removed
    a_poly_hat(drop_alpha_idx, x) {
        let a_poly_hat_val = math.complex(1, 0);
        for (let i = 0; i < this.roots.length; i++) {
            const factor = i != drop_alpha_idx ? math.subtract(x, this.roots[i]) : math.complex(1, 0);
            a_poly_hat_val = math.multiply(a_poly_hat_val, factor);
        }
        return a_poly_hat_val;
    }
    // The polynomial x*a(x), with the drop_alpha_idx-th root removed. -1 indicates to drop the x factor.
    xa_poly_hat(drop_alpha_idx, x) {
        if (drop_alpha_idx == -1) {
            return this.a_poly(x);
        } else {
            return math.multiply(x, this.a_poly_hat(drop_alpha_idx, x));
        }
    }

    // The derivative of a(x)
    a_poly_prime(x) {
        let a_poly_prime_val = math.complex(0, 0);
        // The terms of the derivative are simply excluding a factor from the product that defines a_poly
        for (let i = 0; i < this.roots.length; i++) {
            a_poly_prime_val = math.add(a_poly_prime_val, this.a_poly_hat(i, x));
        }

        return a_poly_prime_val;
    }
    // The derivative of x*a(x) = a(x) + x*a'(x)
    xa_poly_prime(x) {
        return math.add(this.a_poly(x), math.multiply(x, this.a_poly_prime(x)));
    }


    isTopSheet(pt) {
        // I'm a little confused about how the sheet works on the unit circle, just use a simple function for the moment.
        // We know on the unit circle that x^{-2}a(x) >= 0. So the sign of the square root labels 
        // , labelled by the sign of re part of y/x^2. This expression is purely real and pos on unit circle.
        return pt.y.im >= 0;
    }
}


class Spectral_Curve {
    hyp;

    F1;
    F2;
    F1rho;
    F2rho;

    root_paths_top;
    root_paths_bottom;
    unitcircle_path;
    F_paths;
    differentials;

    update() {
        this.hyp = new Hyperelliptic_Curve();
        this.hyp.update();

        // First we define all the path segments

        this.unitcircle_path = new UnitCirclePath(this.hyp);

        this.root_paths_top = new Array(this.hyp.roots.length);
        this.root_paths_bottom = new Array(this.hyp.roots.length);
        for (let i = 0; i < this.hyp.roots.length; i++) {
            this.root_paths_top[i] = new BranchPointPath(this.hyp, i, true);
            this.root_paths_bottom[i] = new BranchPointPath(this.hyp, i, false);
        }

        // We define additionally a degree 2 divisor F, that describes the eigenline bundle. 
        // This will be the future location of a root of the BA function. The rho of this point is a pole of the BA function.
        // The sheets are important.  Start at 1 on the upper sheet. Then go out along a straight line in the x plane to a point in the spectral curve. Call this Fi
        let one = math.complex(1,0);
        let F1x = math.complex({ r: (props.r1/2), phi: (props.phi / 2) })
        let F1_path = new OpenPath(this.hyp, one, F1x, true);
        this.F1 = F1_path.lift_endPoint();
        this.F1rho = rho(this.F1);

        let F2x = math.complex({ r: (props.r2/2), phi: (-props.phi / 2) })
        let F2_path = new OpenPath(this.hyp, one, F2x, true);
        this.F2 = F2_path.lift_endPoint();
        this.F2rho = rho(this.F2);
        this.F_paths = [F1_path, F2_path];

        // Now we define the differentials
        this.differentials = [];
        
        let holo0 = new HoloDiff(this.hyp, 0);
        let holo1 = new HoloDiff(this.hyp, 1);
        this.differentials.push(holo0, holo1);

        let lamtilde0 = new LambdatildeDiff(this.hyp, 0);
        let lamtilde1 = new LambdatildeDiff(this.hyp, 1);
        let laminfty = new LambdaInftyDiff(this.hyp);
        this.differentials.push(lamtilde0, lamtilde1, laminfty);  

        let gamma1 = new GammaDiff(this.hyp, this.F1rho);  
        let gamma2 = new GammaDiff(this.hyp, this.F2rho);  
        this.differentials.push(gamma1, gamma2);

        // This was a test, manually calculated the semi-normalized differentials to the homology
        // const apers0 = this.a_periods(0);
        // const a00 = apers0[0].re;
        // const a01 = apers0[1].re;
        // let aholo0bs = [0, math.complex(a01,-a00), math.complex(a01,a00), 0];
        // let aholo0 = new MeroDiff(this.hyp, aholo0bs);
        // const apers1 = this.a_periods(1);
        // const a10 = apers1[0].re;
        // const a11 = apers1[1].re;
        // let aholo1bs = [0, math.complex(a11,-a10), math.complex(a11,a10), 0];
        // let aholo1 = new MeroDiff(this.hyp, aholo1bs);
        // this.differentials.push(aholo0, aholo1);
    }

    integrate_path_diffs(path, diffs) {
        // This will store the value of the integrals of the diffs along the path.
        let ints = new Array(diffs.length);
        for (let i = 0; i < ints.length; i++) {
            // This is the evals of the ith differential along the path.
            let diff_vals = new Array(path.xs.length);
            for (let j = 0; j < diff_vals.length; j++) {
                diff_vals[j] = diffs[i].eval(path.xs[j], path.ys[j]);
            }
            ints[i] = path.integrate(diff_vals);
        }
        return ints;
    }

    unitcircle_ints() {
        let ints = this.integrate_path_diffs(this.unitcircle_path, this.differentials);

        // The lambda_tilde differentials had an exact form subtracted off to make them manageable.
        // This doesn't change the period calculations, but it does change open path calculations.
        ints[2] = math.chain(ints[2]).add(this.differentials[2].exact(this.unitcircle_path.lift_endPoint())).subtract(this.differentials[2].exact(this.unitcircle_path.lift_startPoint())).done();
        ints[3] = math.chain(ints[3]).add(this.differentials[3].exact(this.unitcircle_path.lift_endPoint())).subtract(this.differentials[3].exact(this.unitcircle_path.lift_startPoint())).done();
        return ints;
    }

    root_paths_top_ints(root_idx) {
        return this.integrate_path_diffs(this.root_paths_top[root_idx], this.differentials);
    }
    root_paths_bottom_ints(root_idx) {
        return this.integrate_path_diffs(this.root_paths_bottom[root_idx], this.differentials);
    }

    root_ints(root_idx) {
        const top = this.root_paths_top_ints(root_idx);
        const bottom = this.root_paths_bottom_ints(root_idx);
        const ints = math.subtract(top,bottom);
        return ints;
    }

    a_periods(root_idx) {
        return math.subtract(this.root_ints(root_idx),this.root_ints(root_idx+2))
    }
    b_periods(root_idx) {
        return math.subtract(this.root_ints(root_idx),this.unitcircle_ints())
    }

    F_ints(idx) {
        let ints = this.integrate_path_diffs(this.F_paths[idx], this.differentials);

        // The lambda_tilde differentials had an exact form subtracted off to make them manageable.
        // This doesn't change the period calculations, but it does change open path calculations.
        ints[2] = math.chain(ints[2]).add(this.differentials[2].exact(this.F_paths[idx].lift_endPoint())).subtract(this.differentials[2].exact(this.F_paths[idx].lift_startPoint())).done();
        ints[3] = math.chain(ints[3]).add(this.differentials[3].exact(this.F_paths[idx].lift_endPoint())).subtract(this.differentials[3].exact(this.F_paths[idx].lift_startPoint())).done();
        // This doesn't change the period calculations, but it does change open path calculations.
        return ints;
    }

    // Solve for Baker-Akhiezer function coefficients
    // Returns [c_holo0, c_holo1, c_lam0, c_lam1, c_laminfty,  c_gamma0, c_gamma1, X_0]
    // where the BA differential is: sum of these coefficients times their respective differentials
    ba_coefficients(idx) {
        // We do a cludge that the index 0, 1 happens to be the value at the sym point
        const X0 = idx;

        // Collect all period integrals
        const a0_pers = this.a_periods(0);
        const a1_pers = this.a_periods(1);
        const b0_pers = this.b_periods(0);
        const b1_pers = this.b_periods(1);

        // Integrals along F-paths (root constraints)
        const f1_ints = this.F_ints(0);
        const f2_ints = this.F_ints(1);

        // mirror sym point path integral (the unit circle integral starting from top sheet)
        const invol_ints = this.unitcircle_ints();

        // Build 7×7 coefficient matrix A and 7×1 RHS vector b
        // Unknowns: [c_holo0, c_holo1, c_lam0, c_lam1, c_laminfty,  c_gamma0, c_gamma1]
        //
        // Constraints:
        // 1. a0 period = 0
        // 2. a1 period = 0
        // 3. b0 period = 0
        // 4. b1 period = 0
        // 5. X at F1 = 0 (integral_to_F1 = - X_0)
        // 6. X at F2 = 0 (integral_to_F2 = - X_0)
        // 7. X at mirror sym point = 1-idx (integral_to_MSP = 1-X0 - X0)
        //
        // Again we use the cludge that idx 0,1 should have 1,0 at the mirror of the sym point.


        const rows = [
            // Row 0: a0 period constraint
            a0_pers,
            // Row 1: a1 period constraint
            a1_pers,
            // Row 2: b0 period constraint
            b0_pers,
            // Row 3: b1 period constraint
            b1_pers,
            // Row 4: F1 root constraint
            f1_ints,
            // Row 5: F2 root constraint
            f2_ints,
            // Row 6: X at mirror sym point
            invol_ints
        ];

        const A = math.matrix(rows);
        const rhs = math.matrix([[0], [0], [0], [0], [-X0], [-X0], [1-(2*idx)]]);

        // Solve the linear system A * x = rhs
        const solution = math.lusolve(A, rhs);

        // Extract the solution vector and convert to array
        return math.flatten(solution).toArray();
    }
}



class Path {
    hyp;
    cc_N;
    path_start = 0;
    path_end = 1;

    xs;
    ys;
    // dx/dt 1/y along this path.
    diff;

    constructor(hyp, cc_N = cc_N_default) {
        this.hyp = hyp;
        this.cc_N = cc_N;
    }

    // computes ys lifted along the discrete path xs (array of points)
    lift_y(startTop) {
        this.ys = new Array(this.xs.length);
        for (let i = 0; i < this.ys.length; i++) {
            let yi = math.sqrt(this.hyp.xa_poly(this.xs[i]));

            if (i == 0) {
                // Initialize the first point. 
                // We have a switch to tell us top or bottom sheet.
                const sheet = this.hyp.isTopSheet(make_pt(this.xs[i], yi));
                const needsFlip = startTop !== sheet;
                yi = needsFlip ? yi.neg() : yi;
            }
            else {
                // Otherwise, i>0, and we must ensure a continuous branch of the sqrt roots. 
                const yi_n = yi.neg();
                const yi_prev = this.ys[i - 1];
                yi = dist_c(yi, yi_prev) > dist_c(yi_n, yi_prev) ? yi_n : yi;
            }

            this.ys[i] = yi;
        }
    }

    lift_startPoint() {
        const x0 = this.xs[0];
        const y0 = this.ys[0];
        return make_pt(x0, y0);
    }
    lift_endPoint() {
        const xn = this.xs[this.xs.length - 1];
        const yn = this.ys[this.xs.length - 1];
        return make_pt(xn, yn);
    }

    integrate(fs) {
        let fs_diff = new Array(fs.length);
        for (let i = 0; i < fs.length; i++) {
            let fs_diff_i = math.multiply(fs[i], this.diff[i]);
            fs_diff_i = Number.isFinite(fs_diff_i.re) && Number.isFinite(fs_diff_i.im) ? fs_diff_i : math.complex(0, 0);
            fs_diff[i] = fs_diff_i;
        }
        return cc_int_sample(fs_diff, this.path_start, this.path_end, this.cc_N)
    }
}


class UnitCirclePath extends Path {
  constructor(hyp, cc_N = cc_N_default) {
    // set the hypcurve and the integration N.
    super(hyp, cc_N);

    this.path_start = 0;
    this.path_end = 2*math.PI;

    // get the points
    let ts = cc_interval_points(this.path_start,this.path_end,this.cc_N);
    
    // next we work out the xs. While we're at it, we do eval the derivative of the path x = exp(i t), since its just x' = i exp(i t) = i x.
    this.xs = new Array(ts.length);
    let dpath = new Array(ts.length);
    for (let i = 0; i < ts.length; i++) {
        const x = math.complex({ r: 1, phi: ts[i] });
        this.xs[i] = x;
        dpath[i] = math.multiply(x,math.complex(0, 1));
    }

    // update the ys.
    this.lift_y(true);

    // Compute dx/y differential
    this.diff = new Array(ts.length);
    for (let i = 0; i < ts.length; i++) {
        this.diff[i] = math.divide(dpath[i], this.ys[i]);
    }
  }
}

class BranchPointPath extends Path {
    constructor(hyp, root_idx, startTop, cc_N = cc_N_default) {
        super(hyp, cc_N);
        this.root_idx = root_idx;
        this.alpha = hyp.roots[root_idx];
        this.one_minus_alpha = math.subtract(math.complex(1, 0), this.alpha);
        this.startTop = startTop;

        this.path_start = 0;
        this.path_end = 1;

        // get the points
        let ts = cc_interval_points(this.path_start,this.path_end,this.cc_N);

        // Now we compute the xs.
        this.xs = new Array(ts.length);
        for (let i = 0; i < ts.length; i++) {
            const x = math.chain(this.one_minus_alpha).multiply((1-ts[i])**2).add(this.alpha).done();
            this.xs[i] = x;
        }

        // update the ys and then the yreds.
        this.lift_y(startTop);
        let y_reds = this.lift_yred();

        // Compute dx/y differential
        this.diff = new Array(ts.length);
        const dscale = math.multiply(-2, this.one_minus_alpha);
        for (let i = 0; i < ts.length; i++) {
            this.diff[i] = math.divide(dscale, y_reds[i]);
        }
    }

    // yred lifted along a discrete path of points. Uses ys to get the correct branch at the start.
    lift_yred() {
        let y_reds = new Array(this.ys.length);
        for (let i = 0; i < this.ys.length; i++) {
            let y_redi = math.sqrt(math.multiply(this.one_minus_alpha,this.hyp.xa_poly_hat(this.root_idx, this.xs[i])));

            // Initialize the first point. We have a switch to tell us top or bottom sheet, labelled by the sign of re part of y/x^2. This expression is purely real and pos on unit circle.
            if (i == 0) {
                // We need to initialize the y_red to match y. 
                // It should be that sqrt((1-α)xa_poly_hat) = sqrt(xa_poly) at x=1, the start of the path. 
                // We don't rely on numerical equality for floats, judge by distance.
                if (dist_c(y_redi, this.ys[0]) > dist_c(y_redi, this.ys[0].neg())) {
                    y_redi = y_redi.neg();
                }
            }
            else {
                // Otherwise, i>0, and we must ensure a continuous branch of the sqrt roots. 
                const y_redi_n = y_redi.neg();
                const y_redi_prev = y_reds[i - 1];
                y_redi = dist_c(y_redi, y_redi_prev) > dist_c(y_redi_n, y_redi_prev) ? y_redi_n : y_redi;
            }

            y_reds[i] = y_redi;
        }
        return y_reds;
    }
}

class OpenPath extends Path {
    constructor(hyp, startPoint, endPoint, startTop, cc_N = cc_N_default) {
        super(hyp, cc_N);
        this.startPoint = startPoint;
        this.endPoint = endPoint;
        this.startTop = true;

        this.path_start = 0;
        this.path_end = 1;

        // get the points
        let ts = cc_interval_points(this.path_start,this.path_end,this.cc_N);

        // Now we compute the xs.
        this.xs = new Array(ts.length);
        for (let i = 0; i < ts.length; i++) {
            const x = math.subtract(math.multiply(1-ts[i], this.startPoint), math.multiply(ts[i], this.endPoint));
            this.xs[i] = x;
        }

        // update the ys
        this.lift_y(startTop);

        // Compute dx/y differential
        this.diff = new Array(ts.length);
        const dscale = math.subtract(this.endPoint, this.startPoint);
        for (let i = 0; i < ts.length; i++) {
            this.diff[i] = math.divide(dscale, this.ys[i]);
        }
    }
}




class Differential {
  hyp;
  constructor(hyp) {
    this.hyp = hyp;
  }
  
  // Subclasses override this
  eval(x, y) {
    throw new Error("Subclass must implement eval(x, y)");
  }
}

class HoloDiff extends Differential {
  idx;
  constructor(hyp, idx) {
    super(hyp);
    this.idx = idx;
  }
  eval(x, y) {
    // TODO: In the future implement the general genus correctly.
    if(this.idx === 0) { 
        return poly_coeff([math.complex(1, 0), math.complex(1, 0)], x);
    } else if(this.idx === 1) {
        return poly_coeff([math.complex(0, 1), math.complex(0, -1)], x);
    }
    throw new Error(`HoloDiff: invalid index ${this.idx}`);
  }
}

// A diff of the form b(x)/x dx/y
class MeroDiff extends Differential {
  bs;
  constructor(hyp, bs) {
    super(hyp);
    this.bs = bs;
  }
  eval(x, y) {
    return math.divide(poly_coeff(this.bs, x),x);
  }
}

// The Lambda diff has a double pole no reside at a branch point.
// We subtract off an exact differential with this same pole (and necessarily some poles at infinity), so be able to integrate through the branch point.
// The exact differential doesn't change its periods, but it does affect open paths.
class LambdatildeDiff extends Differential {
  root_idx;
  constructor(hyp, root_idx) {
    super(hyp);
    this.root_idx = root_idx;
  }
  eval(x, y) {
    const alpha = this.hyp.roots[this.root_idx];
    const factor = math.multiply(alpha, this.hyp.a_poly_prime(alpha));
    return math.chain(this.hyp.a_poly_hat(this.root_idx, x)).subtract(this.hyp.a_poly_prime(x)).divide(factor).done();
  }

  exact(pt) {
    const alpha = this.hyp.roots[this.root_idx];
    const x_diff = math.subtract(pt.x, alpha)
    const factor = math.multiply(alpha, this.hyp.a_poly_prime(alpha));
    return math.chain(pt.y).divide(x_diff).divide(factor).multiply(2).done();
  }
}

// The Lambda diff has a double pole no reside at a branch point.
// This one has it at the branch point x=\infty !
class LambdaInftyDiff extends Differential {
  constructor(hyp) {
    super(hyp);
  }
  eval(x, y) {
    return math.pow(x,2);
  }
}

// The Gamma differential has a double pole no residue at pt_pole
class GammaDiff extends Differential {
  pt_pole;
  constructor(hyp, pt_pole) {
    super(hyp);
    this.pt_pole = pt_pole;
  }
  eval(x, y) {
    const x0 = this.pt_pole.x;
    const y0 = this.pt_pole.y;
    const x_diff = math.subtract(x, x0)
    const term1 = math.divide(math.add(y, y0), math.pow(x_diff, 2))
    const term2 = math.chain(0.5).divide(y0).multiply(this.hyp.a_poly_prime(x0)).divide(x_diff).done();
    return math.add(term1, term2);
  }
}

// var spectralg2 = new Spectral_Curve();
// spectralg2.update();
