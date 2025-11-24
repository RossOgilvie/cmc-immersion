// An implementation of Clenshaw–Curtis quadrature.
// The idea is to integrate via a cos substitution
// \int_{-1}^1 f(x) dx = \int_0^\pi f(cos t) sin t dt
// Suppose we have the Fourier cosine series of f(cos t)
// f(cos t) = 0.5a_0 + sum a_n cos(n t)
// We can integrate analytically
// I = a_0 + sum 2a_{2k} / (1-(2k)^2)
// It remains to approximate the Fourier coefficients 
// a_{2k} = 2/ \pi \int_0^\pi f(cos t) cos(2k t) dt
// which we do just by the trapezoid rule.
// Rearrange so that its a dot product of the function values with certain weights.
// By frequency considerations, if you need to calculate a_k then you need at least k+1 points in the trapezoid rule.

// The number of points sampled is 2*cc_N + 1. Mathematically 2N is the order of the Fourier series
const cc_N_default = 100;

// In the trapezoid rule, and others, you only take half the values at the endpoints.
var cc_corr = function(n,N)
{
  if (n == 0 || n == N) { return 1; }
  else { return 2; }
}

/*
-------------------------------------
WEIGHTS 
-------------------------------------
*/

function computeWeights(cc_N) {
  const cc_pi_div_N = math.PI/cc_N;
  const cc_pi_div_twoN = math.PI/(2*cc_N);
  var cc_weights = new Array(cc_N + 1);
  // calculate the first half of the weights 
  for (n = 0; n <= cc_N; n++) {
    var sum = 0;
    for (k = 0; k <= cc_N; k++) {
      sum += cc_corr(k, cc_N) * (math.cos(k * n * cc_pi_div_N)) / (1 - (2 * k) ** 2)
    }
    cc_weights[n] = cc_corr(n, 2 * cc_N) * sum / (2 * cc_N);
  }
  // weights are symmetric
  for (n = cc_N + 1; n <= 2 * cc_N; n++) {
    cc_weights[n] = cc_weights[2 * cc_N - n];
  }
  // the weights should sum to 0.5

  return cc_weights;
}

function createWeightsManager() {
  // Use Map for insertion order -> simple LRU.
  /** @type {Map<number, number>} */
  const cache = new Map();

  function getFromCacheOrCompute(n) {
    const hit = cache.get(n);
    if (hit) {
      // Return a copy to protect internal cache from mutation
      return hit.slice();
    }
    const computed = computeWeights(n);
    cache.set(n, computed);
    return computed.slice();
  }

  return {
    getWeights(n = cc_N_default) {
      if (!Number.isInteger(n) || n <= 0) throw new Error("n must be a positive integer");
      return getFromCacheOrCompute(n);
    },
    clearCache() { cache.clear(); },
    cacheSize() { return cache.size; }
  };
}

/*
-------------------------------------
INTEGRATOR 
-------------------------------------
*/

const cc_weights = createWeightsManager();

dot = function (a,b) {
  var sum = 0;
  for (i = 0; i < a.length && i < b.length; i++)
  {
    sum += a[i]*b[i]
  }
  return sum
}

// This integrates the function f over the interval [a,b]
function cc_int_ab (f,a,b,cc_N=cc_N_default) {
  fs = new Array(2*cc_N+1);
  ws = cc_weights.getWeights(cc_N);

  let ab_s= 0.5*(b + a)
  let ab_d= 0.5*(b - a)
  let cc_pi_div_N = math.PI/cc_N;
  let cc_pi_div_twoN = math.PI/(2*cc_N);
  
  for (n = 0; n <= 2*cc_N; n++)
  {
    val = f(ab_d*math.cos(n*cc_pi_div_twoN) + ab_s)
    if(Number.isFinite(val)){
    fs[n] = val;
    } else {
    fs[n] = 0;
    }
  }
  return dot(ws,fs) * ab_d
}

// Returns the sample points of the interval [a,b] needed for the cc method
function cc_interval_points(a,b,cc_N=cc_N_default) {
  ts = new Array(2*cc_N+1);

  const ab_s= 0.5*(b + a);
  const ab_d= 0.5*(b - a);
  const cc_pi_div_twoN = math.PI/(2*cc_N);
  
  for (n = 0; n <= 2*cc_N; n++)
  {
    ts[n] = ab_s - ab_d*math.cos(n*cc_pi_div_twoN)
  }
  return ts;
}

// Takes an array of sample values and integrates them together using the cc method weights.
function cc_int_sample(fs,a,b,cc_N=cc_N_default) {
  let total_int = math.complex(0, 0);
  const ws = cc_weights.getWeights(cc_N);

  for (n = 0; n <= 2*cc_N; n++)
  {
    total_int = math.add(total_int, math.multiply(ws[n], fs[n]));
  }

  const ab_d= 0.5*(b - a);
  return math.multiply(total_int, ab_d);
}