import urllib
import math

import numpy as np

from scipy.optimize import minimize
from scipy.optimize import nnls
from scipy.spatial.distance import cosine
# from scipy.optimize import basinhopping  # , differential_evolution
# from scipy.optimize import fmin_cobyla
from scipy.stats import entropy

from mutagene.signatures import get_dummy_signatures_lists

import logging
logger = logging.getLogger(__name__)


def multi_kl(p, q):
    """Kullback-Liebler divergence from multinomial p to multinomial q,
    expressed in bits."""
    # Clip before taking logarithm to avoid NaNs (but still exclude
    # zero-probability mixtures from the calculation)
    return np.sum(p * (np.log2(p.clip(1e-10, 1)) - np.log2(q.clip(1e-10, 1))))


def multi_js(p, q):
    """Jensen-Shannon divergence (symmetric) between two multinomials,
    expressed in bits."""
    # D_{JS}(P\|Q) = (D_{KL}(P\|Q) + D_{KL}(Q\|P)) / 2
    m = 0.5 * (p + q)
    return 0.5 * (multi_kl(p, m) + multi_kl(q, m))
    # return 0.5 * (
    #     (q * (np.log2(q.clip(1e-10, 1)) - np.log2(p.clip(1e-10, 1)))).sum(0) +
    #     (p * (np.log2(p.clip(1e-10, 1)) - np.log2(q.clip(1e-10, 1)))).sum(0))


# Define minimization function
def DivergenceKL(x, A, b):
    b = b + 1e-17
    v = A.dot(x) + 1e-17
    return entropy(v, b)


# Define minimization function
def DivergenceJS(x, A, b):
    b = b + 1e-17
    v = A.dot(x) + 1e-17
    m = (v + b) / 2.0
    return 0.5 * (entropy(v, m) + entropy(b, m))


def Cos(x, A, b):
    """ cosine distance """
    return cosine(A.dot(x), b)


def Elastic(x, A, b):
    """ elastic net """
    alpha = 1.0
    l1_ratio = 0.5

    k = np.count_nonzero(x)  # signatures
    n = b.sum()  # samples

    objective = np.linalg.norm(A.dot(x) - b) / (0.5 * n)
    # l1 reg
    objective += alpha * l1_ratio * k**2
    # l2 reg
    # objective += alpha * (1.0 - l1_ratio) * ||w||^2_2
    return objective


# Define minimization function
def Frobenius(x, A, b):
    return np.linalg.norm(A.dot(x) - b)


# Define minimization function
def FrobeniusZero(x, A, b):
    b_hat = A.dot(x)
    b_hat[b == 0.0] = 0.0
    b_hat /= b_hat.sum()
    error = b - b_hat
    return np.linalg.norm(error)


# def DerNegLogLik(x, A, b):    
#     print(b / A.dot(x))
#     DLL = - np.ma.divide(b, A.dot(x)).filled(0)
#     print(x, DLL)
#     return DLL


def NegLogLik(x, A, b):
    """
    Log Likelihood of a mixture of multinomials
    x = parameters of the mixture
    A = signature x 96 channels
    b = the observed profile (includes counts for 96 channels)
    """

    # if x.sum() - 0.00011 > 1.0:
    #     return 1e6

    if x.sum() > 1.0:
        x /= x.sum()

    LL = np.sum(b * np.ma.log(A.dot(x)).filled(0.0))

    # print("NegLL\t{}\t{}\t{}".format(len(x), "\t".join(map(lambda a: "{:.2f}".format(a), x)), -LL))
    return -LL


def NegLogLikOld(x, A, b):
    """
    Log Likelihood of a mixture of multinomials
    x = parameters of the mixture
    A = signature x 96 channels
    b = the observed profile (includes counts for 96 channels)
    """

    # if x.sum() - 0.00011 > 1.0:
    #     return 1e6

    # if x.sum() > 1.0:
    #     x /= x.sum()

    # Compatibility mode with the website:
    if np.sum(x) > 1.0:
        return 1000

    LL = np.sum(b * np.ma.log(A.dot(x)).filled(0))

    # print("LOG\t{}\t{}\t{}".format(len(x), "\t".join(map(lambda a: "{:.2f}".format(a), x)), LL))
    return -LL


def count_threshold(x, threshold=10e-6):
    """ count the number of above-threshold elements in an array """
    return np.sum(x > threshold)


def AIC(x, A, b):
    """
    AIC https://en.wikipedia.org/wiki/Akaike_information_criterion
    """
    k = count_threshold(x)
    # n = b.sum()
    return 2 * NegLogLik(x, A, b) + 2 * k


def AICc(x, A, b):
    """
    AIC with small sample correction
    """
    k = count_threshold(x)
    n = b.sum()
    if n - k - 1 <= 0:
        return AIC(x, A, b)
    else:
        return AIC(x, A, b) + 2 * k * (k + 1) / (n - k - 1)


def BIC(x, A, b):
    """
    BIC https://en.wikipedia.org/wiki/Bayesian_information_criterion
    """
    k = count_threshold(x)
    n = b.sum()
    # print(2 * NegLogLik(x, A, b) + k * np.log(n))
    return 2 * NegLogLik(x, A, b) + k * np.log(n)


IDENTIFY_MIN_FUNCTIONS = {
    'frobenius': Frobenius,
    'frobeniuszero': FrobeniusZero,
    'cos': Cos,
    'elastic': Elastic,
    'kl': DivergenceKL,
    'divergencekl': DivergenceKL,
    'js': DivergenceJS,
    'divergencejs': DivergenceJS,
    'mle': NegLogLik,  # MLE maximizes LogLik or minimizes NegLogLik
    'compat': NegLogLikOld,  # MLE with added context-independent signatures (old compatibility mode)
    'aicc': AICc,  # AIC corrected for small samples
    'bic': BIC,  # BIC
    'mlez': NegLogLik,  # MLE with added context-independent signatures (different name Z given for benchmarking)
    'aiccz': AICc,  # AIC corrected for small samples with added context-independent signatures (different name Z given for benchmarking)
    'bicz': BIC,  # BIC with added context-independent signatures (different name Z given for benchmarking)
}


def get_fingerprint_url(a):
    data = {"s{}".format(i): v for i, v in enumerate(a)}
    return urllib.parse.urlencode(data)


def decompose_mutational_profile_counts(profile, signatures, func="Frobenius", others_threshold=0.05, global_optimization=None, enable_dummy=None):

    config = {
        'enable_dummy': True,
        'show_dummy': False,
        'show_dummy_aggregate': True,
        'global_optimization': False,
    }

    if global_optimization is not None:
        config['global_optimization'] = global_optimization

    if enable_dummy is not None:
        config['enable_dummy'] = enable_dummy
    else:
        if func.lower().endswith('z'):
            config['enable_dummy'] = True

    W = []
    results = []

    W, signature_names = list(signatures)
    signature_names = signature_names[:]  # make a copy
    for name in signature_names:
        results.append({
            'accession': 0,
            'id': 0,
            'pid': 0,
            'name': name,
            'annotation': '',
            'score': 0.0,
            'mutations': 0
        })

    # add dummy signatures
    if config['enable_dummy']:
        for i, (name, values) in enumerate(get_dummy_signatures_lists()):
            dummy = "d" + str(i)
            signature_names.append(dummy)
            W = np.append(W.T, np.array([values, ]), axis=0).T
            results.append({
                'accession': 0,
                'profile': get_fingerprint_url(values),
                'id': 0,
                'pid': 0,
                'name': dummy,
                'annotation': 'Dummy ' + name,
                'score': 0.0,
                'mutations': 0
            })
    # print(W)
    # print(signature_names)

    v = np.array([profile]).ravel()
    v_freq = v / v.sum()

    # Initial guess with NNLS:
    # DATA, residuals = nnls(W.T, X.T.ravel())
    h0, rnorm = nnls(W, v_freq)
    if h0.sum() > 1.0:
        h0 = h0.ravel() / h0.sum()

    # not sure if we need this failsafe option:
    # if np.isnan(h0.sum()):
    #     h0 = np.ones(h0.shape[0]) / h0.shape[0]

    logger.debug("h0 {}".format(h0))

    min_func = IDENTIFY_MIN_FUNCTIONS.get(func.lower(), Frobenius)

    # if debug:
    #     np.set_printoptions(precision=4)
    #     print("--------------------------------------\n")
    #     print("Signature", signatures)
    #     print("NNLS", h0)
    #     print("FRO", round(Frobenius(h0, W, v), 4), "DIV", round(DivergenceKL(h0, W, v), 4))

    # Define constraints and bounds
    # constraints = {'type': 'eq', 'fun': lambda x: np.sum(x) - 1}
    constraints = {'type': 'ineq', 'fun': lambda x: 1.0 - np.sum(x)}
    # bounds = [[0., None],[0., None],[0., None]]

    # bounds = [[0., 1.0] for _ in range(h0.shape[0])]
    bounds = [[0.0, 1.0] for _ in range(h0.shape[0])]

    # v == counts
    # v_freq == frequency (normalized counts)
    if min_func in (NegLogLik, AIC, AICc, BIC):
        v_target = v
    else:
        v_target = v_freq

    def print_fun(x, f, accepted):
        print("{} at minimum {} accepted {}".format(x, f, accepted))

    class RandomDisplacementBounds(object):
        """random displacement with bounds"""

        def __init__(self, xmin=0.0, xmax=1.0, stepsize=0.3):
            self.xmin = xmin
            self.xmax = xmax
            self.stepsize = stepsize

        def __call__(self, x):
            """take a random step but ensure the new position is within the bounds"""
            if np.any(np.isnan(x)):
                return x
            while True:
                xnew = x + np.random.uniform(-self.stepsize, self.stepsize, x.shape)
                xnew = np.where(xnew < 0.0, 0.0, xnew)
                if xnew.sum() > 1.0:
                    xnew += 1.0 - xnew.sum()
                xnew = np.abs(xnew)
                if np.all(xnew <= self.xmax) and np.all(xnew >= self.xmin):
                    break
            return xnew

    # class RandomDisplacementBounds(object):
    #     """random displacement with bounds"""
    #     def __init__(self, xmin, xmax, stepsize=0.1):
    #         self.xmin = xmin
    #         self.xmax = xmax
    #         self.stepsize = stepsize

    #     def __call__(self, x):
    #         """take a random step but ensure the new position is within the bounds"""
    #         while True:
    #             # xnew = x + np.random.uniform(-self.stepsize, self.stepsize, x.shape)
    #             inew = np.random.randint(0, x.shape[0])
    #             xnew = x.copy()
    #             xnew[inew] += np.random.uniform(-self.stepsize, self.stepsize)
    #             # print("it", np.sum(xnew), xnew)
    #             if np.all(xnew <= self.xmax) and np.all(xnew >= self.xmin) and np.sum(xnew) <= 1.0:
    #                 break
    #         return xnew

    # define the new step taking routine and pass it to basinhopping
    take_step = RandomDisplacementBounds()

    def accept_test(f_new, x_new, f_old, x_old):
        if np.any(np.isnan(x_new)):
            return False
        if x_new.sum() > 1.0:
            return False
        if np.any(x_new > 1.0) or np.any(x_new < 0.0) or np.any(np.isnan(x_new)):
            return False
        return True

    def get_constraint(i, upper=True):
        def a(x):
            return 1.0 - x[i]

        def b(x):
            return x[i] - 0.0

        return a if upper else b

    cobyla_constraints = [
        {'type': 'ineq', 'fun': lambda x: 1.0 - np.sum(x)}  # total upper
    ]
    for i in range(h0.shape[0]):
        cobyla_constraints.append(
            {'type': 'ineq', 'fun': get_constraint(i, True)}  # , 'jac': get_constraint(i, True)}  # upper
        )
        cobyla_constraints.append(
            {'type': 'ineq', 'fun': get_constraint(i, False)}  # , 'jac': get_constraint(i, False)}  # lower
        )

    ##############################################################################
    if config['global_optimization']:
        # L-BFGS-B COBYLA SLSQP
        """
        minout = basinhopping(NegLogLik, h0, minimizer_kwargs={  # 'jac': DerNegLogLik,
            'args': (W, v_target), 'method': 'SLSQP',
            'constraints': cobyla_constraints, 'options': {'disp': False, 'ftol': 0.00001, 'eps': 0.0000001}}, niter=100, T=0.001, take_step=take_step, accept_test=accept_test, callback=print_fun)  #, take_step=take_step , callback=print_fun)

        if debug:
            print(minout, type(minout))
            print(minout.message)
            # print("Status: # failures", minout.res.minimization_failures)

        if np.any(np.isnan(minout.x)):
            print("FAILED!!!!!!!!!!!!!!!!")
            if debug:
                print("MINIMIZATION FAILED:", minout.message, minout.nit)
            # fallback decomposition:
            h = h0.ravel() / h0.sum()
            # h = h0.ravel() * 0.0
        else:
            h = minout.x
            if debug:
                print("MAX LIK", h, round(-NegLogLik(h, W, v_target), 4))
                print("MIN FRO", h, round(min_func(h, W, v_target), 4))
                print("FRO", round(Frobenius(h, W, v_target), 4), "DIV", round(DivergenceKL(h, W, v), 4))
        """
        # Bayesian Optimization
        from bayes_opt import BayesianOptimization
        from functools import partial

        optimizer = BayesianOptimization(
            f=partial(min_func, A=W, b=v_target),
            # pbounds=pbounds,
            verbose=2,  # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
            random_state=1,
        )

    ##############################################################################
    # debug = True
    if not config['global_optimization']:
        minout = minimize(
            min_func, h0, args=(W, v_target),
            method='SLSQP',
            bounds=bounds, constraints=constraints,
            options={'maxiter': 500}
        )

        if minout.success:
            logger.debug("MINIMIZATION: {} {}".format(minout.message, minout.nit))
            h = minout.x
            logger.debug("MAX LIK {} {}".format(h, round(-NegLogLik(h, W, v_target), 4)))
            # print("LIK", round(-NegLogLik(h, W, v), 4), "DIV", round(DivergenceKL(h, W, v), 4))
        else:
            logger.debug("MINIMIZATION FAILED:{} {}".format(minout.message, minout.nit))
            # Minimization did not converge
            # Use our initial guess, but normalize it:
            h = h0.ravel() / h0.sum()

    ##############################################################################
    N_mutations = int(math.ceil(v.sum()))

    for i in range(len(results)):
        results[i]['score'] = h[i]
        results[i]['mutations'] = round(N_mutations * h[i])

    # remove dummy signatures if enabled
    if config['enable_dummy']:
        results = results[:-6]
        signature_names = signature_names[:-6]

    # residuals = min_func(h, W, v)
    reconstructed_profile = W.dot(h)
    ##############################################################################

    below_threshold = []
    above_threshold = []
    # import operator
    # for r in results:
    # print(results)
    j = 0
    for r in sorted(results, key=lambda item: item['score'], reverse=True):
        if others_threshold > 0.0 and round(r['score'], 2) <= others_threshold:
            below_threshold.append(r)
        else:
            j += 1
            r['id'] = j
            above_threshold.append(r)
    results = above_threshold

    # if len(below_threshold) == 1:
    #     # print("ONLY ONE")
    #     results.append(below_threshold[0])
    #     below_threshold = []

    # sum up other signatures
    other_signatures = 0.0
    for r in below_threshold:
        other_signatures += r['score']

    if round(other_signatures, 2) > 0.0:
        results.append({
            'accession': 0,
            'id': 100,
            'pid': 0,
            'name': 'Other signatures',
            'annotation': 'Signatures with individual contrubution &le; 0.05',
            'profile': '',
            'score': other_signatures,
            'mutations': round(N_mutations * other_signatures)
        })
        for j, r in enumerate(below_threshold):
            if round(r['score'], 2) > 0.0:
                r['id'] = 100 + j + 1
                r['pid'] = 100
                results.append(r)

    ll = -NegLogLik(h, W, v_target)
    results.append({
        'accession': 0,
        'id': 201,
        'pid': 0,
        'name': 'LogLik',
        'annotation': '',
        'profile': '',
        'mutations': '',
        'score': ll,
    })

    frobenius = Frobenius(h, W, v_freq)
    results.append({
        'accession': 0,
        'id': 202,
        'pid': 0,
        'name': 'Frobenius',
        'annotation': '',
        'profile': '',
        'mutations': '',
        'score': frobenius,
    })

    frobeniuszero = FrobeniusZero(h, W, v_freq)
    results.append({
        'accession': 0,
        'id': 203,
        'pid': 0,
        'name': 'FrobeniusZero',
        'annotation': '',
        'profile': '',
        'mutations': '',
        'score': frobeniuszero,
    })

    divergencejs = DivergenceJS(h, W, v_freq)
    results.append({
        'accession': 0,
        'id': 204,
        'pid': 0,
        'name': 'DivergenceJS',
        'annotation': '',
        'profile': '',
        'mutations': '',
        'score': divergencejs,
    })

    divergencekl = DivergenceKL(h, W, v_freq)
    results.append({
        'accession': 0,
        'id': 205,
        'pid': 0,
        'name': 'DivergenceKL',
        'annotation': '',
        'profile': '',
        'mutations': '',
        'score': divergencekl,
    })

    summary = []
    summary.append({
        'accession': 0,
        'name': 'Query profile',
        'annotation': 'Query profile derived from mutations uploaded by the user',
        'profile': get_fingerprint_url(profile)
    })

    summary.append({
        'accession': 0,
        'name': 'Reconstructed profile',
        'annotation': 'Profile represented as a combination of signatures in the table below',
        'profile': get_fingerprint_url(reconstructed_profile)
    })

    return h, summary, results
