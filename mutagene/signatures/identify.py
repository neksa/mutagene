import logging
import math
import urllib

import numpy as np
from scipy.optimize import minimize, nnls
from scipy.spatial.distance import cosine

# from scipy.optimize import basinhopping  # , differential_evolution
# from scipy.optimize import fmin_cobyla
from scipy.stats import entropy

from mutagene.signatures import get_dummy_signatures_lists
from mutagene.signatures.likelihood import (
    DEFAULT_BACKGROUND_SELECTION,
    DEFAULT_METHOD,
    fit_likelihood,
    negative_log_likelihood,
    validate_inputs,
)
from mutagene.signatures.sparse import SPARSE_METHODS, fit_sparse

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
    """cosine distance"""
    return cosine(A.dot(x), b)


def Elastic(x, A, b):
    """elastic net"""
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


def NegLogLik(x, A, b):
    """Return count negative log likelihood without implicit weight normalization."""
    return negative_log_likelihood(x, A, b)


def NegLogLikScaled(x, A, b):
    """Return negative log likelihood per mutation."""
    return NegLogLik(x, A, b) / np.sum(b)


def NegLogLikOld(x, A, b):
    """Preserve the compatibility name with the same count likelihood definition."""
    return NegLogLik(x, A, b)


def count_threshold(x, threshold=10e-6):
    """count the number of above-threshold elements in an array"""
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
    "frobenius": Frobenius,
    "frobeniuszero": FrobeniusZero,
    "cos": Cos,
    "elastic": Elastic,
    "kl": DivergenceKL,
    "divergencekl": DivergenceKL,
    "js": DivergenceJS,
    "divergencejs": DivergenceJS,
    "sparse-aicc": NegLogLik,
    "sparse-bic": NegLogLik,
    "mle-scaled": NegLogLikScaled,
    "mle": NegLogLik,  # MLE maximizes LogLik or minimizes NegLogLik
    "compat": NegLogLikOld,  # MLE with added context-independent signatures (old compatibility mode)
    "aicc": AICc,  # AIC corrected for small samples
    "bic": BIC,  # BIC
    "mlez": NegLogLik,  # MLE with added context-independent signatures (different name Z given for benchmarking)
    "aiccz": AICc,  # AIC corrected for small samples with added context-independent signatures (different name Z given for benchmarking)
    "bicz": BIC,  # BIC with added context-independent signatures (different name Z given for benchmarking)
}


def get_fingerprint_url(a):
    data = {f"s{i}": v for i, v in enumerate(a)}
    return urllib.parse.urlencode(data)


def decompose_mutational_profile_counts(
    profile,
    signatures,
    func=DEFAULT_METHOD,
    others_threshold=0.05,
    global_optimization=None,
    enable_dummy=None,
    diagnostics=None,
    background_selection=DEFAULT_BACKGROUND_SELECTION,
):

    config = {
        "enable_dummy": True,
        "show_dummy": False,
        "show_dummy_aggregate": True,
        "global_optimization": False,
    }

    if global_optimization is not None:
        config["global_optimization"] = global_optimization

    if enable_dummy is not None:
        config["enable_dummy"] = enable_dummy
    else:
        if func.lower().endswith("z"):
            config["enable_dummy"] = True

    W = []
    results = []

    W, signature_names = list(signatures)
    signature_names = list(signature_names)
    original_totals = np.asarray(W, dtype=float).sum(axis=0)
    W, v = validate_inputs(W, np.asarray(profile, dtype=float))
    if len(signature_names) != W.shape[1]:
        raise ValueError("Signature names and matrix columns differ")
    if func.lower() not in IDENTIFY_MIN_FUNCTIONS:
        raise ValueError(f"Unknown fitting method: {func}")
    if diagnostics is not None:
        diagnostics.update(
            requested_method=func,
            reference_normalized=True,
            reference_column_sum_min=float(original_totals.min()),
            reference_column_sum_max=float(original_totals.max()),
        )
    for name in signature_names:
        results.append(
            {
                "accession": 0,
                "id": 0,
                "pid": 0,
                "name": name,
                "annotation": "",
                "score": 0.0,
                "mutations": 0,
            }
        )

    # Background components use the fixed SBS96 channel order.
    if config["enable_dummy"] and W.shape[0] != 96:
        raise ValueError("Context-independent background requires 96 channels")
    if config["enable_dummy"]:
        for i, (name, values) in enumerate(get_dummy_signatures_lists()):
            dummy = "d" + str(i)
            signature_names.append(dummy)
            W = np.append(
                W.T,
                np.array(
                    [
                        values,
                    ]
                ),
                axis=0,
            ).T
            results.append(
                {
                    "accession": 0,
                    "profile": get_fingerprint_url(values),
                    "id": 0,
                    "pid": 0,
                    "name": dummy,
                    "annotation": "Dummy " + name,
                    "score": 0.0,
                    "mutations": 0,
                }
            )
    # print(W)
    # print(signature_names)

    v_sum = v.sum()
    v_freq = v / v_sum
    if np.any((W.sum(axis=1) == 0) & (v > 0)):
        raise ValueError("The signature reference cannot explain an observed channel")

    min_func = IDENTIFY_MIN_FUNCTIONS[func.lower()]
    constraints = {"type": "eq", "fun": lambda x: np.sum(x) - 1.0}
    bounds = [[0.0, 1.0] for _ in range(W.shape[1])]

    # v == counts
    # v_freq == frequency (normalized counts)
    if min_func in (NegLogLik, NegLogLikScaled, NegLogLikOld, AIC, AICc, BIC):
        v_target = v
    else:
        v_target = v_freq

    if config["global_optimization"]:
        # Nothing was ever wired up for this branch, so the exposures stayed
        # unassigned and the call died with NameError several screens further
        # down, long after the choice that caused it.
        raise NotImplementedError(
            "global_optimization is not implemented; omit it to use local minimization"
        )

    if func.lower() in SPARSE_METHODS:
        h = fit_sparse(
            W,
            v,
            len(signatures[1]),
            SPARSE_METHODS[func.lower()],
            diagnostics,
            background_selection,
        )
    elif func.lower() in ("mle", "mle-scaled", "mlez", "compat"):
        h = fit_likelihood(W, v, diagnostics)
    else:
        h0, _ = nnls(W, v_freq)
        h0 = h0 / h0.sum()
        minout = minimize(
            min_func,
            h0,
            args=(W, v_target),
            method="SLSQP",
            bounds=bounds,
            constraints=constraints,
            options={"maxiter": 500},
        )

        if diagnostics is not None:
            diagnostics.update(
                optimizer_success=bool(minout.success),
                optimizer_message=str(minout.message),
                optimizer_iterations=int(minout.nit),
                fallback=None if minout.success else "normalized_nnls",
                estimator=func if minout.success else "normalized_nnls",
            )
        if minout.success:
            logger.debug(f"MINIMIZATION: {minout.message} {minout.nit}")
            h = np.maximum(minout.x, 0)
            h /= h.sum()
            logger.debug(f"MAX LIK {h} {round(-NegLogLik(h, W, v_target), 4)}")
        else:
            logger.warning(
                f"Optimization fails: {minout.message}. Using normalized NNLS estimates."
            )
            # Minimization did not converge
            # Use our initial guess, but normalize it:
            h = h0.ravel() / h0.sum()

    ##############################################################################
    N_mutations = int(math.ceil(v.sum()))

    for i in range(len(results)):
        results[i]["score"] = h[i]
        results[i]["mutations"] = round(N_mutations * h[i])

    # remove dummy signatures if enabled
    if config["enable_dummy"]:
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
    for r in sorted(results, key=lambda item: item["score"], reverse=True):
        if others_threshold > 0.0 and round(r["score"], 2) <= others_threshold:
            below_threshold.append(r)
        else:
            j += 1
            r["id"] = j
            above_threshold.append(r)
    results = above_threshold

    # if len(below_threshold) == 1:
    #     # print("ONLY ONE")
    #     results.append(below_threshold[0])
    #     below_threshold = []

    # sum up other signatures
    other_signatures = 0.0
    for r in below_threshold:
        other_signatures += r["score"]

    if round(other_signatures, 2) > 0.0:
        results.append(
            {
                "accession": 0,
                "id": 100,
                "pid": 0,
                "name": "Other signatures",
                "annotation": "Signatures with individual contrubution &le; 0.05",
                "profile": "",
                "score": other_signatures,
                "mutations": round(N_mutations * other_signatures),
            }
        )
        for j, r in enumerate(below_threshold):
            if round(r["score"], 2) > 0.0:
                r["id"] = 100 + j + 1
                r["pid"] = 100
                results.append(r)

    ll = -NegLogLik(h, W, v)
    results.append(
        {
            "accession": 0,
            "id": 201,
            "pid": 0,
            "name": "LogLik",
            "annotation": "",
            "profile": "",
            "mutations": "",
            "score": ll,
        }
    )

    frobenius = Frobenius(h, W, v_freq)
    results.append(
        {
            "accession": 0,
            "id": 202,
            "pid": 0,
            "name": "Frobenius",
            "annotation": "",
            "profile": "",
            "mutations": "",
            "score": frobenius,
        }
    )

    frobeniuszero = FrobeniusZero(h, W, v_freq)
    results.append(
        {
            "accession": 0,
            "id": 203,
            "pid": 0,
            "name": "FrobeniusZero",
            "annotation": "",
            "profile": "",
            "mutations": "",
            "score": frobeniuszero,
        }
    )

    divergencejs = DivergenceJS(h, W, v_freq)
    results.append(
        {
            "accession": 0,
            "id": 204,
            "pid": 0,
            "name": "DivergenceJS",
            "annotation": "",
            "profile": "",
            "mutations": "",
            "score": divergencejs,
        }
    )

    divergencekl = DivergenceKL(h, W, v_freq)
    results.append(
        {
            "accession": 0,
            "id": 205,
            "pid": 0,
            "name": "DivergenceKL",
            "annotation": "",
            "profile": "",
            "mutations": "",
            "score": divergencekl,
        }
    )

    summary = []
    summary.append(
        {
            "accession": 0,
            "name": "Query profile",
            "annotation": "Query profile derived from mutations uploaded by the user",
            "profile": get_fingerprint_url(profile),
        }
    )

    summary.append(
        {
            "accession": 0,
            "name": "Reconstructed profile",
            "annotation": "Profile represented as a combination of signatures in the table below",
            "profile": get_fingerprint_url(reconstructed_profile),
        }
    )

    return h, summary, results
