"""Select signature subsets with constrained likelihood refits and explicit complexity costs."""

import numpy as np

from mutagene.signatures.likelihood import fit_likelihood, validate_inputs

SPARSE_METHODS = {"sparse-aicc": "aicc", "sparse-bic": "bic"}


def information_criterion(nll, dimensions, burden, criterion):
    """Count simplex dimensions; retain nuisance dimensions when null signatures remain available."""
    if criterion == "bic":
        return 2 * nll + dimensions * np.log(burden)
    if criterion == "aicc":
        if burden <= dimensions + 1:
            return np.inf
        return (
            2 * nll + 2 * dimensions + 2 * dimensions * (dimensions + 1) / (burden - dimensions - 1)
        )
    raise ValueError("Sparse selection supports aicc or bic")


def fit_sparse(
    matrix, counts, known_count, criterion="bic", diagnostics=None, background_selection="fixed"
):
    """Select known signatures and apply the requested background selection policy.

    Backward elimination and forward additions use exact constrained likelihood refits.
    This local search does not guarantee a globally optimal subset.
    Candidate optimization failures stop selection rather than silently changing its objective.
    """
    matrix, counts = validate_inputs(matrix, counts)
    if not 0 < known_count <= matrix.shape[1]:
        raise ValueError("Provide a valid known-signature count")
    information_criterion(0, 0, counts.sum(), criterion)
    burden = float(counts.sum())
    observed = counts > 0
    if background_selection not in ("fixed", "joint", "block"):
        raise ValueError("Background selection supports fixed, joint, or block")
    if background_selection == "block":
        return _fit_background_block(matrix, counts, known_count, criterion, diagnostics)
    nuisance = tuple(range(known_count, matrix.shape[1])) if background_selection == "fixed" else ()
    optional_count = known_count if background_selection == "fixed" else matrix.shape[1]
    cache = {}
    iterations = 0
    retries = 0
    refit_gaps = []

    def fit(subset):
        nonlocal iterations, retries
        subset = tuple(sorted(subset))
        if subset in cache:
            return cache[subset]
        columns = subset + nuisance
        if not columns:
            return None
        design = matrix[observed][:, columns]
        if np.any(design.sum(axis=1) == 0):
            return None
        info = {}
        weights = fit_likelihood(matrix[:, columns], counts, info)
        iterations += info["optimizer_iterations"]
        retries += info["optimizer_retries"]
        refit_gaps.append(info["optimality_gap"])
        nll = info["nll_per_mutation"] * burden
        dimensions = len(columns) - 1
        value = information_criterion(nll, dimensions, burden, criterion)
        expanded = np.zeros(matrix.shape[1])
        expanded[list(columns)] = weights
        cache[subset] = (expanded, value, nll)
        return cache[subset]

    selected = tuple(range(optional_count))
    current = fit(selected)
    if current is None:
        raise ValueError("The signature reference cannot explain an observed channel")
    # Remove numerical zeros before searching subsets with information criteria.
    selected = tuple(i for i in selected if current[0][i] > 1e-7)
    current = fit(selected)
    if current is None or not np.isfinite(current[1]):
        raise ValueError("The mutation burden is too small for this information criterion")
    trace = [
        {
            "action": "start",
            "known_count": sum(i < known_count for i in selected),
            "background_count": len(nuisance) + sum(i >= known_count for i in selected),
            "criterion": current[1],
        }
    ]
    while True:
        best = current
        best_subset = selected
        for remove in selected:
            candidate = tuple(i for i in selected if i != remove)
            result = fit(candidate)
            if result is not None and result[1] < best[1] - 1e-7:
                best, best_subset = result, candidate
        if best_subset != selected:
            selected, current = best_subset, best
            trace.append(
                {
                    "action": "remove",
                    "known_count": sum(i < known_count for i in selected),
                    "background_count": len(nuisance) + sum(i >= known_count for i in selected),
                    "criterion": current[1],
                }
            )
            continue
        for add in range(optional_count):
            if add in selected:
                continue
            candidate = tuple(sorted((*selected, add)))
            result = fit(candidate)
            if result is not None and result[1] < best[1] - 1e-7:
                best, best_subset = result, candidate
        if best_subset == selected:
            break
        selected, current = best_subset, best
        trace.append(
            {
                "action": "add",
                "known_count": sum(i < known_count for i in selected),
                "background_count": len(nuisance) + sum(i >= known_count for i in selected),
                "criterion": current[1],
            }
        )
    if diagnostics is not None:
        diagnostics.update(
            optimizer_success=True,
            optimizer_message="Subset search reaches a local optimum",
            optimizer_iterations=int(iterations),
            fallback=None,
            estimator="sparse_multinomial_mle",
            maximum_refit_optimality_gap=float(max(refit_gaps)),
            subset_refits=len(cache),
            optimizer_retries=retries,
            selected_known_indices=[i for i in selected if i < known_count],
            selected_background_indices=[
                i - known_count for i in (*selected, *nuisance) if i >= known_count
            ],
            model_dimensions=len(selected) + len(nuisance) - 1,
            background_selection=background_selection,
            selection_criterion=criterion,
            criterion_value=float(current[1]),
            null_policy=(
                "Retain all null components as nuisance parameters"
                if background_selection == "fixed"
                else "Select null components with known signatures"
            ),
            selection_trace=trace,
        )
    return current[0]


def _fit_background_block(matrix, counts, known_count, criterion, diagnostics):
    """Compare complete information criteria with and without the background block."""
    with_info = {}
    with_weights = fit_sparse(matrix, counts, known_count, criterion, with_info, "fixed")
    candidates = [(with_weights, with_info)]
    # A missing support channel makes the model without background infeasible.
    if np.all(matrix[counts > 0, :known_count].sum(axis=1) > 0):
        without_info = {}
        without_weights = fit_sparse(
            matrix[:, :known_count], counts, known_count, criterion, without_info, "fixed"
        )
        expanded = np.zeros(matrix.shape[1])
        expanded[:known_count] = without_weights
        candidates.append((expanded, without_info))
    weights, chosen = min(candidates, key=lambda candidate: candidate[1]["criterion_value"])
    if diagnostics is not None:
        diagnostics.update(chosen)
        diagnostics.update(
            background_selection="block",
            null_policy="Compare complete models with and without the background block",
            block_criterion_values=[info["criterion_value"] for _, info in candidates],
            optimizer_iterations=sum(info["optimizer_iterations"] for _, info in candidates),
            optimizer_retries=sum(info["optimizer_retries"] for _, info in candidates),
            subset_refits=sum(info["subset_refits"] for _, info in candidates),
            maximum_refit_optimality_gap=max(
                info["maximum_refit_optimality_gap"] for _, info in candidates
            ),
        )
    return weights
