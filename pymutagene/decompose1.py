def decompose_mutational_profile(
        profile, signature_type, func="Frobenius",
        debug=False, others_threshold=0.05):

    assert False, "This function is obsolete. Use decompose_mutational_profile_counts()"

    if isinstance(signature_type, str):
        W = []
        signature_list = db.query(Signature).filter(Signature.signature_type == signature_type).order_by(Signature.id)
        results = []
        for s in signature_list:
            # print(s.name)
            W.append(get_signature_values_list(s))
            results.append({
                'accession': s.id,
                'name': s.name,
                'annotation': s.annotation,
                'score': 0.0
            })

        W = np.array(W).T
    else:
        W = signature_type

    v = np.array([profile]).ravel()

    # print(W)
    # print(v)

    # Initial guess with NNLS:
    # DATA, residuals = nnls(W.T, X.T.ravel())
    h0, rnorm = nnls(W, v)
    if h0.sum() > 1.0:
        h0 = h0.ravel() / h0.sum()

    min_funcs = {
        'frobenius': Frobenius,
        'kl': DivergenceKL,
        'divergencekl': DivergenceKL,
        'js': DivergenceJS,
        'divergencejs': DivergenceJS
    }
    min_func = min_funcs.get(func.lower(), Frobenius)

    if debug:
        np.set_printoptions(precision=4)
        print("--------------------------------------\n")
        print("Signature", signature_type)
        print("NNLS", h0)
        print("FRO", round(Frobenius(h0, W, v), 4), "DIV", round(DivergenceKL(h0, W, v), 4))

    # Define constraints and bounds
    # constraints = {'type': 'eq', 'fun': lambda x: np.sum(x) - 1}
    constraints = {'type': 'ineq', 'fun': lambda x: 1.0 - np.sum(x)}
    # bounds = [[0., None],[0., None],[0., None]]
    bounds = [[0., 1.0] for _ in range(h0.shape[0])]

    # minout = basinhopping(min_func, h0, args=(W, v), method='SLSQP', bounds=bounds, constraints=constraints)
    # if minout.success:
    #     if debug:
    #         print("MINIMIZATION:", minout.message, minout.nit)
    #     h = minout.x
    #     if debug:
    #         print("MIN FRO", h, round(min_func(h, W, v), 4))
    #         print("FRO", round(Frobenius(h, W, v), 4), "DIV", round(DivergenceKL(h, W, v), 4))
    # else:
    #     if debug:
    #         print("MINIMIZATION FAILED:", minout.message, minout.nit)
    #     # Minimization did not converge
    #     # Use our initial guess, but normalize it:
    #     h = h0.ravel() / h0.sum()

    # Call minimization subject to these values
    minout = minimize(min_func, h0, args=(W, v), method='SLSQP', bounds=bounds, constraints=constraints)
    if minout.success:
        if debug:
            print("MINIMIZATION:", minout.message, minout.nit)
        h = minout.x
        if debug:
            print("MIN FRO", h, round(min_func(h, W, v), 4))
            print("FRO", round(Frobenius(h, W, v), 4), "DIV", round(DivergenceKL(h, W, v), 4))
    else:
        if debug:
            print("MINIMIZATION FAILED:", minout.message, minout.nit)
        # Minimization did not converge
        # Use our initial guess, but normalize it:
        h = h0.ravel() / h0.sum()

    """
    # Call minimisation subject to these values
    minout = minimize(Frobenius, h0, args=(W, v), method='SLSQP', bounds=bounds, constraints=constraints)
    if minout.success:
        # print("MINIMIZATION:", minout.message, minout.nit)
        h = minout.x
        print("MIN FRO", h)
        print("FRO", round(Frobenius(h, W, v), 4), "DIV", round(Divergence(h, W, v), 4))
    else:
        # Minimization did not converge
        # Use our initial guess, but normalize it:
        h = h0.ravel() / h0.sum()

    # Call minimisation subject to these values
    minout = minimize(Divergence, h0, args=(W, v), method='SLSQP', bounds=bounds, constraints=constraints)
    if minout.success:
        # print("MINIMIZATION:", minout.message, minout.nit)
        h = minout.x
        print("MIN DIV", h)
        print("FRO", round(Frobenius(h, W, v), 4), "DIV", round(Divergence(h, W, v), 4))
    else:
        # Minimization did not converge
        # Use our initial guess, but normalize it:
        h = h0.ravel() / h0.sum()
    """

    for i in range(len(results)):
        results[i]['score'] = h[i]

    # residuals = min_func(h, W, v)
    residuals = Frobenius(h, W, v)
    reconstructed_profile = W.dot(h)

    below_threshold = []
    above_threshold = []
    # import operator
    # for r in results:
    # print(results)
    for r in sorted(results, key=lambda item: item['score'], reverse=True):
        if round(r['score'], 2) <= others_threshold:
            below_threshold.append(r)
        else:
            above_threshold.append(r)
    results = above_threshold

    if len(below_threshold) == 1:
        # print("ONLY ONE")
        results.append(below_threshold[0])
        below_threshold = []

    # sum up other signatures
    other_signatures = 0.0
    for r in below_threshold:
        other_signatures += r['score']

    if other_signatures > 0.0:
        results.append({
                'accession': 0,
                'name': 'Other signatures',
                'annotation': 'Signatures with individual contrubution &le; 0.05',
                'profile': '',
                'score': other_signatures,
            })

    results.append({
            'accession': 0,
            'name': 'Residuals',
            'annotation': '',
            'profile': '',
            'score': residuals,
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
            'annotation': 'Profile reconstructed as a linear combination of signatures in the table below',
            'profile': get_fingerprint_url(reconstructed_profile)
        })

    # print(summary)
    # print(results)
    return h, summary, results
