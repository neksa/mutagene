

def analyze_all(method='L'):
    "Run all analysis"

    def format_numbers(d):
        d['score'] = "{:6.4g}".format((d['score']))
        return d

    samples = []

    q = db.query(Signature).\
        filter(Signature.active == 1).\
        filter(Signature.mut_type == 'A').\
        filter(Signature.normalization == 'MUT').\
        filter(Signature.signature_type == 'CT')
    for s in q:
        print("CT: ", s.name)
        if s.id == 34904:
            print("SKIP")
            continue
        sample_q = db.query(Signature).\
            join(Sample).\
            options(subqueryload(Signature.data)).\
            filter(Signature.cancer_type == s.cancer_type).\
            filter(Signature.primary_site == s.primary_site).\
            filter(Signature.signature_type == "SPL").\
            filter(Signature.mut_type == "A").\
            filter(Signature.normalization == "MUT").\
            order_by(Signature.id)

        for sample in sample_q:
            # print(sample.name)
            samples.append(sample)

    prefix = "/Users/gonceare/projects/RA"
    for sample in samples:
        # if sample.id > 51:
        #     break

        # if sample.id != 1402:
        #     continue

        TCGA_ID = str(sample.id)
        mutational_profile = get_signature_values_list(sample)
        mutational_profile_counts = get_signature_values_list(sample, counts=True)
        oname = prefix + "/ALL.profiles/{}.profile".format(TCGA_ID)
        print("ONAME:", oname)
        with open(oname, 'w') as o:
            query_formatted = format_profile(mutational_profile)
            o.write(query_formatted)

        oname = prefix + "/ALL.profiles/{}.counts".format(TCGA_ID)
        with open(oname, 'w') as o:
            query_formatted = format_profile(mutational_profile_counts, counts=True)
            o.write(query_formatted)

        oname = prefix + "/ALL.profiles/{}.mutations".format(TCGA_ID)
        # print("ONAME:", oname)
        with open(oname, 'w') as o:
            o.write("{}".format(sample.n_mutations))

        # classification_method = "rf"
        # cancer_type_matches = [format_numbers(d) for d in identify_cancer_type(mutational_profile, classification_method)]
        # primary_site_matches = [format_numbers(d) for d in identify_primary_site(mutational_profile, classification_method)]

        # if method == 'M':
        #     _, contributing_signatures_A = decompose_mutational_profile(mutational_profile, "CLA")
        #     _, contributing_signatures_B = decompose_mutational_profile(mutational_profile, "CLB")
        #     _, contributing_signatures_cosmic = decompose_mutational_profile(mutational_profile, "IS")
        if method == 'L':
            _, contributing_signatures_A = decompose_mutational_profile_counts(mutational_profile_counts, "CLA", 'MLE')
            _, contributing_signatures_B = decompose_mutational_profile_counts(mutational_profile_counts, "CLB", 'MLE')
            _, contributing_signatures_cosmic = decompose_mutational_profile_counts(mutational_profile_counts, "IS", 'MLE')
        if method == 'J':
            _, contributing_signatures_A = decompose_mutational_profile_counts(mutational_profile_counts, "CLA", 'js')
            _, contributing_signatures_B = decompose_mutational_profile_counts(mutational_profile_counts, "CLB", 'js')
            _, contributing_signatures_cosmic = decompose_mutational_profile_counts(mutational_profile_counts, "IS", 'js')
        if method == 'R':
            _, contributing_signatures_A = decompose_mutational_profile_counts(mutational_profile_counts, "CLA", 'frobenius')
            _, contributing_signatures_B = decompose_mutational_profile_counts(mutational_profile_counts, "CLB", 'frobenius')
            _, contributing_signatures_cosmic = decompose_mutational_profile_counts(mutational_profile_counts, "IS", 'frobenius')
        if method == 'Z':
            _, contributing_signatures_A = decompose_mutational_profile_counts(mutational_profile_counts, "CLA", 'frobeniuszero')
            _, contributing_signatures_B = decompose_mutational_profile_counts(mutational_profile_counts, "CLB", 'frobeniuszero')
            _, contributing_signatures_cosmic = decompose_mutational_profile_counts(mutational_profile_counts, "IS", 'frobeniuszero')

        oname = prefix + "/ALL.decompose/{}_{}-5.txt".format(TCGA_ID, method)
        with open(oname, 'w') as o:
            for v in contributing_signatures_A:
                o.write("{}\t{}\n".format(v['name'], v['score']))
            # o.write("residuals\t{}\n".format(residuals_A))

        oname = prefix + "/ALL.decompose/{}_{}-10.txt".format(TCGA_ID, method)
        with open(oname, 'w') as o:
            for v in contributing_signatures_B:
                o.write("{}\t{}\n".format(v['name'], v['score']))
            # o.write("residuals\t{}\n".format(residuals_B))

        oname = prefix + "/ALL.decompose/{}_{}-30.txt".format(TCGA_ID, method)
        with open(oname, 'w') as o:
            for v in contributing_signatures_cosmic:
                o.write("{}\t{}\n".format(v['name'], v['score']))
            # o.write("residuals\t{}\n".format(residuals_cosmic))

