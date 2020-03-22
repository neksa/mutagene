#
# Code copied from mutagene website
#

# from sqlalchemy.sql.expression import bindparam
# from collections import defaultdict

# from .mutations import *
# from ..database.models import *
# from ..database import Session as db


def export_cohorts(PATH):

    ######################################################
    # Export profiles
    ######################################################
    profiles = []
    q = db.query(Signature.id, Signature.name, Signature.mut_type, Signature.n_mutations,
                 Signature.nr_mutations, Signature.n_samples, Signature.signature_type, Signature.cancer_type_id,
                 CancerType.name.label("cancer_type"), PrimarySite.name.label("primary_site"),
                 CancerType.code.label("cancer_code")).\
        join(CancerType).join(PrimarySite).\
        filter(Signature.signature_type == "CT").\
        filter(Signature.mut_type == "A").\
        filter(Signature.normalization == "MUT").\
        filter(Signature.active == 1).\
        order_by(Signature.name)
    profiles.extend(q)

    q = db.query(Signature.id, Signature.name, Signature.mut_type, Signature.n_mutations,
                 Signature.nr_mutations, Signature.n_samples, Signature.signature_type,
                 bindparam("cancer_type", "Pan-cancer"), bindparam("primary_site", "All tissues"),
                 bindparam("cancer_code", ""), bindparam("cancer_projects", "")).\
        filter(Signature.signature_type == "PC").\
        filter(Signature.mut_type == "A").\
        filter(Signature.normalization == "MUT").one()
    profiles.append(q)

    # Bening Pan-tissue
    q = db.query(Signature.id, Signature.name, Signature.mut_type, Signature.n_mutations,
                 Signature.nr_mutations, Signature.n_samples, Signature.signature_type,
                 bindparam("cancer_type", "Benign"), bindparam("primary_site", "All tissues"),
                 bindparam("cancer_code", ""), bindparam("cancer_projects", "")).\
        filter(Signature.signature_type == "BEN").\
        filter(Signature.mut_type == "A").\
        filter(Signature.primary_site == None).one()
    profiles.append(q)

    # dbSNP
    q = db.query(Signature.id, Signature.name, Signature.mut_type, Signature.n_mutations,
                 Signature.nr_mutations, Signature.n_samples, Signature.signature_type,
                 bindparam("cancer_type", "SNP"), bindparam("primary_site", "Germline"),
                 bindparam("cancer_code", ""), bindparam("cancer_projects", "")).\
        filter(Signature.signature_type == "SNP").\
        filter(Signature.mut_type == "A").\
        filter(Signature.primary_site == None).one()
    profiles.append(q)

    # profiles.extend(q)
    # profiles.extend(q)

    # for each cancer type
    all_cancer_types = []
    for profile in profiles:
        # skip the second LUAD dataset
        if profile.id == 34904:
            continue
        all_cancer_types.append(profile.id)
        print("Processing {}".format(profile.name))
        data = {}
        s = db.query(Signature).get(profile.id)
        for d in s.data:
            idx = d.c53 + d.XY
            # data[idx] = d.freq
            data[idx] = d.count

        result = ""
        for p5 in nucleotides:
            for p3 in nucleotides:
                for x in "CT":
                    for y in nucleotides:
                        if x != y:
                            result += "{}[{}>{}]{}\t{}\n".format(p5, x, y, p3, data.get(p5 + p3 + x + y, 0), 8)

        fname = "{}/{}.profile".format(PATH, profile.name.replace(" ", "_"))
        with open(fname, 'w') as o:
            o.write("#NSAMPLES\t{}\n".format(profile.n_samples))
            o.write(result)

        ######################################################
        # Export mutations for each profile
        ######################################################
        if profile.signature_type == "CT":
            mm = db.query(Mutation).join(SignatureMutation).filter(SignatureMutation.signature_id == int(profile.id))
        elif profile.signature_type == "PC":
            mm = db.query(Mutation).join(SignatureMutation).filter(SignatureMutation.signature_id.in_(all_cancer_types))
        aa_mutations = defaultdict(int)
        na_mutations = defaultdict(int)

        for m in mm:
            if len(m.grc_pos) == 0 or len(m.mut_aa) == 0:
                continue
            chrom, pos = m.grc_pos.split('-')[0].split(':')
            # print(m.mut_aa)
            # print(m.mut_cds)
            mut_aa = m.mut_aa.split('.')[1]
            X, Y = m.context.XY
            aa_mutations[(m.gene.name, mut_aa)] += 1
            na_mutations[(chrom, pos, X, Y)] += 1

        if len(aa_mutations):
            with open("{}/{}.aa_mutations.txt".format(PATH, profile.name.replace(" ", "_")), 'w') as aa:
                for key, val in aa_mutations.items():
                    aa.write("{}\t{}\t{}\n".format(*key, val))

        if len(na_mutations):
            with open("{}/{}.dna_mutations.txt".format(PATH, profile.name.replace(" ", "_")), 'w') as na:
                for key, val in na_mutations.items():
                    na.write("{}\t{}\t{}\t{}\t{}\n".format(*key, val))
