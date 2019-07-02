

    # @manager.command
    # def dump_nmf_decomp():
    #     for AB in ('A', 'B'):
    #         deconv = db.query(CachedObjects).get('deconvolution_' + AB).value

    #         # load samples in 5D or 10D
    #         W = deconv['W']
    #         H = deconv['H'].T
    #         ids = deconv['ids']

    #         for i, ID in enumerate(ids):
    #             profile = db.query(Signature).get(int(ID))
    #             v = np.array([get_signature_values_list(profile)]).ravel()
    #             x = np.array(H[i]).squeeze()
    #             residuals = np.linalg.norm(W.dot(x) - v)  # Frobenius norm
    #             N = '5' if AB == 'A' else '10'
    #             oname = "/Users/gonceare/projects/RA/ALL.decompose/{}_NMF{}.txt".format(ID, N)
    #             with open(oname, 'w') as o:
    #                 for j in range(H.shape[1]):
    #                     name = "MUTAGENE {}.{}".format(AB, j + 1)
    #                     o.write("{}\t{}\n".format(name, H[i, j]))
    #                 # print(ID, residuals)
    #                 o.write("Residuals\t{}\n".format(residuals))


