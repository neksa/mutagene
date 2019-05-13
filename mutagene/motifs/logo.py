from .motifs import get_enrichment
from matplotlib import plt

"""
def create_logo(mutations_with_context, motif, pos, ref, alt, range_size=50):
    # TO-DO: sort by nuc height for logo

    x_values = [-1, 0, 1]

    get_enrichment(mutations_with_context, motif, pos, ref, alt, range_size)
    c_pre_count = 0
    g_pre_count = 0
    t_pre_count = 0
    a_pre_count = 0

    c_post_count = 0
    g_post_count = 0
    t_post_count = 0
    a_post_count = 0

    pre_sum = 0
    post_sum = 0

    for char in pre_base:
        if char == "C":
            c_pre_count += 1
        elif char == "A":
            a_pre_count += 1
        elif char == "T":
            t_pre_count += 1
        elif char == "G":
            g_pre_count += 1
        pre_sum += 1

    for char in post_base:
        if char == "C":
            c_post_count += 1
        elif char == "A":
            a_post_count += 1
        elif char == "T":
            t_post_count += 1
        elif char == "G":
            g_post_count += 1
        post_sum += 1

    pre_scale = [c_pre_count * 1.0 / pre_sum, a_pre_count * 1.0 / pre_sum, t_pre_count * 1.0 / pre_sum, g_pre_count * 1.0 / pre_sum]
    post_scale = [c_post_count * 1.0 / post_sum, a_post_count * 1.0 / post_sum, t_post_count * 1.0 / post_sum, g_post_count * 1.0 / post_sum]

    pre_bit = 0
    for scale in pre_scale:
        pre_bit += scale * (math.log(2, scale))
    pre_bit *= -1

    post_bit = 0
    for po_scale in post_scale:
        post_bit += po_scale * (math.log(2, po_scale))
    post_bit *= -1

    ALL_PRE_NUCS = [('C', pre_scale[0]), ('A', pre_scale[1]), ('T', pre_scale[2]), ('G', pre_scale[3])]
    ALL_POST_NUCS = [('C', post_scale[0]), ('A', post_scale[1]), ('T', post_scale[2]), ('G', post_scale[3])]

    fig = plt.figure()
    fig.set_size_inches(len(motif), 2.5)
    fig.set_size_inches(len(ALL_PRE_NUCS),2.5)
    ax = fig.add_subplot(111)
    ax.tick_params(axis='x', which='major', labelsize=8)
    plt.xticks((0.15, 0.5, 0.85), x_values)

    plt.ylabel('bits')

    trans_offset = transforms.offset_copy(ax.transAxes,
                                      fig=fig,
                                      x=0,
                                      y=0,
                                      units='points')

    bit_pre_count = 0
    bit_post_count = 0

    ALL_PRE_NUCS_SORTED = []
    count = 0
    while count < 4:
        index = 0
        min_index = 0
        mini = ALL_PRE_NUCS[0]
        for nuc1, bit1 in ALL_PRE_NUCS:
            if bit1 < mini[1]:
                mini = nuc1, bit1
                min_index = index
            index += 1
        ALL_PRE_NUCS_SORTED.insert(count, mini)
        ALL_PRE_NUCS.pop(min_index)
        count += 1

    ALL_POST_NUCS_SORTED = []
    count = 0
    while count < 4:
        index = 0
        min_index = 0
        mini = ALL_POST_NUCS[0]
        for nuc1, bit1 in ALL_POST_NUCS:
            if bit1 < mini[1]:
                mini = nuc1, bit1
                min_index = index
            index += 1
        ALL_POST_NUCS_SORTED.insert(count, mini)
        ALL_POST_NUCS.pop(min_index)
        count += 1

    for nuc, bit in ALL_PRE_NUCS_SORTED:
        txt = plt.text(0.15, bit_pre_count, nuc, transform=trans_offset, horizontalalignment='center', fontsize=160 * bit, color=COLOR_SCHEME[nuc],
              weight = 'bold')
        bit_pre_count += bit
        txt.set_clip_on(False)

    plt.text(0.5, 0, ref, transform=trans_offset, horizontalalignment='center', fontsize=140, color=COLOR_SCHEME[ref],
             weight = 'bold')

    for nuc, bit in ALL_POST_NUCS_SORTED:
        txt = plt.text(0.85,bit_post_count, nuc, transform = trans_offset, horizontalalignment='center', fontsize=160*bit, color=COLOR_SCHEME[nuc],
              weight = 'bold')
        bit_post_count += bit
        txt.set_clip_on(False)

    plt.show()

    return pre_bit, post_bit

"""
