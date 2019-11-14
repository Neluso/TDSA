import os


dir_path = './output'


def write_data(freq, E_ref_w, E_sam_w):
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)

    wh_ref = open(dir_path + 'ref_spectrum.txt', 'w')
    wh_sam = open(dir_path + 'sam_spectrum.txt', 'w')
    for i in range(freq.size):
        wh_ref.write(str(freq[i]) + ',' + str(E_ref_w[i]))
        wh_sam.write(str(freq[i]) + ',' + str(E_sam_w[i]))
    wh_ref.close()
    wh_sam.close()
    return 0