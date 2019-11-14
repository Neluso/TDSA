import os


def write_data(freq, E_ref_w, freq_sam, E_sam_w, sam_file):
    dir_path = './output'
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)
    dir_path = './output/' + sam_file + '/'
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)

    wh_ref = open(dir_path + 'ref_spectrum.txt', 'w')
    wh_sam = open(dir_path + 'sam_spectrum.txt', 'w')
    for i in range(freq.size):
        wh_ref.write(str(freq[i]) + ',' + str(E_ref_w[i]) + '\n')
        wh_sam.write(str(freq_sam[i]) + ',' + str(E_sam_w[i]) + '\n')
    wh_ref.close()
    wh_sam.close()
    return dir_path