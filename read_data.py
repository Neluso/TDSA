from numpy import array


def read_data(namefile):
    try:
        print('Opening ' + namefile)
        fh = open('./data/' + namefile)
    except:
        print('Error opening ' + namefile)
        quit()
    print(namefile + ' opened successfully')
    data = fh.read()
    data = data.split('\n')
    x = list()
    y = list()
    for item in data:
        item = item.split(',')
        if item == ['']:
            break
        x.append(float(item[0]))
        y.append(float(item[1]))
    x = array(x)
    y = array(y)
    return x, y
