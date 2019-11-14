from numpy import array


def read_data(filename):
    try:
        print('Opening ' + filename.split('/')[-1])
        fh = open(filename)
    except:
        print('Error opening ' + filename.split('/')[-1])
        print('Operation terminated. Please check if files exist.')
        quit()
    print(filename.split('/')[-1] + ' opened successfully')
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
