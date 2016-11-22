import collections
import os


WORK_DIR = os.path.dirname(os.path.realpath(__file__))
attrs = ['mcg', 'gvh', 'lip', 'chg', 'aac', 'alm1', 'alm2']
train_dataset = list()
test_dataset = list()


def get_data():
    global train_dataset, test_dataset
    lines = [line.rstrip('\n') for line in open(WORK_DIR + '/ecoli.data')]
    count = 0
    for line in lines:
        count += 1
        data = []
        for token in line.split('  '):
            if token.strip():
                data.append(token.strip())
        for i in range(1, 8):
            data[i] = float(data[i])
        d = dict(zip(['seq_name', 'mcg', 'gvh', 'lip', 'chg',
                      'aac', 'alm1', 'alm2', 'class'], data))
        if count <= len(lines) * 2 / 3:
            train_dataset.append(d)
        else:
            test_dataset.append(d)
    return (train_dataset, test_dataset)


def find_max_min():
    MAX = {}
    MIN = {}
    dataset = train_dataset + test_dataset
    for attr in attrs:
        MIN[attr] = dataset[0][attr]
        MAX[attr] = dataset[0][attr]

    for item in train_dataset:
        for attr in attrs:
            if float(item[attr]) < MIN[attr]:
                MIN[attr] = item[attr]
            if float(item[attr]) >= MAX[attr]:
                MAX[attr] = item[attr]

    # print('\n-------> MIN: ' + str(MIN))
    # print('\n-------> MAX: ' + str(MAX))
    return (MAX, MIN)


def get_likelihood_function(MAX, MIN):
    Z = {}
    LOW = {}
    MED = {}
    HIG = {}

    for attr in attrs:
        Z[attr] = MAX[attr] - MIN[attr]
        LOW[attr] = [MIN[attr], MIN[attr] + Z[attr] / 3]
        MED[attr] = [MIN[attr] + Z[attr] / 6, MIN[attr] + 5 * Z[attr] / 6]
        HIG[attr] = [MIN[attr] + 2 * Z[attr] / 3, MIN[attr] + Z[attr]]

    # print('\n-------> Z: ' + str(Z))
    # print('\n-------> LOW: ' + str(LOW))
    # print('\n-------> MED: ' + str(MED))
    # print('\n-------> HIG: ' + str(HIG))

    return (Z, LOW, MED, HIG)


def gen_fuzzy_rules(MIN, Z, dataset):
    RULES = dict()
    FIRED_RULES = list()
    for item in dataset:
        rule = dict()
        temp = dict()
        for k, v in item.items():
            temp = rule.copy()
            if k == 'seq_name' or k == 'class':
                continue
            if MIN[k] < v and v < MIN[k] + Z[k] / 6:
                if len(rule) == 0:
                    temp['LOW'] = (3 * MIN[k] + Z[k] - 3 * v) / Z[k]
                else:
                    for _k, _v in rule.items():
                        new_k = _k + '+LOW'
                        temp[new_k] = (
                            3 * MIN[k] + Z[k] - 3 * v) / Z[k] * temp[_k]
                        del temp[_k]
            elif MIN[k] + Z[k] / 6 < v and v < MIN[k] + Z[k] / 3:
                if len(rule) == 0:
                    temp['LOW'] = (3 * MIN[k] + Z[k] - 3 * v) / Z[k]
                    temp['MED'] = (6 * v - 6 * MIN[k] - Z[k]) / (2 * Z[k])
                else:
                    for _k, _v in rule.items():
                        new_k = _k + '+LOW'
                        new_k_2 = _k + '+MED'
                        temp[new_k] = (
                            3 * MIN[k] + Z[k] - 3 * v) / Z[k] * temp[_k]
                        temp[new_k_2] = (3 * MIN[k] + Z[k] -
                                         3 * v) / Z[k] * temp[_k]
                        del temp[_k]
            elif MIN[k] + Z[k] / 3 < v and v < MIN[k] + 2 * Z[k] / 3:
                if len(rule) == 0:
                    if v < MIN[k] + Z[k] / 2:
                        temp['MED'] = (6 * v - 6 * MIN[k] - Z[k]) / (2 * Z[k])
                    else:
                        temp['MED'] = (
                            6 * MIN[k] + 5 * Z[k] - 6 * v) / (2 * Z[k])
                else:
                    for _k, _v in rule.items():
                        new_k = _k + '+MED'
                        if v < MIN[k] + Z[k] / 2:
                            temp[new_k] = (6 * v - 6 * MIN[k] - Z[k]
                                           ) / (2 * Z[k]) * temp[_k]
                        else:
                            temp[new_k] = (6 * MIN[k] + 5 * Z[k] -
                                           6 * v) / (2 * Z[k]) * temp[_k]
                        del temp[_k]
            elif MIN[k] + 2 * Z[k] / 3 < v and v < MIN[k] + 5 * Z[k] / 6:
                if len(rule) == 0:
                    temp['MED'] = (6 * MIN[k] + 5 * Z[k] - 6 * v) / (2 * Z[k])
                    temp['HIG'] = (3 * v - 3 * MIN[k] - 2 * Z[k]) / Z[k]
                else:
                    for _k, _v in rule.items():
                        new_k = _k + '+MED'
                        new_k_2 = _k + '+HIG'
                        temp[new_k] = (6 * MIN[k] + 5 * Z[k] - 6 *
                                       v) / (2 * Z[k]) * temp[_k]
                        temp[new_k_2] = (3 * v - 3 * MIN[k] -
                                         2 * Z[k]) / Z[k] * temp[_k]
                        del temp[_k]
            else:
                if len(rule) == 0:
                    temp['HIG'] = (3 * v - 3 * MIN[k] - 2 * Z[k]) / Z[k]
                else:
                    for _k, _v in rule.items():
                        new_k = _k + '+HIG'
                        temp[new_k] = (3 * v - 3 * MIN[k] - 2 *
                                       Z[k]) / Z[k] * temp[_k]
                        del temp[_k]
            rule = temp.copy()

        if dataset == train_dataset:
            for _k, _v in rule.items():
                if _v <= 0:
                    continue
                if _k not in RULES:
                    RULES[_k] = [_v, item['class']]
                else:
                    if RULES[_k][0] < _v:
                        RULES[_k][0] = _v
                        RULES[_k][1] = item['class']
        else:
            max_v = 0
            max_k = ''
            for _k, _v in rule.items():
                if _v > max_v:
                    max_v = _v
                    max_k = _k
            FIRED_RULES.append({max_k: max_v, 'class': item['class']})
    if dataset == train_dataset:
        RULES = collections.OrderedDict(sorted(RULES.items()))
        print('\n -------> FILTER RULES: {}\n'.format(len(RULES)))
        return RULES
    else:
        print('\n -------> FILTER FIRED SHIT RULES: {}\n'.format(len(FIRED_RULES)))
        return FIRED_RULES


def main():
    get_data()
    MAX, MIN = find_max_min()
    Z, LOW, MED, HIG = get_likelihood_function(MAX, MIN)
    RULES = gen_fuzzy_rules(MIN, Z, train_dataset)
    FIRED_RULES = gen_fuzzy_rules(MIN, Z, test_dataset)
    TMP = []
    for e in FIRED_RULES:
        for _k in e.keys():
            if _k == 'class' or _k in TMP:
                continue
            if _k in RULES.keys():
                TMP.append(_k)
                print(sorted(e.items()))


if __name__ == '__main__':
    main()
