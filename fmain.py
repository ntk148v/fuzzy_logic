import collections
import os
from random import shuffle

WORK_DIR = os.path.dirname(os.path.realpath(__file__))
attrs = ['mcg', 'gvh', 'lip', 'chg', 'aac', 'alm1', 'alm2']
train_dataset = list()
test_dataset = list()


def get_data():
    global train_dataset, test_dataset
    dataset = list()
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
        dataset.append(d)
    shuffle(dataset)
    train_dataset = dataset[:len(dataset)*2/3]
    test_dataset = dataset[len(dataset)*2/3:]
    return (train_dataset, test_dataset)


def find_max_min():
    MAX = {}
    MIN = {}
    dataset = train_dataset + test_dataset
    for attr in attrs:
        MIN[attr] = dataset[0][attr]
        MAX[attr] = dataset[0][attr]

    for item in dataset:
        for attr in attrs:
            if float(item[attr]) < MIN[attr]:
                MIN[attr] = item[attr]
            if float(item[attr]) >= MAX[attr]:
                MAX[attr] = item[attr]

    print('\n-------> MIN: ' + str(MIN))
    print('\n-------> MAX: ' + str(MAX))
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

    print('\n-------> Z: ' + str(Z))
    print('\n-------> LOW: ' + str(LOW))
    print('\n-------> MED: ' + str(MED))
    print('\n-------> HIG: ' + str(HIG))

    return (Z, LOW, MED, HIG)


def gen_fuzzy_rules(MIN, Z, dataset, FILTED_RULES=None):
    count = 0
    RULES = dict()
    FIRED_RULES = list()
    for item in dataset:
        rule = dict()
        temp = dict()
        for k, v in item.items():
            if k == 'seq_name' or k == 'class':
                continue
            if v == MIN[k]:
                if len(rule) == 0:
                    temp['LOW'] = 1
                else:
                    for _k, _v in rule.items():
                        m = _v
                        new_k = _k + '+LOW'
                        temp[new_k] = m
                        del temp[_k]
            elif MIN[k] < v and v <= MIN[k] + Z[k] / 6:
                if len(rule) == 0:
                    temp['LOW'] = (3 * MIN[k] + Z[k] - 3 * v) / Z[k]
                else:
                    for _k, _v in rule.items():
                        m = _v
                        new_k = _k + '+LOW'
                        if ((3 * MIN[k] + Z[k] - 3 * v) / Z[k]) >= m:
                            temp[new_k] = m
                        else: temp[new_k] = (3 * MIN[k] + Z[k] - 3 * v) / Z[k]
                        #temp[new_k] = (
                        #    3 * MIN[k] + Z[k] - 3 * v) / Z[k] * temp[_k]
                        del temp[_k]
            elif MIN[k] + Z[k] / 6 < v and v <= MIN[k] + Z[k] / 3:
                if len(rule) == 0:
                    temp['LOW'] = (3 * MIN[k] + Z[k] - 3 * v) / Z[k]
                    temp['MED'] = (6 * v - 6 * MIN[k] - Z[k]) / (2 * Z[k])
                else:
                    for _k, _v in rule.items():
                        m = _v
                        new_k = _k + '+LOW'
                        new_k_2 = _k + '+MED'
                        if ((3 * MIN[k] + Z[k] - 3 * v) / Z[k]) >= m:
                            temp[new_k] = m
                        else:
                            temp[new_k] = (3 * MIN[k] + Z[k] - 3 * v) / Z[k]

                        if ((6 * v - 6 * MIN[k] - Z[k]) / (2 * Z[k])) >= m:
                            temp[new_k_2] = m
                        else:
                            temp[new_k_2] = (6 * v - 6 * MIN[k] - Z[k]) / (2 * Z[k])

                        #temp[new_k] = (
                        #    3 * MIN[k] + Z[k] - 3 * v) / Z[k] * temp[_k]
                        #temp[new_k_2] = (3 * MIN[k] + Z[k] -
                        #                 3 * v) / Z[k] * temp[_k]
                        del temp[_k]
            elif MIN[k] + Z[k] / 3 < v and v <= MIN[k] + 2 * Z[k] / 3:
                if len(rule) == 0:
                    if v < MIN[k] + Z[k] / 2:
                        temp['MED'] = (6 * v - 6 * MIN[k] - Z[k]) / (2 * Z[k])
                    else:
                        temp['MED'] = (
                            6 * MIN[k] + 5 * Z[k] - 6 * v) / (2 * Z[k])
                else:
                    for _k, _v in rule.items():
                        m = _v
                        new_k = _k + '+MED'
                        if v < MIN[k] + Z[k] / 2:
                            if ((6 * v - 6 * MIN[k] - Z[k]) / (2 * Z[k])) >= m:
                                temp[new_k] = m
                            else:
                                temp[new_k] = (6 * v - 6 * MIN[k] - Z[k]
                                            ) / (2 * Z[k])
                            #temp[new_k] = (6 * v - 6 * MIN[k] - Z[k]
                            #               ) / (2 * Z[k]) * temp[_k]
                        elif v ==  MIN[k] + Z[k] / 2:
                            temp[new_k] = m
                        else:
                            if ((6 * MIN[k] + 5 * Z[k] - 6 * v) / (2 * Z[k])) >= m:
                                temp[new_k] = m
                            else:
                                temp[new_k] = (6 * MIN[k] + 5 * Z[k] -
                                            6 * v) / (2 * Z[k])
                            #temp[new_k] = (6 * MIN[k] + 5 * Z[k] -
                            #               6 * v) / (2 * Z[k]) * temp[_k]
                        del temp[_k]
            elif MIN[k] + 2 * Z[k] / 3 < v and v <= MIN[k] + 5 * Z[k] / 6:
                if len(rule) == 0:
                    temp['MED'] = (6 * MIN[k] + 5 * Z[k] - 6 * v) / (2 * Z[k])
                    temp['HIG'] = (3 * v - 3 * MIN[k] - 2 * Z[k]) / Z[k]
                else:
                    for _k, _v in rule.items():
                        m = _v
                        new_k = _k + '+MED'
                        new_k_2 = _k + '+HIG'
                        if ((6 * MIN[k] + 5 * Z[k] - 6 * v) / (2 * Z[k])) >= m:
                            temp[new_k] = m
                        else:
                            temp[new_k] = (6 * MIN[k] + 5 * Z[k] - 6 * v) / (2 * Z[k])

                        if ((3 * v - 3 * MIN[k] - 2 * Z[k]) / Z[k]) >= m:
                            temp[new_k_2] = m
                        else:
                            temp[new_k_2] = (3 * v - 3 * MIN[k] - 2 * Z[k]) / Z[k]
                        #temp[new_k] = (6 * MIN[k] + 5 * Z[k] - 6 *
                        #               v) / (2 * Z[k]) * temp[_k]
                        #temp[new_k_2] = (3 * v - 3 * MIN[k] -
                        #                 2 * Z[k]) / Z[k] * temp[_k]
                        del temp[_k]
            elif MIN[k] + 5 * Z[k] / 6 < v and v < MIN[k] + Z[k]:
                if len(rule) == 0:
                    temp['HIG'] = (3 * v - 3 * MIN[k] - 2 * Z[k]) / Z[k]
                else:
                    for _k, _v in rule.items():
                        m = _v
                        new_k = _k + '+HIG'
                        if ((3 * v - 3 * MIN[k] - 2 * Z[k]) / Z[k]) >= m:
                            temp[new_k] = m
                        else:
                            temp[new_k] = (3 * v - 3 * MIN[k] - 2 * Z[k]) / Z[k]
                        #temp[new_k] = (3 * v - 3 * MIN[k] - 2 *
                        #               Z[k]) / Z[k] * temp[_k]
                        del temp[_k]
            elif  v == MIN[k]+Z[k]:
                if len(rule) == 0:
                    temp['HIG'] = 1
                else:
                    for _k, _v in rule.items():
                        m = _v
                        new_k = _k + '+HIG'
                        temp[new_k] = m
                        del temp[_k]
            rule = temp.copy()
            count += len(rule)
        if dataset == train_dataset:
            for _k, _v in rule.items():
                if _v <= 0:
                    continue
                if not RULES.has_key(_k):
                    RULES[_k] = [_v, item['class']]
                else:
                    if RULES[_k][0] < _v:
                        RULES[_k][0] = _v
                        RULES[_k][1] = item['class']
        else:
            # max_v = 0
            # max_k = ''
            # for _k, _v in rule.items():
            #     if _v > max_v:
            #         max_v = _v
            #         max_k = _k
            # FIRED_RULES.append({max_k: max_v, 'class': item['class']})
            FIRED_RULES.append({'rules': rule, 'class': item['class']})
    if dataset == train_dataset:
        # RULES = collections.OrderedDict(sorted(RULES.items()))
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
    # Tinh do chinh xac cua suy dien
    num_of_corrected_class = 0
    for e in FIRED_RULES:
        tmp = dict()
        for _k, _v in e['rules'].items():
            if RULES.has_key(_k):
                if tmp.has_key((RULES.get(_k))[1]):
                    tmp[(RULES.get(_k))[1]] += _v
                else:
                    tmp[(RULES.get(_k))[1]] = _v
        max_class = ''
        max_values = 0
        for _k, _v in tmp.items():
            if _v > max_values :
                max_values = _v
                max_class = _k
        if max_class == e['class']:
            num_of_corrected_class += 1
    print("Precision: {}, corrected class: {}, class test: {}".format(str(float(num_of_corrected_class)/len(test_dataset)*100), num_of_corrected_class, len(test_dataset)))
if __name__ == '__main__':
    main()
