from multiset import MultiSet


def vec_to_multiset(table):
    result = [MultiSet() for _ in range(4)]
    
    for row in table:
        for index, multiset in enumerate(result):
            multiset.push(row[index])
    
    return result