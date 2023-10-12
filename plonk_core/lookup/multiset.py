from dataclasses import dataclass
from typing import List
from plonk_core.src import utils
from collections import defaultdict
from field import field
import gmpy2
@dataclass
class MultiSet:
    
    elements:List

    def push(self, element):
        self.elements.append(element)

    def pad(self, n):
        #judge whether n is power of 2
        assert n & (n - 1) == 0  
        #if self is empty, push one 0
        if not self.elements:
            self.push(0)  
        #if size < n, push 0
        while len(self.elements) < n:
            self.elements.append(self.elements[0])  # use the first element to pad
        
    def compress(self, alpha):
        compress_poly = utils.Multiset_lc(self, alpha)
        compress_poly = MultiSet(compress_poly)
        return compress_poly

    def combine_split(self, f_elements:'MultiSet'):
        # create buckets and init
        counters = defaultdict(gmpy2.mpz)
        for element in self.elements:
            counters[element.value] += 1

        # Insert the elements of f into the corresponding bucket and 
        # check whether there is a corresponding element in t
        for element in f_elements.elements:
            if element.value in counters and counters[element.value] > 0:
                counters[element.value] += 1
            else:
                raise ValueError("ElementNotIndexed")

        # Split s into two alternating halves evens and odd

        evens = []
        odds = []
        parity =0
        for key, value in counters.items():
            key = field(key,self.elements[0].params)
            half_count = value//2
            evens.extend([key for _ in range(half_count)])
            odds.extend([key for _ in range(half_count)])
            if value % 2 ==1:
                if parity == 1:
                    odds.append(key)
                    parity = 0
                else:
                    evens.append(key)
                    parity = 1

        return evens, odds