import math

# takes around 90 seconds for 1023 minterms
# takes around 10 minutes for 2047 minterms or 11 vars 
MAX_FUNCTION_VARS = 10 

# converts the terms to binary strings of size corresponding to number of variables
def terms_to_bin_str(terms: list[int], number_of_vars: int):
    res = [] 
    fmt_str = "#0" + str(number_of_vars + 2) + 'b' 
    for num in terms:
        res.append(format(num, fmt_str)[2:]) 
    
    return res



#https://en.wikipedia.org/wiki/Quine%E2%80%93McCluskey_algorithm
def __check_dash_align(min1: str, min2: str):
    for i in range(len(min1)): 
        #If one minterm has dashes and the other does not then the minterms cannot be merged. 
        if min1[i] != '-' and min2[i] == '-': 
            return False 
    return True

def __check_minterm_diff(min1, min2):
    #minterm1 and minterm2 are strings representing all of the currently found prime implicants and merged 
    #minterms. Examples include '01--' and '10-0'
    #integer representation of minterm1 and minterm2 with the dashes removed, these are replaced with 0
    str_m1 = ""
    for c in min1:
        if c == "-":
            str_m1 += '0' 
        else:
            str_m1 += c

    str_m2 = ""
    for c in min2:
        if c == "-":
            str_m2 += '0' 
        else:
            str_m2 += c

    m1 = int(str_m1, 2)
    m2 = int(str_m2, 2)

    res =  m1 ^ m2
    return res != 0 and (res & res - 1) == 0

def __merge_minterms(min1: str, min2: str):
    mergedMinterm = ""
    for i in range(len(min1)): 
        #If the bits vary then replace it with a dash, otherwise the bit remains in the merged minterm.
        if min1[i] != min2[i]:
            mergedMinterm += '-'
        else:
            mergedMinterm += min1[i] 
    return mergedMinterm

def prime_implicants(minterms: list[str]):
    p_implicants = []
    merges = [False] * len(minterms)
    num_merges = 0
    mergedMinterm = ""
    minterm1 = ""
    minterm2 = "" 

    for i in range(len(minterms)):
        for j in range(i+1,len(minterms)):
            minterm1 = minterms[i]
            minterm2 = minterms[j]
            if __check_dash_align(minterm1, minterm2) and __check_minterm_diff(minterm1, minterm2):
                mergedMinterm = __merge_minterms(minterm1, minterm2) 
                if mergedMinterm not in p_implicants:
                    p_implicants.append(mergedMinterm)
                num_merges += 1
                merges[i] = True
                merges[j] = True

    #Filtering all minterms that have not been merged as they are prime implicants. Also removing duplicates
    for j in range(len(minterms)):
        if merges[j] == False and minterms[j] not in p_implicants: 
            p_implicants.append(minterms[j])

    #if no merges have taken place then all of the prime implicants have been found so return, otherwise 
    #keep merging the minterms. 
    if num_merges == 0:
        return p_implicants 
    else: 
        return prime_implicants(p_implicants)


def __match_str(p_implicant: str, minterm: str):
    for i in range(len(p_implicant)):
        if p_implicant[i] == '-' and not minterm[i].isdigit():
            return False
        elif p_implicant[i] != '-' and p_implicant[i] != minterm[i]:
            return False

    return True 


def create_prime_implicant_chart(p_implicants: list[str], minterms: list[str]): 
    p_implicant_chart: dict[str,str] = {}
    #Creating the empty chart with the prime implicants as the key and empty strings as the value. 
    for i in range(len(p_implicants)):
        p_implicant_chart[p_implicants[i]] = ""
    
    for key in p_implicant_chart.keys():
        p_implicant= key
        #Convert the "-" to "\d" which can be used to find the row of ticks above.
        for j in range(len(minterms)):
            #If there is a match between the regular expression and the minterm than append a 1 otherwise 0. 
            if __match_str(p_implicant,minterms[j]):
                p_implicant_chart[p_implicant] += "1"
            else: 
                p_implicant_chart[p_implicant] += "0"

    #The prime implicant chart is complete so return the completed chart. 
    return p_implicant_chart



# https://en.wikipedia.org/wiki/Petrick%27s_method
# pos is represented as set of tuples 
# where each term in the tuple represents a sum
# each tuple is a product
def __to_pos(implicant_chart: dict[str, str]) -> set[tuple[int]]:
    temp_pos = {}
    for implicant_index,(_,value) in enumerate(implicant_chart.items()):
        for col, bit in enumerate(value):
            # ensure all our numbers are multiples of 2
            temp_item = 2**implicant_index

            if bit == '1':
                try:
                    temp_pos[col].append(temp_item)
                except KeyError:
                    temp_pos[col] = [temp_item]

    pos_res = set()
    for _, sum in temp_pos.items():
        pos_res.add(tuple(sum))

    return pos_res


# convert tuple of ints to just ints 
# where each int is a term in the sum
def __to_sop(pos: set[tuple[int]]):
    if len(pos) == 1:
        res = []
        # get the only sum in the set
        sum = pos.pop()
       
        # loop through each term in the sum and add it to the result
        for term in sum:
            res.append(term)
        return res

    sop1: set[int] = set()
    sop2: set[int] = set()
   
    # expand the first two products 
    sum1 = pos.pop()
    sum2 = pos.pop()

    for i in range(len(sum1)):
        for j in range(len(sum2)):
            sop1.add(sum1[i] | sum2[j])

    # expand each product
    # since we are using sets we alternate 
    # between two sets that store the current sop
    # we expanding until all the pos terms have been used up
    while len(pos) > 0:
        temp_sum = pos.pop()
        if len(sop1) == 0:
            while(len(sop2) > 0):
                val = sop2.pop()
                for e in temp_sum:
                    sop1.add(e | val)
        else:
            while(len(sop1) > 0):
                val = sop1.pop()
                for e in temp_sum:
                    sop2.add(e | val)


    sop = list(sop1)
    if len(sop1) == 0:
        sop = list(sop2)



    # we handle absorption
    i = 0 
    j = 0
    while i < len(sop):
        j = i + 1
        while j < len(sop):
            # x ^ y = x, then y can be removed
            x = sop[i]
            y = sop[j]

            if x ^ y == x:
                sop.pop(j)

            j += 1
        i += 1


    # remove any remaing duplicates 
    return list(set(sop))

        




def __p_implicant_to_bool_expr(p_implicant: str, is_fsm=False):
    res = ""
    num_terms = 0
    for index, bit in enumerate(p_implicant): 
        if bit == '-':
            continue
        # convert each value to an uppercase letter
        char = chr(65 + index)
        if is_fsm and index == len(p_implicant) - 1:
            char = "X"

        if bit == "0":
            res += char + "'"
        else:
            res += char
        num_terms += 1

    return (res, num_terms)



def __get_minimal_expr(sop: list[int], implicant_chart: dict[str, str], is_fsm=False):
    fewest_terms = []
    smallest_len = sum([1 for i in range(sop[0].bit_length()) if sop[0] & (1 << i)])


    # getting the product(s) with the fewest terms 
    for product in sop:
        # each term is represented by a one in binary so we count the number of 1's in the product
        term_count= sum([1 for i in range(product.bit_length()) if product & (1 << i)])
        if term_count < smallest_len:
            fewest_terms = [product]
            smallest_len = term_count
        elif term_count == smallest_len:
            fewest_terms.append(product)

   
   
    # current the prime implicant keys to a list 
    keys = list(implicant_chart.keys())

    smallest_term_size = None
    smallest_expr = []
    
    # converting each item into a boolean expression and finding 
    # the smallest
    for product in fewest_terms: 
        current_term_size = 0
        expr = []
        # loop through each bit in the product
        for index in range(product.bit_length()):

            # if the bit is a one
            if product & (1 << index):
                p_implicant = keys[index]
                (b_expr, term_size) = __p_implicant_to_bool_expr(p_implicant, is_fsm)
                expr.append(b_expr)
                current_term_size += term_size
        
        if smallest_term_size == None or current_term_size < smallest_term_size:
            smallest_term_size = current_term_size 
            smallest_expr = expr

    
    return smallest_expr
       

def get_expr(number_of_vars: int, minterms: list[int], dc: list[int] = [], is_fsm=False):
    if number_of_vars > MAX_FUNCTION_VARS:
        print(f"Cannot have functions with more than ${MAX_FUNCTION_VARS}")
        return


    if len(minterms) == 0:
        return "F = 0"

    terms = minterms + dc

    max_minterm = 2**number_of_vars - 1

    for term in terms:
        if term > max_minterm or term < 0:
            raise ValueError(f"Invalid Minterm: {term} max term for {number_of_vars} vars is {max_minterm}")


    # remove any duplicates
    terms = set(terms)

    
    if len(terms) == 2**number_of_vars:
        return "F = 1"

    terms = list(terms)



    bin_terms = terms_to_bin_str(terms, number_of_vars)

    
    p_implicants = prime_implicants(bin_terms)

    # when creating the prime implicant chart we don't use the dont_care terms 
    # so we only care about the min terms 
    chart = create_prime_implicant_chart(p_implicants, bin_terms[:len(minterms)])

    pos = __to_pos(chart)

    sop = __to_sop(pos)

    mininum = __get_minimal_expr(sop, chart, is_fsm)

    expr = ""
    for product in mininum:
        expr += product + " + "

    return expr[:-3]



def get_truth_table(number_of_vars: int, minterms: list[int], dc: list[int] = []):
    if number_of_vars > MAX_FUNCTION_VARS:
        print("Cannot handle functions greater than ", MAX_FUNCTION_VARS, sep="")
        return 

    max_num= 2 ** number_of_vars

    binary_numbers = terms_to_bin_str([x for x in range(max_num + 1)], number_of_vars)

    max_num_len = len(str(max_num)) 

    num_spaces = 0 
    if max_num_len - 3 > 0: 
        num_spaces = max_num_len - 3
    else: 
        num_spaces = 1

   
    print("Dec", " " * num_spaces, "| ", sep = "", end="")
    # print the input variables
    for i in range(65, 65 + number_of_vars):
        print(chr(i), end=" ")

    print(" O/P")
    print("-" * (number_of_vars * 2 + 10 + max_num_len))

    padding = " "

    if max_num_len < 3: 
        padding = (3-max_num_len+1) * padding

    for i in range(max_num):
        fmt_str = "#0" + str(max_num_len) + "d"
        print(format(i, fmt_str),padding,  sep="", end="| ")
        for bit in binary_numbers[i]:
            print(bit, end=" ")

        output = "0"
        if i in minterms: 
            output = "1"
        elif i in dc:
            output = "X"

        print("| ", output, sep="")


class FiniteStateMachine:
    def __init__(self, is_mealy=False):
        self.is_mealy = is_mealy
        self.states = None
        self.largest_state = 0
        self.next_pow_2 = 0 # next power of two after the largest state
        self.largest_minterm_length = 0 # bit length of the largest minterm

   
    def add_state_table(self, states: dict[int,list[str]]):
        self.states = states


    def print_state_table(self):
        if self.states == None:
            print(None)
            return

        for state in self.states.keys(): 
            arr = self.states[state]
            print('S' + str(state), end="    ")
            print('S' + str(arr[0]),  'S' + str(arr[1]),sep="    ", end="    ")


            if self.is_mealy:
                print(arr[2], arr[3], sep="    ")
            else:
                print(arr[2])


    def get_num_flip_flops(self):
        if self.states is None:
            return 0
        
        max = 0 
        for state in self.states.keys():
            if state > max:
                max = state


        self.largest_state = max 
        self.next_pow_2 = 2**(math.ceil(math.log2(max + 1)))
        self.largest_minterm_length = math.ceil(math.log2(self.next_pow_2 << 1))

        return int(math.log2(max)) + 1


    class __Terms:
        def __init__(self) -> None:
           self.minterms = []
           self.dont_care = []


        def add_minterm(self, term: int):
            self.minterms.append(term)


        def add_dont_care(self, dc: int):
            self.dont_care.append(dc)

        def get_terms(self):
            return (self.minterms, self.dont_care)





    # convert the state table to a list of minterms and don't cares
    def __get_minterms(self, num_flip_flops):
        if self.states is None:
            return None

        next_states_minterms = []

        largest_minterm = (self.largest_state << 1) + 1
        for i in range(num_flip_flops):
            # create the next state terms objects
            next_states_minterms.append(FiniteStateMachine.__Terms())

            # add don't cares if there are any
            for dc in range(largest_minterm + 1, 2**self.largest_minterm_length):
                next_states_minterms[i].add_dont_care(dc)



        output = FiniteStateMachine.__Terms()

        if self.is_mealy:
            for dc in range(largest_minterm + 1, 2**self.largest_minterm_length):
                output.add_dont_care(dc)
        else:
            for dc in range(self.largest_state + 1, self.next_pow_2):
                output.add_dont_care(dc)



        for current_state in self.states.keys():

            # next state for x = 0 and x = 1
            x_0 = self.states[current_state][0]
            x_1 = self.states[current_state][1]

            # two min terms for each state
            min_a = current_state << 1
            min_b = (current_state << 1) + 1


            if x_0 == 'x':
                for i in range(num_flip_flops):
                    next_states_minterms[i].add_dont_care(min_a)
            else:
                # convert the next state to binary 
                # loop the binary number and add the minterm to the list if its a one
                temp= terms_to_bin_str([int(x_0)], num_flip_flops)
                x_0_bin_str = temp[0]
                for i in range(num_flip_flops):
                    if x_0_bin_str[i] == '1':
                        next_states_minterms[i].add_minterm(min_a) 

            if x_1 == 'x':
                for i in range(num_flip_flops):
                    next_states_minterms[i].add_dont_care(min_b)
            else:
                temp= terms_to_bin_str([int(x_1)], num_flip_flops)
                x_1_bin_str = temp[0]
                for i in range(num_flip_flops):
                    if x_1_bin_str[i] == '1':
                        next_states_minterms[i].add_minterm(min_b)
            



            # for a moore machine, the output is not dependent on x 
            # so the min terms will be 2 times smaller for the output
            temp_min_a = min_a
            if not self.is_mealy:
                temp_min_a >>= 1
           
            
            z_0 = self.states[current_state][2]
            if z_0 == '1':
                output.add_minterm(temp_min_a)
            elif z_0 != '0':
                output.add_dont_care(temp_min_a)

            if self.is_mealy:
                z_1 = self.states[current_state][3]
                if z_1 == '1':
                    output.add_minterm(min_b)
                elif z_1 != '0':
                    output.add_dont_care(min_b)
            
            
        return (next_states_minterms, output)




    def get_boolean_expressions(self):
        num_flip_flops = self.get_num_flip_flops()
        num_vars= self.largest_minterm_length


        temp = self.__get_minterms(num_flip_flops)
        if temp is None:
            return

        (next_states_minterms, output) = temp


        for i in range(len(next_states_minterms)):
            (mterms, dont_care) = next_states_minterms[i].get_terms()
            flip_flop_var = chr(65 + i) + "+"

            expr = get_expr(num_vars, mterms, dont_care, is_fsm=True)
            print(flip_flop_var, expr, sep=" = ")



        (z_minterms, z_dc)= output.get_terms()
        z_expr = "" 
        if self.is_mealy:
            z_expr = get_expr(num_vars, z_minterms, z_dc, is_fsm=True) 
        else:
            z_expr = get_expr(num_vars - 1, z_minterms, z_dc, is_fsm=True)
 
        print("Z = ", z_expr)



class SequenceDetector(FiniteStateMachine):

    # checks to see which state the fsm should go to if the input is not the next bit in the sequence
    def __find_sub_sequence(self, other_seq: str, each_state_sequences: list[str]):
        while len(other_seq) > 0:
            for state,state_seq in enumerate(each_state_sequences):
                if state_seq == other_seq:
                    return str(state)

            other_seq = other_seq[1:]

        # if no sub sequence was found 
        # we go back to reset
        return '0'


    def __init__(self, sequence: str, allow_overlap=True, is_mealy=False):
        super().__init__(is_mealy)


        #verify the string only contains ones and zeros
        for c in sequence:
            if c != '1' and c != '0':
                raise ValueError("Sequence Must be in Binary")


        state_table: dict[int, list[str]] = {}

 
        # for the first pass through we only care about advancing to the next sequence 
        # we use n as a placeholder
        for cs,bit in enumerate(sequence):
            temp = ['N', 'N', '0', '0']
            temp[int(bit)] = str(cs + 1)

            if cs == len(sequence) - 1 and is_mealy:
                temp[int(bit)+ 2] = '1'
                temp[int(bit)] = 'N'

            if is_mealy:
                state_table[cs] = temp
            else:
                state_table[cs] = temp[:-1]


        if not is_mealy:
            state_table[len(sequence)] = ['N', 'N', '1']

        # holds the sequence for each state 
        # index 0 is reset so there is no sequence
        each_state_sequences = [""]

        for state in state_table.keys():
            temp = state_table[state]
            n_index = 0
            if temp[1] == 'N':
                n_index = 1

            # stay at reset if we don't advance
            if state == 0:
                temp[n_index] = '0'
                continue


            current_seq = sequence[:state]
            each_state_sequences.append(current_seq)

            if not allow_overlap and state == len(state_table) - 1:
                # go back to reset on a moore machine if you don't allow overlap 
                if not is_mealy:
                    temp[0] = '0'
                    temp[1] = '0'
                    break

                # mealy machine only the input that detects 
                # the sequence goes back to reset

                # if the output at for x = 0 is a 1 
                # seq is detected so we go back to reset for next state when x = 0
                # check where to go when x = 1
                if temp[2] == '1':
                    temp[0] = '0'
                    non_advancing_seq = current_seq + str( 1)
                    temp[1] = self.__find_sub_sequence(non_advancing_seq, each_state_sequences)

                # output is 1 when x = 1 
                # next state when x = 1 is then 0 
                # check where to go when x = 0
                else:
                    temp[1] = '0'
                    non_advancing_seq = current_seq + str( 0)
                    temp[0] = self.__find_sub_sequence(non_advancing_seq, each_state_sequences)

                break


            non_advancing_seq = current_seq + str(n_index)
            temp[n_index] = self.__find_sub_sequence(non_advancing_seq, each_state_sequences)

            if state == len(state_table) - 1:
                n_index =  int(not n_index)
                non_advancing_seq = current_seq + str(n_index)
                temp[n_index] = self.__find_sub_sequence(non_advancing_seq, each_state_sequences)
        
        self.add_state_table(state_table)
