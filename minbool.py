# kind of arbitarty right now 
# more testing need to be done 
# this may be even too much for this implementation 
MAX_FUNCTION_VARS = 16

# converts the terms to binary strings of size corresponding to number of variables
def terms_to_bin_str(terms: list[int], number_of_vars: int):
    res = [] 
    fmt_str = "#0" + str(number_of_vars + 2) + 'b' 
    for num in terms:
        res.append(format(num, fmt_str)[2:]) 
    
    return res



#https://en.wikipedia.org/wiki/Quine%E2%80%93McCluskey_algorithm
def check_dash_align(min1: str, min2: str):
    for i in range(len(min1)): 
        #If one minterm has dashes and the other does not then the minterms cannot be merged. 
        if min1[i] != '-' and min2[i] == '-': 
            return False 
    return True

def check_minterm_diff(min1, min2):
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

def merge_minterms(min1: str, min2: str):
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
            if check_dash_align(minterm1, minterm2) and check_minterm_diff(minterm1, minterm2):
                mergedMinterm = merge_minterms(minterm1, minterm2) 
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


def match_str(p_implicant: str, minterm: str):
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
            if match_str(p_implicant,minterms[j]):
                p_implicant_chart[p_implicant] += "1"
            else: 
                p_implicant_chart[p_implicant] += "0"

    #The prime implicant chart is complete so return the completed chart. 
    return p_implicant_chart

# https://en.wikipedia.org/wiki/Petrick%27s_method
def to_pos(implicant_chart: dict[str, str]):
    pos:dict[int, str] = {}
    for implicant_index,(key,value) in enumerate(implicant_chart.items()):
        for col, bit in enumerate(value):
            item = "X" + str(implicant_index)
            if bit == '1':
                try:
                    pos[col] += item + " + "
                except KeyError:
                    pos[col] = item + " + "

    #removing trailing space + space
    # convert it to a list 
    res = []
    for key, value in pos.items():
        res.append(value[:-3])
    return res

# apply absorption law 
# X + XY = X
def absorption(sop: list[str]):
    i = 0
    j = 1
    moves = 0
    while i < len(sop): 
        while j < len(sop): 
            x = sop[i].split("*")
            y = sop[j].split("*")
            can_absorb = True

            # X + XY = X
            for c in x:
                if c not in y:
                    can_absorb = False

            if can_absorb:
                # X + XY = X
                # just remove the all sop[j] because it is equal to XY 
                sop.remove(sop[j])
                moves += 1
                j += 1
                continue
           
            can_absorb = True
            # XY + X = X
            for c in y:
                if c not in x:
                    can_absorb = False

            # XY + X = X
            if can_absorb:
                moves += 1
                sop.remove(sop[i])
                break

            j += 1

        i += 1
        j = i + 1
    
    # if no absorptions were made we break 
    if moves == 0:
        return sop
    else:
        return absorption(sop)

# recursive method that distributes 2 sums at a time (x + y)(w + z)
# puts the distributed product into a new list 
# everything is distributed once the list len = 1
def to_sop(sop: list[str]):
    if len(sop) == 1:
        return sop

    new_sop = []
    index = 0

    for i in range(0,len(sop) - 1, 2):
        # convert to sets to remove duplicates because x + x = x 
        # we also need to apply absoprtion to reduce the size of the lists 
        # or else it takes exponetially longer
        paren1 = list(set(sop[i].split(" + ")))
        paren1 = list(set(absorption(paren1)))
        paren2 = list(set(sop[i + 1].split(" + ")))
        paren2 = list(set(absorption(paren2)))

        new_sop.append("")
        #distributing each value 
        for val1 in paren1:
            for val2 in paren2: 
                #X * X = X
                if val1 == val2: 
                    new_sop[index] += val1 +  " + "
                else:
                    new_sop[index] += val1 + "*" + val2 + " + "
        #remove trailing " + "
        new_sop[index] = new_sop[index][:-3]
        index += 1

    # since we do 2 at a time if we have an odd number 
    # we need to append it to our list
    if len(sop) % 2 == 1:
        new_sop.append(sop[-1])
    return to_sop(new_sop)




def simplify_sop(sop: list[str]):
    # pos to sop puts everything in a list at the first index 
    # split at the + to get each product
    sop = sop[0].split(" + ")

    #remove duplicates in the products  
    # XYX = XY 
    for index, product in enumerate(sop): 
        temp = "".join(sorted(set(product.split("*")))) 
        sop[index] = temp

    sop = list(set(sop))

    #format it better
    #add * between terms
    res = []
    for prod in sop: 
        new_prod = ""
        for i in range(len(prod)):
            new_prod += prod[i]
            if i != len(prod) - 1 and prod[i + 1] == "X" and i != 0: 
                new_prod += "*"
        
        res.append(new_prod)

    sop = res

    absorption(sop)     

    sop = list(set(sop))

    return sop





def p_implicant_to_bool_expr(p_implicant: str):
    res = ""
    num_terms = 0
    for index, bit in enumerate(p_implicant): 
        if bit == '-':
            continue
        # convert each value to an uppercase letter
        char = chr(65 + index)
        if bit == "0":
            res += char + "'"
        else:
            res += char
        num_terms += 1

    return (res, num_terms)



def get_minimal_expr(sop: list[str], implicant_chart: dict[str, str]):
    fewest_terms = []
    smallest_len = len(sop[0])

    # getting the product(s) with the fewest terms 
    for product in sop: 
        term_count = 0
        for c in product: 
            if c == "X":
                term_count += 1
        if term_count < smallest_len:
            fewest_terms = [product]
            smallest_len = term_count
        elif term_count == smallest_len:
            fewest_terms.append(product)

   
    # expand each product to see which one has the least amount of terms 
   
    # current the prime implicant keys to a list 
    keys = list(implicant_chart.keys())

    smallest_term_size = None
    smallest_expr = []
    
    # converting each item into a boolean expression and finding 
    # the smallest
    for product in fewest_terms: 
        terms = product.split("*")
        current_term_size = 0
        expr = []
        for term in terms:
            # convert each term to its prime implicant 
            # each term is Xindex into the prime implicant list 
            key = int(term[1:])
            p_implicant = keys[key]
            (b_expr, term_size) = p_implicant_to_bool_expr(p_implicant)
            expr.append(b_expr)
            current_term_size += term_size
        
        if smallest_term_size == None or current_term_size < smallest_term_size:
            smallest_term_size = current_term_size 
            smallest_expr = expr

    
    return smallest_expr
       

def get_expr(number_of_vars: int, minterms: list[int], dc: list[int] = []):
    if number_of_vars > MAX_FUNCTION_VARS:
        print(f"Cannot have functions with more than ${MAX_FUNCTION_VARS}")
        return

    terms = minterms + dc
    bin_terms = terms_to_bin_str(terms, number_of_vars)

    p_implicants = prime_implicants(bin_terms)

    # when creating the prime implicant chart we don't use the dont_care terms 
    # so we only care about the min terms 
    chart = create_prime_implicant_chart(p_implicants, bin_terms[:len(minterms)])

    pos = to_pos(chart)

    sop = simplify_sop(to_sop(pos)) 

    mininum = get_minimal_expr(sop, chart)

    expr = ""
    for product in mininum:
        expr += product + " + "

    return expr[:-3]



def get_truth_table(number_of_vars: int, minterms: list[int], dc: list[int] = []):
    # need to do some testing
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

