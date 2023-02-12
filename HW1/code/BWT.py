# BWT.py
# HW1, Computational Genomics, Spring 2022
# andrewid:

# WARNING: Do not change the file name, or the function signatures below.
# Autograder expects these names exactly.

def rle(s):
    """Run Length Encoder
    Args: s, string to be encoded
    Returns: RLE(s)
    """
    final_str = ""
    idx = 0
    while idx < len(s):
        curr_char = s[idx]
        final_str = "".join([final_str, curr_char])
        idx2 = 1
        while (idx+idx2) < len(s):
            if s[idx] != s[idx+idx2]:
                break
            idx2 += 1

        if idx2 > 2:
            final_str = "".join([final_str, curr_char, str(idx2)])
            idx += idx2
        else:
            idx += 1

    return final_str

def bwt_encode(s: str):
    """Burrows-Wheeler Transform
    Args: s, string, which must not contain '{' or '}'
    Returns: BWT(s), which contains '{' and '}'
    """
    s = "".join(["{", s, "}"])
    temp_dict = {} # Using this dict to aid with sorting.
    for i in range(len(s)-1, 0, -1): # We're going to index from the end. This will act as the sorting part too.
        key = s[i]
        val = s[i-1]
        if key in temp_dict:
            temp_dict[key] = "".join([temp_dict[key], val])
        else:
            temp_dict[key] = val
    key = s[0]
    val = "}"
    if key in temp_dict:
        temp_dict[key] = "".join([temp_dict[key], val])
    else:
        temp_dict[key] = val
    
    # first_char = temp_dict["}"]
    # del temp_dict["}"]
    sorted_keys = sorted(list(temp_dict.keys()))
    bwt_s = "".join((temp_dict[key] for key in sorted_keys))
    # bwt_s = "".join((first_char, bwt_s))
    return bwt_s

def bwt_decode(bwt):
    """Inverse Burrows-Wheeler Transform
    Args: bwt, BWT'ed string, which should contain '{' and '}'
    Returns: reconstructed original string s, must not contains '{' or '}'
    """
    rank_count_dict = {}
    for x in bwt:
        if x in rank_count_dict:
            rank_count_dict[x][1] += 1
        else:
            rank_count_dict[x] = [0, 1]
    
    # rank_count_dict["}"][0] = 0
    sorted_keys = sorted(list(rank_count_dict.keys()))
    print(sorted_keys)
    rank = 0
    for key in sorted_keys:
        # if key != "}":
        rank_count_dict[key][0] = rank
        rank += rank_count_dict[key][1]
        rank_count_dict[key][1] = 0 # Reset counts to 0, we are reusing this position to reconstruct string
    
    final_string = ""
    next_char = bwt[rank_count_dict["}"][0]]
    while next_char != "{":
        final_string = "".join([next_char, final_string])

        idx = rank_count_dict[next_char][0] + rank_count_dict[next_char][1]
        rank_count_dict[next_char][1] += 1
        next_char = bwt[idx]

    return final_string

def test_string(s):
    compressed = rle(s)
    bwt = bwt_encode(s)
    compressed_bwt = rle(bwt)
    reconstructed = bwt_decode(bwt)
    template = "{:25} ({:3d}) {}"
    print(template.format("original", len(s), s))
    print(template.format("bwt_enc(orig)", len(bwt), bwt))
    print(template.format("bwt_dec(bwt_enc(orig))", len(reconstructed), reconstructed))
    print(template.format("rle(orig)", len(compressed), compressed))
    print(template.format("rle(bwt_enc(orig))", len(compressed_bwt), compressed_bwt))
    print()
    print()

if __name__ == "__main__":
    # Add more of your own strings to explore for question (i)
    test_strings = ["WOOOOOHOOOOHOOOO!",
                    "scottytartanscottytartanscottytartanscottytartan",
                    "BANANA",
                    "steelerssteelerssteelerssteelers"]
    for s in test_strings:
        test_string(s)
