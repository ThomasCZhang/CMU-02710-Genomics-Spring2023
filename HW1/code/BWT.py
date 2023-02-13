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
        idx2 = 1
        while (idx+idx2) < len(s):
            if s[idx] != s[idx+idx2]:
                break
            idx2 += 1
        if idx2 != 1:
            final_str = "".join([final_str, curr_char, curr_char, str(idx2)])
        else:
            final_str = "".join([final_str, curr_char])
        idx += idx2

    return final_str

def bwt_encode(s: str):
    """Burrows-Wheeler Transform
    Args: s, string, which must not contain '{' or '}'
    Returns: BWT(s), which contains '{' and '}'
    """
    s = "".join(["{", s, "}"])
    suffix_array = [s[i:] for i in range(len(s))]
    suffix_array.sort()
    bwt_s = ""
    for suffix in suffix_array:
        idx = len(s)-len(suffix)
        if idx == 0:
            bwt_s = "".join([bwt_s, "}"])
        else:
            bwt_s = "".join([bwt_s, s[idx-1]])
    return bwt_s

def bwt_decode(bwt):
    """Inverse Burrows-Wheeler Transform
    Args: bwt, BWT'ed string, which should contain '{' and '}'
    Returns: reconstructed original string s, must not contains '{' or '}'
    """
    counts = {}
    for key in bwt:
        if key in counts:
            counts[key] += 1
        else:
            counts[key] = 1
    
    rank = {}
    count = 0
    for key in sorted(list(counts.keys())):
        rank[key] = count
        count += counts[key]
        counts[key] = 0

    idx_linked_list = [0 for i in range(len(bwt))]
    for idx0, val in enumerate(bwt):
        idx_linked_list[idx0] = rank[val] + counts[val]
        counts[val] += 1

    curr_idx = rank["}"]
    decoded_string = ""
    while bwt[curr_idx] != "{":
        decoded_string = "".join([bwt[curr_idx], decoded_string])
        curr_idx = idx_linked_list[curr_idx]

    return decoded_string
    
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
    test_strings = [
                    "WOOOOOHOOOOHOOOO!",
                    "scottytartanscottytartanscottytartanscottytartan",
                    "BANANA",
                    "steelerssteelerssteelerssteelers"
    ]
    for s in test_strings:
        test_string(s)
