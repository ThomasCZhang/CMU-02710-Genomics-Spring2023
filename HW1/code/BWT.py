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
    matrix = ["" for i in range(len(bwt))]
    for i in range(len(bwt)):
        for j in range(len(bwt)):
            matrix[j] = "".join([bwt[j], matrix[j]])
        matrix.sort()
    
    for string in matrix:
        if string[len(bwt)-1] == "}":
            return string[1:len(bwt)-1]

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
