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
    raise NotImplementedError

def bwt_encode(s):
    """Burrows-Wheeler Transform
    Args: s, string, which must not contain '{' or '}'
    Returns: BWT(s), which contains '{' and '}'
    """
    raise NotImplementedError

def bwt_decode(bwt):
    """Inverse Burrows-Wheeler Transform
    Args: bwt, BWT'ed string, which should contain '{' and '}'
    Returns: reconstructed original string s, must not contains '{' or '}'
    """
    raise NotImplementedError

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
                    "scottytartanscottytartanscottytartanscottytartan"]
    for s in test_strings:
        test_string(s)
