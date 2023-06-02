from typing import List


def ordinal_str(i: int) -> str:
    """
    Return the integer parameter i as a string with the correct ordinal ending.
    e.g.
    1  -> "1st"
    2  -> "2nd"
    5  -> "5th"
    11 -> "11th"
    21 -> "21st"
    """
    last_dig = i % 10
    second_to_last = (i//10) % 10
    if (second_to_last != 1):
        if (last_dig == 1):
            return f"{i}st"
        elif (last_dig == 2):
            return f"{i}nd"
        elif (last_dig == 3):
            return f"{i}rd"
    return f"{i}th"


def compressed_resid_list(resids: List[int]) -> str:
    """
    Takes a list of ints and returns the compressed ranges. Does not sort the input.
    e.g.
    [1,2,3,4,6,8,9,10] -> " 1:4 6 8:10"
    """
    if (len(resids) == 0):
        return ""
    if (len(resids) == 1):
        return f" {resids[0]}"
    out = ""
    begin_cont = resids[0]
    for iprev, r in enumerate(resids[1:]):
        if (r == resids[iprev]+1):
            continue
        if (resids[iprev] == begin_cont):
            out += f" {begin_cont}"
        else:
            out += f" {begin_cont}:{resids[iprev]}"
        begin_cont = r
    if (r == begin_cont):
        out += f" {begin_cont}"
    else:
        out += f" {begin_cont}:{r}"
    return out
