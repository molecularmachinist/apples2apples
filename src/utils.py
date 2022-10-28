from typing import List


def ordinal_str(i: int) -> str:
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


def compressed_resid_list(resids: List[str]) -> str:
    if (len(resids) == 0):
        return ""
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
