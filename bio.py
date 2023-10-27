"""Bio Hacking"""
import re

genetic_code = {
    "gcu": "A",
    "gcc": "A",
    "gca": "A",
    "gcg": "A",
    "auu": "I",
    "auc": "I",
    "aua": "I",
    "cgu": "R",
    "cgc": "R",
    "cga": "R",
    "cgg": "R",
    "aga": "R",
    "agg": "R",
    "cuu": "L",
    "cuc": "L",
    "cua": "L",
    "cug": "L",
    "uua": "L",
    "uug": "L",
    "aau": "N",
    "aac": "N",
    "aaa": "K",
    "aag": "K",
    "gau": "D",
    "gac": "D",
    "aug": "M",
    "uuu": "F",
    "uuc": "F",
    "ugu": "C",
    "ugc": "C",
    "ccu": "P",
    "ccc": "P",
    "cca": "P",
    "ccg": "P",
    "caa": "Q",
    "cag": "Q",
    "ucu": "S",
    "ucc": "S",
    "uca": "S",
    "ucg": "S",
    "agu": "S",
    "agc": "S",
    "gaa": "E",
    "gag": "E",
    "acu": "T",
    "acc": "T",
    "aca": "T",
    "acg": "T",
    "ugg": "W",
    "ggu": "G",
    "ggc": "G",
    "gga": "G",
    "ggg": "G",
    "uau": "Y",
    "uac": "Y",
    "cau": "H",
    "cac": "H",
    "guu": "V",
    "guc": "V",
    "gua": "V",
    "gug": "V",
    "uaa": "$",
    "uga": "$",
    "uag": "$",
}


def dec(fp: str) -> str:
    """Transcription and Translation.

    :param fp: File path.
    """

    def rf(s: str):
        for i in range(0, len(s), 3):
            yield s[i : i + 3]

    with open(fp, "r") as f:
        s = re.sub(r"\s+", "", f.read()).replace("t", "u")

        res = ""
        codons = rf(s)
        while True:
            try:
                codon = next(codons)
                for k, v in genetic_code.items():
                    if codon == k:
                        res += v

            except StopIteration:
                break

    return res.split("$")


if __name__ == "__main__":
    p = dec("data/corona.txt")

    p = "".join(list(filter(None, p)))
    print(p)
