"""Bio Hacking"""
import re
from collections import OrderedDict

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
    "uaa": "*",
    "uga": "*",
    "uag": "*",
}


def transcription(fp: str) -> str:
    """Transcription. DNA to RNA."""
    with open(fp, "r") as f:
        return re.sub(r"\s+", "", f.read()).replace("t", "u")


def translation(rna: str) -> str:
    """Translation. RNA to Protein.

    :param rna: RNA string.
    :returns: Protein string.
    """

    def rf(s: str):
        """Reading frame. https://en.wikipedia.org/wiki/Reading_frame."""
        for i in range(0, len(s), 3):
            yield s[i : i + 3]

    res = ""
    codons = rf(rna)
    while True:
        try:
            codon = next(codons)
            for k, v in genetic_code.items():
                if codon == k:
                    res += v
        except StopIteration:
            break

    return "".join(list(filter(None, res.split("*"))))


def sars_cov2(fp: str) -> OrderedDict:
    """SARS-CoV-2 genome organization.
        https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2

    :param fp: File path.
    :returns: OrderedDict of all genes coding for protein products.
    """
    sc = OrderedDict()
    rna = transcription(fp)

    # gene="ORF1ab"
    # product="ORF1ab polyprotein"
    orf1a = rna[265:13468]
    orf1b = rna[13467:21555]
    sc["orf1ab"] = translation(orf1a + orf1b)
    # gene="S"
    # gene_synonym="spike glycoprotein"
    # product="surface glycoprotein"
    sc["s"] = translation(rna[21562:25384])
    # gene="ORF3a"
    # product="ORF3a protein"
    sc["orf3a"] = translation(rna[25392:26220])
    # gene="E"
    # product="envelope protein"
    sc["e"] = translation(rna[26244:26472])
    # gene="M"
    # product="membrane glycoprotein"
    sc["m"] = translation(rna[26522:27191])
    # gene="ORF6"
    # product="ORF6 protein"
    sc["orf6"] = translation(rna[27201:27387])
    # gene="ORF7a"
    # product="ORF7a protein"
    sc["orf7a"] = translation(rna[27393:27759])
    # gene="ORF8"
    # product="ORF8 protein"
    sc["orf8"] = translation(rna[27893:28259])
    # gene="N"
    # product="nucleocapsid phosphoprotein"
    sc["n"] = translation(rna[28273:29533])
    # gene="ORF10"
    # product="ORF10 protein"
    sc["orf10"] = translation(rna[29557:29674])

    return sc


if __name__ == "__main__":
    sc = sars_cov2("data/corona.txt")

    print("SARS-CoV-2 genome organization")
    for k, v in sc.items():
        print(f"{k}:\n{v}\n")
