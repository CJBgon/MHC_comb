#!/usr/bin/env python
import json
import sys
import argparse
import itertools
from itertools import chain
from collections import Iterable


def Parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--json", dest="jfile", type=str, default=None,
                        help="Path to .JSON output by arcasHLA.")
    parser.add_argument("-r", "--resolution", dest="res", type=int, default=3,
                        help="At what resolution should the HLA types be returned? default = 3"
                             ": HLA-A:02:01:01")
    parser.add_argument("-o", "--object", dest="type", type=str, default="both",
                        help="which HLA types to return: 'I', 'II' or 'both'")
    parser.add_argument("-s", "--sep", dest="sep", type=str, default=",",
                        help="What separator to use for output. default= ',' ")
    Options = parser.parse_args()
    if not Options.jfile:
        parser.error("no input file provided.")
    if Options.jfile:
        if Options.jfile[-4:] != "json":
            parser.error("Invalid filetype for -j, please provide a .json.")

    return Options

def jsonproc(Options):
    # read jason file into dict:
    with open(Options.jfile) as jsonfile:
        data = json.load(jsonfile)
    dictkeys = data.keys()
    for keys in dictkeys:
        if len(data[keys]) == 2:
            data.update({keys: [':'.join(data[keys][0].split(":")[0:Options.res]),
                                ':'.join(data[keys][1].split(":")[0:Options.res])]})
        elif len(data[keys]) == 1:
            data.update({keys: ':'.join(data[keys][0].split(":")[0:Options.res])})

    return data



def reform(alleles, res):
    all = list(alleles)
    for i, j in enumerate(all):
        all[i] = ':'.join(j.split(":")[0:res])
    ci = all[0:6]
    for i, j in enumerate(ci):
        if j[:4] != "HLA-":
            all[i] = ("HLA-"+j)
    return all


def flatten(lis):
    for item in lis:
        if isinstance(item, Iterable) and not isinstance(item, str):
            for x in flatten(item):
                yield x
        else:
            yield item


def non_recursive_flatten(xs):
    res = []

    def loop(ys):
        for i in ys:
            if isinstance(i, list):
                loop(i)
            else:
                res.append(i)
    loop(xs)
    return res



def gencom(mhcdict):
    # get the MHCII dictionary
    pairings = {'DM': ['DMA', 'DMB'],
                'DO': ['DOA', 'DOB'],
                'DP': ['DPA1', 'DPB1'],
                'DQ': ['DQA1', 'DQB1'],
                'DR': ['DRA', 'DRB1', 'DRB3']}
    ls = []
    for i, j in pairings.items():
        A = [j[0]]
        B = list(chain.from_iterable([j[1:]]))
        comb1 = [mhcdict[x] for x in A]
        comb2 = list(non_recursive_flatten([mhcdict[x] for x in B]))
        list(itertools.chain(comb2))
        for m in comb1:
            ls.append(m)
        for m in comb2:
            ls.append(m)
        # dont use, we need combinations not permutations.
        # combslist = list(list(zip(r, p)) for (r, p) in list(zip(repeat(comb1), itertools.permutations(comb2))))
        combslist = [(x,y) for x in comb1 for y in comb2]
        for allele in combslist:

            a = list(non_recursive_flatten(allele))
            ls.append("-".join([x for x in a]))

    return ls


def main():
    Options = Parser()
    allele_dict = jsonproc(Options)
    if Options.type == "both":
        mhcir = non_recursive_flatten(list(allele_dict.values())[0:3])
        mhci = reform(alleles=mhcir, res = Options.res)
        mhcii = gencom(allele_dict)
        out = mhci + mhcii
    elif Options.type == "I":
        mhcir = non_recursive_flatten(list(allele_dict.values())[0:3])
        out = reform(alleles=mhcir, res=Options.res)
    elif Options.type == "II":
        out = gencom(allele_dict)
    else:
        error("MHC allele option not recognized: 'I', 'II' or 'both.")


    returnStr = ''
    for item in out:
        returnStr += str(item) + Options.sep
    print(returnStr)


if __name__ == "__main__":
    main()
