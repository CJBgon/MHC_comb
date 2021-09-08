#!/usr/bin/env python
import json
import sys
import argparse
import itertools
import pandas as pd
import numpy as np
import warnings
from collections.abc import Iterable
from collections import defaultdict
from itertools import chain, cycle
from functools import reduce


def Parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--arcashla", dest="jfile", type=str, default=None,
                        help="Path to .JSON output by arcasHLA.")
    parser.add_argument("-x", "--xhla", dest="xhla", type=str, default=None,
                        help="Path to xHLA tsv output.")
    parser.add_argument("-s", "--seq2hlai", dest="seqi", type=str, default=None,
                        help="path to MHC-I seq2HLA output: ClassI.HLAgenotype4digits.")
    parser.add_argument("-q", "--seq2hlaii", dest="seqii", type=str, default=None,
                        help="path to MHC-II seq2HLA output: ClassII.HLAgenotype4digits.")
    parser.add_argument("-p", "--polysolver", dest="polys", type=str, default=None,
                        help="path to polysolver winners.txt file.")
    parser.add_argument("-r", "--resolution", dest="res", type=int, default=3,
                        help="At what resolution should the HLA types be returned? default = 3"
                             ": HLA-A:02:01:01")
    parser.add_argument('-w', '--weight', nargs='*', help='Set weights to determine consensus',
                        default='arcasHLA Polysolver seq2HLA xHLA')
    parser.add_argument("-o", "--object", dest="type", type=str, default="both",
                        help="which HLA types to return: 'I', 'II' or 'both'")
    parser.add_argument("-d", "--sep", dest="sep", type=str, default=",",
                        help="What separator to use for output. default= ',' ")
    parser.add_argument("-l", "--log", dest="log", default=False, action='store_true',
                        help="boolean if a log file should be created. DEFAULT = False")
    parser.add_argument("--lohhla", dest="lohhla", default=False, action='store_true',
                        help="boolean if a log file should be created. DEFAULT = False")
    parser.add_argument("-f", "--logpath", dest="logpath", type=str, default=None,
                        help="directory to output logs and results in.")


    Options=parser.parse_args()
    if (
            not Options.jfile and
            not Options.xhla  and
            not Options.seqi  and
            not Options.polys and
            not Options.seqii
        ):
        parser.error("no input file provided.")

    if Options.jfile:
        if Options.jfile[-4:] != "json":
            parser.error("Invalid filetype for -a, please provide a .json.")

    if Options.log or Options.lohhla:
        if not Options.logpath:
            parser.error('No log file location provided with -f or --lohhla')

    return Options


def custom_formatwarning(msg, *a, **b):
    return str(msg) + '\n'


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
            data.update({keys: [':'.join(data[keys][0].split(":")[0:Options.res])]})

    return data


def seqproc(seqone=None, seqtwo=None):
    # as either one or two files can be provided by seq2HLA, we set the variables to None by default.
    ridict = {}
    if seqone:
        ri = pd.read_csv(seqone, sep='\t', header=0, usecols=[0, 1, 3])
        for allele in ri.iloc[:, 0]:
            ridict[allele] = non_recursive_flatten(
                ri.loc[ri['#Locus'] == allele, ['Allele 1', 'Allele 2']].values.tolist()
            )
    if seqtwo:
        rii = pd.read_csv(seqtwo, sep='\t', header=0, usecols=[0, 1, 3])
        for allele in rii.iloc[:, 0]:
            ridict[allele] = non_recursive_flatten(
                rii.loc[rii['#Locus'] == allele, ['Allele 1', 'Allele 2']].values.tolist()
            )

    return ridict


def xhlaproc(tsvfile):
    data = pd.read_csv(tsvfile, sep="\t")
    alleles = data.loc[:, 'full']
    split_type = [y.split('*')[0] for y in alleles]
    keys = set([y.split('*')[0] for y in alleles])

    xhladict = {}
    for x in keys:
        l = []
        for i, val in enumerate(split_type):
            if x == val:
                l.append(alleles[i])
        xhladict[x] = l

    return xhladict


def polyproc(polyfile):
    data = pd.read_csv(polyfile, sep="\t", header=None)
    a = data.loc[0].tolist()[1:]
    b = data.loc[1].tolist()[1:]
    c = data.loc[2].tolist()[1:]
    a1 = []
    b1 = []
    c1 = []
    polydict = {}
    for x in a:
        a1.append("A*" + ":".join(x.split('_')[2:]))
        polydict["A"] = a1
    for x in b:
        b1.append("B*" + ":".join(x.split('_')[2:]))
        polydict["B"] = b1
    for x in c:
        c1.append("C*" + ":".join(x.split('_')[2:]))
        polydict["C"] = c1

    return polydict


def breaktie(tie, haplo, dictlist, we_filt):
    # for the non-consencus haplotype:
    # iterate through the weights, if it's in the matching dictlist add the score based on the iteration.
    resdict = {}
    for t in tie:
        score = 0
        for a, b in we_filt.items():
            if t in dictlist[a][haplo]:
                score += 1/b
        resdict[t] = score
    # this does require python 3.6+
    resdict_sorted = {k: v for k, v in sorted(resdict.items(), key=lambda item: item[1], reverse=True)}
    return resdict_sorted


def top_n(d, n):
    # return the top value from our list
    dct = defaultdict(list)
    for k, v in d.items():
        dct[v].append(k)
    # return sorted(dct.values())[-n:][::-1]  
    return [dct[s] for s in sorted(dct)[::-1][:n]]


def consensus(di, we, log, logf):
    # weight = dict(zip(['1', '0.75', '0.5', '0.25'], we))
    # what haplotypes / keys do we have available?
    # keys = [k.keys() for k in dictlist]
    # a = [x + ['a'] for x in a]

    # remove typers with no input.
    u = []
    for d in di:
        if di[d] is None:
            u.append(d)
    for n in u:
        di.pop(n)
        if n in we:
            raise ValueError('Not all assigned weights have matched input.')
            input()
            sys.exit(1)
        else:  # weights are assigned by the user. as such they aren't pre-defined and dont need pruning.
            val = range(1, len(we)+1, 1)
            weight = dict(zip(we, val))
    if len(u) == 0:
        val = range(1, len(we) + 1, 1)
        weight = dict(zip(we, val))

    # List available haplotypes:
    keys = [list(di[d].keys()) for d in di]
    setofkeys = set(non_recursive_flatten(keys))
    con = {}

    # if either the typer with the highest assigned weight, or two or more typers only have one allele
    # we are most likely looking at a homozygous patient. therefore, pick only one allele.
    # Likely, this requires some testing and optimisation.


    for d in setofkeys:
        usedict = {}
        usewe = {}
        noweight = []
        for a in di:
            if d in di[a].keys():
                usedict[a] = di[a]
            else:
                noweight.append(a)
        for k, value in weight.items():
            if k not in noweight:
                usewe[k] = value

        check = []
        if d in ['A', 'B', 'C', 'DOA', 'DMA', 'DPA1', 'DPB1', 'DOB', 'DMB',
                 'DRA', 'DRA1', 'DRA2', 'DQA', 'DQB', 'DRB', 'DQA1', 'DQA2', 'DQB1', 'DQB2', 'DRB1', 'DRB3', 'DRB5']:
            for a in usedict:
                if isinstance(usedict[a][d], list):
                    check.append([x for x in usedict[a][d]])
                else:
                    check.append([usedict[a][d]])
            check2 = set(non_recursive_flatten(check))

            # determine homozygousity or allele options:
            hom = 0
            hom2_check = False
            hom1_check = False
            hom2_allele = []
            homcount = []
            # xHLA en arcasHLA return one in the case of homo. poly and seq2HLA return two. adjust at the conform.
            for key, value in usewe.items():
                if not isinstance(usedict[key][d], list):
                    break  # in case only a single value was given.
                if len(usedict[key][d]) < 2:
                    break
                if value == 1:
                    if usedict[key][d][0] == usedict[key][d][1]:
                        hom1_check = True
                        hom_allele = usedict[key][d][0]
                if usedict[key][d][0] == usedict[key][d][1]:  # cant use set() because that would count 1 input as homo.
                    hom += 1
                    hom2_allele.append(usedict[key][d][0])
            if hom > 1 and len(set(hom2_allele)) == 1:  # in the case of 2 homozygous calls, 1 allele
                hom2_check = True
                tophom = list(set(hom2_allele))
            elif hom > 1 and len(set(hom2_allele)) > 1:  # in the case of 2 homozygous calls, different alleles.
                hom2_check = True
                for i in set(hom2_allele):
                    homcount.append(non_recursive_flatten(check).count(i))
                hcounts = dict(zip(set(hom2_allele), homcount))
                tiehom = top_n(hcounts, 1)[0]
                tophom = list(breaktie(tiehom, d, usedict, usewe).keys())[:1]

            # not all typers return all HLA-II haplotypes. if there is only one option in all typers dont just
            # merge them, as we can't be sure they are aimed at the maternal / paternal genes.
            # determine how many top hits to return based on the input:
            ret_int = 1
            s_check = []
            # using a roundabout way just in case one typer does not return a list of options, rather one string.
            # which would return a ret_int the the length of that string, not the max amount of options returned.
            for c in usedict:
                # had to add the len(check2) as homozygous calls in typers would return {[x1, x1]} instead of
                # {x1}, which would lead to an list index out of range.
                if isinstance(usedict[c][d], list) and len(check2) > 1:
                    s_check.append(len(usedict[c][d]))
                    ret_int = max(s_check)

            # s_check = [len(usedict[c][d]) for c in usedict]  # count how many items are in each dict

            if ret_int > 2:
                print(
                    'Too many options per haplotype, perhaps the typer has grouped HLA-DRB1 and HLA-DRB2 under HLA-DRB?'
                )
                raise ValueError('more than 2 options per HLA haplotype.')
                input()
                sys.exit(1)

            if hom1_check:
                con[d] = [hom_allele]
            elif hom2_check and not hom1_check:
                con[d] = [tophom]
            elif len(check2) == ret_int:  # check if there are only 2 haplotypes given ( complete consensus)
                con[d] = list(check2)
            else:  # count what the top detected haplotypes are:
                cnts = []
                for i in check2:
                    cnts.append(non_recursive_flatten(check).count(i))
                counts = dict(zip(check2, cnts))
                topcount = top_n(counts, ret_int)   # top candidates based on counts
                # check if there is a consensus based on counts:
                if len(topcount[0]) == ret_int:  # a tie for the top two counts
                    con[d] = topcount[0]
                elif len(topcount[0]) == 1 and len(topcount[1]) == 1 and ret_int == 2:  # no tie in the top two counts
                    con[d] = [topcount[0][0], topcount[1][0]]
                elif len(topcount[0]) > ret_int:  # break tie based on assigned weights.
                    tie = topcount[0]
                    con[d] = list(breaktie(tie, d, usedict, usewe).keys())[:ret_int]
                elif len(topcount[1]) > 1:  # this will only occur if ret_int > 1.
                    tie = topcount[1]
                    con[d] = [topcount[0][0], list(breaktie(tie, d, usedict, usewe).keys())[0]]
                else:
                    print("Could not reach consensus on haplotypes.")

            if log:
                # check for each typer if the result corresponds with the input.
                # there are three instances we are interested in:
                # concensus: both typer and con are exactly the same (regardless of order)
                # type I error : either the typer or the con has an homozygous call while the other does not.
                # type II error: the typer or con has got one or two mismatches, but the amount of alleles is the same.

                for k in usedict.keys():
                    cons = True
                    typei = 0
                    typeii = 0
                    if isinstance(usedict[key][d], list):
                        haplolist = usedict[key][d]
                    else:
                        haplolist = [usedict[key][d]]

                    if sorted(haplolist) == sorted(con[d]):
                        cons = True
                    else:
                        cons = False
                        if  len(haplolist) != len(con[d]):
                            typei = len(haplolist) - len(con[d])

                        for i in haplolist:
                            if i not in con[d]:
                                typeii +=1

                    with open(logf, 'a+') as lo:
                        lo.write('\t'.join(
                            [str(k), str(d), str(cons), str(typei), str(typeii)]
                        ) + '\n')

        else:
            warnings.formatwarning = custom_formatwarning
            warnings.warn(str(d + ': Not recognized key-type for haplotype.'))

    return con


def reform(dicto, res):
    loci = ['A', 'B', 'C', 'E', 'G', 'F', 'H', 'K', 'L',
            'DRA', 'DMA', 'DMB', 'DOA', 'DOB',
            'DRB1', 'DRB3', 'DRB4', 'DRB5', 'DQA1',
            'DQB1', 'DQB2',  'DPA1', 'DPB1']
    copydict = dicto
    iterkey = list(copydict.keys())
    for x, y in copydict.items():
        for pos, val in enumerate(y):
            clean = ':'.join(val.split(":")[0:res])
            y[pos] = clean.strip("\'")  # remove trailing " ' "

    # remove 'HLA' or 'MHC' from the keys:
    for key in iterkey:
        if key.startswith('HLA-') or key.startswith('MHC-'):
            newkey = key[4:]
            copydict[newkey] = copydict[key]
            del copydict[key]

    # check if the keys are in our key list:
    # we need to unify the keys as we want to combine different typers per haplotype.
    for key in iterkey:
        if key not in loci:
            num = [x.split('*')[0] for x in copydict[key]]
            for loc, name in enumerate(num):
                if name in copydict:
                    copydict[name].append(copydict[key][loc])
                else:
                    copydict[name] = [copydict[key][loc]]
            del copydict[key]

    iterkey = list(copydict.keys())
    for key in iterkey:
        if key in ['A', 'B', 'C'] and len(copydict[key]) == 1:
            copydict[key].append(copydict[key][0])

    return copydict


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
                'DQ': ['DQA1', 'DQB1', 'DQB2'],
                'DR': ['DRA', 'DRB1', 'DRB3']}
    ls = []
    for i, j in pairings.items():
        A = [j[0]]
        B = list(chain.from_iterable([j[1:]]))
        comb1 = list(non_recursive_flatten([mhcdict[x] for x in A if x in mhcdict]))
        comb2 = list(non_recursive_flatten([mhcdict[x] for x in B if x in mhcdict]))
        list(itertools.chain(comb2))
        for m in comb1:
            ls.append(m)
        for m in comb2:
            ls.append(m)
        # dont use, we need combinations not permutations.
        # combslist = list(list(zip(r, p)) for (r, p) in list(zip(repeat(comb1), itertools.permutations(comb2))))
        combslist = [(x, y) for x in comb1 for y in comb2]
        for allele in combslist:
            a = list(non_recursive_flatten(allele))
            ls.append("-".join([x for x in a]))

    return ls


def fileout(mhcdict, path):
    ndict = {}
    for k in ['A', 'B', 'C']:
        usedict = mhcdict[k]  # change to k here
        if len(usedict) < 2:  # output requires 2 HLA alleles even in case of homozygosity
            iterc = cycle(usedict)
            for _ in range(1):
                usedict.append(next(iterc))
        # now bring it back to the format we need:
        repls = ('*', '_'), (':', '_')
        check = []
        for x in usedict:
            check.append(reduce(lambda a, kv: a.replace(*kv), repls, x))
        nkey='HLA-'+k
        retcheck = ['hla_'+r.lower() for r in check]
        ndict[nkey] = retcheck
    mhcdf = pd.DataFrame(dict([(r, pd.Series(i)) for r, i in ndict.items()]))
    mhcdft = mhcdf.T
    mhcdft.to_csv(path_or_buf=path, sep='\t', index=True, header=False)




def main():
    Options = Parser()
    # get all dicts of available options.
    if Options.jfile:
        arcasraw = jsonproc(Options)
        arcasdict = reform(dicto=arcasraw, res=Options.res)
    else:
        arcasdict = None

    if Options.xhla:
        xhlaraw = xhlaproc(tsvfile=Options.xhla)
        xhladict = reform(dicto=xhlaraw, res=Options.res)
    else:
        xhladict = None

    if Options.seqi or Options.seqii:
        seq2raw = seqproc(seqone=Options.seqi, seqtwo=Options.seqii)
        seq2dict = reform(dicto=seq2raw, res=Options.res)
    else:
        seq2dict = None

    if Options.polys:
        polyraw = polyproc(polyfile=Options.polys)
        polydict = reform(dicto=polyraw, res=Options.res)
    else:
        polydict = None

    dictlist = {'arcasHLA': arcasdict, 'xHLA': xhladict, 'Polysolver': polydict, 'seq2HLA': seq2dict}
    if Options.log:
        if Options.logpath[-1] != '/':
            loguse = Options.logpath + '/MHC_comb_log.txt'
        else:
            loguse = Options.logpath + 'MHC_comb_log.txt'

        condict = consensus(di=dictlist, we=Options.weight, log=Options.log, logf=loguse)
    else:
        condict = consensus(di=dictlist, we=Options.weight, log=False, logf=None)

    mii = {'DOA', 'DMA', 'DPA1', 'DPB1', 'DOB', 'DMB',
           'DRA', 'DQA1', 'DQB1', 'DQB2', 'DRB1', 'DRB3', 'DRB5'}

    mi = {'A', 'B', 'C'}

    if Options.type == "both":
        mhci = {k: condict[k] for k in condict.keys() & mi}
        mhcfi = non_recursive_flatten(list(mhci.values()))
        mhcfin = ['HLA-'+t for t in mhcfi]
        mhciir = {k: condict[k] for k in condict.keys() & mii}
        mhcboth = non_recursive_flatten(list(
            gencom(mhciir)
        )) + mhcfin
        out = mhcboth

    elif Options.type == 'I':
        mhci = {k: condict[k] for k in condict.keys() & mi}
        mhcfi = non_recursive_flatten(list(mhci.values()))
        out = ['HLA-'+t for t in mhcfi]

    elif Options.type == 'II':
        mhciir = {k: condict[k] for k in condict.keys() & mii}
        out = non_recursive_flatten(list(
            gencom(mhciir)
        ))

    else:
        raise ValueError("MHC allele option not recognized: 'I', 'II' or 'both.")

    # we need condict for this, therefore we need to run it later.
    if Options.lohhla:
        if Options.logpath[-1] != '/':
            loghla = Options.logpath + '/mhccomb.winners.hla.txt'
        else:
            loghla = Options.logpath + 'mhccomb.winners.hla.txt'

        fileout(mhcdict=condict, path=loghla)


    returnStr = ''
    for item in out:
        returnStr += str(item) + Options.sep
    print(returnStr)


if __name__ == "__main__":
    main()
