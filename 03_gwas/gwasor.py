#!/usr/bin/env python3
import pandas as pd
import argparse
from subprocess import Popen, PIPE
from os.path import basename
import os

# ---------------------------------------------------------------------------
# External tool paths (override via environment variables, see config.env)
# ---------------------------------------------------------------------------
PLINK     = os.environ.get('PLINK',     'plink')
EMMAX     = os.environ.get('EMMAX',     'emmax-intel64')
EMMAX_KIN = os.environ.get('EMMAX_KIN', 'emmax-kin-intel64')
FLASHPCA  = os.environ.get('FLASHPCA',  'flashpca2')
QQNORM_R  = os.environ.get('QQNORM_R',  'qqnorm_pheno.R')
QQMAN_R   = os.environ.get('QQMAN_R',
             os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'qqman4Pop_lecttuce.R'))

def makedir(args):
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    print(f"make output directory: {args.outdir}")


def keep(args):
    print( "keep pheno that in pop family list")
    famfile = args.bfile + ".tfam"
    fam = pd.read_csv(famfile,sep=" ",header=None,index_col=0,dtype=str)
    raw_pheno = pd.read_csv(args.pheno,header=None,sep="\t",index_col=0,dtype=str)
    raw_pheno.columns = ["trait"]
    pheno = pd.concat([fam,raw_pheno],axis=1).dropna()
    outPheno = args.outdir + "/" +   "keep_" + basename(args.pheno)
    pheno[[1, "trait"]].to_csv(outPheno,header=None,sep="\t",index=None)
    outkeep =  args.outdir + "/" +  "keep_" + basename(args.pheno) + ".list"
    pheno[1].to_csv(outkeep,header=None,sep="\t",index=1)
    print(f"result file {outPheno, outkeep}")
    return outPheno, outkeep

def qqnorm(args,outPheno):
    if args.subcommand == 'fine_gwas':
        outqqnorm = args.pheno
    else:
        outqqnorm =  args.outdir + "/" + basename(outPheno) + ".qqnorm"
        cmd = QQNORM_R + ' ' + outPheno + " " +  outqqnorm
        print(cmd)
        qqnorm=Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
        qq_stdout,qq_stderr= qqnorm.communicate()
        if qqnorm.returncode==0:
            print ('*************** qqnorm of the phenottype done  ***************')
        else:
            print (qq_stderr)
    return outqqnorm

def final_emmax_pheno(args,outPheno,outnorm):
    if args.bfile:
        famfile = args.bfile + ".tfam"
    else:
        famfile = args.pheno + ".tfam"
    if args.binary == "yes":
        pheno = outPheno
    else:
        pheno = outnorm
    fam = pd.read_csv(famfile, sep=" ", header=None,dtype=str)
    fam.index=fam[0]
    trait = pd.read_csv(pheno, header=None, sep="\t",dtype=str)
    trait.columns = ["Taxa", "P"]
    trait.index = trait["Taxa"]
    outdf = pd.concat([fam, trait], axis=1).dropna()
    emmax_pheno =  args.outdir + "/" + basename(args.pheno) + ".pheno"
    outdf[[0,1,"P"]].to_csv(emmax_pheno,sep=" ",na_rep="NA",header=None, index=None)
    print("convert qqnorm file to final emmax pheno: ", emmax_pheno)
    return emmax_pheno


def genotype_out(args):
## out genotype ##
    if hasattr(args, 'keeplst'):
        keep_sample_list = args.keeplst
    else:
        keep_sample_list =  args.outdir + "/" +  "keep_" + basename(args.pheno) + ".list"
    if hasattr(args, 'outgt'):
        outgt =  args.outdir + "/" + args.outgt
    else:
        outgt =  args.outdir + "/" + basename(args.pheno)
    cmd= PLINK + ' --bfile ' + args.bfile + ' --keep ' + keep_sample_list + " --maf 0.05 --out " +  outgt + " --chr-set 27 " + " --recode12 --output-missing-genotype 0 --transpose && " + PLINK + ' --tfile ' + outgt + " --chr-set 27 --make-bed --out " + outgt
    print(cmd)
    genotype =Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
    gt_stdout,gt_stderr= genotype.communicate()
    if genotype.returncode==0:
        print('*************** genotype file done ***************')
    else:
        print(gt_stderr)
    return outgt

def makekin(args):
    if args.subcommand == 'pri_gwas' and args.kin == "BN":
        gt = args.outdir + "/" + basename(args.pheno)
        cmd = EMMAX_KIN + "  -v -d 10 "   + gt  + " -o "  + gt  + ".aBN.kinf"
    elif args.subcommand == 'pri_gwas' and args.kin == "IBS":
        cmd = EMMAX_KIN + " -s -v -d 10 " + gt  + " -o "  + gt + ".aIBS.kinf"
    elif args.subcommand == 'emmax-kin' and args.kin == "BN":
        gt = args.bfile
        cmd = EMMAX_KIN + "  -v -d 10 " + gt  +  " -o "   + args.outdir +"/" + gt  + ".aBN.kinf" 
    else:
        cmd = EMMAX_KIN + " -s -v -d 10 " + gt + " -o "  + args.outdir +"/" + gt + ".aIBS.kinf"
    print(cmd)
    kin = Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
    kin_stdout, kin_stderr= kin.communicate()
    if kin.returncode==0:
        return '*************** genotype kinship down  ***************'
    else:
        return kin_stderr

def flashpca(args):
    cmd = FLASHPCA + " --bfile " + args.bfile + " --suffix " + args.bfile + " -n 10 -m 50000"
    print(cmd)
    flashpca = Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
    pca_stdout, pca_stderr= flashpca.communicate()
    pcsFile = "pcs" + args.bfile
    pcs = pd.read_csv(pcsFile, header=None, sep="\t")
    colnames = pcs.columns.tolist()
    colnames.insert(2, "third")
    pcs = pcs.reindex(columns=colnames)
    pcs["third"] = [1 for i in range(pcs.shape[0])]
    pcs.to_csv( args.outdir + "/" + "total_pca_" + args.bfile, sep="\t",index=None, header=None)

    for i in range(args.pca,11):
        pca = i + 3
        pcas = pcs.iloc[1:, 0:pca]
        outpca = args.outdir + "/" +  args.bfile + ".pca" + str(i)
        pcas.to_csv(outpca, sep="\t", index=None, header=None)
        print(f"write {outpca} pcas")
    if flashpca.returncode==0:
        return '*************** genotype kinship down  ***************'
    else:
        return pca_stderr

def emmaxWithPca(args):
    if args.subcommand == 'pri_gwas':
        gt =   args.outdir + "/" + basename(args.pheno)
    else:
        gt = args.bfile
    if args.kin == "BN":
        kinship = gt  + ".aBN.kinf"
    elif args.kin == "IBS":
           kinship = gt  + ".aIBS.kinf"
    out =  args.outdir + "/" + basename(args.pheno)
    pheno = args.outdir + "/" + basename(args.pheno) + ".pheno"
    pca = args.pca
    cmd= EMMAX + ' -E 22000000 -v -d 10 -t '  + gt +  " -p " + pheno + " -c " + pca + " -k " + kinship + " -o " + out
    print(cmd)
    emmax=Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
    emmax_stdout,emmax_stderr= emmax.communicate()
    if emmax.returncode==0:
        return '*************** Association  Done with pca ***************'
    else:
        return emmax_stderr

def emmaxNoPca(args):
    if args.subcommand == "pri_gwas":
        gt =  args.outdir + "/" + basename(args.pheno)
    else:
        gt = args.bfile
    if args.kin == "BN":
        kinship = gt  + ".aBN.kinf"
    elif args.kin == "IBS":
        kinship = gt  + ".aIBS.kinf"
    
    out =  args.outdir + "/" + basename(args.pheno)
    pheno =  args.outdir + "/"  + basename(args.pheno) + ".pheno"
    cmd= EMMAX + ' -E 22000000 -v -d 10 -t ' + gt +  " -p " + pheno + " -k " + kinship + " -o " + out
    print(cmd)
    emmax=Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
    emmax_stdout,emmax_stderr= emmax.communicate()
    if emmax.returncode==0:
        return '*************** Association  Done with no pca ***************'
    else:
        return emmax_stderr

def emmax2qqman(args):
    """
    convert emmax out .ps file to qqman plot data
    """
    if hasattr(args, 'pheno'):
        out =  args.outdir + "/" + basename(args.pheno)
        ps = out  + ".ps"
        qqman = out + ".qqman" 
    elif args.ps:
        ps = args.ps
        qqman = ps.split(".ps")[0] + ".qqman" 
    with open(ps , 'r') as fin:
        lines = fin.readlines()
    fin.close()    
    print(f"converting emmax_out.ps to emmax_out.qqman for plot qq and manhattan figures: {qqman}")
    fileOut = open(qqman, 'w')
    for line in lines:
        newline = line.split()[0] + "\t" +  line.split()[0][2:3] + "\t" + line.split()[0][4:] + "\t" + line.split()[3]
        print(newline, file=fileOut)
    fileOut.close()
def isoform_emmax2qqman(args):
    import pandas as pd
    if hasattr(args, 'pheno'):
        out =  args.outdir + "/" + basename(args.pheno)
        ps = out  + ".ps"
        qqmanout = out + ".qqman"
    elif args.ps:
        ps = args.ps
        qqmanout = ps.split(".ps")[0] + ".qqman"
    # Use the -map argument when provided (tab-delimited isoform genetic map),
    # otherwise fall back to <outdir>/<pheno>.map next to the .ps file.
    if getattr(args, 'map', None):
        map_path = args.map
    else:
        map_path = out + ".map"
    isoform_map = pd.read_csv(map_path, header=None,sep="\t",index_col=None)
    isoform_map.columns = ["chr", "markerID", "genetic_pos", "physical_pos"]
    isoform_map.index = isoform_map["markerID"]
    psdf = pd.read_csv(ps, header=None,sep="\t",index_col=None)
    psdf.columns = ["snp", "beta", "se", "p-value"]
    psdf.index = psdf["snp"]

    outman = isoform_map.loc[psdf.index]
    outman["p-value"] = psdf["p-value"]
    outman[["chr", "physical_pos", "p-value"]].to_csv(qqmanout,header=None,sep="\t",index=True)

def qqmanplot(args):
    sig = args.sig
    if hasattr(args, 'pheno'):
        out = args.outdir + "/" +  basename(args.pheno)
        qqman =  out + ".qqman"
    elif args.qqman:
        qqman = args.qqman
    lambda_value =  qqman + ".lambda"
    cmd= QQMAN_R + ' ' +  qqman +  " " +   qqman + " " + sig +  " > " + lambda_value 
    print(cmd)
    man=Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
    man_stdout, man_stderr= man.communicate()
    if man.returncode==0:
        return '*************** manhatan and qq plot  Done ***************'
    else:
        return man_stderr


def main():
    parser = argparse.ArgumentParser(description='Perform multiple tasks in gwas writed by Derek: 13414960404@163.com')
    subparsers = parser.add_subparsers(title='subcommands', description='valid subcommands', dest='subcommand')

    task1_parser = subparsers.add_parser('pri_gwas', help='from raw phentype file to emmax gwas result(phenotype will be normlized,subgroup genotype,kinship file will be generated)')
    task1_parser.add_argument('-p',dest='pheno', type=str, help='phenotype file no columns name, fam ID at the first column, trait at the second', required=True)
    task1_parser.add_argument('-b',dest='bfile', type=str, help='big population genotype file in plink fromat to extract target subpopulation genotype file', required=True)
    task1_parser.add_argument("-sig", dest= "sig", help="significance line height which equal to -log10(1 or 0.05 /(effective snps)), caculate by gec", required=True)
    task1_parser.add_argument("-kin", dest= "kin", help="BN or IBS", default="BN")
    task1_parser.add_argument("-c", dest="pca", help="pca file format for emmax")
    task1_parser.add_argument("-o", dest="outdir",help="the output directory default is ./ ", default="./")
    task1_parser.add_argument("-t", dest="binary",help="if binary trait value 0/1 fill yes, default is no, setting to yes the input pheno will not be qqnorm", default="no")
    task1_parser.add_argument("-map", dest="map",help="the map file 4 pav gwas ")

    task2_parser = subparsers.add_parser('fine_gwas', help='from qqnorm pheno file to emmax gwas result,kinship;genotype already done before')
    task2_parser.add_argument('-p',dest='pheno',type=str, help='qqnorm of phenotype file', required=True)
    task2_parser.add_argument('-b',dest= "bfile", type=str, help='genotype file in plink format', required=True)
    task2_parser.add_argument("-c", dest="pca", help="pca file format for emmax")
    task2_parser.add_argument("-sig", dest= "sig", help="significance line height which equal to -log10(1 or 0.05 /(effective snps))")
    task2_parser.add_argument("-kin", dest= "kin", help="BN or IBS kin-ship in emmax association analysis", default="BN")
    task2_parser.add_argument("-o", dest="outdir",help="the output directory default is ./ ", default="./")
    task2_parser.add_argument("-t", dest="binary",help="if binary trait 0/1 fill 'yes',else,default is no", default="no")
    task2_parser.add_argument("-map", dest="map",help="the map file 4 pav gwas ")

    task3_parser = subparsers.add_parser('emmax-kin', help='generate emmax-kin ship matrix')
    task3_parser.add_argument('-b',dest= "bfile", type=str, help='genotype file in plink format', required=True)
    task3_parser.add_argument("-kin", dest= "kin", help="BN or IBS", default="BN")
    task3_parser.add_argument("-o", dest="outdir",help="the output directory default is ./ ", default="./")

    task6_parser = subparsers.add_parser('emmax_genotype', help='generate emmax genotyoe file')
    task6_parser.add_argument('-b',dest= "bfile", type=str, help='genotype file in plink format', required=True)
    task6_parser.add_argument('-keep', dest='keeplst',help='a fam ID list has two columns, first col is the same as second col ')
    task6_parser.add_argument('-gt', dest='outgt',help='the prefix of the output genotype file')
    task6_parser.add_argument("-o", dest="outdir",help="the output directory default is ./ ", default="./")

    task7_parser = subparsers.add_parser('pca', help='flashpca2 to perform pca analysis of genotype file')
    task7_parser.add_argument('-b',dest= "bfile", type=str, help='genotype file in plink format', required=True)
    task7_parser.add_argument('-mp', dest='pca',help='at least pca number to select for emmax analysis',default=2,type=int)
    task7_parser.add_argument("-o", dest="outdir",help="the output directory default is ./ ", default="./")

    task4_parser = subparsers.add_parser('ps2plotdata', help='to only perform qq plot and manhatan plot task we should convert emmax_result.ps file to snpID chr position pvalue format')
    task4_parser.add_argument('-ps',dest= "ps", type=str, help='emmax out result')
    task4_parser.add_argument("-o", dest="outdir",help="the output directory default is ./ ", default="./")
    task4_parser.add_argument("-map", dest="map",help="the map file 4 pav gwas ")

    task5_parser = subparsers.add_parser('plot', help='perform qq plot and manhatan plot use my.qqman file')
    task5_parser.add_argument('-qqman',dest= "qqman", type=str, help='emmax out result in my.qqman format')
    task5_parser.add_argument("-sig", dest= "sig", help="significance line height which equal to -log10(1 or 0.05 /(effective snps))",required=True)
    task5_parser.add_argument("-o", dest="outdir",help="the output directory default is ./ ", default="./")

    args = parser.parse_args()
    if args.subcommand == 'pri_gwas':
        if args.pca:
            print("gwas with pca", args.pca)
            makedir(args)
            outPheno, outkeep = keep(args)
            outnorm = qqnorm(args,outPheno)
            outgt = genotype_out(args)
            emmax_pheno = final_emmax_pheno(args,outPheno,outnorm)
            makekin(args)
            emmaxWithPca(args)
            if not args.map:
                emmax2qqman(args)
            else:
                isoform_emmax2qqman(args)

            qqmanplot(args)
        if not args.pca:
            print("gwas without pca")
            makedir(args)
            outPheno, outkeep = keep(args)
            outnorm = qqnorm(args,outPheno)
            outgt = genotype_out(args)
            emmax_pheno = final_emmax_pheno(args,outPheno,outnorm)
            makekin(args)
            emmaxNoPca(args)
            if not args.map:
                emmax2qqman(args)
            else:
                isoform_emmax2qqman(args)
            qqmanplot(args)
    elif args.subcommand == 'fine_gwas':
        if args.pca:
            print("gwas with pca", args.pca)
            makedir(args)
            outPheno, outkeep = keep(args)
            outnorm = qqnorm(args,outPheno)
            emmax_pheno = final_emmax_pheno(args,outPheno,outnorm)
            emmaxWithPca(args)
            if not args.map:
                emmax2qqman(args)
            else:
                isoform_emmax2qqman(args)
            qqmanplot(args)
        if not args.pca:
            print("gwas without pca")
            makedir(args)
            outPheno, outkeep = keep(args)
            outnorm = qqnorm(args,outPheno)
            emmax_pheno = final_emmax_pheno(args,outPheno,outnorm)
            emmaxNoPca(args)
            if not args.map:
                emmax2qqman(args)
            else:
                isoform_emmax2qqman(args)
            qqmanplot(args)
    elif args.subcommand == 'emmax-kin':
        makedir(args)
        makekin(args)
    elif args.subcommand == 'plot':
        makedir(args)
        qqmanplot(args)
    elif args.subcommand == 'ps2plotdata':
        makedir(args)
        if not args.map:
            emmax2qqman(args)
        else:
            isoform_emmax2qqman(args)

    elif args.subcommand == 'emmax_genotype':
        makedir(args)
        genotype_out(args)
    elif args.subcommand == 'pca':
        makedir(args)
        flashpca(args)
    else:
        print("not a valid subcommand")

if __name__ == '__main__':
    main()
