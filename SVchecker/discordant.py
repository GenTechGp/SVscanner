# Takes in two files containing lists ID SUPPORT READS ... as the first three columns and outputs ID SUPPORT_1 SUPPORT_2 DISCORDANT
# python3 discordant -c $supportingReads -i $inversionSupport
import sys
import pandas as pd
import argparse


def getDiscordant(row):
    inversion_reads = row['inversionREADS'].split(',') if isinstance(row['inversionREADS'], str) else []
    checker_reads = row['checkerREADS'].split(',') if isinstance(row['checkerREADS'], str) else []
    
    discordant_reads = set(inversion_reads) ^ set(checker_reads)
    
    discordant_count = len(discordant_reads)  # Count the discordant reads
    
    if discordant_reads:
        return discordant_count, ','.join(discordant_reads)
    else:
        return 0, 'NA'


if __name__ == "__main__":
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser(description="inversion",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-o",
        help="output file name [stdout]")
    parser.add_argument("-c",
        help="checker file")
    parser.add_argument("-i",
        help="SV caller file")


    args = parser.parse_args()
    inversionSupport = args.i
    checkerSupport = args.c


    # Read files into pandas dataframes
    inversionSupportDF = pd.read_csv(inversionSupport, sep='\t', usecols=[0, 1, 2, 3], header=None)
    checkerSupportDF = pd.read_csv(checkerSupport, sep='\t', usecols=[0, 1, 2], header=None)
    # Merge dataframes on ID column

    mergedDF = pd.merge(inversionSupportDF, checkerSupportDF, on=0)

    mergedDF.columns = ['ID', 'inversionSUPPORT', 'inversionREADS', 'FILTER', 'checkerSUPPORT', 'checkerREADS']

    # Calculate discordant reads
    # mergedDF['discordantREADS'] = mergedDF.apply(getDiscordant, axis=1)
    mergedDF[['discordantCount', 'discordantREADS']] = mergedDF.apply(lambda row: pd.Series(getDiscordant(row)), axis=1)
    if args.o:    
        mergedDF.to_csv(args.o, sep='\t', columns=['ID', 'inversionSUPPORT', 'checkerSUPPORT', 'discordantCount', 'discordantREADS', 'FILTER'], index=False)
    else:
        print(mergedDF.to_string(index=False))
